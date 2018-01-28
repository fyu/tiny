#include "tiny_bundle.h"

#include <unordered_map>

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <glog/logging.h>

#include "furry/common/str.h"
// #include "base/timer.h"
// #include "geo/lightfield/sfm/settings.pb.h"
// #include "geo/lightfield/sfm/solvers/project.h"
// #include "third_party/eigen3/Eigen/Eigen"
// #include "third_party/ceres/include/ceres/ceres.h"

#include "numeric.h"
#include "tiny_bundle_util.h"
// #include "util.h"

// using namespace lightfield_sfm;
using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

using ceres::CostFunction;
using ceres::AutoDiffCostFunction;
using ceres::LossFunction;
using ceres::SoftLOneLoss;
using ceres::Problem;
using ceres::Solver;
using ceres::SubsetParameterization;
using ceres::LocalParameterization;
using ceres::AutoDiffLocalParameterization;
// using ceres::Covariance;

using namespace furry;
using namespace furry::tiny;

namespace {

class AffinePoint;
class HomogeneousPoint;

// typedef AffinePoint PointUtil;
typedef HomogeneousPoint PointUtil;

}

namespace {

PointUtil *point_util; //  = new AffinePoint();

}

// class PointUtil {
//  public:
//   // template <typename T>
//   // virtual void TransformPoint(const T* const angles,
//   //                             const T* const center,
//   //                             const T* const in_point,
//   //                             T* out_point) const = 0;

//   // // return false if the point is not finite
//   // template <typename T>
//   // virtual bool ConvertTo3d(const T* const in_p, T* out_p) const = 0;

//   virtual int point_size() const = 0;

//   virtual LocalParameterization* local_parameterization() const = 0;

//   virtual void MakePoint(double x, double y, double z, double *point) const = 0;

//   virtual ~PointUtil() = 0;
// };

namespace {

class AffinePoint {
 public:

  static const int kPointSize = 3;

  template <typename T>
  void TransformPoint(const T* const angles,
                      const T* const center,
                      const T* const in_point,
                      T* out_point) const {
    T R[3][3];
    EulerAnglesToRotationMatrix(angles, R);

    T t_point[3];
    for (int i = 0; i < 3; ++i) {
      t_point[i] = in_point[i] - center[i];
    }

    for (int i = 0; i < 3; ++i) {
      out_point[i] = T(0);
      for (int j = 0; j < 3; ++j) {
        out_point[i] += R[i][j] * t_point[j];
      }
    }
  }

  // return false if the point is not finite
  template <typename T>
  bool ConvertTo3d(const T* const in_p, T* out_p) const {
    Copy<3>(in_p, out_p);
    return true;
  }

  LocalParameterization* local_parameterization() const {
    return nullptr;
  }

  void MakePoint(double x, double y, double z, double *p) const {
    p[0] = x;
    p[1] = y;
    p[2] = z;
  }

 private:
};

template <typename T>
void Normalize(const T* const in_p, T *out_p) {
  T z = in_p[0] * in_p[0] + in_p[1] * in_p[1] + in_p[2] * in_p[2] + in_p[3] * in_p[3];
  if (z > T(0)) {
    z = sqrt(z);
    out_p[0] = in_p[0] / z;
    out_p[1] = in_p[1] / z;
    out_p[2] = in_p[2] / z;
    out_p[3] = in_p[3] / z;
  } else {
    out_p[0] = T(0);
    out_p[1] = T(0);
    out_p[2] = T(0);
    out_p[3] = T(0);
  }
}

struct HomogeneousPointPlus {
  template <typename T>
  bool operator () (const T* x_, const T* delta_, T* x_plus_delta_) const {
    Matrix<T, 4, 1> x;
    Matrix<T, 4, 1> u;
    Normalize(x_, x.data());
    u[0] = T(0);
    Copy<3>(delta_, u.data() + 1);
    Matrix<T, 4, 1> e;
    e = x;
    e[0] -= T(1);
    T beta = T(2) / e.dot(e);
    Matrix<T, 4, 1> x_plus_u = u - beta * x * x.transpose() * u + x;
    x_plus_u[0] += beta * x.dot(u);
    // Copy<4>(x_plus_u.data(), x_plus_delta_);
    Normalize(x_plus_u.data(), x_plus_delta_);
    return true;
  }

  // bool operator () (const double* x_, const double* delta_, double* x_plus_delta_) const {
  //   typedef double T;
  //   Matrix<T, 4, 1> x;
  //   Matrix<T, 4, 1> u;
  //   Normalize(x_, x.data());
  //   u[0] = T(0);
  //   Copy<3>(delta_, u.data() + 1);
  //   Matrix<T, 4, 1> e;
  //   e = x;
  //   e[0] -= T(1);
  //   T beta = T(2) / e.dot(e);
  //   Matrix<T, 4, 1> x_plus_u = u - beta * x * x.transpose() * u + x;
  //   x_plus_u[0] += beta * x.dot(u);
  //   // Copy<4>(x_plus_u.data(), x_plus_delta_);
  //   Normalize(x_plus_u.data(), x_plus_delta_);
  //   LOG(INFO) << x.transpose() << ' ' << u.transpose() << ' '
  //             << x_plus_u.transpose();
  //   return true;
  // }
};

class HomogeneousPoint {
 public:

  static const int kPointSize = 4;

  HomogeneousPoint () {
    local_parameterization_ = new AutoDiffLocalParameterization<
      HomogeneousPointPlus, 4, 3>;
  }

  ~HomogeneousPoint() {
    // delete local_parameterization_;
  }

  template <typename T>
  void TransformPoint(const T* const angles,
                      const T* const center,
                      const T* const in_point,
                      T* out_point) const {
    T R[3][3];
    EulerAnglesToRotationMatrix(angles, R);

    T t_point[3];
    for (int i = 0; i < 3; ++i) {
      t_point[i] = in_point[i] - center[i] * in_point[3];
    }

    for (int i = 0; i < 3; ++i) {
      out_point[i] = T(0);
      for (int j = 0; j < 3; ++j) {
        out_point[i] += R[i][j] * t_point[j];
      }
    }
  }

  // return false if the point is not finite
  template <typename T>
  bool ConvertTo3d(const T* const in_p, T* out_p) const {
    if ((T(0) < in_p[3] && in_p[3] < T(1e-15)) ||
        (in_p[3] < T(0) && in_p[3] > T(-1e-15)))
       return false;
    out_p[0] = in_p[0] / in_p[3];
    out_p[1] = in_p[1] / in_p[3];
    out_p[2] = in_p[2] / in_p[3];
    return true;
  }

  LocalParameterization* local_parameterization() const {
    return local_parameterization_;
  }

  void MakePoint(double x, double y, double z, double *p) const {
    double ap[4]; // afine point
    ap[0] = x;
    ap[1] = y;
    ap[2] = z;
    ap[3] = 1.0;
    Normalize(ap, p);
    // LOG(INFO) << "Making " << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << p[3];
  }

 private:

  LocalParameterization *local_parameterization_;
};


template <typename T>
void TransformPoint(const T* const angles,
                    const T* const center,
                    const T* const in_point,
                    T* out_point) {
  T R[3][3];
  EulerAnglesToRotationMatrix(angles, R);

  T t_point[3];
  for (int i = 0; i < 3; ++i) {
    t_point[i] = in_point[i] - center[i];
  }

  for (int i = 0; i < 3; ++i) {
    out_point[i] = T(0);
    for (int j = 0; j < 3; ++j) {
      out_point[i] += R[i][j] * t_point[j];
    }
  }
}


struct SimpleCameraResidual {
  // u and v are the coordinates of the detected feature relative to the
  // principal point (i.e., after subtracting the principal point).
  SimpleCameraResidual(double u, double v): u_(u), v_(v) {
    CHECK(!std::isnan(u_));
    CHECK(!std::isnan(v_));
  }

  template <typename A>
  bool ApplyInternals(A const z[3],  // camera-centric 3D point.
                             A const f[2],  // focal length (x/y).
                             A const p[2],  // principal point.
                             A const k,     // lens distortion.
                             A x[2]) const {
    // Lens distortion model. Let
    //   u = z0 / z2
    //   v = z1 / z2
    // be the dehomogenized point of the calibrated camera (i.e. before
    // any intrinsics have been applied).
    A const u = z[0] / z[2];
    A const v = z[1] / z[2];

    // Lens distortion factor: lens_factor * [u, v] is a distorted point.

    // Final result, including focal length.
    x[0] = f[0] * u + p[0];
    x[1] = f[1] * v + p[1];


    if (z[2] > 0.0) {
      return true;
    } else {
      return false;
    }
  }

  // externals is 6 parameters:
  //    a) index 0:2 = euler angles
  //    b) index 3:5 = translation
  // internals is 2 parameters f and k
  // point is 3 parameters (X,Y,Z).
  // residuals are 2 values.
  template <typename T>
  bool operator()(const T* const internals,
                  const T* const externals,
                  const T* const point,
                  T* residuals) const {
    const T& focal = internals[0];
    const T& k = internals[1];

    T p[3];
    point_util->TransformPoint(externals, externals + 3, point, p);
    // camera projection.
    T f[2];
    f[0] = focal;
    f[1] = focal;
    T principal_pt[2];
    principal_pt[0] = T(0.0);
    principal_pt[1] = T(0.0);

    ApplyInternals(p, f, principal_pt, k, residuals);

    residuals[0] = residuals[0] - T(u_);
    residuals[1] = residuals[1] - T(v_);
    return true;
  }

  double u_;
  double v_;
};

class SmoothDepthPrior {
 public:
  SmoothDepthPrior(double weight) : weight_(weight) {
    // static_assert(PointUtil::kPointSize == 3, "SmoothDepthPrior only uses 3d points");
  }

  template <typename T>
  void CosAngleResidual(const T* const p0,
                        const T* const p1,
                        const T* const depth_dir,
                        T* residual) const {
    T dist = T(0);
    dist += (p0[0] - p1[0]) * (p0[0] - p1[0]);
    dist += (p0[1] - p1[1]) * (p0[1] - p1[1]);
    dist += (p0[2] - p1[2]) * (p0[2] - p1[2]);

    if (dist < T(1e-15)) {
      *residual = T(0);
    } else {
      dist = sqrt(dist);
      T point_dir[3];
      point_dir[0] = (p0[0] - p1[0]) / dist;
      point_dir[1] = (p0[1] - p1[1]) / dist;
      point_dir[2] = (p0[2] - p1[2]) / dist;
      *residual = T(0);
      *residual += point_dir[0] * depth_dir[0];
      *residual += point_dir[1] * depth_dir[1];
      *residual += point_dir[2] * depth_dir[2];
      // *residual = acos(*residual);
    }
  }

  template <typename T>
  bool operator () (const T* const p0,
                    const T* const p1,
                    T* residual) const {
    static const T depth_dir[] = {T(0), T(0), T(1)};
    T p03d[3], p13d[3];
    if (!point_util->ConvertTo3d(p0, p03d)) {
      *residual = T(0);
      return true;
    }
    if (!point_util->ConvertTo3d(p1, p13d)) {
      *residual = T(0);
      return true;
    }
    CosAngleResidual(p03d, p13d, depth_dir, residual);
    *residual *= T(weight_);
    return true;
  }
 private:
  double weight_;
};

void BuildSmoothDepthPriorResidualBlocks(
    const vector<TrackPair> &track_neighbors,
    unordered_map<track_id_t, int> &track_index,
    double smooth_depth_weight,
    double *points,
    ResidualBlocks *residual_blocks) {
  CostFunction *cost_function;
  vector<double*> param_blocks(2);
  for (auto neighbor : track_neighbors) {
    param_blocks[0] = points + PointUtil::kPointSize * track_index[neighbor.track(0)];
    param_blocks[1] = points + PointUtil::kPointSize * track_index[neighbor.track(1)];
    if (param_blocks[0] == param_blocks[1]) {
      LOG(FATAL) << "Duplicated param blocks: "
                 << neighbor.track(0) << ' ' << neighbor.track(1) << ' '
                 << track_index[neighbor.track(0)] << ' '
                 << track_index[neighbor.track(1)] << ' '
                 << track_index.size();
    }
    cost_function =
        new AutoDiffCostFunction<SmoothDepthPrior, 1, PointUtil::kPointSize, PointUtil::kPointSize>(
            new SmoothDepthPrior(neighbor.score() * smooth_depth_weight));
    residual_blocks->Add(cost_function, NULL, param_blocks);
  }
}

} // namespace

// bool CompCameras(const pair<docid_hash_t, CameraMapValue>& c1,
//                  const pair<docid_hash_t, CameraMapValue>& c2) {
//   return global_docid_hash_to_filepath[c1.first] <
//       global_docid_hash_to_filepath[c2.first];
// }

// vector<pair<docid_hash_t, CameraMapValue>>
// GetSortedCameras(const Model& model) {
//   vector<pair<docid_hash_t, CameraMapValue>> cameras;
//   auto& camera_map = model.GetCameras();
//   for (auto& camera : camera_map) {
//     cameras.push_back(camera);
//   }
//   sort(cameras.begin(), cameras.end(), CompCameras);
//   return cameras;
// }

namespace furry {
namespace tiny {

bool Bundle(const tiny::BundleSettings& settings,
            const vector<TrackPair> *track_neighbors,
            Model* model) {
  point_util = new PointUtil;
  const int num_cameras = model->NumCameras();
  const int num_tracks = model->NumTracks();

  const int intrinsics_size = 2;
  double intrinsics[intrinsics_size]; //f and k
  intrinsics[0] = settings.f_mean();
  intrinsics[1] = settings.k_mean();

  // hash_map<docid_hash_t, int> camera_index;
  unordered_map<track_id_t, int> track_index;
  vector<ResidualBlocks> all_residual_blocks;

  const int camera_block_size = 6;
  const int point_block_size = PointUtil::kPointSize;
  const int num_parameters =
      num_cameras * camera_block_size + num_tracks * point_block_size;
  vector<double> initial_parameters(num_parameters);
  double* x = initial_parameters.data();

  Problem problem;
  CostFunction* cost_function;

  cost_function = new AutoDiffCostFunction<CameraPrior, 2, 2>(
      new CameraPrior(settings));
  static const int fix_r_constants[] = {0, 1, 2};
  static const int fix_z_constants[] = {5};
  static const int fix_r_z_constants[] = {0, 1, 2, 5};
  static const int fix_t_constants[] = {3, 4, 5};
  SubsetParameterization* camera_parameterization = nullptr;

  if (settings.zero_rotation() && settings.zero_depth()) {
    camera_parameterization =
        new SubsetParameterization(6, vector<int>(fix_r_z_constants,
                                                  fix_r_z_constants + 4));
  } else if (settings.zero_depth()) {
    camera_parameterization =
        new SubsetParameterization(6, vector<int>(fix_z_constants,
                                                  fix_z_constants + 1));
  } else if (settings.zero_rotation()) {
    camera_parameterization =
        new SubsetParameterization(6, vector<int>(fix_r_constants,
                                                  fix_r_constants + 3));
  } else if (settings.zero_translation()) {
    camera_parameterization =
        new SubsetParameterization(6, vector<int>(fix_t_constants,
                                                  fix_t_constants + 3));
  }

  problem.AddResidualBlock(cost_function, NULL, intrinsics);

  // Initialize camera parameters
  // const CameraMap& cameras = model->GetCameras();
  // auto cameras = GetSortedCameras(*model);
  int camera_count = 0;
  int fixed_camera_count = 0;
  int movable_camera_count = 0;
  problem.AddParameterBlock(intrinsics, intrinsics_size);
  if (settings.fix_internals()) {
    problem.SetParameterBlockConstant(intrinsics);
  }
  ResidualBlocks angle_residual_blocks("Angles");
  // for (auto it = cameras.begin(); it != cameras.end(); ++it) {
  for (int i_camera = 0; i_camera < model->NumCameras(); ++i_camera) {
    // LOG(INFO) << global_docid_hash_to_filepath[it->first];
    // const docid_hash_t docid_hash = it->first;
    auto camera = model->GetCamera(i_camera);
    double* camera_params = x + camera_block_size * camera_count;
    // camera_index[docid_hash] = camera_count++;
    ++camera_count;

    // Initialize camera parameters.
    InitializeCameraBlock(*camera, camera_params);
    problem.AddParameterBlock(camera_params, camera_block_size);

    // CostFunction* cost_function =
    //     CreateCameraPrior(it->second.camera, settings);
    // problem.AddResidualBlock(cost_function, NULL, camera_params);
    if (settings.angle_std() > 0) {
      cost_function =
          new AutoDiffCostFunction<VectorPrior<3>, 3, 6>(
              new VectorPrior<3>(0, settings.angle_mean(),
                                 settings.angle_std()));
      angle_residual_blocks.Add(cost_function, NULL,
                                      vector<double*>(1, camera_params));
      // problem.AddResidualBlock(cost_function, NULL, camera_params);
    }

    // if (!what.OptimiseCamera(docid_hash)) {
    if (settings.fix_externals()) {
      ++fixed_camera_count;
      problem.SetParameterBlockConstant(camera_params);
    } else {
      if (camera_parameterization) {
        problem.SetParameterization(camera_params, camera_parameterization);
      }
      ++movable_camera_count;
    }
  }
  // LOG(INFO) << "Fixed: " << fixed_camera_count;
  all_residual_blocks.push_back(angle_residual_blocks);

  if (settings.smooth_camera_trajectory()) {
    // double* p0;
    // double* p1;
    // double* p2;
    vector<double*> p(3);
    ResidualBlocks residual_blocks("Smooth Camera Trajectory");
    for (int i = 1; i < model->NumCameras() - 1; ++i) {
      // p0 = CameraCenterPtr(x + camera_block_size * (i - 1));
      // p1 = CameraCenterPtr(x + camera_block_size * i);
      // p2 = CameraCenterPtr(x + camera_block_size * (i + 1));
      p[0] = x + camera_block_size * (i - 1);
      p[1] = x + camera_block_size * i;
      p[2] = x + camera_block_size * (i + 1);
      cost_function =
          new AutoDiffCostFunction<SmoothCurvePrior, 3, 6, 6, 6>(
              new SmoothCurvePrior(
                  3,settings.smooth_camera_trajectory_weight()));
      residual_blocks.Add(cost_function, NULL, p);
      // problem.AddResidualBlock(cost_function, NULL, p0, p1, p2);
    }
    all_residual_blocks.push_back(residual_blocks);
  }

  if (settings.angle_adj_std() > 0) {
    // double* r0;
    // double* r1;
    vector<double*> r(2);
    ResidualBlocks residual_blocks("Adjacent Angles");
    for (int i = 0; i < model->NumCameras() - 1; ++i) {
      r[0] = CameraRotationPtr(x + camera_block_size * i);
      r[1] = CameraRotationPtr(x + camera_block_size * (i + 1));
      cost_function =
          new AutoDiffCostFunction<
            VectorDiffPrior<3>, 3, camera_block_size, camera_block_size>(
                new VectorDiffPrior<3>(0, settings.angle_adj_mean(),
                                       settings.angle_adj_std()));
      // problem.AddResidualBlock(cost_function, NULL, r0, r1);
      residual_blocks.Add(cost_function, NULL, r);
    }
    all_residual_blocks.push_back(residual_blocks);
  }

  // CHECK_EQ(camera_count, num_cameras);
  // CHECK_LE(movable_camera_count, what.num_cameras());

  vector<double*> parameters(3);
  parameters[0] = intrinsics;

  LossFunction *loss_function = nullptr;
  if (settings.use_robustifier()) {
    loss_function = new SoftLOneLoss(1.0);
  }

  // const TrackMap& tracks = model->GetTracks();
  int track_count = 0;
  int fixed_track_count = 0;
  // Iterate over tracks.
  ResidualBlocks projection_residual_blocks("Projection");
  // for (auto it = tracks.begin(); it != tracks.end(); ++it) {
  for (int i_track = 0; i_track < model->NumTracks(); ++i_track) {
    auto track = model->GetTrack(i_track);
    auto &point_3d = track->GetPoint();
    // const Track& track = it->second.track;
    // const Point& point = it->second.point;
    // rushmore::numeric::VectorD3 point_3d = point.position();

    track_index[i_track] = track_count;

    parameters[2] = x + camera_block_size * num_cameras +
        point_block_size * track_count;
    point_util->MakePoint(point_3d[0], point_3d[1], point_3d[2], parameters[2]);

    problem.AddParameterBlock(parameters[2], point_block_size,
                              point_util->local_parameterization());

    // Iterate over each camera projection of the track.
    for (int i = 0; i < track->NumCameras(); ++i) {
      // const docid_hash_t docid_hash = track[i].docid_hash;
      // hash_map<docid_hash_t, int>::const_iterator cam_id =
      //     camera_index.find(docid_hash);
      // CHECK(cam_id != camera_index.end());
      auto camera = model->GetCameraById(track->GetCameraId(i));
      auto &detection = camera->GetDetectionOfTrack(track->GetId());

      parameters[1] = x + camera_block_size * track->GetCameraId(i);

      double detection_u, detection_v;
      CostFunction* cost_function;

      // rushmore::numeric::VectorD2 pp =
      //     model->GetCamera(docid_hash).principalPoint();
      // TODO(fyu) set principal point position mode carefully
      auto pp = camera->GetPrincipalPoint();
      detection_u = detection.x - pp[0];
      detection_v = detection.y - pp[1];
      // detection_u = detection.x;
      // detection_v = detection.y;
      cost_function =
          new AutoDiffCostFunction <SimpleCameraResidual, 2, 2,
                                    camera_block_size, point_block_size>(
              new SimpleCameraResidual(detection_u, detection_v));
      projection_residual_blocks.Add(cost_function, loss_function, parameters);
    }

    // if (!what.OptimisePoint(it->first)) {
    if (settings.fix_tracks()) {
      problem.SetParameterBlockConstant(parameters[2]);
      ++fixed_track_count;
    }

    ++track_count;
  }
  all_residual_blocks.push_back(projection_residual_blocks);

  if (settings.with_smooth_depth()) {
    if (track_neighbors == nullptr) {
      LOG(ERROR) << "Don't have track neighbors for smooth depth prior";
    }
    ResidualBlocks residual_blocks("Smooth Depth Prior");
    BuildSmoothDepthPriorResidualBlocks(
        *track_neighbors,
        track_index,
        settings.smooth_depth().weight(),
        x + num_cameras * camera_block_size,
        &residual_blocks);
    all_residual_blocks.push_back(residual_blocks);
  }


  for (auto& blocks : all_residual_blocks) {
    blocks.AddTo(&problem);
  }

  CHECK_EQ(track_count, num_tracks);
  LOG(INFO) << "cameras: " << num_cameras << " fixed: " << fixed_camera_count
            << " tracks: " << num_tracks << " fixed: " << fixed_track_count;
  const int delta_cameras = num_cameras - fixed_camera_count;

  auto solver_options = settings.solver_options();
  Solver::Options options;
  options.max_num_iterations = solver_options.max_num_iterations();
  options.num_threads = solver_options.num_threads();
  options.num_linear_solver_threads = solver_options.num_threads();
  // options.linear_solver_type = ceres::SPARSE_SCHUR;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout =
      solver_options.minimizer_progress_to_stdout();
  options.use_inner_iterations = solver_options.use_inner_iterations();

  // Too large a trust region radius makes the linear solver
  // cranky.
  options.max_trust_region_radius = 1e7;

  if (delta_cameras > 800) {
    options.function_tolerance = 1e-6;
    options.gradient_tolerance = 1e-6;
    options.parameter_tolerance = 1e-6;
    options.linear_solver_type = ceres::ITERATIVE_SCHUR;
    options.preconditioner_type = ceres::SCHUR_JACOBI;
    LOG(INFO) << "Switching to ITERATIVE_SCHUR: SCHUR_JACOBI";
  }

  // WallTimer timer;
  // timer.Start();

  Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  // timer.Stop();
  if (solver_options.print_summary()) {
    cout << summary.FullReport() << '\n';
  }
  // LOG(INFO) << "Ceres solver time: " << timer.Get()
  //           << ", termination type = "  << summary.termination_type;

  // if (settings.calc_covariance()) {
  //   Covariance::Options cov_options;
  //   Covariance covariance(cov_options);
  //   vector<pair<const double*, const double*> > covariance_blocks;

  //   camera_count = 0;
  //   for (auto it = cameras.begin(); it != cameras.end(); ++it) {
  //     // const docid_hash_t docid_hash = it->first;
  //     if (what.OptimiseCamera(it->first)) {
  //       double* camera_params = x + camera_block_size * camera_count;
  //       covariance_blocks.push_back(make_pair(camera_params, camera_params));
  //     }
  //     ++camera_count;
  //   }
  //   track_count = 0;
  //   for (auto it = tracks.begin(); it != tracks.end(); ++it) {
  //     if (what.OptimisePoint(it->first)) {
  //       double* track_params = x + camera_block_size * num_cameras +
  //           point_block_size * track_count;
  //       covariance_blocks.push_back(make_pair(track_params, track_params));
  //     }
  //     ++track_count;
  //   }
  //   // CHECK(covariance.Compute(covariance_blocks, &problem));
  //   covariance.Compute(covariance_blocks, &problem);
  //   vector<Matrix<double, 6, 6, RowMajor>> cov_cameras(cameras.size());
  //   vector<Matrix<double, 3, 3, RowMajor>> cov_points(tracks.size());
  //   camera_count = 0;
  //   for (auto it = cameras.begin(); it != cameras.end(); ++it) {
  //     // const docid_hash_t docid_hash = it->first;
  //     double* camera_params = x + camera_block_size * camera_count;
  //     // covariance_blocks.push_back(make_pair(camera_params, camera_params));
  //     covariance.GetCovarianceBlock(camera_params, camera_params,
  //                                   cov_cameras[camera_count].data());
  //     ++camera_count;
  //   }
  //   track_count = 0;
  //   for (auto it = tracks.begin(); it != tracks.end(); ++it) {
  //     double* track_params = x + camera_block_size * num_cameras +
  //         point_block_size * track_count;
  //     covariance.GetCovarianceBlock(track_params, track_params,
  //     cov_points[track_count].data());
  //     ++track_count;
  //   }

  //   LOG(INFO) << "Angle Uncertainty";
  //   camera_count = 0;
  //   for (auto it = cameras.begin(); it != cameras.end(); ++it) {
  //     if (what.OptimiseCamera(it->first)) {
  //       LOG(INFO) << cov_cameras[camera_count].diagonal().transpose();
  //     }
  //     ++camera_count;
  //   }
  //   LOG(INFO) << "Point Uncertainty";
  //   track_count = 0;
  //   for (auto it = tracks.begin(); it != tracks.end(); ++it) {
  //     if (what.OptimisePoint(it->first)) {
  //       LOG(INFO) << cov_points[track_count].diagonal().transpose();
  //     }
  //     ++track_count;
  //   }
  // }

  // Copy solved camera pose back to the model.
  camera_count = 0;
  double quat[4];
  double *center;
  // double f = intrinsics[0];
  // double k = intrinsics[1];
  // for (auto it = cameras.begin(); it != cameras.end(); ++it) {
  for (int i = 0; i < model->NumCameras(); ++i) {
    // docid_hash_t docid_hash = it->first;
    auto camera = model->GetCamera(i);

    double* camera_params = x + camera_block_size * camera_count;
    // LOG(INFO) << GetFilename(camera.GetDebugUrl()) << " angles: "
    //           << camera_params[0] << ' '
    //           << camera_params[1] << ' '
    //           << camera_params[2];
    center = camera_params + 3;
    EulerAnglesToQuaternion(camera_params, quat);
    // rushmore::numeric::VectorD2 pp = camera.principalPoint();
    // camera.SetFromParams(camera.type(), f, pp.data, k, quat, center);
    camera->SetCenter(center);
    camera->SetQuaternion(quat);
    camera_count++;
  }

  for (auto& blocks : all_residual_blocks) {
    LOG(INFO) << blocks.Info();
  }

  CHECK_EQ(camera_count, num_cameras);

  // Copy solved 3d locations of tracks back to the model.
  track_count = 0;
  // for (TrackMap::const_iterator it = tracks.begin(); it != tracks.end();
  //      ++it) {
  for (int i = 0; i < model->NumTracks(); ++i) {
    auto track = model->GetTrack(i);
    double* raw_point = x + camera_block_size * num_cameras +
        point_block_size * track_count;
    double p[3];
    if (point_util->ConvertTo3d(raw_point, p)) {
      // point.SetElements(p);
      track->SetPoint(p);
    } else {
      LOG(INFO) << "Infinite point " << raw_point[0] << ' ' << raw_point[1] << ' ' << raw_point[2] << ' ' << raw_point[3];;
    }
    ++track_count;
  }
  CHECK_EQ(track_count, num_tracks);

  delete point_util;

  return true;
}


} // tiny
} // furry

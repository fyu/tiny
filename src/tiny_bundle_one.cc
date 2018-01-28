#include "tiny_bundle_one.h"

#include <unordered_map>

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <glog/logging.h>

#include "furry/common/str.h"

#include "numeric.h"
#include "tiny_bundle_util.h"

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
using ceres::Covariance;

namespace {

class AffinePointOne;
class HomogeneousPointOne;

// typedef AffinePointOne PointParameterization;
typedef HomogeneousPointOne PointParameterization;

class HomogeneousPointOne {
 public:

  static const int kPointSize = 1;

  static HomogeneousPointOne* NewInstance(const tiny::BundleSettings& settings,
                               const Model &model,
                               int track_id) {
    auto detection = model.GetCamera(settings.ref_camera_index())
        ->GetDetectionInRetinaCoord(track_id);
    return new HomogeneousPointOne(detection.x, detection.y);
  }

  HomogeneousPointOne(double x, double y)
      : ref_x_(x), ref_y_(y) {}

  template <typename T>
  void TransformPoint(const T* const angles,
                      const T* const center,
                      const T* const in_point,
                      T* out_point) const {
    T R[3][3];
    EulerAnglesToRotationMatrix(angles, R);

    T t_point[3];
    t_point[0] = T(ref_x_) - center[0] * *in_point;
    t_point[1] = T(ref_y_) - center[1] * *in_point;
    t_point[2] = T(1)      - center[2] * *in_point;

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
    if ((T(0) < in_p[0] && in_p[0] < T(1e-15)) ||
        (in_p[0] < T(0) && in_p[0] > T(-1e-15)))
       return false;
    out_p[0] = T(ref_x_) / *in_p;
    out_p[1] = T(ref_y_) / *in_p;
    out_p[2] = T(1)      / *in_p;
    return true;
  }

  LocalParameterization* local_parameterization() const {
    return nullptr;
  }

  void MakePoint(double x, double y, double z, double *p) const {
    *p = 1 / z;
  }

  int size() const {
    return kPointSize;
  }

 private:
  double ref_x_;
  double ref_y_;
};

class AffinePointOne {
 public:

  static const int kPointSize = 1;

  static AffinePointOne* NewInstance(const tiny::BundleSettings& settings,
                                     const Model &model,
                                     int track_id) {
    auto detection = model.GetCamera(settings.ref_camera_index())
        ->GetDetectionInRetinaCoord(track_id);
    // VLOG(1) << detection;
    return new AffinePointOne(detection.x, detection.y);
  }

  AffinePointOne(double x, double y)
      : ref_x_(x), ref_y_(y) {}

  template <typename T>
  void TransformPoint(const T* const angles,
                      const T* const center,
                      const T* const in_point,
                      T* out_point) const {
    T R[3][3];
    EulerAnglesToRotationMatrix(angles, R);

    T t_point[3];
    t_point[0] = T(ref_x_) * *in_point - center[0];
    t_point[1] = T(ref_y_) * *in_point - center[1];
    t_point[2] = T(1)      * *in_point - center[2];

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
    out_p[0] = T(ref_x_) * *in_p;
    out_p[1] = T(ref_y_) * *in_p;
    out_p[2] = T(1)      * *in_p;
    return true;
  }

  LocalParameterization* local_parameterization() const {
    return nullptr;
  }

  void MakePoint(double x, double y, double z, double *p) const {
    *p = z;
  }

  int size() const {
    return kPointSize;
  }

 private:
  double ref_x_;
  double ref_y_;
};

struct SimpleCameraResidual {
  // u and v are the coordinates of the detected feature relative to the
  // principal point (i.e., after subtracting the principal point).
  SimpleCameraResidual(double u, double v,
                       const PointParameterization *point_p)
      : u_(u), v_(v), point_p_(point_p) {
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
    point_p_->TransformPoint(externals, externals + 3, point, p);
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
  const PointParameterization *point_p_;
};

class SmoothDepthPrior {
 public:
  SmoothDepthPrior(double weight,
                   PointParameterization *point0_p,
                   PointParameterization *point1_p)
      : weight_(weight), point0_p_(point0_p), point1_p_(point1_p) {
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
    if (!point0_p_->ConvertTo3d(p0, p03d)) {
      *residual = T(0);
      return true;
    }
    if (!point1_p_->ConvertTo3d(p1, p13d)) {
      *residual = T(0);
      return true;
    }
    CosAngleResidual(p03d, p13d, depth_dir, residual);
    *residual *= T(weight_);
    return true;
  }
 private:
  double weight_;
  PointParameterization *point0_p_;
  PointParameterization *point1_p_;
};

void BuildSmoothDepthPriorResidualBlocks(
    const vector<TrackPair> &track_neighbors,
    const vector<PointParameterization*> &point_ps,
    unordered_map<track_id_t, int> &track_index,
    double smooth_depth_weight,
    double *points,
    ResidualBlocks *residual_blocks) {
  CostFunction *cost_function;
  vector<double*> param_blocks(2);
  for (auto neighbor : track_neighbors) {
    auto point0_p = point_ps[track_index[neighbor.track(0)]];
    auto point1_p = point_ps[track_index[neighbor.track(1)]];
    param_blocks[0] = points + point0_p->size() * track_index[neighbor.track(0)];
    param_blocks[1] = points + point1_p->size() * track_index[neighbor.track(1)];
    if (param_blocks[0] == param_blocks[1]) {
      LOG(FATAL) << "Duplicated param blocks: "
                 << neighbor.track(0) << ' ' << neighbor.track(1) << ' '
                 << track_index[neighbor.track(0)] << ' '
                 << track_index[neighbor.track(1)] << ' '
                 << track_index.size();
    }
    cost_function =
        new AutoDiffCostFunction<SmoothDepthPrior, 1,
                                 PointParameterization::kPointSize,
                                 PointParameterization::kPointSize>(
            new SmoothDepthPrior(neighbor.score() * smooth_depth_weight,
                                 point0_p, point1_p));
    residual_blocks->Add(cost_function, NULL, param_blocks);
  }
}

} // namespace

namespace furry {
namespace tiny {

bool BundleOne(const tiny::BundleSettings& settings,
               const vector<TrackPair> *track_neighbors,
               Model* model) {

  // for (int i = 0; i < model->NumTracks(); ++i) {
  //   cout << model->GetTrack(i)->GetPoint().transpose() << '\n';
  // }
  const int num_cameras = model->NumCameras();
  const int num_tracks = model->NumTracks();

  const int intrinsics_size = 2;
  double intrinsics[intrinsics_size]; //f and k
  intrinsics[0] = settings.f_mean();
  intrinsics[1] = settings.k_mean();

  unordered_map<camera_id_t, int> camera_index;
  unordered_map<track_id_t, int> track_index;
  vector<ResidualBlocks> all_residual_blocks;

  const int camera_block_size = 6;
  const int point_block_size = PointParameterization::kPointSize;
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
    camera_index[camera->GetId()] = camera_count++;
    // ++camera_count;

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
    if (settings.fix_externals() || i_camera == settings.ref_camera_index()) {
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

  LossFunction* loss_function = nullptr;
  if (settings.use_robustifier()) {
    loss_function = new SoftLOneLoss(1.0);
  }
  // LossFunction* loss_function = NULL;

  int track_count = 0;
  int fixed_track_count = 0;

  // Initialize point parameterization
  vector<PointParameterization*> point_ps(model->NumTracks());
  for (int i = 0; i < model->NumTracks(); ++i) {
    point_ps[i] = PointParameterization::NewInstance(settings, *model, i);
  }

  // Iterate over tracks.
  ResidualBlocks projection_residual_blocks("Projection");
  // for (auto it = tracks.begin(); it != tracks.end(); ++it) {
  for (int i_track = 0; i_track < model->NumTracks(); ++i_track) {
    auto track = model->GetTrack(i_track);
    auto &point_3d = track->GetPoint();

    auto point_p = point_ps[i_track];

    track_index[i_track] = track_count;

    parameters[2] = x + camera_block_size * num_cameras +
        point_block_size * track_count;
    point_p->MakePoint(point_3d[0], point_3d[1], point_3d[2], parameters[2]);

    problem.AddParameterBlock(parameters[2], point_block_size,
                              point_p->local_parameterization());

    // Iterate over each camera projection of the track.
    for (int i = 0; i < track->NumCameras(); ++i) {
      // const docid_hash_t docid_hash = track[i].docid_hash;
      // hash_map<docid_hash_t, int>::const_iterator cam_id =
      //     camera_index.find(docid_hash);
      // CHECK(cam_id != camera_index.end());
      auto camera = model->GetCameraById(track->GetCameraId(i));
      auto &detection = camera->GetDetectionOfTrack(track->GetId());

      parameters[1] = x + camera_block_size * camera_index[camera->GetId()];

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
          new AutoDiffCostFunction <
            SimpleCameraResidual, 2, 2, camera_block_size, point_block_size>(
                new SimpleCameraResidual(detection_u, detection_v, point_p));
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
        point_ps,
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
  // options.function_tolerance = 1e-15;
  // options.gradient_tolerance = 1e-15;
  // options.parameter_tolerance = 1e-15;

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

  if (settings.calc_covariance()) {
    Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::SPARSE_QR;
    cov_options.num_threads = solver_options.num_threads();
    cov_options.apply_loss_function = false;
    vector<pair<const double*, const double*>> covariance_blocks;
    double *track_parameters = x + camera_block_size * num_cameras;
    for (int i = 0; i < model->NumTracks(); ++i) {
      // auto track = model->GetTrack(i);
      double *p = track_parameters + point_block_size * track_index.at(i);
      covariance_blocks.push_back(make_pair(p, p));
    }

    // Hold camera pose parameters
    for (int i = 0; i < model->NumCameras(); ++i) {
      problem.SetParameterBlockConstant(x + camera_block_size * i);
    }

    Covariance covariance(cov_options);
    if (!covariance.Compute(covariance_blocks, &problem)) {
      LOG(WARNING) << "Fail to compute covariance";
    } else {
      vector<Matrix<double, PointParameterization::kPointSize, PointParameterization::kPointSize>> track_covariance(model->NumTracks());
      for (int  i = 0; i < model->NumTracks(); ++i) {
        double *p = track_parameters + point_block_size * track_index.at(i);
        covariance.GetCovarianceBlock(p, p, track_covariance[i].data());
      }
      LOG(INFO) << "Writing covariance of points to model";
      for (int i = 0; i < model->NumTracks(); ++i) {
        model->GetTrack(i)->SetCovariance(
            vector<double>(track_covariance[i].data(),
                           track_covariance[i].data() +
                           track_covariance[i].rows() *
                           track_covariance[i].cols()));
      }
    }
  }

  // Copy solved camera pose back to the model.
  camera_count = 0;
  double quat[4];
  double *center;

  double f[2];
  double k[2] = {0, 0};
  f[0] = f[1] = intrinsics[0];
  k[0] = intrinsics[1];

  for (int i = 0; i < model->NumCameras(); ++i) {
    // docid_hash_t docid_hash = it->first;
    auto camera = model->GetCamera(i);

    double* camera_params = x + camera_block_size * camera_count;

    center = camera_params + 3;
    EulerAnglesToQuaternion(camera_params, quat);

    camera->SetCenter(center);
    camera->SetQuaternion(quat);

    camera->SetInternals(f, k);

    camera_count++;
  }

  for (auto& blocks : all_residual_blocks) {
    LOG(INFO) << blocks.Info();
  }

  CHECK_EQ(camera_count, num_cameras);

  // Copy solved 3d locations of tracks back to the model.
  track_count = 0;

  for (int i = 0; i < model->NumTracks(); ++i) {
    auto track = model->GetTrack(i);
    double* raw_point = x + camera_block_size * num_cameras +
        point_block_size * track_count;
    double p[3];
    if (point_ps[i]->ConvertTo3d(raw_point, p)) {
      // point.SetElements(p);
      track->SetPoint(p);
    } else {
      LOG(INFO) << "Infinite point " << raw_point[0] << ' ' << raw_point[1] << ' ' << raw_point[2] << ' ' << raw_point[3];;
    }
    ++track_count;
  }
  CHECK_EQ(track_count, num_tracks);

  // delete point_util;

  return true;
}

} // tiny
} // furry

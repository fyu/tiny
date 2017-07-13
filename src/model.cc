#include "model.h"

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <dirent.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <regex>
#include <sstream>

#include <glog/logging.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/common/cv.h"
#include "furry/common/path.h"
#include "furry/common/str.h"

using namespace Eigen;
using namespace cv;
using namespace std;

namespace {

// // Using the description in http://en.wikipedia.org/wiki/Laguerre's_method
// double QuinticLaguerreSolver(const double a_0, const double a_1,
//                              const double a_3, const double a_5,
//                              const double x_0, const double epsilon,
//                              const int maxIter) {
//   double x = x_0;
//   for (int k = 0; k < maxIter; ++k) {
//     const double x2 = x * x;
//     const double pk = a_0 + x * (a_1 + x2 * (a_3 + a_5 * x2));

//     if (fabs(pk) < 1e-10)
//       return x;

//     const double dpk = a_1 + x2 * (3.0 * a_3 + 5.0 * a_5 * x2);
//     const double ddpk = x * (6 * a_3 + 20 * a_5 * x2);
//     const double g = dpk / pk;
//     const double h = g * g - ddpk / pk;
//     const double j = sqrt(fabs(4. * (5. * h - g * g)));

//     const double deltax = (g > 0) ? 5.0 / (g + j) : 5.0 / (g - j);
//     if (fabs(deltax) < epsilon)
//       break;
//     x = x - deltax;
//   }
//   return x;
// }

// void RadialDistort(const double u,   const double v,
//                    const double c_u, const double c_v,
//                    const double f_u, const double f_v,
//                    const double k_1, const double k_2,
//                    double* up,       double* vp) {
//   const double x  = (u - c_u) / f_u;
//   const double y  = (v - c_v) / f_v;
//   const double r2 = x * x + y * y;
//   const double d  = 1 + r2 * (k_1 + k_2 * r2);
//   (*up) = c_u + (u-c_u) * d;
//   (*vp) = c_v + (v-c_v) * d;
// }

// void RadialUndistort(const double up,  const double vp,
//                      const double c_u, const double c_v,
//                      const double f_u, const double f_v,
//                      const double k_1, const double k_2,
//                      double* u, double* v) {
//   const double kMinDistanceFromCenter = 1e-6;
//   const double kMaxIter = 10;
//   const double kEpsilon = 1e-6;

//   const double u1 = (up - c_u) / f_u;
//   const double v1 = (vp - c_v) / f_v;

//   if (fmax(fabs(u1), fabs(v1)) < kMinDistanceFromCenter) {
//     (*u) = up;
//     (*v) = vp;
//     return;
//   }

//   double x = u1;
//   double y = v1;
//   if (fabs(u1) > fabs(v1)) {
//     const double alpha = v1 / u1;
//     const double beta = 1 + alpha * alpha;
//     const double a_0 = - u1;
//     const double a_1 = 1;
//     const double a_3 = k_1 * beta;
//     const double a_5 = k_2 * beta * beta;
//     x = QuinticLaguerreSolver(a_0, a_1, a_3, a_5, u1, kEpsilon, kMaxIter);
//     y = alpha * x;
//   } else {
//     const double alpha = u1 / v1;
//     const double beta = 1 + alpha * alpha;
//     const double a_0 = - v1;
//     const double a_1 = 1;
//     const double a_3 = k_1 * beta;
//     const double a_5 = k_2 * beta * beta;
//     y = QuinticLaguerreSolver(a_0, a_1, a_3, a_5, v1, kEpsilon, kMaxIter);
//     x = alpha * y;
//   }

//   (*u) = f_u * x + c_u;
//   (*v) = f_v * y + c_v;
// }

// string AddDropboxPath(const string &path) {
//   string fullpath = "/Users/fy/Dropbox/" + path;
//   if (furry::Exists(fullpath)) {
//     return fullpath;
//   }
//   fullpath = "/n/fs/graphicslab/fy/Dropbox/" + path;
//   if (furry::Exists(fullpath)) {
//     return fullpath;
//   }
//   return path;
// }

string AddNfsPath(const string &path) {
  string fullpath = "/Volumes/" + path;
  if (furry::Exists(fullpath)) {
    return fullpath;
  }
  fullpath = "/n/fs/" + path;
  if (furry::Exists(fullpath)) {
    return fullpath;
  }
  return path;
}

string RemoveNfsPath(const string &path) {
  size_t it;
  if ((it = path.find("/n/fs/")) != string::npos && it == 0) {
    return path.substr(it + strlen("/n/fs/"));
  }
  if ((it = path.find("/Volumes/graphicslab")) != string::npos && it == 0) {
    return path.substr(it + strlen("/Volumes/"));
  }
  return path;
}

}

namespace {

Eigen::Vector3d TriangulateMidpoint(const std::vector<Eigen::Matrix3d> &rs,
                                    const std::vector<Eigen::Vector3d> &cs,
                                    const std::vector<Eigen::Matrix3d> &ks,
                                    const std::vector<Eigen::Vector2d> &ps) {
  int num_images = rs.size();
  std::vector<Eigen::Vector3d> us;
  us.reserve(num_images);
  for (int i = 0; i < num_images; ++i)
    us.push_back(ks[i].inverse() * ps[i].homogeneous());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_images * 3, num_images);
  Eigen::Matrix<double, 3, Eigen::Dynamic> diag_A(3, num_images);
  Eigen::VectorXd b(num_images * 3);
  Eigen::Vector3d sum_cs = Eigen::Vector3d::Zero();
  for (int i = 0; i < num_images; ++i) {
    A.block<3, 1>(i * 3, i) = rs[i] * us[i] / us[i][2];
    int j = (i + 1) % num_images;
    A.block<3, 1>(i * 3, j) = -rs[j] * us[j] / us[j][2];
    b.segment<3>(i * 3) = cs[j] - cs[i];
    sum_cs += cs[i];
    diag_A.col(i) = A.block<3, 1>(i * 3, i);
  }
  Eigen::VectorXd a = A.fullPivLu().solve(b);
  return (diag_A * a + sum_cs) / num_images;
}

}

namespace furry {
namespace tiny {

Matrix3d Camera::Internals::K() const {
  Matrix3d k = Matrix3d::Identity();
  k(0, 0) = f[0];
  k(1, 1) = f[1];
  k(0, 2) = p[0];
  k(1, 2) = p[1];
  k(0, 1) = skew;
  return k;
}

Eigen::Quaterniond Camera::GetQuaternion() const {
  return Quaterniond(externals_.q[0], externals_.q[1],
                     externals_.q[2], externals_.q[3]);
}

void Camera::SetQuaternion(const Eigen::Quaterniond &q) {
  externals_.q[0] = q.w();
  externals_.q[1] = q.x();
  externals_.q[2] = q.y();
  externals_.q[3] = q.z();
}


void Camera::Internals::SetParams(double* f_, double* k_, double* p_, double s_) {
  f[0] = f_[0];
  f[1] = f_[1];
  k[0] = k_[0];
  k[1] = k_[1];
  p[0] = p_[0];
  p[1] = p_[1];
  skew = s_;
}

void Camera::Internals::fprint(FILE* file) const {
  fprintf(file, "%0.6lf %0.6lf %0.6lf %0.6lf %0.6lf %0.6lf %0.6lf",
          k[0], k[1], f[0], f[1], p[0], p[1], skew);
}

Matrix3d Camera::Externals::R() const {
  return Quaterniond(q[0], q[1], q[2], q[3]).matrix();
}

Vector3d Camera::Externals::C() const {
  return c;
}

void Camera::Externals::SetParams(double* q_, double* c_) {
  q[0] = q_[0];
  q[1] = q_[1];
  q[2] = q_[2];
  q[3] = q_[3];
  c[0] = c_[0];
  c[1] = c_[1];
  c[2] = c_[2];
}

void Camera::Externals::fprint(FILE* file) const {
  fprintf(file, "%0.6lf %0.6lf %0.6lf %0.6lf %0.6lf %0.6lf %0.6lf",
          c[0], c[1], c[2], q[0], q[1], q[2], q[3]);
}

Matrix<double, 3, 4> Camera::P() const {
  return internals_.K() * externals_.R() *
      (Matrix<double, 3, 4>() <<
       Matrix3d::Identity(), -externals_.C()).finished();
}

Eigen::Matrix3d Camera::R() const {
  return externals_.R();
}

Eigen::Vector3d Camera::C() const {
  return externals_.C();
}

Eigen::Matrix3d Camera::K() const {
  return internals_.K();
}

const Camera::Internals& Camera::internals() const {
  return internals_;
}

const Camera::Externals& Camera::GetExternals() const {
  return externals_;
}

const Eigen::Vector3d& Camera::GetCenter() const {
  return externals_.c;
}

Eigen::Vector3d Camera::GetUp() const {
  return R() * Vector3d(0, -1, 0);
}

Eigen::Vector3d Camera::GetTowards() const {
  return R() * Vector3d(0, 0, 1);
}

double Camera::GetFovX() const {
  return atan(internals_.p[0] / internals_.f[0]) * 2;
}

double Camera::GetFovY() const {
  return atan(internals_.p[1] / internals_.f[1]) * 2;
}

Eigen::Vector3d Camera::GetEulerAngles() const {
  auto &q = externals_.q;
  return Quaterniond(q[0], q[1], q[2], q[3]).matrix().eulerAngles(0, 1, 2);
}

Eigen::Vector2d Camera::GetPrincipalPoint() const {
  return Eigen::Vector2d(internals_.p[0], internals_.p[1]);
}

void Camera::SetInternals(double *f, double *k) {
  internals_.SetParams(f, k, internals_.p, internals_.skew);
}

void Camera::SetCenter(double *c) {
  SetCenter(Vector3d(c[0], c[1], c[2]));
}

void Camera::SetCenter(const Vector3d &c) {
  externals_.c = c;
}

void Camera::SetCenter(double x, double y, double z) {
  SetCenter(Vector3d(x, y, z));
}

void Camera::SetQuaternion(double *q) {
  SetQuaternion(Vector4d(q[0], q[1], q[2], q[3]));
}

void Camera::SetQuaternion(const Vector4d &q) {
  externals_.q = q;
}

void Camera::SetQuaternion(double w, double x, double y, double z) {
  SetQuaternion(Vector4d(w, x, y, z));
}

void Camera::SetRotationMatrix(const Eigen::Matrix3d &R) {
  SetQuaternion(Eigen::Quaterniond(R));
}

void Camera::SetEulerAngles(double rx, double ry, double rz) {
  Matrix3d m;
  m = AngleAxisd(rx, Vector3d::UnitX()) *
      AngleAxisd(ry, Vector3d::UnitY()) *
      AngleAxisd(rz, Vector3d::UnitZ());
  SetRotationMatrix(m);
}

void Camera::SetParams(double* f, double* k, double* p, double s,
                       double* q, double* c) {
  internals_.SetParams(f, k, p, s);
  externals_.SetParams(q, c);
}

void Camera::SetFilename(const std::string& filename) {
  filename_ = filename;
}

const std::string& Camera::GetFilename() const {
  return filename_;
}

void Camera::SetDirpath(const std::string& dirpath) {
  dirpath_ = dirpath;
}

const std::string& Camera::GetDirpath() const {
  return dirpath_;
}

string Camera::GetImagePath() const {
  return StringPrintf("%s/%s", dirpath_.c_str(), filename_.c_str());
}

void Camera::fprint(FILE* file) const {
  internals_.fprint(file);
  fprintf(file, " ");
  externals_.fprint(file);
  fprintf(file, " ");
  fprintf(file, "%s", filename_.c_str());
}

Eigen::Matrix3d Camera::HomographyFrom(const Camera& c, double d) const {
  auto K1 = c.K();
  auto R1 = c.R();
  auto C1 = c.C();
  auto K2 = K();
  auto R2 = R();
  auto C2 = C();

  Matrix3d H  = d * K2 * R2 * R1.transpose() * K1.inverse();
  H.col(2) += K2 * R2 * (C1 - C2);
  // H.col(2) += K2 * (C2 - R2 * R1.transpose() * C1);
  // Matrix3d H = d * K2 * R2;
  // H.col(2) = K2 * C2;
  // H *= K1.inverse();
  return H;
}

Eigen::Matrix3d Camera::RotationHomography() const {
  return K() * R().transpose() * K().inverse();
}

Eigen::Vector3d Camera::ProjectPoint(const Eigen::Vector3d& p) const {
  return K() * R() * (p - C());
}

cv::Mat Camera::ReadImage() const {
  return const_cast<Camera*>(this)->ReadImage();
}

cv::Mat Camera::ReadImage() {
  if (image_.empty()) {
    string image_name = StringPrintf(
        "%s/%s", AddNfsPath(dirpath_).c_str(), filename_.c_str());
    image_ = imread(image_name);
    CHECK(!image_.empty()) << "Failed to read " << image_name;
  }
  return image_;
}

cv::Mat Camera::ReadUndistortedImage() {
  return this->Undistort(ReadImage());
}

cv::Mat Camera::ShowDetections() const {
  return DrawPoints(ReadImage(), detections_, 5);
}

cv::Mat Camera::Undistort(const Mat &image) const {
  cv::Mat undistorted_image(image.size(), image.type());
  cv::Mat camera_matrix(3, 3, CV_32F);
  auto K = internals_.K();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      camera_matrix.at<float>(i, j) = K(i, j);
    }
  }
  vector<float> dist_coeffs(4, 0);
  dist_coeffs[0] = internals_.k[0];
  dist_coeffs[1] = internals_.k[1];
  cv::undistort(image, undistorted_image, camera_matrix, dist_coeffs);
  return undistorted_image;
}

template <typename T>
cv::Mat MakeMask(const cv::Mat image, std::function<uint8_t (const T&)> func) {
  cv::Mat mask(image.size(), CV_8UC1);
  for (int r = 0; r < image.rows; ++r) {
    for (int c = 0; c < image.cols; ++c) {
      mask.at<uint8_t>(r, c) = func(image.at<T>(r, c));
    }
  }
  return mask;
}

cv::Mat Camera::RemoveRotation(const cv::Mat &image, Mat *mask) const {
  Mat h = To<Mat>(this->RotationHomography());
  Mat rect_image;
  cv::warpPerspective(image, rect_image, h, image.size(), INTER_LINEAR);
  if (mask != nullptr) {
    assert(image.type() == CV_8UC3);
    *mask = MakeMask<Vec3b>(rect_image,
                     [](const Vec3b &p) -> uint8_t {
                       if (p == Vec3b(0, 0, 0)) {
                         return 0;
                       } else {
                         return 255;
                       }
                     });
  }
  return rect_image;
}

int Camera::NumDetections() const {
  return detections_.size();
}

void Camera::AddDetection(const Point2 &detection, int track_id) {
  detections_.push_back(detection);
  track_id_index_.insert(make_pair(track_id, track_ids_.size()));
  track_ids_.push_back(track_id);
}

const Point2& Camera::GetDetection(int i) const {
  return detections_[i];
}

Point2 Camera::GetDetectionInRetinaCoord(int i) const {
  Eigen::Vector2d norm_coord = Eigen::Vector3d(K().colPivHouseholderQr().solve(Eigen::Vector3d(detections_[i].x, detections_[i].y, 1))).hnormalized();
  // VLOG(1) << detections_[i] << ' ' << norm_coord.transpose();
  return Point2(norm_coord[0], norm_coord[1]);
}

track_id_t Camera::GetTrackId(int i) const {
  return track_ids_[i];
}

const Point2& Camera::GetDetectionOfTrack(track_id_t id) const {
  // VLOG(1) << track_id_index_.size() << ' ' << id << ' ' << GetId() << '\n';
  return detections_[track_id_index_.at(id)];
}

Point2 Camera::GetDetectionOfTrackInRetinaCoord(track_id_t id) const {
  return GetDetectionInRetinaCoord(track_id_index_.at(id));
}

cv::Vec3b Camera::GetColorOfTrack(track_id_t id) const {
  return const_cast<Camera*>(this)->ReadImage().at<cv::Vec3b>(GetDetectionOfTrack(id));
}

void Camera::RemoveTrack(int index) {
  auto track_id = track_ids_[index];
  track_id_index_.at(track_ids_.back()) = index;
  detections_[index] = detections_.back();
  track_ids_[index] = track_ids_.back();
  track_id_index_.erase(track_id);
  detections_.resize(detections_.size() - 1);
  track_ids_.resize(track_ids_.size() - 1);
}

void Camera::RemoveTrackById(track_id_t track_id) {
  auto index = track_id_index_.at(track_id);
  RemoveTrack(index);
}

void Camera::ResetTrackIdIndex() {
  track_id_index_.clear();
  for (size_t i = 0; i < track_ids_.size(); ++i) {
    track_id_index_[track_ids_[i]] = i;
  }
}

void Track::SetPoint(double *p) {
  SetPoint(p[0], p[1], p[2]);
}

void Track::SetPoint(double x, double y, double z) {
  point_[0] = x;
  point_[1] = y;
  point_[2] = z;
}

void Track::SetPoint(const Eigen::Vector3d &point) {
  point_ = point;
}

const Eigen::Vector3d& Track::GetPoint() const {
  return point_;
}

const Eigen::Vector3ub& Track::GetColor() const {
  return color_;
}

void Track::SetColor(const Eigen::Vector3i &color) {
  color_[0] = color[0];
  color_[1] = color[1];
  color_[2] = color[2];
}

void Track::SetColor(uint8_t r, uint8_t g, uint8_t b) {
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
}

void Track::AddCameraId(camera_id_t camera_id) {
  camera_id_index_[camera_id] = camera_ids_.size();
  camera_ids_.push_back(camera_id);
}

int Track::NumCameras() const {
  return camera_ids_.size();
}

camera_id_t Track::GetCameraId(int i) const {
  return camera_ids_[i];
}

void Track::RemoveCamera(int index) {
  auto camera_id = camera_ids_[index];
  camera_id_index_[camera_ids_.back()] = index;
  camera_ids_[index] = camera_ids_.back();
  camera_id_index_.erase(camera_id);
  camera_ids_.resize(camera_ids_.size() - 1);
}

void Track::RemoveCameraById(camera_id_t camera_id) {
  auto index = camera_id_index_.at(camera_id);
  RemoveCamera(index);
}

bool Model::ReadFile(const std::string &prefix) {
  ReadCameras(StringPrintf("%s_cameras.txt", prefix.c_str()), &cameras_);
  ReadTracks(StringPrintf("%s_tracks.txt", prefix.c_str()), &tracks_);
  image_dir_ = cameras_[0]->GetDirpath();
  next_camera_id_ = cameras_.size();
  next_track_id_ = tracks_.size();
  for (int i = 0; i < NumCameras(); ++i) {
    camera_id_index_[i] = i;
  }
  for (int i = 0; i < NumTracks(); ++i) {
    track_id_index_[i] = i;
  }
  return true;
}

bool Model::ReadImageDir(const std::string &image_dir,
                         const std::string &image_name_pattern) {
  vector<string> all_image_names;
  if (!ListDir(image_dir, "", &all_image_names)) {
    LOG(WARNING) << "Failed to read dir " << image_dir;
    return false;
  }
  vector<string> image_names;
  regex e(image_name_pattern);
  for (auto &name : all_image_names) {
    if (regex_match(name, e)) {
      image_names.push_back(name);
    }
  }
  if (image_names.empty()) {
    LOG(WARNING) << "Can't find valid images in dir " << image_dir;
    return false;
  }
  double f[2] = {1781, 1781};
  double k[2] = {0, 0};
  double p[2] = {959.5, 539.5};
  double s = 0;
  double q[4] = {1, 0, 0, 0};
  double c[3] = {0, 0, 0};
  cameras_.resize(image_names.size());
  for (size_t i = 0; i < cameras_.size(); ++i) {
    auto camera = new Camera;
    camera->SetFilename(image_names[i]);
    camera->SetDirpath(image_dir);
    cameras_[i] = camera;
    camera_id_index_[i] = i;
  }
  Mat image = cameras_[0]->ReadImage();
  p[0] = (image.cols - 1) / 2;
  p[1] = (image.rows - 1) / 2;
  static const double kNexusFocusLenghMagicNumber = 0.9382;
  f[0] = f[1] = kNexusFocusLenghMagicNumber * image.cols;
  for (size_t i = 0; i < cameras_.size(); ++i) {
    cameras_[i]->SetParams(f, k, p, s, q, c);
  }
  sort(cameras_.begin(), cameras_.end(),
       [](const Camera* const c0, const Camera* const c1) {
         return c0->GetFilename() < c1->GetFilename();
       });
  for (size_t i = 0; i < cameras_.size(); ++i) {
    cameras_[i]->SetId(i);
  }
  next_camera_id_ = cameras_.size();
  image_dir_ = image_dir;

  return true;
}

void Model::WriteFile(const std::string &prefix, bool write_covariance) {
  VLOG(1) << "Writing cameras\n";
  WriteCameras(StringPrintf("%s_cameras.txt", prefix.c_str()),
               cameras_, RemoveNfsPath(image_dir_));
  VLOG(1) << "Writing tracks\n";
  WriteTracks(StringPrintf("%s_tracks.txt", prefix.c_str()), tracks_);
  if (write_covariance) {
    WriteCovariance(StringPrintf("%s_tracks_cov.txt", prefix.c_str()), tracks_);
  }
}

int Model::NumTracks() const {
  return tracks_.size();
}

int Model::NumCameras() const {
  return cameras_.size();
}

const Camera* Model::GetCamera(int i) const {
  return cameras_[i];
}

Camera* Model::GetCamera(int i) {
  return cameras_[i];
}

camera_id_t Model::GetCameraId(int i) const {
  return cameras_[i]->GetId();
}

const Track* Model::GetTrack(int i) const {
  return tracks_[i];
}

Track* Model::GetTrack(int i) {
  return tracks_[i];
}

track_id_t Model::GetTrackId(int i) const {
  return tracks_[i]->GetId();
}

Eigen::Vector3d Model::GetMeanPoint() const {
  Eigen::Vector3d mean(0, 0, 0);
  for (auto track : tracks_) {
    mean += track->GetPoint();
  }
  return mean / tracks_.size();
}

void Model::GetTrackBoundingBox(Eigen::Vector3d *min_coord,
                                Eigen::Vector3d *max_coord) const {
  *min_coord = tracks_[0]->GetPoint();
  *max_coord = tracks_[0]->GetPoint();
  for (auto track : tracks_) {
    auto &p = track->GetPoint();
    if (p[0] < (*min_coord)[0]) (*min_coord)[0] = p[0];
    if (p[1] < (*min_coord)[1]) (*min_coord)[1] = p[1];
    if (p[2] < (*min_coord)[2]) (*min_coord)[2] = p[2];
    if (p[0] > (*max_coord)[0]) (*max_coord)[0] = p[0];
    if (p[1] > (*max_coord)[1]) (*max_coord)[1] = p[1];
    if (p[2] > (*max_coord)[2]) (*max_coord)[2] = p[2];
  }
}


void Model::AddTrack(Track *track,
                     const std::vector<Point2> &detections) {
  CHECK_EQ(track->NumCameras(), detections.size());
  track_id_t track_id = GetAvailableTrackId();
  track->SetId(track_id);
  track_id_index_[track_id] = tracks_.size();
  tracks_.push_back(track);
  for (int i = 0; i < track->NumCameras(); ++i) {
    GetCameraById(track->GetCameraId(i))->AddDetection(detections[i], track_id);
  }
}

void Model::AddCamera(Camera *camera) {
  camera_id_t camera_id = GetAvailableCameraId();
  camera->SetId(camera_id);
  camera_id_index_[camera_id] = cameras_.size();
  cameras_.push_back(camera);
  for (int i = 0; i < camera->NumDetections(); ++i) {
    GetTrackById(camera->GetTrackId(i))->AddCameraId(camera_id);
  }
}

void Model::RemoveCameras(int begin, int end) {
  CHECK_LE(end, NumCameras());
  if (end < 0) end = NumCameras();
  CHECK_LE(begin, end);
  if (begin == end) return;
  vector<track_id_t> track_ids;
  for (int i = begin; i < end; ++i) {
    auto camera = GetCamera(i);
    track_ids.resize(camera->NumDetections());
    for (int j = 0; j < camera->NumDetections(); ++j) {
      track_ids[j] = camera->GetTrackId(j);
    }
    for (auto id : track_ids) {
      RemoveDetection(id, camera->GetId());
    }
    camera_id_index_.erase(camera->GetId());
    delete camera;
  }

  int num_rest = cameras_.size() - end;
  for (int i = 0; i < num_rest; ++i) {
    cameras_[i + begin] = cameras_[end + i];
    camera_id_index_[cameras_[i + begin]->GetId()] = i + begin;
  }

  cameras_.resize(cameras_.size() - (end - begin));
}

void Model::RemoveTrack(int index) {
  auto track = tracks_[index];
  auto track_id = track->GetId();
  track_id_index_[tracks_.back()->GetId()] = index;
  tracks_[index] = tracks_.back();
  track_id_index_.erase(track_id);
  tracks_.resize(tracks_.size() - 1);

  for (int i = 0; i < track->NumCameras(); ++i) {
    GetCameraById(track->GetCameraId(i))->RemoveTrackById(track_id);
  }

  delete track;
}

void Model::RemoveTrackById(track_id_t track_id) {
  auto index = track_id_index_[track_id];
  RemoveTrack(index);
}

void Model::RemoveTracksById(const std::vector<track_id_t> &ids) {
  for (auto id : ids) {
    RemoveTrackById(id);
  }
}

void Model::RemoveTracks() {
  for (auto track : tracks_) {
    delete track;
  }
  tracks_.resize(0);
  for (auto camera : cameras_) {
    camera->detections_.resize(0);
    camera->track_ids_.resize(0);
    camera->track_id_index_.clear();
  }
  track_id_index_.clear();
  next_track_id_ = 0;
}

void Model::RemoveEmptyTracks() {
  vector<track_id_t> tracks_to_remove;
  for (int i = 0; i < NumTracks(); ++i) {
    auto track = GetTrack(i);
    if (track->NumCameras() == 0) {
      tracks_to_remove.push_back(track->GetId());
    }
  }
  for (auto id : tracks_to_remove) {
    RemoveTrackById(id);
  }
}

void Model::ColorTracks() {
  for (int i_track = 0; i_track < NumTracks(); ++i_track) {
    Vector3i total_color(0, 0, 0);
    auto track = GetTrack(i_track);
    for (int i = 0; i < track->NumCameras(); ++i) {
      auto color = GetCameraById(track->GetCameraId(i))->GetColorOfTrack(i_track);
      for (int j = 0; j < 3; ++j) {
        total_color[j] += color[2 - j];
      }
    }
    total_color /= track->NumCameras();
    track->SetColor(total_color);
  }
}

void Model::RemoveDetection(track_id_t track_id, camera_id_t camera_id) {
  GetTrackById(track_id)->RemoveCameraById(camera_id);
  GetCameraById(camera_id)->RemoveTrackById(track_id);
}

const Point2& Model::GetDetection(
    track_id_t track_id, camera_id_t camera_id) const {
  return GetCameraById(camera_id)->GetDetectionOfTrack(track_id);
}

Point2 Model::Project(track_id_t track_id, camera_id_t camera_id) {
  Vector3d p = GetCameraById(camera_id)->P() *
      GetTrackById(track_id)->GetPoint().homogeneous();
  return Point2(p[0] / p[2], p[1] / p[2]);
}

double Model::GetProjectionError(track_id_t track_id, camera_id_t camera_id) {
  auto camera = GetCameraById(camera_id);
  auto track = GetTrackById(track_id);
  Vector3d p3 = camera->P() * track->GetPoint().homogeneous();
  Point2 p2(p3[0] / p3[2], p3[1] / p3[2]);
  Point2 d = camera->GetDetectionOfTrack(track_id);
  return sqrt((p2.x - d.x) * (p2.x - d.x) + (p2.y - d.y) * (p2.y - d.y));
}

void Model::TriangulateTrack(int index) {
  auto track = GetTrack(index);
  auto track_id = track->GetId();
  vector<Matrix3d> rs(track->NumCameras());
  vector<Vector3d> cs(track->NumCameras());
  vector<Matrix3d> ks(track->NumCameras());
  vector<Vector2d> ps(track->NumCameras());
  for (int i = 0; i < track->NumCameras(); ++i) {
    auto camera = GetCameraById(track->GetCameraId(i));
    rs[i] = camera->R().transpose();
    cs[i] = camera->C();
    ks[i] = camera->K();
    auto d = camera->GetDetectionOfTrack(track_id);
    ps[i] = Vector2d(d.x, d.y);
  }
  track->SetPoint(TriangulateMidpoint(rs, cs,ks, ps));
}

void Model::TriangulateTrackById(track_id_t track_id) {
  TriangulateTrack(track_id_index_.at(track_id));
}

void Model::RelabelTracks() {
  unordered_map<track_id_t, track_id_t> dict;
  for (size_t i = 0; i < tracks_.size(); ++i) {
    auto track = tracks_[i];
    dict[track->GetId()] = i;
    track->SetId(i);
  }
  for (int i = 0; i < cameras_.size(); ++i) {
    auto camera = cameras_[i];
    for (size_t j = 0; j < camera->track_ids_.size(); ++j) {
      camera->track_ids_[j] = dict.at(camera->track_ids_[j]);
    }
    camera->ResetTrackIdIndex();
  }
}

void Model::ReorderCameras(const std::vector<int> &camera_indexes) {
  CHECK_EQ(camera_indexes.size(), cameras_.size());
  std::vector<Camera*> old_cameras_(cameras_.begin(), cameras_.end());
  for (size_t i = 0; i < cameras_.size(); ++i) {
    cameras_[i] = old_cameras_[camera_indexes[i]];
    camera_id_index_[cameras_[i]->GetId()] = i;
  }
}

Model::~Model() {
  for (auto camera : cameras_)
    delete camera;
  for (auto track : tracks_)
    delete track;
}

void Model2Pmvs(const std::string &out_dir, const Model &model) {
  string cmd;
  cmd = StringPrintf("mkdir %s/visualize/", out_dir.c_str());
  LOG(INFO) << "Running: " << cmd;
  system(cmd.c_str());
  cmd = StringPrintf("mkdir %s/txt/", out_dir.c_str());
  LOG(INFO) << "Running: " << cmd;
  system(cmd.c_str());
  cmd = StringPrintf("mkdir %s/models/", out_dir.c_str());
  LOG(INFO) << "Running: " << cmd;
  system(cmd.c_str());

  string filename;
  for (int i = 0; i < model.NumCameras(); ++i) {
    auto camera = model.GetCamera(i);
    filename = StringPrintf("%s/visualize/%08d.jpg", out_dir.c_str(), i);
    LOG(INFO) << "Writing: " << filename;
    imwrite(filename, camera->ReadImage());
    // system(cmd.c_str());
  }

  ofstream ofs;
  for (int i = 0; i < model.NumCameras(); ++i) {
    filename = StringPrintf("%s/txt/%08d.txt", out_dir.c_str(), i);
    ofs.open(filename);
    ofs << "CONTOUR\n";
    ofs << model.GetCamera(i)->P();
    ofs.close();
  }

  filename = StringPrintf("%s/options.txt", out_dir.c_str());
  ofs.open(filename);
  ofs << "timages -1 0 " << model.NumCameras() << '\n'
      << "oimages 0\n"
      << "level 1\n"
      << "csize 2\n"
      << "threshold 0.7\n"
      << "CPU 12\n";
  ofs.close();
}

void Model2Ply(const std::string &out_path, const Model &model) {
  ofstream ss(out_path);
  ss << "ply\nformat ascii 1.0\n"
      << "element vertex " << model.NumTracks() << '\n'
      << "property float x\n"
      << "property float y\n"
      << "property float z\n"
      << "property uchar red\n"
      << "property uchar green\n"
      << "property uchar blue\n"
      << "end_header\n";
  for (int i = 0; i < model.NumTracks(); ++i) {
    Eigen::Vector3d p = model.GetTrack(i)->GetPoint();
    Eigen::Vector3ub c = model.GetTrack(i)->GetColor();
    ss << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
        << (int)c[0] << ' ' << (int)c[1] << ' ' << (int)c[2] << '\n';
  }
}

void ReadCameras(const std::string& filename,
                 std::vector<Camera*>* cameras) {
  int num_cameras;
  char strbuf[1024];
  double k[2];
  double f[2];
  double p[2];
  double s;
  double c[3];
  double q[4];

  FILE* file = fopen(filename.c_str(), "r");
  CHECK(file != NULL) << "Can't read " << filename;
  fscanf(file, "%d %s", &num_cameras, strbuf);
  cameras->resize(num_cameras);
  for (size_t i = 0; i < cameras->size(); ++i) {
    (*cameras)[i] = new Camera;
    (*cameras)[i]->SetId(i);
  }
  string dirpath = strbuf;
  // dirpath = AddDropboxPath(dirpath);

  for (int i = 0; i < num_cameras; ++i) {
    fscanf(file, "%lf %lf %lf %lf %lf %lf %lf",
           k, k + 1, f, f + 1, p, p + 1, &s);
    fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %s",
           c, c + 1, c + 2,
           q, q + 1, q + 2, q + 3,
           strbuf);
    (*cameras)[i]->SetParams(f, k, p, s, q, c);
    (*cameras)[i]->SetFilename(strbuf);
    (*cameras)[i]->SetDirpath(dirpath);
  }

  vector<Point2> detections;
  vector<int> track_ids;
  Point2 detection;
  int track_id;
  int num_points;
  for (int i = 0; i < num_cameras; ++i) {
    detections.resize(0);
    track_ids.resize(0);
    fscanf(file, "%d", &num_points);
    for (int j = 0; j < num_points; ++j) {
      fscanf(file, "%f %f", &detection.x, &detection.y);
      detections.push_back(detection);
    }
    for (int j = 0; j < num_points; ++j) {
      fscanf(file, "%d", &track_id);
      track_ids.push_back(track_id);
    }
    for (int j = 0; j < num_points; ++j) {
      (*cameras)[i]->AddDetection(detections[j], track_ids[j]);
    }
  }
  fclose(file);
}

void WriteCameras(const std::string& filename,
                  const std::vector<Camera*>& cameras,
                  const std::string& dirpath) {
  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "%lu %s\n", cameras.size(), dirpath.c_str());
  for (auto& camera : cameras) {
    camera->fprint(file);
    fprintf(file, "\n");
  }
  for (auto &camera : cameras) {
    fprintf(file, "%d ", camera->NumDetections());
    for (int i = 0; i < camera->NumDetections(); ++i) {
      fprintf(file, "%0.6f %0.6f ",
              camera->GetDetection(i).x,
              camera->GetDetection(i).y);
    }
    for (int i = 0; i < camera->NumDetections(); ++i) {
      fprintf(file, "%lld ", camera->GetTrackId(i));
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

void ReadTracks(const std::string& filename,
                std::vector<Track*>* tracks) {
  FILE* file = fopen(filename.c_str(), "r");
  CHECK(file != NULL) << "Can't read " << filename;
  int num_tracks;
  fscanf(file, "%d", &num_tracks);
  tracks->resize(num_tracks);
  Vector3d p;
  Vector3i c;
  for (int i = 0; i < num_tracks; ++i) {
    (*tracks)[i] = new Track;
    (*tracks)[i]->SetId(i);
  }
  for (int i = 0; i < num_tracks; ++i) {
    fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
    (*tracks)[i]->SetPoint(p);
  }
  for (int i = 0; i < num_tracks; ++i) {
    fscanf(file, "%d %d %d", &c[0], &c[1], &c[2]);
    (*tracks)[i]->SetColor(c);
  }
  int num_cameras, camera_id;
  for (auto &track : *tracks) {
    fscanf(file, "%d", &num_cameras);
    for (int j = 0; j < num_cameras; ++j) {
      fscanf(file, "%d", &camera_id);
      track->AddCameraId(camera_id);
    }
  }
  fclose(file);
}

void WriteTracks(const std::string& filename,
                 const std::vector<Track*>& tracks) {
  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "%lu\n", tracks.size());
  for (auto &track : tracks) {
    auto &p = track->GetPoint();
    fprintf(file, "%0.6f %0.6f %0.6f\n", p[0], p[1], p[2]);
  }
  for (auto &track : tracks) {
    auto &c = track->GetColor();
    fprintf(file, "%d %d %d\n", c[0], c[1], c[2]);
  }
  for (auto &track : tracks) {
    fprintf(file, "%d ", track->NumCameras());
    for (int i = 0; i < track->NumCameras(); ++i) {
      fprintf(file, "%lld ", track->GetCameraId(i));
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

void WriteCovariance(const std::string &filename,
                     const std::vector<Track*> &tracks) {
  FILE* file = fopen(filename.c_str(), "w");
  for (auto track : tracks) {
    auto p = track->GetPoint();
    auto &cov = track->GetCovariance();
    fprintf(file, "%lf %lf %lf ", p[0], p[1], p[2]);
    for (auto v : cov) {
      fprintf(file, "%lf ", v);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

} // tiny
} // furry

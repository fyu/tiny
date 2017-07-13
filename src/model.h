#ifndef FURRY_TINY_MODEL_H
#define FURRY_TINY_MODEL_H

#include <string>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <unordered_map>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include <glog/logging.h>

#include "furry/common/eigen.h"

namespace furry {
namespace tiny {

// typedef cv::Point2d Point2;
typedef cv::Point2f Point2;

class Model;

class Identifiable {
 public:
  typedef int64_t IdType;
  void SetId(IdType id);
  IdType GetId() const;

 private:
  IdType id_;
};

typedef Identifiable::IdType track_id_t;
typedef Identifiable::IdType camera_id_t;

class Camera : public Identifiable {

  // Projection matrix:
  // K * R * [I -C]
 public:

  struct Internals {
    Eigen::Matrix3d K() const;
    void fprint(FILE* file) const;
    void SetParams(double* f, double* k, double* p, double s);
    double f[2];
    double k[2];
    double p[2];
    double skew;
  };

  struct Externals {
    Eigen::Matrix3d R() const;
    Eigen::Vector3d C() const;
    void fprint(FILE* file) const;
    void SetParams(double* q, double* c);
    Eigen::Vector4d q;
    Eigen::Vector3d c;
  };

  Eigen::Matrix<double, 3, 4> P() const;
  Eigen::Matrix3d R() const;
  Eigen::Vector3d C() const;
  Eigen::Matrix3d K() const;
  Eigen::Quaterniond GetQuaternion() const;

  void fprint(FILE* file) const;
  void SetParams(double* f, double* k, double* p, double s,
                 double* q, double* c);
  void SetFilename(const std::string& filename);
  const std::string& GetFilename() const;
  void SetDirpath(const std::string& dirpath);
  const std::string& GetDirpath() const;
  std::string GetImagePath() const;
  Eigen::Matrix3d HomographyFrom(const Camera& c, double d) const;
  // K * R^-1 * K^-1
  Eigen::Matrix3d RotationHomography() const;
  Eigen::Vector3d ProjectPoint(const Eigen::Vector3d& p) const;

  const Internals& internals() const;

  const Externals& GetExternals() const;

  const Eigen::Vector3d& GetCenter() const;
  Eigen::Vector3d GetUp() const;
  Eigen::Vector3d GetTowards() const;
  double GetFovX() const;
  double GetFovY() const;
  Eigen::Vector3d GetEulerAngles() const;
  Eigen::Vector2d GetPrincipalPoint() const;

  void SetInternals(double *f, double *k);
  void SetCenter(double *c);
  void SetCenter(const Eigen::Vector3d &c);
  void SetCenter(double x, double y, double z);
  void SetQuaternion(double *q);
  void SetQuaternion(const Eigen::Vector4d &q);
  void SetQuaternion(double w, double x, double y, double z);
  void SetQuaternion(const Eigen::Quaterniond &q);
  void SetRotationMatrix(const Eigen::Matrix3d &R);
  void SetEulerAngles(double rx, double ry, double rz);

  cv::Mat ReadImage();
  cv::Mat ReadImage() const;
  cv::Mat ReadUndistortedImage();
  cv::Mat ShowDetections() const;
  cv::Mat Undistort(const cv::Mat &image) const;
  cv::Mat RemoveRotation(const cv::Mat &image, cv::Mat *mask = nullptr) const;

  int NumDetections() const;
  void AddDetection(const Point2 &detection, int track_id);
  const Point2& GetDetection(int i) const;
  Point2 GetDetectionInRetinaCoord(int i) const;
  track_id_t GetTrackId(int i) const;
  const Point2& GetDetectionOfTrack(track_id_t id) const;
  Point2 GetDetectionOfTrackInRetinaCoord(track_id_t id) const;
  cv::Vec3b GetColorOfTrack(track_id_t id) const;
  void RemoveTrack(int index);
  void RemoveTrackById(track_id_t track_id);

 private:
  void ResetTrackIdIndex();
 private:
  friend class Model;

  Internals internals_;
  Externals externals_;
  std::string filename_;
  std::string dirpath_;
  cv::Mat image_;

  std::vector<Point2> detections_;
  std::vector<track_id_t> track_ids_;
  std::unordered_map<track_id_t, int> track_id_index_;
};

class Track : public Identifiable {
 public:
  const Eigen::Vector3d& GetPoint() const;
  void SetPoint(double *p);
  void SetPoint(double x, double y, double z);
  void SetPoint(const Eigen::Vector3d &point);
  const Eigen::Vector3ub& GetColor() const;
  void SetColor(uint8_t r, uint8_t g, uint8_t b);
  void SetColor(const Eigen::Vector3i &color);
  void SetColor(const cv::Vec3b &color);
  void AddCameraId(camera_id_t camera_id);
  int NumCameras() const;
  camera_id_t GetCameraId(int i) const;
  void RemoveCamera(int index);
  void RemoveCameraById(camera_id_t camera_id);

  void SetCovariance(const std::vector<double> &covariance);
  const std::vector<double>& GetCovariance() const;
 private:
  Eigen::Vector3d point_;
  Eigen::Vector3ub color_;
  std::vector<int> camera_ids_;
  std::vector<double> covariance_;
  std::unordered_map<camera_id_t, int> camera_id_index_;
};

class Model {
 public:
  bool ReadFile(const std::string &prefix);
  bool ReadImageDir(const std::string &image_dir,
                    const std::string &image_name_pattern);
  void WriteFile(const std::string &prefix, bool write_covariance = false);
  int NumTracks() const;
  int NumCameras() const;
  const Camera* GetCamera(int i) const;
  Camera* GetCamera(int i);
  camera_id_t GetCameraId(int i) const;
  const Track* GetTrack(int i) const;
  Track* GetTrack(int i);
  track_id_t GetTrackId(int i) const;
  Eigen::Vector3d GetMeanPoint() const;

  const Camera* GetCameraById(camera_id_t camera_id) const;
  Camera* GetCameraById(camera_id_t camera_id);
  const Track* GetTrackById(track_id_t track_id) const;
  Track* GetTrackById(track_id_t track_id);
  void GetTrackBoundingBox(Eigen::Vector3d *max_coord,
                           Eigen::Vector3d *min_coord) const;

  void AddTrack(Track *track, const std::vector<Point2> &detections);
  void AddCamera(Camera *camera);
  void RemoveCameras(int begin, int end = -1);
  void RemoveTrack(int index);
  void RemoveTrackById(track_id_t track_id);
  void RemoveTracksById(const std::vector<track_id_t> &ids);
  void RemoveTracks();
  void RemoveEmptyTracks();
  void ColorTracks();

  void RemoveDetection(track_id_t track_id, camera_id_t camera_id);
  const Point2& GetDetection(track_id_t track_id, camera_id_t camera_id) const;
  Point2 Project(track_id_t track_id, camera_id_t camera_id);
  double GetProjectionError(track_id_t track_id, camera_id_t camera_id);

  void TriangulateTrack(int index);
  void TriangulateTrackById(track_id_t track_id);

  void RelabelTracks();
  void ReorderCameras(const std::vector<int> &camera_indexes);

  ~Model();

 private:
  camera_id_t GetAvailableCameraId();
  track_id_t GetAvailableTrackId();

 private:
  std::vector<Camera*> cameras_;
  std::vector<Track*> tracks_;
  std::unordered_map<camera_id_t, int> camera_id_index_;
  std::unordered_map<track_id_t, int> track_id_index_;
  camera_id_t next_camera_id_ = 0;
  track_id_t next_track_id_ = 0;
  std::string image_dir_;
};

void Model2Pmvs(const std::string &out_dir, const Model &model);

void Model2Ply(const std::string &out_path, const Model &model);

void ReadCameras(const std::string& filename,
                 std::vector<Camera*>* cameras);

void WriteCameras(const std::string& filename,
                  const std::vector<Camera*>& cameras,
                  const std::string& dirpath);

void ReadTracks(const std::string& filename,
                std::vector<Track*>* tracks);

void WriteTracks(const std::string& filename,
                 const std::vector<Track*> &tracks);

void WriteCovariance(const std::string& filename,
                      const std::vector<Track*> &tracks);

inline void Identifiable::SetId(IdType id) {
  id_ = id;
}

inline Identifiable::IdType Identifiable::GetId() const {
  return id_;
}

inline void Track::SetCovariance(const std::vector<double> &covariance) {
  covariance_ = covariance;
}

inline const std::vector<double>& Track::GetCovariance() const {
  return covariance_;
}

inline camera_id_t Model::GetAvailableCameraId() {
  return next_camera_id_++;
}

inline track_id_t Model::GetAvailableTrackId() {
  return next_track_id_++;
}

inline const Camera* Model::GetCameraById(camera_id_t camera_id) const {
  return const_cast<Model*>(this)->GetCameraById(camera_id);
}

inline Camera* Model::GetCameraById(camera_id_t camera_id) {
  // VLOG(1) << camera_id_index_.size();
  return GetCamera(camera_id_index_.at(camera_id));
}

inline const Track* Model::GetTrackById(track_id_t track_id) const {
  return const_cast<const Model*>(this)->GetTrackById(track_id);
}

inline Track* Model::GetTrackById(track_id_t track_id) {
  return GetTrack(track_id_index_.at(track_id));
}


} // tiny
} // furry

#endif // FURRY_TINY_MODEL_H

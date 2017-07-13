#ifndef FURRY_COMMON_POSE_H_
#define FURRY_COMMON_POSE_H_

#include <istream>
#include <ostream>

#include "opencv2/core/core.hpp"

#include <Eigen/Dense>

namespace furry
{

class Pose
{
 public:
  static Pose Identity();

  Pose();
  Pose(const Pose& p);
  Pose(const Eigen::Vector3d& p, const Eigen::Quaterniond &q);
  Pose(double x, double y, double z, double qx, double qy, double qz, double qw);
  // double pos_x;
  // double pos_y;
  // double pos_z;
  // double quat_x;
  // double quat_y;
  // double quat_z;
  // double quat_w;

  Eigen::Quaterniond quaternion() const;
  Eigen::Vector4d quaternion_coeffs() const;

  Eigen::Matrix<double, 3, 4> toProjectiveMatrix() const;
  Eigen::Vector3d towards() const;
  Eigen::Vector3d up() const;
  Eigen::Vector3d right() const;
  Eigen::Vector3d position() const;
  Eigen::Matrix3d rotation_matrix() const;
  Eigen::Matrix<double, 3, 4> projective_matrix() const;
  Eigen::Matrix4d full_projective_matrix() const;
  //Eigen::Matrix4d inverse_projective_matrix() const;
  Eigen::Matrix3d RotationMatrixFrom(const Pose& p) const;
  Eigen::Quaterniond RotationQuaterionFrom(const Pose& p) const;
  Eigen::Vector3d RotationAnglesFrom(const Pose& p);
  Eigen::Vector3d TranslationFrom(const Pose& p) const;
  Pose TransformIn(const Pose &p) const;
  double DistanceTo(const Pose& p) const;
  double Distance2To(const Pose &p) const;

  Pose& set_position(const Eigen::Vector3d &p);
  Pose& set_rotation(const Eigen::Quaterniond &q);

  Pose& operator = (const Pose &p);

  friend Pose LinearInterpolate(const Pose& p0, const Pose& p1, double t);

  friend std::istream& operator >> (std::istream & is, Pose& pose);

  friend std::ostream& operator << (std::ostream & os, const Pose & pose);

 protected:
  Eigen::Vector3d pos_;
  Eigen::Quaterniond quat_;

};

// linear interpolation between two poses
// conceptually p = (1-t) * p0 + t * p1
Pose
LinearInterpolate(const Pose& p0, const Pose& p1, double t);

std::istream&
operator >> (std::istream & is, Pose& pose);

std::ostream&
operator << (std::ostream & os, const Pose & pose);

////////////////////////////////////////////////////////////
// Pose Inline
////////////////////////////////////////////////////////////

inline Eigen::Matrix3d Pose::rotation_matrix() const {
  return quat_.toRotationMatrix();
}

inline Eigen::Quaterniond Pose::quaternion() const {
  return quat_;
}
inline Eigen::Vector4d Pose::quaternion_coeffs() const {
  return quat_.coeffs();
}

inline double Pose::DistanceTo(const Pose &p) const {
  return sqrt(Distance2To(p));
}

inline double Pose::Distance2To(const Pose &p) const {
  return (pos_ - p.pos_).array().square().sum();
}

inline Pose Pose::TransformIn(const Pose &p) const {
  return Pose(TranslationFrom(p), RotationQuaterionFrom(p));
}

inline Pose&
Pose::operator = (const Pose &p) {
  this->pos_ = p.pos_;
  this->quat_ = p.quat_;
  return *this;
}

} // furry

#endif //_FURRY_POSE_H_

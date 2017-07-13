#include "furry/common/pose.h"
#include "furry/common/geometry.h"

#include <iostream>

using cv::Matx;
using cv::Matx33d;
using cv::Vec3d;

//using namespace std;

namespace furry
{

Pose Pose::Identity() {
  return Pose(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity());
}

Pose::Pose() : pos_(0, 0, 0), quat_(0, 0, 0, 0)
{
}

Pose::Pose(const Pose& p) : pos_(p.pos_), quat_(p.quat_)
{
}

Pose::Pose(const Eigen::Vector3d& p, const Eigen::Quaterniond& q) :
  pos_(p), quat_(q)
{
}

Pose::Pose(double x, double y, double z, double qx, double qy, double qz, double qw) :
  pos_(x, y, z), quat_(qw, qx, qy, qz)
{
}

Pose
LinearInterpolate(const Pose& p0, const Pose& p1, double t)
{
  return Pose(p0.pos_ * (1-t) + p1.pos_ * t, p0.quat_.slerp(t, p1.quat_));
}

// Pose
// operator + (const Pose& p0, const Pose& p1)
// {
//   return Pose(p0.pos_x + p1.pos_x,
//               p0.pos_y + p1.pos_y,
//               p0.pos_z + p1.pos_z,
//               p0.quat_x + p1.quat_x,
//               p0.quat_y + p1.quat_y,
//               p0.quat_z + p1.quat_z,
//               p0.quat_w + p1.quat_w);
// }

// Pose
// operator - (const Pose& p0, const Pose& p1)
// {
//   return Pose(p0.pos_x - p1.pos_x,
//               p0.pos_y - p1.pos_y,
//               p0.pos_z - p1.pos_z,
//               p0.quat_x - p1.quat_x,
//               p0.quat_y - p1.quat_y,
//               p0.quat_z - p1.quat_z,
//               p0.quat_w - p1.quat_w);
// }

std::istream&
operator >> (std::istream& is, Pose& p)
{
  is >> p.pos_.x() >> p.pos_.y() >> p.pos_.z()
     >> p.quat_.x() >> p.quat_.y() >> p.quat_.z() >> p.quat_.w();
  return is;
}

std::ostream&
operator << (std::ostream& os, const Pose& p)
{
    os << p.pos_.x() << " "
       << p.pos_.y() << " "
       << p.pos_.z() << " "
       << p.quat_.x() << " "
       << p.quat_.y() << " "
       << p.quat_.z() << " "
       << p.quat_.w();
    return os;
}

Eigen::Matrix<double, 3, 4> Pose::toProjectiveMatrix() const
{
  Eigen::Matrix<double, 3, 4> t = Eigen::Matrix<double, 3, 4>::Identity();
  t(0, 3) = -pos_.x();
  t(1, 3) = -pos_.y();
  t(2, 3) = -pos_.z();

  Eigen::Matrix3d r = quat2mat<double, 3, 3>(quat_.x(), quat_.y(), quat_.z(), quat_.w());

  //cout << t << endl << r << endl;

  return r.transpose() * t;
}

Eigen::Matrix<double, 3, 4> Pose::projective_matrix() const {
  // Eigen::Matrix<double, 3, 4> t = Eigen::Matrix<double, 3, 4>::Identity();
  // t(0, 3) = -pos_.x();
  // t(1, 3) = -pos_.y();
  // t(2, 3) = -pos_.z();

  // Eigen::Matrix3d r = quat2mat<double, 3, 3>(quat_.x(), quat_.y(), quat_.z(), quat_.w());

  // //cout << t << endl << r << endl;
  // return r.transpose() * t;
  auto rt = rotation_matrix().transpose();
  return (Eigen::Matrix<double, 3, 4>() << rt, -rt * position()).finished();
}

Eigen::Matrix4d Pose::full_projective_matrix() const {
  return (Eigen::Matrix4d() << projective_matrix(),
          0, 0, 0, 1).finished();
}

Eigen::Vector3d Pose::towards() const
{
  //  return quat2mat<double, 3, 3>(quat_.x(), quat_.y(), quat_.z(), quat_.w()) *
  //Eigen::Vector3d(0, 0, 1);
  return quat_ * Eigen::Vector3d(0, 0, -1);
}

Eigen::Vector3d Pose::up() const {
  return quat_ * Eigen::Vector3d::UnitY();
}

Eigen::Vector3d Pose::right() const {
  return quat_ * Eigen::Vector3d(1, 0, 0);
}

Eigen::Vector3d Pose::position() const
{
  return Eigen::Vector3d(pos_.x(), pos_.y(), pos_.z());
}

Pose& Pose::set_position(const Eigen::Vector3d &p) {
  pos_ = p;
  return *this;
}

Pose& Pose::set_rotation(const Eigen::Quaterniond &q) {
  quat_ = q;
  return *this;
}

Eigen::Matrix3d
Pose::RotationMatrixFrom(const Pose& p) const {
  return RotationQuaterionFrom(p).matrix();
  //return q.matrix();
}

Eigen::Quaterniond
Pose::RotationQuaterionFrom(const Pose &p) const {
  // Eigen::Quaterniond q = this->quat_ * p.quat_.conjugate();
  // Eigen::Quaterniond q = this->quat_.conjugate() * p.quat_;
  Eigen::Quaterniond q = p.quat_.conjugate() * this->quat_;
  return q.normalized();
  //  return q;
}

Eigen::Vector3d
Pose::TranslationFrom(const Pose& p) const {
  return p.rotation_matrix().transpose() * (this->pos_ - p.pos_);
}

} // furry

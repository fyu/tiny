#ifndef FURRY_COMMON_GEOMETRY
#define FURRY_COMMON_GEOMETRY

#include <cmath>
#include <limits>

#include <Eigen/Dense>

#include "furry/common/cv.h"

namespace furry
{

// inline Eigen::Vector3d
// normalize_line(Eigen::Vector3d line)
// {
//   double norm = std::sqrt(line[0] * line[0] + line[1] * line[1]);
//   if (norm > std::numeric_limits<double>::epsilon())
//     line /= norm;
//   else
//     line = Eigen::Vector3d(0, 0, 1);
//   return line;
// }

template <typename T, int s, int t>
cv::Vec<T, t> project(const cv::Matx<T, t, s>& m, const cv::Vec<T, s>& v)
{
  return m * v;
}

// project homogeneous coordinates
template <typename T, int s, int t>
cv::Vec<T, t> project(const cv::Matx<T, t+1, s+1>& m, const cv::Vec<T, s>& v)
{
  cv::Mat hv(v);
  hv.push_back(1.);
  cv::Mat r = cv::Mat(m) * cv::Mat(hv);
  for (int i = 0; i < t; ++i)
  {
    r.at<T>(i, 0) /= r.at<T>(t, 0);
  }
  r.pop_back(1);
  return r;
}

template <typename T, int s, int t>
Eigen::Matrix<T, t, 1> project(const Eigen::Matrix<T, t, s>& m, const Eigen::Matrix<T, s, 1>& v)
{
  return m * v;
}

// project homogeneous coordinates
template <typename T, int s, int t>
Eigen::Matrix<T, t, 1> project(const Eigen::Matrix<T, t+1, s+1>& m,
                               const Eigen::Matrix<T, s, 1>& v)
{
  Eigen::Matrix<T, t+1, 1> pv = m * v.homogeneous();
  return pv.hnormalized();
}

template <typename T, int r, int c>
Eigen::Matrix<T, r, c> quat2mat(T x, T y, T z, T w)
{
  static_assert(r >= 3 && c>= 3,
                "Rotation matrix must be at least 3 by 3");

  Eigen::Matrix<T, r, c> m = Eigen::Matrix<T, r, c>::Identity();
  //std::cout << w << ' ' << x << ' ' << y << ' ' << z << '\n';
  m(0, 0) = 1 - 2 * (y * y + z * z);
  m(0, 1) = 2 * (x * y - z * w);
  m(0, 2) = 2 * (x * z + y * w);
  m(1, 0) = 2 * (x * y + z * w);
  m(1, 1) = 1 - 2 * (x * x + z * z);
  m(1, 2) = 2 * (y * z - x * w);
  m(2, 0) = 2 * (x * z - y * w);
  m(2, 1) = 2 * (y * z + x * w);
  m(2, 2) = 1 - 2 * (x * x + y * y);

  return m;
}

template <typename Scalar>
Eigen::Matrix<Scalar, 3, 3> quat2rm(const Eigen::Quaternion<Scalar> &quat) {
  return quat2mat<Scalar, 3, 3>(quat.x(), quat.y(), quat.z(), quat.w());
}

// template <typename T, int d> Eigen::Matrix<T, d-1, 1>
// homo2eu(const Eigen::Matrix<T, d, 1>& p)
// {
//   return p.head<d-1>() / p[d-1];
// }

// Return the line end points crossed by the frame of the given
// width and height
bool LineEndPoints(const Eigen::Vector3d &line, int width, int height,
                   Eigen::Vector2d *p0, Eigen::Vector2d *p1);

Eigen::Matrix3d Quaternion2RotationMatrix(double x, double y, double z, double w);

template <typename Scalar> Eigen::Matrix<Scalar, 3, 1>
    QuaternionToAngles(const Eigen::Quaternion<Scalar> &quat) {
  auto np = quat.normalized();
  Scalar q0 = np.w();
  Scalar q1 = np.x();
  Scalar q2 = np.y();
  Scalar q3 = np.z();
  Eigen::Matrix<Scalar, 3, 1> angles;
  angles.x() = - atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
  angles.y() = - asin(2 * (q0 * q2 - q3 * q1));
  angles.z() = - atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
  return angles;
}

// rotation angles to rotation matrix
template <typename Scalar> Eigen::Matrix<Scalar, 3, 3>
    AnglesToRotationMatrix(const Eigen::Matrix<Scalar, 3, 1> &angles) {
  Eigen::Matrix<Scalar, 3, 3> rx, ry, rz;
  Scalar x = angles.x();
  Scalar y = angles.y();
  Scalar z = angles.z();
  rx << 1, 0, 0, 0, cos(x), sin(x), 0, -sin(x), cos(x);
  ry << cos(y), 0, -sin(y), 0, 1, 0, sin(y), 0, cos(y);
  rz << cos(z), sin(z), 0, -sin(z), cos(z), 0, 0, 0, 1;
  return rz * ry * rx;
}

template <typename Scalar>
void DecomposeQuaternion(const Eigen::Quaternion<Scalar> &quat,
                         Scalar &angle,
                         Eigen::Matrix<Scalar, 3, 1> &axis) {
  auto unit_quat = quat.normalized();
  angle = acos(unit_quat.w()) * 2;
  axis = quat.vec() / sin(angle / 2);
}

} // furry

#endif

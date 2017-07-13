#include "numeric.h"

// #include "image/rushmore/numeric/arraytypes.h"
// #include "image/rushmore/numeric/quaternion.h"

using namespace Eigen;

namespace furry {
namespace tiny {

// Eigen::Vector3d ToEigen(const lightfield_sfm::Point& p) {
//   return Eigen::Vector3d(p[0], p[1], p[2]);
// }

// lightfield_sfm::Point ToLfPoint(const Eigen::Vector3d& v) {
//   lightfield_sfm::Point point;
//   point.set_type(0);
//   point.SetElements(v.data());
//   return point;
// }

// void RotationToQuaternion(const Matrix3d& m, rushmore::numeric::VectorD4 *q) {
//   double d[4];
//   double R[3][3];
//   for (int i = 0; i < 3; ++i) {
//     for (int j = 0; j < 3; ++j) {
//       R[i][j] = m(i, j);
//     }
//   }
//   rushmore::numeric::RotationToQuaternion(R, d);
//   q->CopyFrom(d);
// }

// void QuaternionToEulerAngles(const rushmore::numeric::VectorD4& q,
//                              double* angles) {
//   Eigen::Matrix3d mat = ToEigen(rushmore::numeric::QuaternionToRotation(q));
//   Eigen::Vector3d v = mat.eulerAngles(0, 1, 2);
//   angles[0] = v[0];
//   angles[1] = v[1];
//   angles[2] = v[2];
// }

double uniform_rand() {
  return (double)std::rand() / RAND_MAX;
}

Eigen::Vector3d QuaternionToEulerAngles(const Eigen::Vector4d &q) {
  Quaterniond quat(q[0], q[1], q[2], q[3]);
  return quat.matrix().eulerAngles(0, 1, 2);
}

void EulerAnglesToQuaternion(const double* const angles,
                             double* q) {
  Quaterniond quat = AngleAxisd(angles[0], Vector3d::UnitX()) *
      AngleAxisd(angles[1], Vector3d::UnitY()) *
      AngleAxisd(angles[2], Vector3d::UnitZ());
  q[0] = quat.w();
  q[1] = quat.x();
  q[2] = quat.y();
  q[3] = quat.z();
}

} // tiny
} // furry

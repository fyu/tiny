#ifndef TINY_NUMERIC_H
#define TINY_NUMERIC_H

#include <Eigen/Dense>

namespace furry {
namespace tiny {

double rand();

typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix3d;


double uniform_rand();

// template <class T, int n>
// Eigen::Matrix<T, n, 1> ToEigen(const rushmore::numeric::Array1D<T, n>& array) {
//   Eigen::Matrix<T, n, 1> eigen_vec;
//   array.CopyTo(eigen_vec.data());
//   return eigen_vec;
// }

// template <class T, int m, int n, int lay>
// Eigen::Matrix<T, m, n, Eigen::RowMajor> ToEigen(
//     const rushmore::numeric::Array2D<T, m, n, lay>& matrix) {
//   Eigen::Matrix<T, m, n, Eigen::RowMajor> eigen_matrix;
//   matrix.CopyToRowMajor(eigen_matrix.data());
//   return eigen_matrix;
// }

// Eigen::Vector3d ToEigen(const lightfield_sfm::Point& p);

// template <class T, int n>
// rushmore::numeric::Array1D<T, n> ToRmArray(
//     const Eigen::Matrix<T, n, 1>& array) {
//   rushmore::numeric::Array1D<T, n> result;
//   result.CopyFrom(array.data());
//   return result;
// }

// lightfield_sfm::Point ToLfPoint(const Eigen::Vector3d& v);

// void RotationToQuaternion(const Matrix3d& m, rushmore::numeric::VectorD4 *q);

// void QuaternionToEulerAngles(const rushmore::numeric::VectorD4& q,
//                              double* angles);

Eigen::Vector3d QuaternionToEulerAngles(const Eigen::Vector4d &q);

void EulerAnglesToQuaternion(const double* const angles,
                             double* q);

template <typename T>
void EulerAnglesToRotationMatrix(const T* const angles, T R[][3]) {
  T x = angles[0];
  T y = angles[1];
  T z = angles[2];
  R[0][0] = cos(y) * cos(z);
  R[0][1] = - cos(y) * sin(z);
  R[0][2] = sin(y);
  R[1][0] = cos(z) * sin(x) * sin(y) + cos(x) * sin(z);
  R[1][1] = cos(x) * cos(z) - sin(x) * sin(y) * sin(z);
  R[1][2] = - cos(y) * sin(x);
  R[2][0] = - cos(x) * cos(z) * sin(y) + sin(x) * sin(z);
  R[2][1] = cos(z) * sin(x) + cos(x) * sin(y) * sin(z);
  R[2][2] = cos(x) * cos(y);
}

template <typename T>
void EulerAnglesToAppRotationMatrix(const T* const angles, T R[][3]) {
  R[0][0] = T(1);
  R[0][1] = - angles[2];
  R[0][2] = angles[1];
  R[1][0] = angles[2];
  R[1][1] = T(1);
  R[1][2] = - angles[0];
  R[2][0] = - angles[1];
  R[2][1] = angles[0];
  R[2][2] = T(1);
}

} // tiny
} // furry

#endif

#ifndef FURRY_COMMON_MATH
#define FURRY_COMMON_MATH

#include <limits>
#include <opencv2/core/core.hpp>

namespace furry
{

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const double NaN = std::numeric_limits<double>::quiet_NaN();
const cv::Point2f nan_point = cv::Point2f(NaN, NaN);
const double infinity = std::numeric_limits<double>::infinity();

template <typename T> bool
equal(const T& v1, const T& v2)
{
  return v1 == v2;
}

template <> inline bool
equal<double>(const double& v1, const double& v2)
{
  return fabs(v1 - v2) < std::numeric_limits<double>::epsilon();
}

template <> inline bool
equal<float>(const float& v1, const float& v2)
{
  return fabs(v1 - v2) < std::numeric_limits<float>::epsilon();
}

template <typename T> bool
IsZero(const T& x) {
  static const T pos_epsilon = T(1e-12);
  static const T neg_epsilon = - pos_epsilon;
  if (x > neg_epsilon && x < pos_epsilon) {
    return true;
  } else {
    return false;
  }
}

template <> inline bool
IsZero(const double& v) {
  return equal<double>(v, 0);
}

template <> inline bool
IsZero(const float& v) {
  return equal<float>(v, 0);
}

template <typename T> bool
IsNonZero(const T& v) {
  return !IsZero(v);
}

// template <typename T, int r, int c>
// cv::Matx<T, r, c> quat2mat(T x, T y, T z, T w)
// {
//   static_assert(r >= 3 && c>= 3,
//                 "Rotation matrix must be at least 3 by 3");

//   cv::Matx<T, r, c> m = cv::Matx<T, r, c>::eye();
//   m(0, 0) = 1 - 2 * (y * y + z * z);
//   m(0, 1) = 2 * (x * y - z * w);
//   m(0, 2) = 2 * (x * z + y * w);
//   m(1, 0) = 2 * (x * y + z * w);
//   m(1, 1) = 1 - 2 * (x * x + z * z);
//   m(1, 2) = 2 * (y * z - x * w);
//   m(2, 0) = 2 * (x * z - y * w);
//   m(2, 1) = 2 * (y * z + x * w);
//   m(2, 2) = 1 - 2 * (x * x + y * y);

//   return m;
// }

/**
 * find the real roots of the cubic polynomial
 * coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 */
std::vector<double>
cubic_roots(const double coeff[4]);

float
cube_root(float x);

double
cube_root(double x);

template <typename T> T
square(const T &v) {
  return v * v;
}

template <int Size, typename T> void
Normalize(T *v) {
  T scale = T(0);
  for (int i = 0; i < Size; ++i)
    scale += v[i] * v[i];
  scale = 1 / sqrt(scale);
  for (int i = 0; i < Size; ++i)
    v[i] = v[i] * scale;
}

} // furry

#endif // FURRY_COMMON_MATH

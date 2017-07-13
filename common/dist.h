#ifndef FURRY_COMMON_DIST_H_
#define FURRY_COMMON_DIST_H_

#include <vector>
#include <cmath>

#include <Eigen/Dense>

namespace furry {

template <typename T>
double Distance2(const std::vector<T> &v0, const std::vector<T> &v1) {
  double dist = 0;
  assert(v0.size() == v1.size());
  for (size_t i = 0; i < v0.size(); ++i) {
    dist += (v0[i] - v1[i]) * (v0[i] - v1[i]);
  }
  return dist;
}

template <typename T>
double Distance(const std::vector<T> &v0, const std::vector<T> &v1) {
  double dist = 0;
  assert(v0.size() == v1.size());
  for (size_t i = 0; i < v0.size(); ++i) {
    dist += (v0[i] - v1[i]) * (v0[i] - v1[i]);
  }
  return std::sqrt(dist);
}

template <typename Derived> typename Derived::Scalar
Distance2(const Eigen::MatrixBase<Derived> &p0,
          const Eigen::MatrixBase<Derived> &p1) {
  auto diff = p0 - p1;
  return diff.cwiseProduct(diff).sum();
}

template <typename Derived> double
Distance(const Eigen::MatrixBase<Derived> &p0,
         const Eigen::MatrixBase<Derived> &p1) {
  return sqrt((double)Distance2(p0, p1));
}

} // furry

#endif

#ifndef FURRY_ALG_EUCLIDEAN_H_
#define FURRY_ALG_EUCLIDEAN_H_

#include <cmath>

#include <Eigen/Dense>

namespace furry {

template <typename Derived> Derived
LinearInterpolate(const Eigen::MatrixBase<Derived> &p0,
                  const Eigen::MatrixBase<Derived> &p1,
                  double t) {
  typedef typename Derived::Scalar Scalar;
  return Scalar(1 - t) * p0 + Scalar(t) * p1;
}

} // furry

#endif // FURRY_ALG_EUCLIDEAN_H_

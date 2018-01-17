#ifndef FURRY_COMMON_EIGEN_H_
#define FURRY_COMMON_EIGEN_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include <boost/functional/hash.hpp>

#include "furry/common/cast.h"

namespace Eigen {

// typedef Matrix<uint8_t, 3, 1> Vector3ub;
typedef Matrix<uint8_t, 3, 1> Vector3ub;

} // Eigen

namespace furry {

template <typename Scalar>
std::ostream& operator << (std::ostream &os,
                           const Eigen::Quaternion<Scalar>& q) {
  os << q.w() << ' ' << q.x() << ' ' << q.y() << ' ' << q.z();
  return os;
}

template <typename T0, typename T1>
struct ToImpl<std::true_type,
              Eigen::Matrix<T0, 2, Eigen::Dynamic>,
              std::vector<Eigen::Matrix<T1, 2, 1>>> {
  static Eigen::Matrix<T0, 2, Eigen::Dynamic> cast(
      const std::vector<Eigen::Matrix<T1, 2, 1>> &ps) {
    Eigen::Matrix<T0, 2, Eigen::Dynamic> m(2, ps.size());
    for (size_t i = 0; i < ps.size(); ++i) {
      m.col(i) = ps[i];
    }
    return m;
  }
};

} // furry

namespace std {
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
class numeric_limits<Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> > {
  typedef Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> DataType;
  typedef typename DataType::Scalar Scalar;
 public:
  static DataType max() {
    return DataType::Constant(std::numeric_limits<Scalar>::max());
  }
  static DataType lowest() {
    return DataType::Constant(std::numeric_limits<Scalar>::lowest());
  }
};
}

namespace std {

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct hash<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> > {
  size_t operator () (const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> &m) const {
    size_t seed = 0;
    std::hash<_Scalar> hasher;
    for (int i = 0; i < m.rows(); ++i) {
      for (int j = 0; j < m.cols(); ++j) {
        boost::hash_combine(seed, hasher(m(i, j)));
      }
    }
    return seed;
  }
};

} // std

namespace boost {
namespace serialization {

template <class Archive, class Derived>
void serialize(Archive &ar, const Eigen::MatrixBase<Derived> &m_,
               const unsigned int version) {
  Eigen::MatrixBase<Derived> &m = const_cast<Eigen::MatrixBase<Derived>>(m_);
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j) {
      ar & m.coeff(i, j);
    }
  }
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(Archive &ar,
               Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> &m,
               const unsigned int version) {
  // _Scalar *data = m.data();
  // int size = m.rows() * m.cols();
  // std::for_each(data, data + size,
  //               [&ar] (_Scalar &v) {
  //                 ar & v;
  //               });
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j) {
      ar & m(i, j);
    }
  }
}

} // serialization
} // boost

#endif

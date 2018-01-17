#ifndef FURRY_COMMON_POINT_H_
#define FURRY_COMMON_POINT_H_

#include <stdexcept>
#include <type_traits>
#include <iostream>
#include <limits>

#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>

// #include "R2Shapes/R2Shapes.h"
// #include "R3Shapes/R3Shapes.h"

#include "furry/common/define.h"

namespace furry {

typedef float Scalar2d;
// typedef float Scalar3d;
typedef double Scalar3d;
typedef cv::Point_<Scalar2d> CvPoint2;
typedef CvPoint2 Point2;
typedef Eigen::Matrix<Scalar2d, 2, 1> Vector2;
typedef Eigen::Matrix<Scalar3d, 3, 1> Vector3;
typedef Vector3 Point3;

std::istream& operator >> (std::istream& is, Point3 &p);

#if 0

struct Point3 {
  Point3() {}

  Point3(double x, double y, double z) {
    d[0] = x;
    d[1] = y;
    d[2] = z;
  }

  Point3(const R3Point &p) {
    d[0] = p.X();
    d[1] = p.Y();
    d[2] = p.Z();
  }

  double& x() {
    return d[0];
  }

  double x() const {
    return d[0];
  }

  double& y() {
    return d[1];
  }

  double y() const {
    return d[1];
  }

  double& z() {
    return d[2];
  }

  double z() const {
    return d[2];
  }

  void*& data() {
    return data_ptr;
  }

  const void* data() const {
    return data_ptr;
  }

  void set_data(void* data_ptr_) {
    data_ptr = data_ptr_;
  }

  double& operator [] (int i) {
    return d[i];
  }

  double operator [] (int i) const {
    return d[i];
  }

  double d[3] = {-1, -1, -1};
  void *data_ptr = nullptr;
};

std::ostream& operator << (std::ostream& os, const Point3 &p);
bool operator == (const Point3 &p0, const Point3 &p1);

struct Point2 {

  typedef double Scalar;

  Point2() {}

  Point2(Scalar x, Scalar y) {
    d[0] = x;
    d[1] = y;
  }

  Point2(const R2Point &p) {
    d[0] = p.X();
    d[1] = p.Y();
  }

  Point2(const Point2 &p) {
    d[0] = p.d[0];
    d[1] = p.d[1];
  }

  Scalar& x() {
    return d[0];
  }

  Scalar x() const {
    return d[0];
  }

  Scalar& y() {
    return d[1];
  }

  Scalar y() const {
    return d[1];
  }

  void*& data() {
    return data_ptr;
  }

  const void* data() const {
    return data_ptr;
  }

  void set_data(void* data_ptr_) {
    data_ptr = data_ptr_;
  }

  operator cv::Point () const {
    return cv::Point(x() + 0.5, y() + 0.5);
  }

  Scalar& operator [] (int i) {
    return d[i];
  }

  Scalar operator [] (int i) const {
    return d[i];
  }

  Scalar d[2] = {-1, -1};
  void *data_ptr = nullptr;
};



double Distance2(const Point2 &p0, const Point2 &p1);
std::ostream& operator << (std::ostream& os, const Point2 &p);
std::istream& operator >> (std::istream& is, Point2 &p);

#endif

template <typename PointType>
struct PointTraits;

template <typename Scalar, int Rows, int Flags, int MaxRows, int MaxCols>
struct PointTraits<Eigen::Matrix<Scalar, Rows, 1, Flags, MaxRows, MaxCols>> {

  typedef Eigen::Matrix<Scalar, Rows, 1> PointType;
  typedef Scalar ValueType;
  static const size_t kNumDimensions = Rows;

  struct X {
    Scalar operator () (const PointType &p) {
      return p.x();
    }
  };
  struct Y {
    Scalar operator () (const PointType &p) {
      return p.y();
    }
  };
  struct Z {
    Scalar operator () (const PointType &p) {
      return p.z();
    }
  };
  struct W {
    Scalar operator () (const PointType &p) {
      return p.w();
    }
  };

  struct Get {
    Scalar operator () (const PointType &p, size_t index) {
      return p[index];
    }
  };

  struct SetData {
    void operator () (PointType &p, const ValueType *values) {
      p = PointType(values);
    }
  };

  struct Data {
    const ValueType* operator () (const PointType *p) {
      return p->data();
    }
  };
};

template <typename Scalar>
struct PointTraits<cv::Point_<Scalar>> {

  typedef cv::Point_<Scalar> PointType;
  typedef Scalar ValueType;

  static const size_t kNumDimensions = 2;

  struct X {
    Scalar operator () (const PointType &p) {
      return p.x;
    }
  };
  struct Y {
    Scalar operator () (const PointType &p) {
      return p.y;
    }
  };

  struct Get {
    Scalar operator () (const PointType &p, size_t index) {
      // switch (index) {
      //   case 0:
      //     return p.x;
      //   case 1:
      //     return p.y;
      //   default:
      //     throw std::out_of_range();
      // }
      PointType p_ = p;
      return reinterpret_cast<Scalar*>(&p_)[index];
    }
  };

  struct SetData {
    void operator () (PointType &p, const ValueType *values) {
      p.x = values[0];
      p.y = values[1];
    }
  };

  struct Data {
    const ValueType* operator () (const PointType *p) {
      return reinterpret_cast<const ValueType*>(p);
    }
  };
};

template <typename PointType, typename OutputIterator> int
ReadPoints(const std::string &file_name,
           OutputIterator result,
           void *buffer = NULL,
           size_t num_bytes = 0);

// On successful write, return the number of written points
// On failed write, return -1
template <typename InputIterator> int
WritePoints(const std::string &file_name,
            InputIterator first,
            InputIterator last,
            void *buffer = NULL,
            size_t num_bytes = 0);

int WritePoints(const std::string &filename,
                const std::vector<CvPoint2> &points);

int ReadPoints(const std::string &filename,
               std::vector<CvPoint2> *points);

int WritePoints(const std::string &filename,
                const std::vector<Vector2> &points);

int ReadPoints(const std::string &filename,
               std::vector<Vector2> *points);

int WritePoints(const std::string &filename,
                const std::vector<Vector3> &points);

int ReadPoints(const std::string &filename,
               std::vector<Vector3> *points);

int ReadPoints(const std::string &filename,
               std::vector<cv::KeyPoint> *points);

int WritePoints(const std::string &filename,
                const std::vector<cv::KeyPoint> &points);


int WritePoints(const std::vector<CvPoint2> &points, FILE *file);

int ReadPoints(FILE *file, std::vector<CvPoint2> *points);

int WritePoints(const std::vector<Vector2> &points, FILE *file);

int ReadPoints(FILE *file, std::vector<Vector2> *points);

int WritePoints(const std::vector<Vector3> &points, FILE *file);

int ReadPoints(FILE *file, std::vector<Vector3> *points);

int WritePoints(const std::vector<cv::KeyPoint> &points, FILE *file);

int ReadPoints(FILE *file, std::vector<cv::KeyPoint> *points);

template <typename PointType0, typename PointType1>
typename std::enable_if<PointTraits<PointType0>::kNumDimensions ==
                        PointTraits<PointType1>::kNumDimensions,
                        typename PointTraits<PointType0>::ValueType>::type
EuclideanDistance(const PointType0 &p0, const PointType1 &p1);

// Square of Euclidean distance
template <typename PointType0, typename PointType1>
typename std::enable_if<PointTraits<PointType0>::kNumDimensions ==
                        PointTraits<PointType1>::kNumDimensions,
                        typename PointTraits<PointType0>::ValueType>::type
EuclideanDistance2(const PointType0 &p0, const PointType1 &p1);

template <typename InputIterator, typename PointType> inline InputIterator
NearestNeighbor(InputIterator first, InputIterator last, const PointType &p);

} // furry

namespace std {

template <>
class numeric_limits<furry::Point2> {
typedef furry::Scalar2d Scalar;
 public:
  static furry::Point2 max() {
    return furry::Point2(numeric_limits<Scalar>::max(),
                  numeric_limits<Scalar>::max());
  }

  static furry::Point2 lowest() {
    return furry::Point2(numeric_limits<Scalar>::lowest(),
                  numeric_limits<Scalar>::lowest());
  }
};

template <>
struct hash<furry::Point3> {
  size_t operator () (const furry::Point3 &p) const {
    size_t seed = 0;
    std::hash<double> hash_double;
    boost::hash_combine(seed, hash_double(p.x()));
    boost::hash_combine(seed, hash_double(p.y()));
    boost::hash_combine(seed, hash_double(p.z()));
    return seed;
  }
};

} // std

#include "furry/common/point-inl.h"

#endif // FURRY_COMMON_POINT_H_

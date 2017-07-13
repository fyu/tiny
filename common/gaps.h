#ifndef FURRY_COMMON_GAPS_H_
#define FURRY_COMMON_GAPS_H_

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "opencv2/opencv.hpp"
#include "R2Shapes/R2Shapes.h"
#include "R3Shapes/R3Shapes.h"
#include "RNBasics/RNBasics.h"
#include "furry/common/cast.h"
#include "furry/common/point.h"

namespace furry {

cv::Mat R2Image2CvImage(const R2Image *image);

template <typename PtrType>
struct ToImpl<std::true_type, RNArray<PtrType>, std::vector<PtrType>> {
  static RNArray<PtrType> cast(const std::vector<PtrType> &values) {
    RNArray<PtrType> array;
    for (auto &value : values) {
      array.InsertTail(value);
    }
    return array;
  }
};

template <>
struct ToImpl<std::true_type, cv::Mat, R2Image*> {
  static cv::Mat cast(const R2Image *image) {
    return R2Image2CvImage(image);
  }
};

template <>
struct ToImpl<std::true_type, cv::Mat, R2Image> {
  static cv::Mat cast(const R2Image &image) {
    return R2Image2CvImage(&image);
  }
};

// template <>
// struct ToImpl<std::true_type, R2Point, Point2> {
//   static R2Point cast(const Point2 &p) {
//     return R2Point(p.x(), p.y());
//   }
// };

template <>
struct ToImpl<std::true_type, Point2, R2Point> {
  static Point2 cast(const R2Point &p) {
    return Point2(p.X(), p.Y());
  }
};

template <typename Scalar>
struct ToImpl<std::true_type, R2Point, cv::Point_<Scalar>> {
  static R2Point cast(const cv::Point_<Scalar> &p) {
    return R2Point(p.x, p.y);
  }
};

template <typename Scalar>
struct ToImpl<std::true_type, cv::Point_<Scalar>, R2Point> {
  static cv::Point_<Scalar> cast(const R2Point &p) {
    return cv::Point_<Scalar>(To<Scalar>(p.X()), To<Scalar>(p.Y()));
  }
};

template <>
struct ToImpl<std::true_type, R2Point, cv::KeyPoint> {
  static R2Point cast(const cv::KeyPoint &p) {
    return R2Point(p.pt.x, p.pt.y);
  }
};

template <>
struct ToImpl<std::true_type, Eigen::Vector2d, R2Point> {
  static Eigen::Vector2d cast(const R2Point &p) {
    return Eigen::Vector2d(p.X(), p.Y());
  }
};

template <>
struct ToImpl<std::true_type, Eigen::Vector2f, R2Point> {
  static Eigen::Vector2f cast(const R2Point &p) {
    return Eigen::Vector2f(p.X(), p.Y());
  }
};

template <typename T>
struct ToImpl<std::true_type, Eigen::Matrix<T, 2, 1>, R2Vector> {
  static Eigen::Matrix<T, 2, 1> cast(const R2Vector &p) {
    return Eigen::Matrix<T, 2, 1>(p.X(), p.Y());
  }
};

template <typename T>
struct ToImpl<std::true_type, R2Point, Eigen::Matrix<T, 2, 1>> {
  static R2Point cast(const Eigen::Matrix<T, 2, 1> &p) {
    return R2Point(p.x(), p.y());
  }
};

template <typename T>
struct ToImpl<std::true_type, R2Vector, Eigen::Matrix<T, 2, 1>> {
  static R2Vector cast(const Eigen::Matrix<T, 2, 1> &p) {
    return R2Vector(p.x(), p.y());
  }
};

////////////////////////////////////////////////////////////
// R3Shapes
////////////////////////////////////////////////////////////

template <typename T>
struct ToImpl<std::true_type, Eigen::Matrix<T, 3, 1>, R3Point> {
  static Eigen::Matrix<T, 3, 1> cast(const R3Point &p) {
    return Eigen::Matrix<T, 3, 1>(p.X(), p.Y(), p.Z());
  }
};

template <typename T>
struct ToImpl<std::true_type, R3Point, Eigen::Matrix<T, 3, 1>> {
  static R3Point cast(const Eigen::Matrix<T, 3, 1> &p) {
    return R3Point(p.x(), p.y(), p.z());
  }
};

template <typename T>
struct ToImpl<std::true_type, Eigen::Matrix<T, 3, 1>, R3Vector> {
  static Eigen::Matrix<T, 3, 1> cast(const R3Vector &p) {
    return Eigen::Matrix<T, 3, 1>(p.X(), p.Y(), p.Z());
  }
};

template <typename T>
struct ToImpl<std::true_type, R3Vector, Eigen::Matrix<T, 3, 1>> {
  static R3Vector cast(const Eigen::Matrix<T, 3, 1> &p) {
    return R3Vector(p[0], p[1], p[2]);
  }
};

template <>
struct ToImpl<std::true_type, Eigen::Vector4d, R3Quaternion> {
  static Eigen::Vector4d cast(const R3Quaternion &q) {
    return Eigen::Vector4d(q.A(), q.B(), q.C(), q.D());
  }
};

template <>
struct ToImpl<std::true_type, R3Quaternion, Eigen::Vector4d> {
  static R3Quaternion cast(const Eigen::Vector4d &q) {
    return R3Quaternion(q[0], q[1], q[2], q[3]);
  }
};

template <>
struct ToImpl<std::true_type, R3Point, R3Vector> {
  static R3Point cast(const R3Vector &v) {
    return R3Point(v.X(), v.Y(), v.Z());
  }
};

template <>
struct ToImpl<std::true_type, Eigen::Quaterniond, R3Quaternion> {
  static Eigen::Quaterniond cast(const R3Quaternion &q) {
    return Eigen::Quaterniond(q.A(), q.B(), q.C(), q.D());
  }
};

R3Quaternion
Conjugate(const R3Quaternion &q);

R3Quaternion
RotateBetween(const R3Quaternion &from,
              const R3Quaternion &to);

R3Vector
Rotate(const R3Quaternion &q, const R3Vector &v);

R3Point
Rotate(const R3Quaternion &q, const R3Point &p);

Eigen::Vector3d R3QuaternionToAngleAxis(const R3Quaternion &q);

double CalcDistance2(const R3Point &p0, const R3Point &p1);

double CalcDistance(const R3Point &p0, const R3Point &p1);

} // furry

////////////////////////////////////////////////////////////
// IO
////////////////////////////////////////////////////////////

std::ostream& operator << (std::ostream &os, const R3Plane &p);
std::ostream& operator << (std::ostream &os, const R3Point &p);
std::ostream& operator << (std::ostream &os, const R2Point &p);
std::ostream& operator << (std::ostream &os, const R3Quaternion &q);
std::ostream& operator << (std::ostream &os, const R3Vector &v);

#endif // FURRY_COMMON_GAPS_H_

#ifndef FURRY_COMMON_CV_H_
#define FURRY_COMMON_CV_H_

#include <cmath>

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include <functional>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/function_types/function_arity.hpp>

#include "furry/common/cast.h"
#include "furry/common/factory.h"
#include "furry/common/numeric.h"
#include "furry/common/point.h"

namespace furry
{

extern const cv::Scalar kCvWhite;
extern const cv::Scalar kCvBlack;
extern const cv::Scalar kCvBlue;
extern const cv::Scalar kCvRed;
extern const cv::Scalar kCvGreen;
extern const cv::Scalar kCvMagenta;
extern const cv::Scalar kCvCyan;
extern const cv::Scalar kCvYellow;
extern const cv::Scalar kCvGray;
extern const cv::Scalar kCvLightGray;

extern const int kCvColorPlatteSize;
extern const cv::Scalar kCvColorPlatte[6];

////////////////////////////////////////////////////////////
// Cast
////////////////////////////////////////////////////////////

template <> inline cv::Point
to<cv::Point, Eigen::Vector2d>(const Eigen::Vector2d& p) {
  return cv::Point(p.x() + 0.5, p.y() + 0.5);
}

template <> inline cv::Point2d
to<cv::Point2d, Eigen::Vector2d>(const Eigen::Vector2d &p) {
  return cv::Point2d(p.x(), p.y());
}

template <> inline cv::Point2f
to<cv::Point2f, Eigen::Vector2f>(const Eigen::Vector2f &p) {
  return cv::Point2f(p.x(), p.y());
}

template <> inline cv::Point2f
to<cv::Point2f, Eigen::Vector2d>(const Eigen::Vector2d &p) {
  return cv::Point2f(p.x(), p.y());
}

template <> inline Eigen::Vector2d
to<Eigen::Vector2d, cv::Point2d>(const cv::Point2d &p) {
  return Eigen::Vector2d(p.x, p.y);
}

template <> inline Eigen::Vector2f
to<Eigen::Vector2f, cv::Point2f>(const cv::Point2f &p) {
  return Eigen::Vector2f(p.x, p.y);
}

template <> inline Eigen::Vector2i
to<Eigen::Vector2i, cv::Point>(const cv::Point &p) {
  return Eigen::Vector2i(p.x, p.y);
}

template <> inline Eigen::Vector2d
to<Eigen::Vector2d, cv::Point2f>(const cv::Point2f &p) {
  return Eigen::Vector2d(p.x, p.y);
}

template <>
struct ToImpl<std::true_type, cv::Vec3b, cv::Scalar> {
  static cv::Vec3b cast(const cv::Scalar &c) {
    return cv::Vec3b(c[0], c[1], c[2]);
  }
};

template <typename TargetScalar, typename SourceScalar>
struct ToImpl<typename std::enable_if<!std::is_same<TargetScalar, SourceScalar>::value, std::true_type>::type, cv::Point_<TargetScalar>, cv::Point_<SourceScalar>> {
// struct ToImpl<cv::Point_<TargetScalar>, typename std::enable_if<!std::is_same<TargetScalar, SourceScalar>::value, cv::Point_>::type <SourceScalar>> {
  // static std::enable_if<!std::is_same<TargetScalar, SourceScalar>, cv::Point_<TargetScalar> >
  static cv::Point_<TargetScalar>
      cast(const cv::Point_<SourceScalar> &p) {
    return cv::Point_<TargetScalar>(To<TargetScalar>(p.x),
                                    To<TargetScalar>(p.y));
  }
};

// template <typename Scalar>
// struct ToImpl<std::true_type, Point2, cv::Point_<Scalar>> {
//   static Point2 cast(const cv::Point_<Scalar> &p) {
//     return Point2(p.x, p.y);
//   }
// };

template <typename Scalar>
struct ToImpl<std::true_type, cv::KeyPoint, Point2, Scalar> {
  static cv::KeyPoint cast(const Point2 &p, Scalar diameter) {
    return cv::KeyPoint(To<cv::Point2f>(p), diameter);
  }
};

template <typename Scalar>
struct ToImpl<std::true_type, cv::Point_<Scalar>, cv::KeyPoint> {
  static cv::Point_<Scalar> cast(const cv::KeyPoint &p) {
    return To<cv::Point_<Scalar>>(p.pt);
  }
};

template <>
struct ToImpl<std::true_type, Point2, cv::KeyPoint> {
  static Point2 cast(const cv::KeyPoint &p) {
    return Point2(p.pt.x, p.pt.y);
  }
};

template <typename T0, typename T1>
struct ToImpl<std::true_type, Eigen::Matrix<T0, 2, 1>, cv::Point_<T1>> {
  static Eigen::Matrix<T0, 2, 1> cast(const cv::Point_<T1> &p) {
    return Eigen::Matrix<T0, 2, 1>(To<T0>(p.x), To<T0>(p.y));
  }
};

template <typename T0, typename T1>
struct ToImpl<std::true_type, cv::Point_<T0>, Eigen::Matrix<T1, 2, 1>> {
  static cv::Point_<T0> cast (const Eigen::Matrix<T1, 2, 1> &p) {
    return cv::Point_<T0>(To<T0>(p.x()), To<T0>(p.y()));
  }
};

template <typename T0, typename T1>
struct ToImpl<std::true_type, Eigen::Matrix<T0, 2, Eigen::Dynamic>,
              std::vector<cv::Point_<T1>>> {
  static Eigen::Matrix<T0, 2, Eigen::Dynamic> cast(const std::vector<cv::Point_<T1>> &ps) {
    Eigen::Matrix<T0, 2, Eigen::Dynamic> m(2, ps.size());
    for (size_t i = 0; i < ps.size(); ++i) {
      m.col(i) = To<Eigen::Matrix<T0, 2, 1>>(ps[i]);
    }
    return m;
  }
};

template <typename T, int n>
struct ToImpl<std::true_type, Eigen::Matrix<T, n, 1>, cv::Vec<T, n>> {
  static Eigen::Matrix<T, n, 1> cast (const cv::Vec<T, n> &v) {
    // return cv::Point_<T0>(To<T0>(p.x()), To<T0>(p.y()));
    Eigen::Matrix<T, n, 1> m;
    for (int i = 0; i < n; ++i) {
      m[i] = v[i];
    }
    return m;
  }
};

// template <typename T, typename Derived>
// struct ToImpl<
//   typename std::enable_if<(Derived::rows == 2 && Derived::cols == 1) ||
//     (Derived::cols == 2 && Derived::rows == 1), std::true_type>::type, cv::Point_<T>, Eigen::MatrixBase<Derived>> {
//   static cv::Point_<T> cast (const Eigen::MatrixBase<Derived> &p) {
//     return cv::Point_<T>(To<T>(p.x()), To<T>(p.y()));
//   }
// };

template <typename Scalar, typename SizeType>
struct MakeImpl<typename std::enable_if<
                std::is_arithmetic<SizeType>::value, cv::KeyPoint>::type,
              cv::Point_<Scalar>,
              SizeType> {
  static cv::KeyPoint Do(const cv::Point_<Scalar> &p, SizeType size) {
    return cv::KeyPoint(p.x, p.y, size);
  }
};

template <typename Scalar>
cv::Point_<Scalar> FlipY(const cv::Point_<Scalar> &p, double height) {
  return cv::Point_<Scalar>(p.x, height - p.y);
}

cv::KeyPoint FlipY(const cv::KeyPoint &p, double height);

////////////////////////////////////////////////////////////
//  Feature 2D
////////////////////////////////////////////////////////////

std::ostream& operator << (std::ostream &os, const cv::KeyPoint &p);

////////////////////////////////////////////////////////////
// Rect
////////////////////////////////////////////////////////////

template<typename _Tp> bool
operator <= (const cv::Rect_<_Tp>& r1, const cv::Rect_<_Tp>& r2)
{
  return (r1 & r2) == r1;
}

template <typename DataType> std::ostream&
operator << (std::ostream& os, const cv::Rect_<DataType>& r)
{
  os << r.x << " " << r.y << " " << r.width << " " << r.height;
  return os;
}

///////////////////////////////////////////////////////////
// Vec
////////////////////////////////////////////////////////////

template<typename DataType, int n> std::ostream&
operator << (std::ostream& os, const cv::Vec<DataType, n>& v)
{
  os << v[0];
  for (int i = 1; i < n; ++i)
  {
    os << ' ' << v[i];
  }
  return os;
}

double RgbToLuminosity(const cv::Vec3b &c);

double BgrToLuminosity(const cv::Vec3b &c);

////////////////////////////////////////////////////////////
// Point_
////////////////////////////////////////////////////////////

template <> inline cv::Point
to<cv::Point>(const cv::Point2d& p) {
  return cv::Point(p.x + 0.5, p.y + 0.5);
}

template <> inline cv::Point
to<cv::Point>(const cv::Point2f &p) {
  return cv::Point(p.x + 0.5, p.y + 0.5);
}

template <typename T> cv::Point_<T>
operator / (const cv::Point_<T>& p, double d)
{
  return cv::Point_<T>(p.x / d, p.y / d);
}

template <typename T>
cv::Point_<T> operator + (const cv::Point_<T> &p, const cv::Vec<T, 2> &v) {
  return cv::Point_<T>(p.x + v[0], p.y + v[1]);
}

template <typename T>
cv::Point_<T> operator + (const cv::Vec<T, 2> &v, const cv::Point_<T> &p) {
  return cv::Point_<T>(p.x + v[0], p.y + v[1]);
}

inline cv::Point2f operator + (const cv::Point2f &p, const cv::Vec2f &v) {
  return cv::Point2f(p.x + v[0], p.y + v[1]);
}

inline cv::Point2f operator - (const cv::Point2f &p, const cv::Vec2f &v) {
  return cv::Point2f(p.x - v[0], p.y - v[1]);
}


template <> bool
equal<cv::Point2d>(const cv::Point2d& p0, const cv::Point2d& p1);

std::ostream& operator << (std::ostream& os, const cv::Point2d& p);
std::ostream& operator << (std::ostream& os, const cv::Point3d& p);

template <typename Tp> Tp
Distance2(const cv::Point_<Tp> &p0, const cv::Point_<Tp> &p1) {
  return (p0.x - p1.x) * (p0.x - p1.x) + (p0.y - p1.y) * (p0.y - p1.y);
}

template <typename Tp> Tp
Distance(const cv::Point_<Tp> &p0, const cv::Point_<Tp> &p1) {
  return std::sqrt(Distance2<Tp>(p0, p1));
}

template <typename Tp> Tp
CalcDistance(const cv::Point_<Tp> &p0, const cv::Point_<Tp> &p1) {
  return std::sqrt(Distance2<Tp>(p0, p1));
}

template <typename T> Eigen::Matrix<T, 2, 1>
Cv2Eigen(const cv::Point_<T> &p) {
  return Eigen::Matrix<T, 2, 1>(p.x, p.y);
}

template <typename T> cv::Point_<T>
Eigen2Cv(const Eigen::Matrix<T, 2, 1> &p) {
  return cv::Point_<T>(p.x(), p.y());
}

template <typename T> std::vector<Eigen::Matrix<T, 2, 1> >
Cv2Eigen(const std::vector<cv::Point_<T>> &points) {
  std::vector<Eigen::Matrix<T, 2, 1>> result;
  for (auto it = points.begin(); it != points.end(); ++it)
    result.push_back(Cv2Eigen(*it));
  return result;
}

template <typename T> std::vector<cv::Point_<T>>
Eigen2Cv(const std::vector<Eigen::Matrix<T, 2, 1>> &points) {
  std::vector<cv::Point_<T>> result;
  for (auto it = points.begin(); it != points.end(); ++it)
    result.push_back(Eigen2Cv(*it));
  return result;
}

template <typename T>
struct ToImpl<std::true_type, Eigen::Matrix<T, 2, 1>, cv::Point_<T>> {
  static Eigen::Matrix<T, 2, 1> cast (const cv::Point_<T> &p) {
    return Eigen::Matrix<T, 2, 1>(p.x, p.y);
  }
};

////////////////////////////////////////////////////////////
// Scalar
////////////////////////////////////////////////////////////

cv::Scalar
operator / (const cv::Scalar& s0, const cv::Scalar& s1);

cv::Scalar
operator / (const cv::Scalar &s, double d);

cv::Scalar
operator * (const cv::Scalar& s0, const cv::Scalar& s1);

cv::Scalar
operator + (const cv::Scalar& s0, const cv::Scalar& s1);

cv::Scalar
operator - (const cv::Scalar& s0, const cv::Scalar& s1);

cv::Scalar
abs(const cv::Scalar& s);

double
sum(const cv::Scalar& s);

////////////////////////////////////////////////////////////
// Size
////////////////////////////////////////////////////////////

template <typename T> cv::Size_<T>
operator / (const cv::Size_<T>& s, double c) {
  return cv::Size_<T>(s.width / c, s.height / c);
}

template <typename T> cv::Size_<T>
operator / (const cv::Size_<T>& s, int c) {
  return cv::Size_<T>(s.width / c, s.height / c);
}

template <typename T> cv::Size_<T>
operator * (const cv::Size_<T> &s, double c) {
  return cv::Size_<T>(s.width * c, s.height * c);
}


////////////////////////////////////////////////////////////
// Mat
////////////////////////////////////////////////////////////

bool WriteImage(const cv::Mat &image, FILE *file);
bool ReadImage(FILE *file, cv::Mat *image);

template<typename _Tp, int _rows, int _cols, int _options, int _maxRows, int _maxCols>
struct ToImpl<std::true_type,
              Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols>,
              cv::Mat> {
  static Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols> cast(
      const cv::Mat& src) {
    Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols> dst;
    CV_DbgAssert(src.rows == _rows && src.cols == _cols);
    if( !(dst.Flags & Eigen::RowMajorBit) )
    {
      cv::Mat _dst(src.cols, src.rows, cv::DataType<_Tp>::type,
                   dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
      if( src.type() == _dst.type() )
        transpose(src, _dst);
      else if( src.cols == src.rows )
      {
        src.convertTo(_dst, _dst.type());
        transpose(_dst, _dst);
      }
      else
        cv::Mat(src.t()).convertTo(_dst, _dst.type());
      CV_DbgAssert(_dst.data == (uchar*)dst.data());
    }
    else
    {
      cv::Mat _dst(src.rows, src.cols, cv::DataType<_Tp>::type,
                   dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
      src.convertTo(_dst, _dst.type());
      CV_DbgAssert(_dst.data == (uchar*)dst.data());
    }
    return dst;
  }
};

template<typename _Tp, int _rows, int _cols, int _options, int _maxRows, int _maxCols>
struct ToImpl<std::true_type,
              cv::Mat,
              Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols>> {
  static cv::Mat cast (const Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols> &src) {
    cv::Mat dst;
    if( !(src.Flags & Eigen::RowMajorBit) )
    {
      cv::Mat _src(src.cols(), src.rows(), cv::DataType<_Tp>::type,
                   (void*)src.data(), src.stride()*sizeof(_Tp));
      transpose(_src, dst);
    }
    else
    {
      cv::Mat _src(src.rows(), src.cols(), cv::DataType<_Tp>::type,
                   (void*)src.data(), src.stride()*sizeof(_Tp));
      _src.copyTo(dst);
    }
    return dst;
  }
};

cv::Rect
imframe(const cv::Mat& image);

inline bool out_of_image(cv::Size image_size, cv::Point2d p)
{
	if (std::isnan(p.x) || std::isnan(p.y) ||
		p.x < 0 || p.x >= image_size.width || p.y < 0 || p.y >= image_size.height)
			return true;
		else
			return false;
}

template <typename T>
inline T m_at_bilinear(const cv::Mat & m, const cv::Point2d & p)
{
	if (out_of_image(m.size(), p)) return 0;
	double x1 = std::floor(p.x);
	double x2 = x1 + 1;
	double y1 = std::floor(p.y);
	double y2 = y1 + 1;
	T v1 = (x2 - p.x) * m.at<T>(cv::Point(x1, y1)) + (p.x - x1) * m.at<T>(cv::Point(x2, y1));
	T v2 = (x2 - p.x) * m.at<T>(cv::Point(x1, y2)) + (p.x - x1) * m.at<T>(cv::Point(x2, y2));
	return (y2 - p.y) * v1 + (p.y - y1) * v2;
}

template <typename ElemType, typename Function>
typename std::enable_if<boost::function_types::function_arity<Function>::value == 1>::type
Map(const cv::Mat_<ElemType> &src,
    Function func,
    cv::Mat_<typename boost::function_types::result_type<Function>::type> &dst) {
  for (int r = 0; r < src.rows; ++r) {
    for (int c = 0; c < src.cols; ++c) {
      dst(r, c) = func(src(r, c));
    }
  }

  return;
}

// template <typename ElemType, typename ReturnType, typename ParamType> void
// Map(const cv::Mat_<ElemType> &src,
//     ReturnType (*func) (ParamType),
//     cv::Mat_<ReturnType> &dst) {
//   for (int r = 0; r < src.rows; ++r) {
//     for (int c = 0; c < src.cols; ++c) {
//       dst(r, c) = func(src(r, c));
//     }
//   }

//   return;
// }

template <typename ReturnType, typename ParamType, int Channels> void
Map(cv::InputArray src,
    ReturnType (*func) (const cv::Vec<ParamType, Channels>&),
    cv::OutputArray dst) {
  typedef const cv::Vec<ParamType, Channels>& VecType;
  cv::Mat src_ = src.getMat();
  dst.create(src.size(), cv::DataType<ReturnType>::type);
  cv::Mat dst_ = dst.getMat();
  assert(src_.channels() == Channels);
  switch (src_.depth()) {
    case CV_8U:
      Map<cv::Vec<unsigned char, Channels>>(
          (const cv::Mat_<cv::Vec<unsigned char, Channels>>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_16S:
      Map<cv::Vec<short, Channels>>(
          (const cv::Mat_<cv::Vec<short, Channels>>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_32S:
      Map<cv::Vec<int, Channels>>(
          (const cv::Mat_<cv::Vec<int, Channels>>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_32F:
      Map<cv::Vec<float, Channels>>(
          (const cv::Mat_<cv::Vec<float, Channels>>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_64F:
      Map<cv::Vec<double, Channels>>(
          (const cv::Mat_<cv::Vec<double, Channels>>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    default:
      ;
  }
}

template <typename ReturnType, typename ParamType> void
Map(cv::InputArray src,
    ReturnType (*func) (ParamType),
    cv::OutputArray dst) {
  cv::Mat src_ = src.getMat();
  dst.create(src.size(), cv::DataType<ReturnType>::type);
  cv::Mat dst_ = dst.getMat();
  assert(src_.channels() == 1);
  switch (src_.depth()) {
    case CV_8U:
      Map<unsigned char>(
          (const cv::Mat_<unsigned char>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_16S:
      Map<short>(
          (const cv::Mat_<short>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_32S:
      Map<int>(
          (const cv::Mat_<int>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_32F:
      Map<float>(
          (const cv::Mat_<float>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    case CV_64F:
      Map<double>(
          (const cv::Mat_<double>&)src_,
          func,
          (cv::Mat_<ReturnType>&)dst_);
      break;
    default:
      return;
  }
  // switch (src_.channels()) {
  //   case 1:
  //     switch (src_.depth()) {
  //       case CV_8U:
  //         Map<unsigned char, ReturnType, ParamType>(
  //             (const cv::Mat_<unsigned char>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_16S:
  //         Map<short, ReturnType, ParamType>(
  //             (const cv::Mat_<short>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32S:
  //         Map<int, ReturnType, ParamType>(
  //             (const cv::Mat_<int>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32F:
  //         Map<float, ReturnType, ParamType>(
  //             (const cv::Mat_<float>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_64F:
  //         Map<double, ReturnType, ParamType>(
  //             (const cv::Mat_<double>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       default:
  //         return;
  //     }
  //     break;
  //   case 2:
  //     switch (src_.depth()) {
  //       case CV_8U:
  //         Map<cv::Vec<unsigned char, 2>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<unsigned char, 2>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_16S:
  //         Map<cv::Vec<short, 2>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<short, 2>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32S:
  //         Map<cv::Vec<int, 2>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<int, 2>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32F:
  //         Map<cv::Vec<float, 2>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<float, 2>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_64F:
  //         Map<cv::Vec<double, 2>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<double, 2>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       default:
  //         ;
  //     }
  //     break;
  //   case 3:
  //     switch (src_.depth()) {
  //       case CV_8U:
  //         Map<cv::Vec<unsigned char, 3>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<unsigned char, 3>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_16S:
  //         Map<cv::Vec<short, 3>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<short, 3>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32S:
  //         Map<cv::Vec<int, 3>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<int, 3>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_32F:
  //         Map<cv::Vec<float, 3>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<float, 3>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       case CV_64F:
  //         Map<cv::Vec<double, 3>, ReturnType, ParamType>(
  //             (const cv::Mat_<cv::Vec<double, 3>>&)src_,
  //             func,
  //             (cv::Mat_<ReturnType>&)dst_);
  //         break;
  //       default:
  //         ;
  //     }
  //     break;
  //   default:
  //     ;
  // }
  return;
}

////////////////////////////////////////////////////////////
// CV and Eigen conversion
// This code is from OpenCV
// The original copyright notice
/*M///////////////////////////////////////////////////////////////////////////////////////
  //
  //  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
  //
  //  By downloading, copying, installing or using the software you agree to this license.
  //  If you do not agree to this license, do not download, install,
  //  copy or use the software.
  //
  //
  //                          License Agreement
  //                For Open Source Computer Vision Library
  //
  // Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
  // Copyright (C) 2009, Willow Garage Inc., all rights reserved.
  // Third party copyrights are property of their respective owners.
  //
  // Redistribution and use in source and binary forms, with or without modification,
  // are permitted provided that the following conditions are met:
  //
  //   * Redistribution's of source code must retain the above copyright notice,
  //     this list of conditions and the following disclaimer.
  //
  //   * Redistribution's in binary form must reproduce the above copyright notice,
  //     this list of conditions and the following disclaimer in the documentation
  //     and/or other materials provided with the distribution.
  //
  //   * The name of the copyright holders may not be used to endorse or promote products
  //     derived from this software without specific prior written permission.
  //
  // This software is provided by the copyright holders and contributors "as is" and
  // any express or implied warranties, including, but not limited to, the implied
  // warranties of merchantability and fitness for a particular purpose are disclaimed.
  // In no event shall the Intel Corporation or contributors be liable for any direct,
  // indirect, incidental, special, exemplary, or consequential damages
  // (including, but not limited to, procurement of substitute goods or services;
  // loss of use, data, or profits; or business interruption) however caused
  // and on any theory of liability, whether in contract, strict liability,
  // or tort (including negligence or otherwise) arising in any way out of
  // the use of this software, even if advised of the possibility of such damage.
  //
  //M*/
template<typename _Tp, int _rows, int _cols, int _options, int _maxRows, int _maxCols>
    void eigen2cv( const Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols>& src, cv::Mat& dst )
{
  if( !(src.Flags & Eigen::RowMajorBit) )
  {
    cv::Mat _src(src.cols(), src.rows(), cv::DataType<_Tp>::type,
             (void*)src.data(), src.stride()*sizeof(_Tp));
    transpose(_src, dst);
  }
  else
  {
    cv::Mat _src(src.rows(), src.cols(), cv::DataType<_Tp>::type,
             (void*)src.data(), src.stride()*sizeof(_Tp));
    _src.copyTo(dst);
  }
}

template<typename _Tp, int _rows, int _cols, int _options, int _maxRows, int _maxCols>
    void cv2eigen( const cv::Mat& src,
                   Eigen::Matrix<_Tp, _rows, _cols, _options, _maxRows, _maxCols>& dst )
{
  CV_DbgAssert(src.rows == _rows && src.cols == _cols);
  if( !(dst.Flags & Eigen::RowMajorBit) )
  {
    cv::Mat _dst(src.cols, src.rows, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    if( src.type() == _dst.type() )
      transpose(src, _dst);
    else if( src.cols == src.rows )
    {
      src.convertTo(_dst, _dst.type());
      transpose(_dst, _dst);
    }
    else
      cv::Mat(src.t()).convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
  else
  {
    cv::Mat _dst(src.rows, src.cols, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    src.convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
}

template<typename _Tp>
void cv2eigen( const cv::Mat& src,
                   Eigen::Matrix<_Tp, Eigen::Dynamic, Eigen::Dynamic>& dst )
{
  dst.resize(src.rows, src.cols);
  if( !(dst.Flags & Eigen::RowMajorBit) )
  {
    cv::Mat _dst(src.cols, src.rows, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    if( src.type() == _dst.type() )
      transpose(src, _dst);
    else if( src.cols == src.rows )
    {
      src.convertTo(_dst, _dst.type());
      transpose(_dst, _dst);
    }
    else
      cv::Mat(src.t()).convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
  else
  {
    cv::Mat _dst(src.rows, src.cols, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    src.convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
}


template<typename _Tp>
void cv2eigen( const cv::Mat& src,
                   Eigen::Matrix<_Tp, Eigen::Dynamic, 1>& dst )
{
  CV_Assert(src.cols == 1);
  dst.resize(src.rows);

  if( !(dst.Flags & Eigen::RowMajorBit) )
  {
    cv::Mat _dst(src.cols, src.rows, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    if( src.type() == _dst.type() )
      transpose(src, _dst);
    else
      cv::Mat(src.t()).convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
  else
  {
    cv::Mat _dst(src.rows, src.cols, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    src.convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
}


template<typename _Tp>
void cv2eigen( const cv::Mat& src,
                   Eigen::Matrix<_Tp, 1, Eigen::Dynamic>& dst )
{
  CV_Assert(src.rows == 1);
  dst.resize(src.cols);
  if( !(dst.Flags & Eigen::RowMajorBit) )
  {
    cv::Mat _dst(src.cols, src.rows, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    if( src.type() == _dst.type() )
      transpose(src, _dst);
    else
      cv::Mat(src.t()).convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
  else
  {
    cv::Mat _dst(src.rows, src.cols, cv::DataType<_Tp>::type,
             dst.data(), (size_t)(dst.stride()*sizeof(_Tp)));
    src.convertTo(_dst, _dst.type());
    CV_DbgAssert(_dst.data == (uchar*)dst.data());
  }
}

////////////////////////////////////////////////////////////
// Drawing function
////////////////////////////////////////////////////////////

std::vector<cv::Scalar> RandomColors(int n);

cv::Scalar RandomColor();

template <typename PointType> cv::Mat
DrawPoints(cv::InputArray image,
           const std::vector<PointType> &points,
           const std::vector<cv::Scalar> &colors,
           int point_radius,
           bool flip_y = false) {
  assert(colors.size() == 1 || colors.size() == points.size());
  std::function<cv::Scalar (int)> get_color;
  if (colors.size() == 1)
    get_color = [&colors](int i) { return colors[0]; };
  else
    get_color = [&colors](int i) { return colors[i]; };
  cv::Mat copy;
  image.getMat().copyTo(copy);
  for (int i = 0; i < points.size(); ++i) {
    cv::Point p = To<cv::Point>(points[i]);
    if (flip_y)
      p.y = copy.size().height - p.y;
    circle(copy,
           p,
           point_radius,
           get_color(i),
           -1);
  }
  return copy;
}

template <typename PointType> cv::Mat
DrawPoints(cv::InputArray image,
           const std::vector<PointType> &points,
           int point_radius,
           cv::Scalar color = cv::Scalar(0, 0, 255),
           bool flip_y = false) {
  cv::Mat copy;
  image.getMat().copyTo(copy);
  for (int i = 0; i < points.size(); ++i) {
    cv::Point p = To<cv::Point>(points[i]);
    if (flip_y)
      p.y = copy.size().height - p.y;
    circle(copy,
           p,
           point_radius,
           color,
           -1);
  }
  return copy;
}

cv::Mat
DrawImages(cv::InputArray image0,
           cv::InputArray image1,
           cv::Size image_size = cv::Size(-1, -1));

void DrawMatches(const std::vector<cv::Point2f> &points0,
                 const std::vector<cv::Point2f> &points1,
                 const std::vector<cv::Scalar> &colors0,
                 const std::vector<cv::Scalar> &colors1,
                 int point_radius,
                 int line_width,
                 cv::Mat *image);

cv::Mat DrawImagesAndPoints(cv::InputArray image0,
                            cv::InputArray image1,
                            const std::vector<cv::Point2f> &points0,
                            const std::vector<cv::Point2f> &points1,
                            int point_radius,
                            int line_width);

template <typename PointType> cv::Mat
DrawImagesAndPoints(cv::InputArray image0,
                    cv::InputArray image1,
                    const std::vector<PointType> &points0,
                    const std::vector<PointType> &points1,
                    const std::vector<cv::Scalar> &colors0,
                    const std::vector<cv::Scalar> &colors1,
                    int point_radius,
                    int line_width = -1) {
  assert(image0.size() == image1.size() && image0.type() == image1.type());
  if (line_width > 0)
    assert(points0.size() == points1.size());
  int width = image0.size().width;
  int height = image0.size().height;
  std::function<cv::Scalar (int)> pick_color0, pick_color1;;

  if (colors0.size() == 1)
    pick_color0 = [&colors0](int i) { return colors0[0]; };
  else
    pick_color0 = [&colors0](int i) { return colors0[i]; };
  if (colors1.size() == 1)
    pick_color1 = [&colors1](int i) { return colors1[0]; };
  else
    pick_color1 = [&colors1](int i) { return colors1[i]; };

  cv::Mat out(height, width * 2, image0.type());
  image0.getMat().copyTo(out(cv::Rect(0, 0, width, height)));
  image1.getMat().copyTo(out(cv::Rect(width, 0, width, height)));
  for (size_t i = 0; i < points0.size(); ++i)
    circle(out,
           // cv::Point(points0[i].x + 0.5, points0[i].y + 0.5),
           To<cv::Point>(points0[i]),
           point_radius,
           pick_color0(i),
           -1);
  for (size_t i = 0; i < points1.size(); ++i)
    circle(out,
           // cv::Point(points1[i].x + width + 0.5, points1[i].y + 0.5),
           To<cv::Point>(points1[i]) + cv::Point(width, 0),
           point_radius,
           pick_color1(i),
           -1);

  if (line_width > 0)
    for (size_t i = 0; i < points0.size(); ++i)
      line(out,
           // cv::Point(points0[i].x + 0.5, points0[i].y + 0.5),
           // cv::Point(points1[i].x + width + 0.5, points1[i].y + 0.5),
           To<cv::Point>(points0[i]),
           To<cv::Point>(points1[i]) + cv::Point(width, 0),
           pick_color0(i),
           line_width);

  return out;
}

template <typename PointType> void
ShowImagesAndPoints(const std::string &window_name,
                    cv::InputArray image0,
                    cv::InputArray image1,
                    const std::vector<PointType> &points0,
                    const std::vector<PointType> &points1,
                    const std::vector<cv::Scalar> &colors0,
                    const std::vector<cv::Scalar> &colors1,
                    int point_radius,
                    int line_width = -1) {
  cv::namedWindow(window_name, CV_WINDOW_NORMAL);
  cv::imshow(window_name,
             DrawImagesAndPoints(image0,
                                 image1,
                                 points0,
                                 points1,
                                 colors0,
                                 colors1,
                                 point_radius,
                                 line_width));
  cv::waitKey();
}

#define WriteImageOrDie(filename, image)              \
  LOG(INFO) << "Writing " << filename;                \
  CHECK(cv::imwrite(filename, image))                 \
  << "Failed to write image \"" << filename << "\"";

#define write_image_or_dir WriteImageOrDie

cv::Mat ReadImageOrDie(const std::string &filename);

} // furry

namespace std {

template <typename Scalar>
struct hash<cv::Point_<Scalar>> {
  size_t operator () (const cv::Point_<Scalar> &p) const {
    size_t seed = 0;
    boost::hash_combine(seed, hash<Scalar>()(p.x));
    boost::hash_combine(seed, hash<Scalar>()(p.y));
    return seed;
  }
};

template <>
struct hash<cv::KeyPoint> {
  size_t operator () (const cv::KeyPoint &p) const {
    return hash<cv::Point2f>()(p.pt);
  }
};

} // std

namespace boost {
namespace serialization {

template <class Archive, class T>
void serialize(Archive &ar, cv::Point_<T> &p, const unsigned int version) {
  ar & p.x & p.y;
}

template <class Archive, class T>
void serialize(Archive &ar, cv::Point3_<T> &p, const unsigned int version) {
  ar & p.x & p.y & p.z;
}


// template <class Archive>
// void serialize(Archive &ar, cv::Point2f &p, const unsigned int version) {
//   ar & p.x & p.y;
// }

} // serialization
} // boost

#endif // FURRY_COMMON_CV_H_

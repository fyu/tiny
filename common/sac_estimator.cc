#include "sac_estimator.hpp"
#include <cmath>
#include "furry/common/geometry.hpp"

#include "opencv2/opencv.hpp"
#include "opencv2/core/eigen.hpp"

namespace furry
{

////////////////////////////////////////////////////////////
// Line2DSacModel
////////////////////////////////////////////////////////////

Line2DSacModel::Line2DSacModel(ParamType param)
{
}

void
Line2DSacModel::fit(const std::vector<DataType>& points)
{
  Eigen::Vector3d h0 = points[0].homogeneous();
  Eigen::Vector3d h1 = points[1].homogeneous();
  data_ = h0.cross(h1);
  // double norm = sqrt(data_[0] * data_[0] + data_[1] * data_[1]);
  // if (norm > std::numeric_limits<double>::epsilon())
  //   data_ /= norm;
  data_ = normalize_line(data_);
}

double
Line2DSacModel::residual(const DataType& point) const
{
  return std::fabs(data_.dot(point.homogeneous()));
}

std::ostream& operator << (std::ostream& os, const Line2DSacModel& model)
{
  os << model.data();
  return os;
}

////////////////////////////////////////////////////////////
// FundamentalMatrixSacModel
////////////////////////////////////////////////////////////

FundmentalMatrixSacModel::FundmentalMatrixSacModel(ParamType)
{
}

void
FundmentalMatrixSacModel::fit(const std::vector<DataType>& points)
{
  furry_debug_assert(points.size() >= kNumConstraints);
  std::vector<cv::Point2d> points1(kNumConstraints);
  std::vector<cv::Point2d> points2(kNumConstraints);
  for (int i = 0; i < kNumConstraints; ++i)
  {
    points1[i].x = points[i].first[0];
    points1[i].y = points[i].first[1];
    points2[i].x = points[i].second[0];
    points2[i].y = points[i].second[1];
  }
  cv::Mat f = cv::findFundamentalMat(points1, points2, CV_FM_8POINT);
  cv::cv2eigen(f, data_);
}

double
FundmentalMatrixSacModel::residual(const DataType& point) const
{
  Eigen::Vector3d epiline = data_ * point.first.homogeneous();
  epiline = normalize_line(epiline);
  double d1 = std::fabs(epiline.dot(point.second.homogeneous()));
  epiline = data_.adjoint() * point.second.homogeneous();
  epiline = normalize_line(epiline);
  double d2 = std::fabs(epiline.dot(point.first.homogeneous()));
  return std::max(d1, d2);
}

} // furry

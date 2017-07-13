#ifndef FURRY_TINY_DEPTH_MAP_H
#define FURRY_TINY_DEPTH_MAP_H

#include <string>
#include <limits>

#include <Eigen/Dense>
#include <glog/logging.h>

#include <opencv2/core/core.hpp>

#include "plane_sweeping.h"

namespace furry {
namespace tiny {

class DepthMap {
 public:
  // typedef uint8_t DepthType;
  typedef float DepthType;
  constexpr static int kCvDepthType = CV_32F; // cv::DataType<DepthType>::type;

  constexpr static DepthType kUnknownDepth = -1.0f; // numeric_limits<DepthType>::max();

  void Allocate(int rows, int cols) {
    LOG(INFO) << "Allocate depth map: " << rows << ' ' << cols;
    depth_.create(rows, cols, kCvDepthType);
    depth_ = kUnknownDepth;
    // depth_ = cv::Mat::ones(rows, cols, kCvDepthType);
  }

  int NumRows() const {
    return depth_.rows;
  }

  int NumCols() const {
    return depth_.cols;
  }

  DepthType GetUnknownDepth() const {
    return kUnknownDepth;
  }

  DepthType GetDepth(int r, int c) const {
    return depth_.at<DepthType>(r, c);
  }

  const PlaneSweepingData& GetPlaneSweepingData() const {
    return psd_;
  }

  void SetPlaneSweepingData(const PlaneSweepingData &psd) {
    psd_ = psd;
  }

  // void SetImage(cv::Mat image) {
  //   image_ = image;
  //   CHECK_EQ(image_.type(), CV_8UC3);
  // }

  void SetDepth(int r, int c, DepthType d) {
    depth_.at<DepthType>(r, c) = d;
  }

  void SetDepthOfPixelIndex(int pixel, DepthType d) {
    *(depth_.ptr<DepthType>() + pixel) = d;
  }

  int NumKnownDepths() const;

  // void SetMinDepth(double min_depth) {
  //   LOG(INFO) << "setting min: " << min_depth;
  //   min_depth_ = min_depth;
  // }

  // void SetMaxDepth(double max_depth) {
  //   LOG(INFO) << "setting max: " << max_depth;
  //   max_depth_ = max_depth;
  // }

  // void SetNumSamples(double num_samples) {
  //   num_samples_ = num_samples;
  // }

  // void SetIntrinsicMatrix(const Eigen::Matrix3d &m) {
  //   intrinsic_matrix_ = m;
  // }

  bool WriteBin(const std::string &filename) const;
  bool WriteDepthImage(const std::string &filename) const;
  bool WritePly(const std::string &filename) const;
  bool WriteInfo(const std::string &filename) const;

 private:
  cv::Mat depth_;
  // cv::Mat image_;
  // double min_depth_;
  // double max_depth_;
  // double num_samples_;
  // Eigen::Matrix3d intrinsic_matrix_;

  PlaneSweepingData psd_;
};

} // tiny
} // furry

#endif // FURRY_TINY_DEPTH_MAP_h

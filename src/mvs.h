#ifndef FURRY_TINY_MVS_H
#define FURRY_TINY_MVS_H

#include <vector>

#include <opencv2/core/core.hpp>
#include <Eigen/Dense>

#include "model.h"
#include "cost_volume.h"
#include "depth_map.h"

namespace furry {
namespace tiny {

class MvsModel {
 public:
  void SetCameras(const std::vector<Camera*>& all_cameras,
                  const std::vector<int>& camera_indexes,
                  int ref_camera_index);
  int num_cameras() const;
  const Camera* camera(int i) const;
  const Camera* ref_camera() const;
  int ref_camera_index() const;
  const cv::Mat& image(int i) const;
  const std::vector<cv::Mat>& images() const;
  std::vector<cv::Mat> MarkRefPointOnImages() const;
  std::vector<cv::Point2f> ReprojectRefPointOnCameras() const;
  std::vector<Eigen::Matrix3d> Homographies(double depth) const;
  std::vector<cv::Mat> ExtractPatches() const;
  std::vector<cv::Mat> ExtractGrayPatches() const;
  std::vector<float> CalcPatchNcc() const;
  void SetCostVolume(CostVolume &&cost_volume);
  void GenerateCostVolume(CostVolume *cost_volume, DepthMap *depth_map);
  bool UseCostVolume() const;
  std::vector<std::vector<cv::Point2f>> GetPhotoConsistency() const;

  const cv::Mat& GetRefImage() const;
  const cv::Point2f& GetRefPoint() const;

  void SetRefPoint(const cv::Point2f& p);
  void set_ref_depth(double d);
  void SetPatchSize(const cv::Size size);
  void SetPatchRadius(int radius);
  void set_patch_stride(double stride);
  void SetDepthRange(float lower, float upper);
  void SetDepthStep(float depth_step);

  double ref_depth() const;
  Eigen::Vector2f depth_range() const;
  float depth_step() const;

 private:
  std::vector<std::vector<cv::Point2f>> CalcNccOnDepth() const;
  std::vector<std::vector<cv::Point2f>> TakeCostVolumeSamples() const;

 private:
  std::vector<Camera*> cameras_;
  int ref_camera_index_;
  std::vector<cv::Mat> images_;
  std::vector<cv::Mat> gray_images_;
  cv::Mat ref_image_;

  cv::Point2f ref_point_;
  double ref_depth_;
  std::vector<Eigen::Matrix3d> homographies_;
  cv::Size patch_size_ = cv::Size(15, 15);
  int patch_stride_ = 1;
  Eigen::Vector2f depth_range_;
  float depth_step_;

  CostVolume cost_volume_;
  bool use_cost_volume_ = false;
};

} // tiny
} // furry

#endif

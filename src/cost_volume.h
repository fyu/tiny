#ifndef FURRY_APPS_TINY_COST_VOLUME_H
#define FURRY_APPS_TINY_COST_VOLUME_H

#include <string>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "grid.h"

#include "plane_sweeping.h"

namespace furry {
namespace tiny {

class CostVolume : public Grid3f {
  typedef float T;
 public:
  CostVolume() {}
  CostVolume(CostVolume &&cost_volume);

  CostVolume& operator = (CostVolume &&cost_volume);
  using Grid3f::operator = ;

  bool Allocate(int width, int height, int num_samples) {
    return Grid3<T>::Allocate(height, width, num_samples);
  }

  int GetWidth() const {
    return GetSize(1);
  }

  int GetHeight() const {
    return GetSize(0);
  }

  int NumSamples() const {
    return GetSize(2);
  }

  // T& GetSample(int x, int y, int s) {
  //   return Grid3<T>::At(y, x, s);
  // }

  // const T& GetSample(int x, int y, int s) const {
  //   return Grid3<T>::At(y, x, s);
  // }

  // const T* GetSamples(int x, int y) const {
  //   return At(y, x);
  // }

  const T& GetSampleOfPixelIndex(int pix, int sample) const {
    return (GetData() + pix * NumSamples())[sample];
  }

  const cv::Vec3b& GetColorOfPixelIndex(int pix) const {
    return *(psd_.image.ptr<cv::Vec3b>() + pix);
  }

  const T* GetSamples(const cv::Point2i &p) const {
    return At(p.y, p.x);
  }

  double GetMinDepth() const {
    return psd_.min_depth;
  }

  double GetMaxDepth() const {
    return psd_.max_depth;
  }

  double GetDepth(int s) const;

  // double GetPatchRadius() const {
  //   return patch_radius_;
  // }

  cv::Mat GetImage() const {
    return psd_.image;
  }

  // const Eigen::Matrix3d& GetIntrinsicMatrix() const {
  //   return intrinsic_matrix_;
  // }

  const PlaneSweepingData& GetPlaneSweepingData() const {
    return psd_;
  }

  void SetPlaneSweepingData(const PlaneSweepingData &psd) {
    psd_ = psd;
  }

  void SetMinDepth(double min_depth) {
    psd_.min_depth = min_depth;
  }

  void SetMaxDepth(double max_depth) {
    psd_.max_depth = max_depth;
  }

  // void SetPatchRadius(int radius) {
  //   patch_radius_ = radius;
  // }

  void SetImage(cv::Mat image) {
    psd_.image = image;
  }

  // void SetIntrinsicMatrix(const Eigen::Matrix3d &m) {
  //   intrinsic_matrix_ = m;
  // }

  void Modulate();

  void GetMinLabels(Eigen::MatrixXi *labels,
                    Eigen::MatrixXf *costs = nullptr) const;

  void RemoveEndLabels();

  bool ReadFile(const std::string &filename);
  bool WriteFile(const std::string &filename) const;

 protected:
  bool ReadFile(FILE *file);
  bool WriteFile(FILE *file) const;

  // bool ReadImage(FILE *file);
  // bool WriteImage(FILE *file) const;

 private:
  // double min_depth_ = -1;
  // double max_depth_ = -1;
  // int patch_radius_ = -1;
  // Eigen::Matrix3d intrinsic_matrix_;
  // cv::Mat image_;

  PlaneSweepingData psd_;
};

} // tiny
} // furry

#endif // FURRY_APPS_TINY_COST_VOLUME_H

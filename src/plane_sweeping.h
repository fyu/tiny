#ifndef FURRY_APPS_TINY_PLANE_SWEEPING_H
#define FURRY_APPS_TINY_PLANE_SWEEPING_H

#include <cstdio>
#include <string>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

namespace furry {
namespace tiny {

struct PlaneSweepingData {
  double min_depth;
  double max_depth;
  int patch_radius;
  int scale;
  int num_samples;
  Eigen::Matrix3d intrinsic_matrix;
  cv::Mat image;

  bool WriteBin(FILE *file) const;
  bool ReadBin(FILE *file);

  bool WriteTxt(const std::string &filename) const;
};

} // tiny
} // furry

#endif

#include "smooth.h"

#include <glog/logging.h>
#include <opencv2/imgproc/imgproc.hpp>

namespace {

void MeanShift(const furry::tiny::MeanShiftSettings &settings,
               const cv::Mat &src,
               cv::Mat *dst) {
  pyrMeanShiftFiltering(src, *dst,
                        settings.sp(),
                        settings.sr(),
                        settings.max_level());
}

}

namespace furry {
namespace tiny {

void SmoothImage(const SmoothSettings &settings,
                 const cv::Mat &src,
                 cv::Mat *dst) {
  if (settings.use_mean_shift()) {
    LOG(INFO) << "Smoothing image with mean shift";
    MeanShift(settings.mean_shift_settings(), src, dst);
  } else if (settings.use_slic()) {
    LOG(ERROR) << "Slic is not implemented";
    *dst = src;
  } else {
    LOG(INFO) << "Image is not smoothed";
    *dst = src;
  }
}

} // tiny
} // furry

#ifndef FURRY_APPS_TINY_SMOOTH_H
#define FURRY_APPS_TINY_SMOOTH_H

#include <opencv2/core/core.hpp>

#include "mrf_settings.pb.h"

namespace furry {
namespace tiny {

void SmoothImage(const SmoothSettings &settings,
                 const cv::Mat &src,
                 cv::Mat *dst);

} // tiny
} // furry

#endif // FURRY_APPS_TINY_SMOOTH_H

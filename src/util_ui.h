#ifndef FURRY_TINY_MVS_UTIL_UI_H
#define FURRY_TINY_MVS_UTIL_UI_H

#include <vector>

#include <opencv2/core/core.hpp>

namespace furry {
namespace tiny {

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

cv::Mat VListImages(const std::vector<cv::Mat> &images);
cv::Mat HListImages(const std::vector<cv::Mat> &images);
cv::Mat GridImages(const std::vector<cv::Mat>& images);

} // tiny
} // furry

#endif // FURRY_TINY_MVS_UTIL_UI_H

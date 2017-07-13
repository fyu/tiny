#include "util_ui.h"

#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

namespace furry {
namespace tiny {

const cv::Scalar kCvWhite = CV_RGB(255, 255, 255);
const cv::Scalar kCvBlack = CV_RGB(0, 0, 0);
const cv::Scalar kCvBlue = CV_RGB(0, 0, 255);
const cv::Scalar kCvRed = CV_RGB(255, 0, 0);
const cv::Scalar kCvGreen = CV_RGB(0, 255, 0);
const cv::Scalar kCvMagenta = CV_RGB(255, 0, 255);
const cv::Scalar kCvCyan = CV_RGB(0, 255, 255);
const cv::Scalar kCvYellow = CV_RGB(255, 255, 0);
const cv::Scalar kCvGray = CV_RGB(128, 128, 128);
const cv::Scalar kCvLightGray = CV_RGB(211, 211, 211);

cv::Mat VListImages(const std::vector<cv::Mat> &images) {
  int width = images[0].cols;
  int height = images[0].rows;
  int num_images = images.size();
  cv::Mat show_image(height * num_images, width, images[0].type());
  for (int i = 0; i < num_images; ++i) {
    images[i].copyTo(show_image(cv::Range(i * height, (i + 1) * height),
                                cv::Range::all()));
  }
  return show_image;
}

cv::Mat HListImages(const std::vector<cv::Mat> &images) {
  int width = images[0].cols;
  int height = images[0].rows;
  int num_images = images.size();
  cv::Mat show_image(height, width * num_images, images[0].type());
  for (int i = 0; i < num_images; ++i) {
    images[i].copyTo(show_image(cv::Range::all(),
                                cv::Range(i * width, (i + 1) * width)));
  }
  return show_image;
}

cv::Mat GridImages(const vector<cv::Mat>& images) {
  int width = images[0].cols;
  int height = images[0].rows;
  int num_images = images.size();
  if (num_images < 4) {
    if (width > height) {
      return VListImages(images);
    } else {
      return HListImages(images);
    }
  } else {
    int long_side = ceil(sqrt((double)num_images));
    int short_side = ceil((double)num_images / long_side);
    int grid_width, grid_height;
    if (width > height) {
      grid_width = short_side;
      grid_height = long_side;
    } else {
      grid_width = long_side;
      grid_height = short_side;
    }
    cv::Mat show_image(grid_height * height,
                       grid_width * width, images[0].type(),
                       kCvWhite);
    int i_image = 0;
    for (int r = 0; r < grid_height && i_image < num_images; ++r) {
      for (int c = 0; c < grid_width && i_image < num_images; ++c) {
        images[i_image].copyTo(show_image(
            cv::Range(r * height, (r + 1) * height),
            cv::Range(c * width, (c + 1) * width)));
        ++i_image;
      }
    }
    return show_image;
  }
}

} // tiny
} // furry

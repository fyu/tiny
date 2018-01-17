#ifndef FURRY_COMMON_COLOR_MAP_H_
#define FURRY_COMMON_COLOR_MAP_H_

#include <initializer_list>

#include <opencv2/core/core.hpp>

#include "furry/common/numeric.h"

namespace furry {

class ColorMap {
 public:
  struct Color {
    Color(std::initializer_list<unsigned char> l) {
      unsigned char *data = reinterpret_cast<unsigned char*>(this);
      int index = 0;
      for (auto it = l.begin(); it != l.end() && index < 3; ++it, ++index) {
        data[index] = *it;
      }
    }
    Color(unsigned char r_, unsigned char g_, unsigned char b_);
    Color();
    unsigned char r;
    unsigned char g;
    unsigned char b;
  };
  virtual Color color(double index) const = 0;
  cv::Scalar cvscalar(double index) const;
  cv::Vec3b cvvec3b(double index) const;
  cv::Mat DrawMap(const cv::Mat &index_mat, int sidebar_width = 100) const;
  template <typename T>
  void AssignMap(const cv::Mat &index_mat,
                 double max_value,
                 double min_value,
                 cv::Mat *cm) const {
    if (cm->rows < index_mat.rows ||
        cm->cols < index_mat.cols ||
        cm->type() != CV_8UC3) {
      *cm = cv::Mat(index_mat.size(), CV_8UC3);
    }
    double range = max_value - min_value;
    if (equal(range, 0.)) {
      (*cm)(cv::Rect(0, 0, index_mat.cols, index_mat.rows)) =
          cv::Mat::zeros(index_mat.size(), CV_8UC3);
      return;
    }
    for (int r = 0; r < index_mat.rows; ++r) {
      for (int c = 0; c < index_mat.cols; ++c) {
        cm->at<cv::Vec3b>(r, c) =
            cvvec3b((index_mat.at<T>(r, c) - min_value) / range);
      }
    }
  }
};

class DivergingColorMap : public ColorMap {
 public:
  // DivergingColorMap();
  Color color(double index) const;
 private:
static const Color colors_[257];
  };

}

#endif // FURRY_COMMON_COLOR_MAP_H_

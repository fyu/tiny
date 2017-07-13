#include "densecrf.h"

#include <glog/logging.h>

#include "furry/common/clock.h"
#include "furry/common/parallel.h"

namespace furry {
namespace tiny {

void GetUnary(const CostVolume &cost_volume,
              Eigen::MatrixXf *unary) {
  int width =  cost_volume.GetWidth();
  int height =  cost_volume.GetHeight();
  int labels = cost_volume.NumSamples();
  unary->resize(labels, width * height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      for (int s = 0; s < labels; ++s) {
        (*unary)(s, i * width + j) = cost_volume.At(i, j, s);
      }
    }
  }
}

TruncatedLinearCompatibility::TruncatedLinearCompatibility(
    float weight, int max_panelty)
    : weight_(weight), max_panelty_(max_panelty) {
}

void TruncatedLinearCompatibility::apply(
    MatrixXf & out_values, const MatrixXf & in_values ) const {
  out_values.resize(in_values.rows(), out_values.cols());
  int chunk_size = 200;
  int num_jobs = in_values.cols() / chunk_size;
  int num_rows = in_values.rows();
  int num_cols = in_values.cols();
  int step = (double)num_cols / num_jobs + 0.5;
  ParFor(0, num_jobs, [&](int b) {
      int i = 0;
      int j = step * b;
      int p = num_rows;
      int q = step;
      if (q + j > num_cols) q = num_cols - j;
      MatrixXf out_block;
      ApplyBlock(in_values.block(i, j, p, q), &out_block);
      out_values.block(i, j, p, q) = out_block;
    });
}

VectorXf TruncatedLinearCompatibility::parameters() const {
  VectorXf r(1);
  r[0] = weight_;
  // r[1] = max_panelty_;
  return r;
}

void TruncatedLinearCompatibility::setParameters( const VectorXf & v ) {
  weight_ = v[0];
  // max_panelty_ = v[1];
}

VectorXf TruncatedLinearCompatibility::gradient(
    const MatrixXf & b, const MatrixXf & Q ) const {
  LOG(FATAL) << "Unimplemented";
  VectorXf r(1);
  return r;
}

void TruncatedLinearCompatibility::ApplyBlock(
    const Eigen::MatrixXf & in_values, Eigen::MatrixXf *out_values) const {
  MatrixXf box;
  // LOG(INFO) << "in values: " << in_values.col(0).transpose();
  int box_radius = (max_panelty_ + 1) / 2;
  // auto t = tic();
  BoxFilter(in_values, box_radius, &box);
  BoxFilter(box, box_radius, out_values);
  // PrintEventTime("Box filter", t);
  *out_values = - weight_ * *out_values;
  // LOG(INFO) << "out values: " << out_values.col(0).transpose();
}

void TruncatedLinearCompatibility::BoxFilter(
    const Eigen::MatrixXf &src, int radius, Eigen::MatrixXf *dst) const {
  dst->resize(src.rows(), src.cols());
  dst->row(0) = src.topRows(radius + 1).colwise().sum();
  for (int r = 1; r < src.rows(); ++r) {
    dst->row(r) = dst->row(r - 1);
    if (r > radius)
      dst->row(r) -= src.row(r - radius - 1);
    if (r + radius < src.rows())
      dst->row(r) += src.row(r + radius);
  }
}

} // tiny
} // furry

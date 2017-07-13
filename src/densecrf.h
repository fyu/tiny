#ifndef FURRY_APPS_TINY_DENSECRF_H
#define FURRY_APPS_TINY_DENSECRF_H

#include <Eigen/Dense>

#include "furry/3rdparty/densecrf/densecrf.h"

#include "cost_volume.h"

namespace furry {
namespace tiny {

void GetUnary(const CostVolume &cost_volume,
              Eigen::MatrixXf *unary);

// weight_ * min(|i - j|, max_panelty_) * v
class TruncatedLinearCompatibility: public LabelCompatibility {
 protected:
  float weight_;
  int max_panelty_;
 public:
  TruncatedLinearCompatibility( float weight, int max_panelty );
  virtual void apply(Eigen::MatrixXf & out_values,
                     const Eigen::MatrixXf & in_values) const;

  // Training and parameters
  virtual Eigen::VectorXf parameters() const;
  virtual void setParameters( const Eigen::VectorXf & v );
  virtual Eigen::VectorXf gradient( const Eigen::MatrixXf & b,
                                    const Eigen::MatrixXf & Q ) const;

 private:
  void BoxFilter(
      const Eigen::MatrixXf &src, int radius, Eigen::MatrixXf *dst) const;
  void ApplyBlock(const Eigen::MatrixXf & in_values,
                  Eigen::MatrixXf *out_values) const;
};

} // tiny
} // furry

#endif

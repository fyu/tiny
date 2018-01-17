#ifndef FURRY_COMMON_CERES_H_
#define FURRY_COMMON_CERES_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <ceres/ceres.h>
#include <ceres/jet.h>
#include <ceres/rotation.h>

#include "furry/common/str.h"

// Ceres-solver autodiff style template functions

namespace furry {

template <int Dim, typename T> void
LinearInterpolate(const T* const p0,
                  const T* const p1,
                  const double t,
                  T* result) {
  for (int i = 0; i < Dim; ++i) {
    result[i] = T(1 - t) * p0[i] + T(t) * p1[i];
  }
}

// Assume matrix is in row major
template <int Rows, int Cols, typename MatrixType, typename VectorType> void
MatrixXVector(const MatrixType* const m,
              const VectorType* const v,
              VectorType* result) {
  for (int r = 0; r < Rows; ++r) {
    result[r] = VectorType(0);
    for (int c = 0; c < Cols; ++c) {
      result[r] += VectorType(m[r * Cols + c]) * v[c];
    }
  }
}

template <int Dim, typename T> void
VectorMinus(const T* const v0, const T* const v1, T* result) {
  for (int i = 0; i < Dim; ++i) {
    result[i] = v0[i] - v1[i];
  }
}

template <int Dim, typename T> void
VectorAddVector(const T* const v0, const T* const v1, T* result) {
  for (int i = 0; i < Dim; ++i) {
    result[i] = v0[i] + v1[i];
  }
}

template <int Dim, typename T> void
ScaleVector(const T* const v, double scale, T* result) {
  for (int i = 0; i < Dim; ++i) {
    result[i] = v[i] * scale;
  }
}

template <int Dim, typename T0, typename T1> void
Assign(const T0 * const v0, T1 *v1) {
  for (int i = 0; i < Dim; ++i) {
    v1[i] = T1(v0[i]);
  }
}

////////////////////////////////////////////////////////////
// Quaternion
////////////////////////////////////////////////////////////


template <int Dim, typename T> void
VectorNorm2(const T* const v, T* result) {
  *result = T(0);
  for (int i = 0; i < Dim; ++i) {
    *result += v[i] * v[i];
  }
}

template <int Dim, typename T> void
VectorNorm(const T* const v, T* result) {
  VectorNorm2<Dim, T>(v, result);
  *result = ceres::sqrt(*result);
}

// result can be v0 or v1
template <int Dim, typename T> void
VectorDotProduct(const T* const v0, const T* const v1, T* result) {
  for (int i = 0; i < Dim; ++i) {
    result[i] = v0[i] * v1[i];
  }
}

template <typename T> void
QuaternionNormalize(const T* const q, T* result) {
  T scale = T(1) / sqrt(q[0] * q[0] +
                        q[1] * q[1] +
                        q[2] * q[2] +
                        q[3] * q[3]);
  result[0] = q[0] * scale;
  result[1] = q[1] * scale;
  result[2] = q[2] * scale;
  result[3] = q[3] * scale;
}

template <typename T> void
QuaternionExponential(const T* const q, T* result) {
  T vector_norm;
  VectorNorm<3, T>(q + 1, &vector_norm);
  T vector_ratio = ceres::sin(vector_norm) / vector_norm;
  T exp_w = ceres::exp(q[0]);
  result[0] = exp_w * ceres::cos(vector_norm);
  result[1] = q[1] * vector_ratio;
  result[2] = q[2] * vector_ratio;
  result[3] = q[3] * vector_ratio;
}

template <typename T> bool
CeresIsZero(const T& x) {
  static const T pos_epsilon = T(1e-12);
  static const T neg_epsilon = - pos_epsilon;
  if (x > neg_epsilon && x < pos_epsilon) {
    return true;
  } else {
    return false;
  }
}

// template <typename T> void
// QuaternionProduct(const T* const q, const T* const p, T* result) {
//   // result[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
//   // result[1] = p[1] * q[0] + p[0] * q[1] + p[3] * q[2] - p[2] * q[3];
//   // result[2] = p[2] * q[0] - p[3] * q[1] + p[0] * q[2] + p[1] * q[3];
//   // result[3] = p[3] * q[0] + p[2] * q[1] - p[1] * q[2] + p[0] * q[3];
//   result[0] = q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3];
//   result[1] = q[1] * p[0] + q[0] * p[1] - q[3] * p[2] + q[2] * p[3];
//   result[2] = q[2] * p[0] + q[3] * p[1] + q[0] * p[2] - q[1] * p[3];
//   result[3] = q[3] * p[0] - q[2] * p[1] + q[1] * p[2] + q[0] * p[3];
// }

// template <typename T> void
// QuaternionConjugateProduct(const T* const q, const T* const p, T* result) {
//   result[0] = q[0] * p[0] + q[1] * p[1] + q[2] * p[2] + q[3] * p[3];
//   result[1] = q[1] * p[0] - q[0] * p[1] + q[3] * p[2] - q[2] * p[3];
//   result[2] = q[2] * p[0] - q[3] * p[1] - q[0] * p[2] + q[1] * p[3];
//   result[3] = q[3] * p[0] + q[2] * p[1] - q[1] * p[2] - q[0] * p[3];
// }

// z w*
template<typename T>
void QuaternionConjugateProduct(const T z[4], const T w[4], T zwc[4]) {
  zwc[0] =   z[0] * w[0] + z[1] * w[1] + z[2] * w[2] + z[3] * w[3];
  zwc[1] = - z[0] * w[1] + z[1] * w[0] - z[2] * w[3] + z[3] * w[2];
  zwc[2] = - z[0] * w[2] + z[1] * w[3] + z[2] * w[0] - z[3] * w[1];
  zwc[3] = - z[0] * w[3] - z[1] * w[2] + z[2] * w[1] + z[3] * w[0];
}


template <typename T> void
UnitQuaternionConjugate(const T* const q, T* result) {
  result[0] = q[0];
  result[1] = -q[1];
  result[2] = -q[2];
  result[3] = -q[3];
}

template <typename T>
void QuaternionConjugate(const T* const q, T* result) {
  result[0] = q[0];
  result[1] = -q[1];
  result[2] = -q[2];
  result[3] = -q[3];
}

template <typename T>
void UnitQuaternionAngle(const T* const quaternion, T *angle) {
  const T& q1 = quaternion[1];
  const T& q2 = quaternion[2];
  const T& q3 = quaternion[3];
  const T sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;
  if (sin_squared_theta > T(0.0)) {
    const T sin_theta = sqrt(sin_squared_theta);
    const T& cos_theta = quaternion[0];

    // If cos_theta is negative, theta is greater than pi/2, which
    // means that angle for the angle_axis vector which is 2 * theta
    // would be greater than pi.
    //
    // While this will result in the correct rotation, it does not
    // result in a normalized angle-axis vector.
    //
    // In that case we observe that 2 * theta ~ 2 * theta - 2 * pi,
    // which is equivalent saying
    //
    //   theta - pi = atan(sin(theta - pi), cos(theta - pi))
    //              = atan(-sin(theta), -cos(theta))
    //
    *angle =
        T(2.0) * ((cos_theta < 0.0)
                  ? atan2(-sin_theta, -cos_theta)
                  : atan2(sin_theta, cos_theta));
  } else {
    // For zero rotation, sqrt() will produce NaN in the derivative since
    // the argument is zero.  By approximating with a Taylor series,
    // and truncating at one term, the value and first derivatives will be
    // computed correctly when Jets are used.
    *angle = T(0.0);
  }
}

template <typename T> void
UnitQuaternionPower(const T* const q, double p, T* result) {
  T angle;
  UnitQuaternionAngle(q, &angle);
  angle /= T(2.0);

  T t;
  if (CeresIsZero(angle)) {
    t = T(p);
  } else {
    t = sin(T(p) * angle) / sin(angle);
  }
  result[0] = ceres::cos(T(p) * angle);
  result[1] = q[1] * t;
  result[2] = q[2] * t;
  result[3] = q[3] * t;
}

template <typename T> void
QuaternionPower(const T* const q, double p, T* result) {
  T unit_q[4];
  QuaternionNormalize(q, unit_q);
  UnitQuaternionPower(unit_q, p, result);
}

#if 0

// naive slerp implementation

template <typename T> void
UnitQuaternionSlerp(const T* const q0, const T* const q1, double t, T* result) {
  T move[4];
  QuaternionConjugateProduct<T>(q1, q0, move);
  QuaternionNormalize<T>(move, move);
  UnitQuaternionPower<T>(move, t, move);
  ceres::QuaternionProduct<T>(move, q0, result);
  QuaternionNormalize<T>(result, result);
}

#endif

#if 1

template <typename T> void
UnitQuaternionSlerp(const T* const q0, const T* const q1, double t, T* result) {
  static const T one = T(1.0 - 1e-12);
  T d = q1[0]*q0[0] + q1[1]*q0[1] + q1[2]*q0[2] + q1[3]*q0[3];
  T absD = ::ceres::abs(d);

  T scale0;
  T scale1;

  if(absD>=one)
  {
    scale0 = T(1 - t);
    scale1 = T(t);
  }
  else
  {
    // theta is the angle between the 2 quaternions
    T theta = acos(absD);
    // std::cout << "theta: " << theta << " d: " << absD << " one: " << one  << '\n';
    T sinTheta = sin(theta);

    scale0 = sin( ( T(1 - t) ) * theta) / sinTheta;
    scale1 = sin( ( T(t) * theta) ) / sinTheta;
  }
  if(d < T(0)) scale1 = -scale1;

  result[0] = scale0 * q0[0] + scale1 * q1[0];
  result[1] = scale0 * q0[1] + scale1 * q1[1];
  result[2] = scale0 * q0[2] + scale1 * q1[2];
  result[3] = scale0 * q0[3] + scale1 * q1[3];

  QuaternionNormalize(result, result);
}

#endif

template <typename T> void
QuaternionSlerp(const T* const q0, const T* const q1, double t, T* result) {
  T unit_q0[4];
  T unit_q1[4];
  QuaternionNormalize(q0, unit_q0);
  QuaternionNormalize(q1, unit_q1);
  UnitQuaternionSlerp(unit_q0, unit_q1, t, result);
}

// template <typename T> void
// QuaternionSlerp(const T* const q0, const T* const q1, double t,  T* result) {
//   T unit_q0[4];
//   T unit_q1[4];
//   QuaternionNormalize(q0, unit_q0);
//   QuaternionNormalize(q1, unit_q1);
//   UnitQuaternionSlerp(unit_q0, unit_q1, t, result);
// }

class ResidualBlock {
 public:
  ResidualBlock(ceres::CostFunction* cost_func,
                ceres::LossFunction* loss_func,
                const std::vector<double*> &params)
      : cost_func_(cost_func),
        loss_func_(loss_func),
        params_(params) {
  }

  void AddTo(ceres::Problem* problem) const {
    problem->AddResidualBlock(cost_func_, loss_func_, params_);
  }

  void Residuals(std::vector<double>* residuals) const {
    residuals->resize(cost_func_->num_residuals());
    cost_func_->Evaluate(params_.data(), residuals->data(), NULL);
  }

  double Cost() const {
    std::vector<double> residuals;
    this->Residuals(&residuals);
    double cost = 0;
    for (auto r : residuals) {
      cost += r * r;
    }
    return cost;
  }

  double ActualCost() const {
    double s[3];
    loss_func_->Evaluate(this->Cost(), s);
    return s[0];
  }

 private:
  ceres::CostFunction* cost_func_;
  ceres::LossFunction* loss_func_;
  std::vector<double*> params_;
};

class ResidualBlocks {
 public:
  ResidualBlocks(const std::string& name) : name_(name) {}

  void Add(ResidualBlock* block) {
    residual_blocks_.push_back(block);
  }

  void Add(ceres::CostFunction* cost_func,
           ceres::LossFunction* loss_func,
           const std::vector<double*> &params) {
    residual_blocks_.push_back(new ResidualBlock(cost_func, loss_func, params));
  }

  void AddTo(ceres::Problem* problem) const {
    for (auto block : residual_blocks_) {
      block->AddTo(problem);
    }
  }

  void Release() {
    for (auto block : residual_blocks_) {
      delete block;
    }
    residual_blocks_.resize(0);
  }

  std::string Info() const {
    double total = 0;
    for (auto block : residual_blocks_) {
      total += block->Cost();
    }
    return StringPrintf("%s: %lf (%lu * %lf ^ 2 * 0.5)",
                        name_.c_str(), total * 0.5,
                        residual_blocks_.size(),
                        sqrt(total / residual_blocks_.size()));
  }

 private:
  std::string name_;
  std::vector<ResidualBlock*> residual_blocks_;
};

} // furry

#endif

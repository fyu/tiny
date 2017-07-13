#ifndef TINY_BUNDLE_UTIL_H
#define TINY_BUNDLE_UTIL_H

#include <string>
#include <vector>

#include <ceres/ceres.h>

#include "furry/common/ceres.h"
#include "furry/common/str.h"

#include "model.h"
#include "ba_settings.pb.h"

namespace furry {
namespace tiny {

void ChangeCoord(const Camera::Externals &externals, Model *model);

template <int dim, typename T1, typename T2>
void Copy(const T1* const src, T2 *dst) {
  for (int i = 0; i < dim; ++i) {
    dst[i] = T2(src[i]);
  }
}

inline double* CameraCenterPtr(double* camera_param) {
  return camera_param + 3;
}

inline double* CameraRotationPtr(double* camera_param) {
  return camera_param;
}

class CameraPrior {
 public:
  CameraPrior(const BundleSettings& settings)
      : f_mean_(settings.f_mean()),
        f_weight_(1 / settings.f_std()),
        k_mean_(settings.k_mean()),
        k_weight_(1 / settings.k_std()) {
  }

  template <typename T>
  bool operator () (const T* const internals,
                    T* residuals) const {
    const T& f = internals[0];
    const T& k = internals[1];
    residuals[0] = f_weight_ * (f - f_mean_);
    residuals[1] = k_weight_ * (k - k_mean_);
    return true;
  }

 private:
  const double f_mean_;
  const double f_weight_;
  const double k_mean_;
  const double k_weight_;
};

template <int Dim>
struct VectorPrior {
  VectorPrior(int start_, double mean_, double std_) {
    start = start_;
    for (int i = 0; i < Dim; ++i) {
      mean[i] = mean_;
      weight[i] = 1 / std_;
    }
  }

  template <typename T>
  bool operator () (const T* const vec, T* residuals) const {
    for (int i = 0; i < Dim; ++i) {
      residuals[i] = T(weight[i]) * (vec[start + i] - mean[i]);
    }
    return true;
  }

  int start;
  double mean[Dim];
  double weight[Dim];
};

template <int Dim>
struct VectorDiffPrior {
  VectorDiffPrior(int start_, double mean_, double std_) {
    start = start_;
    for (int i = 0; i < Dim; ++i) {
      mean[i] = mean_;
      weight[i] = 1 / std_;
    }
  }

  template <typename T>
  bool operator () (const T* const vec0,
                    const T* const vec1,
                    T* residuals) const {
    for (int i = 0; i < Dim; ++i) {
      residuals[i] = T(weight[i]) * (abs(vec0[start + i] - vec1[start + i])
                                     - mean[i]);
    }
    return true;
  }

  int start;
  double mean[Dim];
  double weight[Dim];
};

struct SmoothCurvePrior {
 public:
  SmoothCurvePrior(int padding, double weight) :
      padding_(padding),weight_(weight) {}

  template <typename T>
  bool operator () (const T* const p0_,
                    const T* const p1_,
                    const T* const p2_,
                    T* residual) const {
    const T* p0 = p0_ + padding_;
    const T* p1 = p1_ + padding_;
    const T* p2 = p2_ + padding_;
    T v0[3], v1[3];
    v0[0] = p1[0] - p0[0];
    v0[1] = p1[1] - p0[1];
    v0[2] = p1[2] - p0[2];
    v1[0] = p2[0] - p1[0];
    v1[1] = p2[1] - p1[1];
    v1[2] = p2[2] - p1[2];
    T n0 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
    T n1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    if (n0 > T(1e-10) && n1 > T(1e-10)) {
      n0 = sqrt(n0);
      n1 = sqrt(n1);
      // T l = (n0 + n1) / T(2);
      // residual[0] = T(weight_) * (v0[0] / n0 - v1[0] / n1) / l;
      // residual[1] = T(weight_) * (v0[1] / n0 - v1[1] / n1) / l;
      // residual[2] = T(weight_) * (v0[2] / n0 - v1[2] / n1) / l;
      residual[0] = T(weight_) * (v0[0] / n0 - v1[0] / n1);
      residual[1] = T(weight_) * (v0[1] / n0 - v1[1] / n1);
      residual[2] = T(weight_) * (v0[2] / n0 - v1[2] / n1);
    } else {
      residual[0] = T(0);
      residual[1] = T(0);
      residual[2] = T(0);
    }
    return true;
  }

 private:
  int padding_;
  double weight_;
};

inline void InitializeCameraBlock(const Camera& camera,
                           double* camera_params) {
  auto c = camera.GetCenter();
  auto a = camera.GetEulerAngles();
  camera_params[0] = a[0];
  camera_params[1] = a[1];
  camera_params[2] = a[2];
  camera_params[3] = c[0];
  camera_params[4] = c[1];
  camera_params[5] = c[2];
}

} // tiny
} // furry

#endif

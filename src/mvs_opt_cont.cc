#include "mvs_opt.h"

#include <functional>
#include <vector>

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <glog/logging.h>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/common/ceres.h"
#include "furry/common/parallel.h"

using namespace Eigen;
using namespace std;

namespace {

class Rbf1D {
 public:

  class Multiquadric {
   public:
    Multiquadric(double x0, double r) : x0_(x0), r2_(r * r) {}

    bool operator () (const double *x, double *v, double **j = nullptr) const {
      double diff = (*x - x0_);
      double dist2 = diff * diff;
      *v = sqrt(r2_ + dist2);
      if (j != nullptr) {
        **j = diff / *v;
      }
      return true;
    }

   private:
    double x0_;
    double r2_;
  };

  class Gaussian {
   public:
    Gaussian(double x0, double spacing)
        : x0_(x0), r2_(1.0f / (spacing * spacing)) {}

    bool operator () (const double *x, double *v, double **j = nullptr) const {
      double diff = (*x - x0_);
      // LOG(INFO) << "x: " << *x << " x0_: " << x0_ << " Diff: " << diff;
      double dist2 = diff * diff;
      *v = exp(-dist2 * r2_);
      if (j != nullptr && j[0] != nullptr) {
        **j = -2 * *v * diff * r2_;
        // CHECK(!std::isnan(**j));
      }
      return true;
    }

   private:
    double x0_;
    double r2_;
  };

  typedef std::function<bool (const double*, double*, double**)> KernelType;

  Rbf1D(double spacing) : spacing_(spacing) {}

  bool operator () (const double *x, double *v, double **jac = nullptr) const;

  bool Interpolate(const std::vector<Eigen::Vector2d> &points);

 private:
  std::vector<KernelType> kernels_;
  std::vector<double> coeffs_;
  double max_;
  double spacing_;
};

class Rbf1DCostFunc : public ceres::SizedCostFunction<1, 1> {
 public:
  Rbf1DCostFunc(Rbf1D *rbf) : rbf_(rbf) {}
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const;
 private:
  Rbf1D *rbf_;
};

class Robustifier {
 public:
  virtual bool operator () (double x, double *v, double *jac) const = 0;
};

class CauchyRobustifier : public Robustifier {
 public:
  bool operator () (double x, double *v, double *jac) const;
};

class TrivialRobustifier : public Robustifier {
 public:
  bool operator () (double x, double *v, double *jac) const;
};

class FirstOrderSmoothnessCostFunc : public ceres::SizedCostFunction<1, 1, 1> {
 public:
  FirstOrderSmoothnessCostFunc(double weight, double scale, Robustifier *robust)
      : weight_(weight), scale_(1 / scale), robust_(robust) {}

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const;
 private:
  double scale_;
  double weight_;
  Robustifier *robust_;
};

class SecondOrderSmoothnessCostFunc
    : public ceres::SizedCostFunction<1, 1, 1, 1> {
 public:
  SecondOrderSmoothnessCostFunc(
      double weight, double scale, Robustifier *robust)
    : weight_(weight), scale_(1 / scale), robust_(robust) {}

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const;

 private:
  double scale_;
  double weight_;
  Robustifier *robust_;
};

class ColorConnection {
 public:

  void SetImage(const cv::Mat &image);
  void SetParameters(double weight, double sigma);

  double GetWeight(int i0, int i1) const;

 private:

  template <typename T>
  float Dist(cv::Vec<T, 3> c0, cv::Vec<T, 3> c1) const {
    return sqrt((c0[0] - c1[0]) * (c0[0] - c1[0]) +
                (c0[1] - c1[1]) * (c0[1] - c1[1]) +
                (c0[2] - c1[2]) * (c0[2] - c1[2]));
  }

  cv::Mat lab_image_;
  double weight_;
  double sigma_;
};

bool Rbf1D::operator () (
    const double *x, double *v, double **jac) const {
  *v = 0;
  double kernel_value;
  if (jac == nullptr) {
    for (size_t i = 0; i < kernels_.size(); ++i) {
      kernels_[i](x, &kernel_value, nullptr);
      *v += kernel_value * coeffs_[i];
    }
  } else {
    double jac_value = 0;
    for (size_t i = 0; i < kernels_.size(); ++i) {
      kernels_[i](x, &kernel_value, jac);
      *v += kernel_value * coeffs_[i];
      jac_value += **jac * coeffs_[i];
      // LOG(INFO) << "v: " << *v << " jac: " << jac_value;
    }
    **jac = - jac_value;
  }
  *v = max_ - *v;
  return true;
}


bool Rbf1D::Interpolate(const vector<Vector2d> &points) {
  int num_points = points.size();
  kernels_.resize(0);
  max_ = 0;
  for (int i = 0; i < num_points; ++i) {
    Gaussian kernel(points[i].x(), spacing_);
    kernels_.push_back(kernel);
    if (points[i].y() > max_) {
      max_ = points[i].y();
    }
  }
  MatrixXd A(num_points, num_points);
  VectorXd b(num_points);
  for (int r = 0; r < num_points; ++r) {
    if (points[r].y() > 0) {
      b(r) = max_ - points[r].y();
    } else  {
      b(r) = 0;
    }
    for (int c = 0; c < num_points; ++c) {
      kernels_[c](&points[r].x(), &A(r, c), nullptr);
    }
  }

  VectorXd x = A.ldlt().solve(b);
  // VectorXd x = A.fullPivLu().solve(b);
  // VectorXd c = A * x - b;
  // LOG(INFO) << "Interpolation error: " << c.norm();
  coeffs_.resize(num_points);
  memcpy(coeffs_.data(), x.data(), num_points * sizeof(double));

  return true;
}

bool Rbf1DCostFunc::Evaluate(double const* const* parameters,
                             double* residuals,
                             double** jacobians) const {
  if (parameters[0][0] < - 0 || parameters[0][0] > 63) return false;
  (*rbf_)(parameters[0], residuals, jacobians);
  residuals[0] = sqrt(residuals[0]);
  if (jacobians != nullptr && jacobians[0] != nullptr) {
    // CHECK(residuals[0] != 0);
    if (residuals[0] == 0) {
      jacobians[0][0] = 0;
    } else {
      jacobians[0][0] = 0.5 * jacobians[0][0] / residuals[0];
    }
    // LOG(INFO) << "R: " << residuals[0] << " J: " << jacobians[0][0];
  }
  return true;
}

bool CauchyRobustifier::operator () (double x, double *v, double *jac) const {
  *v = log(1 + x);
  *jac = 1 / (1 + x);
  return true;
}

bool TrivialRobustifier::operator () (double x, double *v, double *jac) const {
  *v = x;
  *jac = 1;
  return true;
}

bool FirstOrderSmoothnessCostFunc::Evaluate(double const* const* parameters,
                                            double* residuals,
                                            double** jacobians) const {
  double diff = parameters[0][0] - parameters[1][0];
  double dist = sqrt(diff * diff);
  double r_v, r_j;
  (*robust_)(dist * scale_, &r_v, &r_j);
  residuals[0] = sqrt(r_v * weight_);
  if (jacobians != nullptr) {
    double jac;
    if (residuals[0] == 0) {
      jac = 0;
    } else {
      jac = scale_ * weight_ * diff * r_j / (2 * dist * residuals[0]);
    }
    if (jacobians[0] != nullptr) jacobians[0][0] = jac;
    if (jacobians[1] != nullptr) jacobians[1][0] = -jac;
  }
  return true;
}

bool SecondOrderSmoothnessCostFunc::Evaluate(double const* const* parameters,
                                            double* residuals,
                                            double** jacobians) const {
  double diff = parameters[0][0] - 2 * parameters[1][0] + parameters[2][0];
  double dist = sqrt(diff * diff);
  double r_v, r_j;
  (*robust_)(dist * scale_, &r_v, &r_j);
  residuals[0] = sqrt(r_v * weight_);
  if (jacobians != nullptr) {
    double jac;
    if (residuals[0] == 0) {
      jac = 0;
    } else {
      jac = scale_ * weight_ * diff * r_j / (2 * dist * residuals[0]);
    }
    if (jacobians[0] != nullptr) jacobians[0][0] = jac;
    if (jacobians[1] != nullptr) jacobians[1][0] = -2 * jac;
    if (jacobians[2] != nullptr) jacobians[2][0] = jac;
  }
  return true;
}

void ColorConnection::SetImage(const cv::Mat &image) {
  cv::cvtColor(image, lab_image_, CV_BGR2Lab);
}

void ColorConnection::SetParameters(double weight, double sigma) {
  weight_ = weight;
  sigma_ = sigma;
}

double ColorConnection::GetWeight(int i0, int i1) const {
  return weight_ * exp(- Dist(*(lab_image_.ptr<cv::Vec3b>() + i0),
                              *(lab_image_.ptr<cv::Vec3b>() + i1)) /
                       sigma_);
}


} // anonymous

namespace furry {
namespace tiny {

void SolveContinuousDepth(
    const ContinuousMrfSettings &settings,
    const CostVolume &cost_volume,
    DepthMap *depthmap) {
  int rows = cost_volume.NumRows();
  int cols = cost_volume.NumCols();
  CHECK_EQ(depthmap->NumRows(), rows);
  CHECK_EQ(depthmap->NumCols(), cols);
  vector<vector<Rbf1D>> continuous_cost_volume(rows);

  LOG(INFO) << "Initializing RBF";
  for (int r = 0; r < rows; ++r) {
    continuous_cost_volume[r].reserve(cols);
    for (int c = 0; c < cols; ++c) {
      continuous_cost_volume[r].push_back(Rbf1D(settings.rbf_spacing()));
    }
  }
  // for (int r = 0; r < rows; ++r) {
  ParFor(0, rows, [&](int r) {
      vector<Vector2d> points(cost_volume.NumSamples());
      for (int c = 0; c < cols; ++c) {
        auto samples = cost_volume.At(r, c);
        for (int s = 0; s < cost_volume.NumSamples(); ++s) {
          points[s].x() = s;
          points[s].y() = samples[s];
        }
        continuous_cost_volume[r][c].Interpolate(points);
      }
    });

  ColorConnection color_connection;
  color_connection.SetImage(cost_volume.GetImage());

  double *depth = new double[rows * cols];

  // Initialize depth
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      depth[r * cols + c] = depthmap->GetDepth(r, c);
    }
  }

  ceres::Problem *problem = new ceres::Problem;

  vector<double*> parameters;

  vector<ResidualBlocks> all_residual_blocks;

  LOG(INFO) << "Initializing data term";
  ResidualBlocks data_blocks("Data");
  parameters.resize(1);
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      parameters[0] = depth + r * cols + c;
      data_blocks.Add(new Rbf1DCostFunc(&continuous_cost_volume[r][c]),
                      nullptr,
                      parameters);
    }
  }
  all_residual_blocks.push_back(data_blocks);
  // data_blocks.AddTo(problem);

  // CauchyRobustifier robustifier;
  TrivialRobustifier robustifier;

  if (settings.use_fos()) {
    color_connection.SetParameters(settings.fos().color_max_weight(),
                                   settings.fos().color_sigma());
    double max_panelty = 0.15 * cost_volume.NumSamples();
    ResidualBlocks fos_blocks("First Order Smoothness");
    parameters.resize(2);
    for (int r = 0; r < rows; ++r) {
      for (int c = 1; c < cols; ++c) {
        int i0 = r * cols + c;
        int i1 = r * cols + c - 1;
        parameters[0] = depth + i0;
        parameters[1] = depth + i1;
        fos_blocks.Add(
            new FirstOrderSmoothnessCostFunc(
                settings.fos().weight() * color_connection.GetWeight(i0, i1),
                max_panelty, &robustifier),
            nullptr, parameters);
      }
    }
    for (int r = 1; r < rows; ++r) {
      for (int c = 0; c < cols; ++c) {
        int i0 = r * cols + c;
        int i1 = (r - 1) * cols + c;
        parameters[0] = depth + i0;
        parameters[1] = depth + i1;
        fos_blocks.Add(
            new FirstOrderSmoothnessCostFunc(
                settings.fos().weight() * color_connection.GetWeight(i0, i1),
                max_panelty, &robustifier),
            nullptr, parameters);
      }
    }
    all_residual_blocks.push_back(fos_blocks);
    // fos_blocks.AddTo(problem);
  }

  if (settings.use_sos()) {
    LOG(INFO) << "Initializing SOS";
    color_connection.SetParameters(settings.sos().color_max_weight(),
                                   settings.sos().color_sigma());
    double max_panelty = 0.15 * cost_volume.NumSamples();
    ResidualBlocks sos_blocks("Second Order Smoothness");
    parameters.resize(3);
    for (int r = 0; r < rows; ++r) {
      for (int c = 1; c < cols - 1; ++c) {
        int i0 = r * cols + c - 1;
        int i1 = r * cols + c;
        int i2 = r * cols + c + 1;
        parameters[0] = depth + i0;
        parameters[1] = depth + i1;
        parameters[2] = depth + i2;
        sos_blocks.Add(
            new SecondOrderSmoothnessCostFunc(
                settings.sos().weight() * color_connection.GetWeight(i0, i2),
                max_panelty, &robustifier),
            nullptr, parameters);
      }
    }
    for (int r = 1; r < rows - 1; ++r) {
      for (int c = 0; c < cols; ++c) {
        int i0 = (r - 1) * cols + c;
        int i1 = r * cols + c;
        int i2 = (r + 1) * cols + c;
        parameters[0] = depth + i0;
        parameters[1] = depth + i1;
        parameters[2] = depth + i2;
        sos_blocks.Add(
            new SecondOrderSmoothnessCostFunc(
                settings.sos().weight() * color_connection.GetWeight(i0, i2),
                max_panelty, &robustifier),
            nullptr, parameters);
      }
    }
    all_residual_blocks.push_back(sos_blocks);
    // sos_blocks.AddTo(problem);
  }

  for (auto &blocks : all_residual_blocks) {
    blocks.AddTo(problem);
  }

  auto solver_options = settings.solver_options();
  ceres::Solver::Options options;
  options.max_num_iterations = solver_options.max_num_iterations();
  options.num_threads = solver_options.num_threads();
  options.num_linear_solver_threads = solver_options.num_threads();
  // options.linear_solver_type = ceres::SPARSE_SCHUR;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout =
      solver_options.minimizer_progress_to_stdout();
  options.use_inner_iterations = solver_options.use_inner_iterations();

  ceres::Solver::Summary summary;
  ceres::Solve(options, problem, &summary);

  if (solver_options.print_summary()) {
    cout << summary.FullReport() << '\n';
  }

  for (auto& blocks : all_residual_blocks) {
    LOG(INFO) << blocks.Info();
  }

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      int index = r * cols + c;
      depthmap->SetDepthOfPixelIndex(index, depth[index]);
    }
  }

  delete problem;
  delete depth;
}

} // tiny
} // furry

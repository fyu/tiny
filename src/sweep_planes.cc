const char *kUsageMessage = "Use plane sweeping algorithm to produce cost volume";

#include <algorithm>
#include <cmath>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/common/clock.h"
#include "furry/common/cv.h"
#include "furry/common/init.h"
#include "furry/common/memory.h"
#include "furry/common/parallel.h"
#include "furry/common/str.h"

#include "model.h"
#include "cost_volume.h"
#include "plane_sweeping.h"

using namespace Eigen;
using namespace furry;
using namespace furry::tiny;
using namespace std;
using namespace cv;

DEFINE_string(in_model, "", "Input Model");
DEFINE_int32(ref_index, 0, "Reference image index");
DEFINE_int32(patch_radius, 3, "");
DEFINE_int32(scale, 0,
             "The effective image resolution is 2 ^ (-scale) of original one");
DEFINE_int32(num_samples, 0, "");
DEFINE_double(min_depth, -1, "");
DEFINE_double(max_depth, -1, "");
DEFINE_string(out_cv, "", "");
DEFINE_int32(p, -1, "");

void GetDepthRange(const furry::tiny::Model &model, int ref_image,
                   double *min_depth, double *max_depth);
void InversePerspectiveSampling(double min, double max, int num_samples,
                                vector<double> *depth_samples);
double InversePerspectiveSampling(const double min, const double max,
                                  const double value);

double CalcNcc(const Mat &p0, const Mat &p1);
double CalcSos(const Mat &p0, const Mat &p1);
double MergeScores(vector<double> &scores);

class MatchFunc {
 public:
  virtual ~MatchFunc() {}
  virtual double operator () (const Mat &p0, const Mat &p1) const = 0;
};

class Sad : public MatchFunc {
 public:
  virtual double operator () (const Mat &p0, const Mat &p1) const;
};

class Ssd : public MatchFunc {
 public:
  virtual double operator () (const Mat &p0, const Mat &p1) const;
};

class Ncc : public MatchFunc {
 public:
  virtual double operator () (const Mat &p0, const Mat &p1) const;
};

int main(int argc, char *argv[]) {
  furry::Init(&argc, &argv, kUsageMessage);
  tbb::task_scheduler_init tbb_init(FLAGS_p);
  CHECK(!FLAGS_in_model.empty());
  furry::tiny::Model model;
  model.ReadFile(FLAGS_in_model);
  double min_depth = FLAGS_min_depth;
  double max_depth = FLAGS_max_depth;
  if (min_depth < 0 || max_depth < 0) {
    GetDepthRange(model, FLAGS_ref_index,
                  (min_depth < 0) ? &min_depth : nullptr,
                  (max_depth < 0) ? &max_depth : nullptr);
  }
  LOG(INFO) << "Depth range: " << min_depth << ' ' << max_depth;

  vector<double> depth_samples;
  InversePerspectiveSampling(min_depth, max_depth, FLAGS_num_samples,
                             &depth_samples);

  vector<Mat> scaled_gray_images(model.NumCameras());
  for (int i = 0; i < model.NumCameras(); ++i) {
    scaled_gray_images[i] = model.GetCamera(i)->ReadImage();
    scaled_gray_images[i].convertTo(scaled_gray_images[i], CV_32F);
    cvtColor(scaled_gray_images[i], scaled_gray_images[i], CV_BGR2GRAY);
    for (int s = 0; s < FLAGS_scale; ++s) {
      pyrDown(scaled_gray_images[i], scaled_gray_images[i]);
      // resize(scaled_gray_images[i], scaled_gray_images[i], Size(0, 0), 0.5, 0.5);
    }
    // GaussianBlur(scaled_gray_images[i], scaled_gray_images[i], cv::Size(3, 3), 1);
  }

  int width = scaled_gray_images[FLAGS_ref_index].cols;
  int height = scaled_gray_images[FLAGS_ref_index].rows;
  int num_images = model.NumCameras();
  int patch_radius = FLAGS_patch_radius;
  furry::tiny::CostVolume cost_volume;
  cost_volume.Allocate(width, height, FLAGS_num_samples);
  cost_volume = 0;

  auto ref_camera = model.GetCamera(FLAGS_ref_index);
  Mat ref_image = scaled_gray_images[FLAGS_ref_index];
  auto match_func = ScopedPtr<MatchFunc>(new Sad);

  furry::ParFor(0, FLAGS_num_samples, [&](int i_sample) {
    double depth = depth_samples[i_sample];
    LOG(INFO) << "Sweeping Plane " << i_sample << " depth: " << depth;
    vector<Mat> homographies(num_images);
    vector<Mat> warped_images;
    auto ref_camra = model.GetCamera(FLAGS_ref_index);
    vector<furry::tiny::Camera> cameras;
    for (size_t i = 0; i < num_images; ++i) {
      Matrix3d h = model.GetCamera(i)->HomographyFrom(*ref_camra, depth);
      double actual_scale = pow(2, FLAGS_scale);
      // LOG(INFO) << "Actual scale: " << actual_scale;
      h.leftCols<2>() *= actual_scale;
      h.row(2) *= actual_scale;
      // h /= h(2, 2);
      homographies[i] = furry::To<Mat>(h);
    }
    for (size_t i = 0; i < num_images; ++i) {
      if (i == FLAGS_ref_index) continue;
      Mat image;
      warpPerspective(scaled_gray_images[i], image, homographies[i],
                      ref_image.size(), INTER_LINEAR + WARP_INVERSE_MAP);
      warped_images.push_back(image);
    }
    // if (i_sample == 4) {
    // for (size_t i = 0; i < warped_images.size(); ++i) {
    //   imwrite(StringPrintf("warp_%04d.jpg", i), warped_images[i]);
    // }
    // }
    vector<double> scores(warped_images.size());
    for (int x = patch_radius; x < width - patch_radius; ++x) {
      for (int y = patch_radius; y < height - patch_radius; ++y) {
        Range row_range(y - patch_radius, y + patch_radius + 1);
        Range col_range(x - patch_radius, x + patch_radius + 1);
        auto p2 = ref_image(row_range, col_range);
        for (size_t i = 0; i < warped_images.size(); ++i) {
          auto p1 = warped_images[i](row_range, col_range);
          // scores[i] = 1 - CalcNcc(p1, p2);
          // scores[i] = CalcSos(p1, p2);
          scores[i] = (*match_func)(p1, p2);
        }
        double score = MergeScores(scores);
        CHECK(!std::isnan(score)) << x << ' ' << y << " gives a nan score";
        cost_volume.Set(y, x, i_sample, score);
      }
    }
    });

  Mat ref_color_image = model.GetCamera(FLAGS_ref_index)->ReadImage();
  for (int i = 0; i < FLAGS_scale; ++i) {
    pyrDown(ref_color_image, ref_color_image);
  }
  LOG(INFO) << "Modulating cost volume";
  CHECK(!cost_volume.HasNaN()) << "has nan before modulation";
  cost_volume.Modulate();
  CHECK(!cost_volume.HasNaN()) << "has nan after modulation";
  LOG(INFO) << "Setting cost volume properties";
  PlaneSweepingData psd;
  psd.min_depth = min_depth;
  psd.max_depth = max_depth;
  psd.patch_radius = patch_radius;
  psd.scale = FLAGS_scale;
  psd.num_samples = FLAGS_num_samples;
  psd.intrinsic_matrix = ref_camera->internals().K();
  psd.image = ref_color_image;
  cost_volume.SetPlaneSweepingData(psd);
  // cost_volume.SetImage(ref_color_image);
  // cost_volume.SetMinDepth(min_depth);
  // cost_volume.SetMaxDepth(max_depth);
  // cost_volume.SetPatchRadius(patch_radius);
  // cost_volume.SetIntrinsicMatrix(ref_camera->internals().K());
  LOG(INFO) << "Writing cost volume to file";
  cost_volume.WriteFile(FLAGS_out_cv);

  return 0;
}

void GetDepthRange(const furry::tiny::Model &model, int ref_index,
                   double *min_depth, double *max_depth) {
  auto ref_camera = model.GetCamera(ref_index);
  Matrix3d R = ref_camera->R();
  Vector3d C = ref_camera->C();
  vector<double> depths(model.NumTracks());
  for (int i = 0; i < model.NumTracks(); ++i) {
    Vector3d p = R * model.GetTrack(i)->GetPoint() + C;
    depths[i] = p[2];
  }
  sort(depths.begin(), depths.end());
  if (min_depth != nullptr) *min_depth = depths[depths.size() * 0.05] * 0.80;
  if (max_depth != nullptr) *max_depth = depths[depths.size() * 0.95] * 1.25;
}

double InversePerspectiveSampling(const double min, const double max,
                                  const double value) {

  return (max * min) / (max - (max - min) * value);
}


void InversePerspectiveSampling(double min, double max, int num_samples,
                                vector<double> *depth_samples) {
  depth_samples->resize(num_samples);
  double step = 1.0 / (num_samples - 1);
  for (int i = 0; i < num_samples; ++i) {
    (*depth_samples)[i] = InversePerspectiveSampling(min, max, i * step);
  }
  // double step = (max - min) / (num_samples - 1);
  // for (int i = 0; i < num_samples; ++i) {
  //   (*depth_samples)[i] = min + step * i;
  // }
}

double CalcNcc(const Mat &p0, const Mat &p1) {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  // assert(p0.type() == CV_8UC1);
  assert(p0.type() == CV_32FC1);
  double mu0, mu1, sigma0, sigma1, cov;
  mu0 = mu1 = sigma0 = sigma1 = cov = 0;
  int num_pixels = p0.cols * p0.rows;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      mu0 += p0.at<float>(r, c);
      mu1 += p1.at<float>(r, c);
    }
  }
  mu0 /= num_pixels;
  mu1 /= num_pixels;

  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      sigma0 += (p0.at<float>(r, c) - mu0) * (p0.at<float>(r, c) - mu0);
      sigma1 += (p1.at<float>(r, c) - mu1) * (p1.at<float>(r, c) - mu1);
      cov += (p0.at<float>(r, c) - mu0) * (p1.at<float>(r, c) - mu1);
    }
  }
  return cov / sqrt(sigma0) / sqrt(sigma1);
}

double CalcSos(const Mat &p0, const Mat &p1) {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  assert(p0.type() == CV_32FC1);
  double sos = 0;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      double diff = p0.at<float>(r, c) - p1.at<float>(r, c);
      sos += diff * diff;
    }
  }
  return sos;
}

double Sad::operator () (const Mat &p0, const Mat &p1) const {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  assert(p0.type() == CV_32FC1 && p1.type() == CV_32FC1);
  double sos = 0;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      double diff = p0.at<float>(r, c) - p1.at<float>(r, c);
      sos += fabs(diff);
    }
  }
  return sos;
}

double Ssd::operator () (const Mat &p0, const Mat &p1) const {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  assert(p0.type() == CV_32FC1);
  double sos = 0;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      double diff = p0.at<float>(r, c) - p1.at<float>(r, c);
      sos += diff * diff;
    }
  }
  return sos;
}

double Ncc::operator () (const Mat &p0, const Mat &p1) const {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  // assert(p0.type() == CV_8UC1);
  assert(p0.type() == CV_32FC1);
  double mu0, mu1, sigma0, sigma1, cov;
  mu0 = mu1 = sigma0 = sigma1 = cov = 0;
  int num_pixels = p0.cols * p0.rows;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      mu0 += p0.at<float>(r, c);
      mu1 += p1.at<float>(r, c);
    }
  }
  mu0 /= num_pixels;
  mu1 /= num_pixels;

  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      sigma0 += (p0.at<float>(r, c) - mu0) * (p0.at<float>(r, c) - mu0);
      sigma1 += (p1.at<float>(r, c) - mu1) * (p1.at<float>(r, c) - mu1);
      cov += (p0.at<float>(r, c) - mu0) * (p1.at<float>(r, c) - mu1);
    }
  }
  return 1 - cov / sqrt(sigma0) / sqrt(sigma1);
}

double MergeScores(vector<double> &scores) {
  double score = 0;
  double valid_ratio = 0.5;
  int num_valid_scores = scores.size() * valid_ratio;
  sort(scores.begin(), scores.end());
  for (int i = 0; i < num_valid_scores; ++i) {
    score += scores[i];
  }
  score /= num_valid_scores;
  return score;
}

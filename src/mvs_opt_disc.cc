#include "mvs_opt.h"

#include <algorithm>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/3rdparty/mrf/GCoptimization.h"
#include "furry/3rdparty/mrf/TRW-S.h"
#include "furry/3rdparty/mrf/MaxProdBP.h"
#include "furry/3rdparty/mrf/ICM.h"

#include "densecrf.h"

using namespace cv;
using namespace Eigen;
using namespace std;

namespace furry {
namespace tiny {

void WtaDepth(const CostVolume &cost_volume, DepthMap *depthmap) {
  LOG(INFO) << "Solving by wta";
  depthmap->Allocate(cost_volume.NumRows(), cost_volume.NumCols());
  for (int r = 0; r < cost_volume.NumRows(); ++r) {
    for (int c = 0; c < cost_volume.NumCols(); ++c) {
      auto samples = cost_volume.At(r, c);
      uint8_t best_label = depthmap->GetUnknownDepth();
      float best_cost = numeric_limits<float>::max();
      // float best_cost = numeric_limits<float>::min();
      for(int s = cost_volume.NumSamples() - 1; s > -1; --s) {
        if (samples[s] >= 0 && samples[s] < best_cost) {
          best_label = s;
          best_cost = samples[s];
        }
      }
      depthmap->SetDepth(r, c, best_label);
    }
  }
}

class UnaryCost {
 public:
  UnaryCost(const CostVolume *cost_volume) : cost_volume_(cost_volume) {}

  double operator () (int pix, int label) {
    return cost_volume_->GetSampleOfPixelIndex(pix, label);
  }

 private:
  const CostVolume *cost_volume_;
};

class BinaryCost {
 public:
  struct Options {
    float panelty;
    float max_panelty;
    float color_max_weight;
    float color_sigma;
  };

  BinaryCost(const CostVolume *cost_volume, const Options &options)
      : cost_volume_(cost_volume), options_(options) {
    // lab_image_.create(cost_volume_->GetImage().size(), CV_32FC3);
    cv::cvtColor(cost_volume_->GetImage(), lab_image_, CV_BGR2Lab);
  }

  double operator () (int pix0, int pix1, int label0, int label1) {
    // return min(options_.panelty * abs(label0 - label1),
    //            options_.max_panelty) *
    //     options_.color_max_weight * exp(
    //         - Dist(*(lab_image_.ptr<cv::Vec3f>() + pix0),
    //                *(lab_image_.ptr<cv::Vec3f>() + pix1))
    //         / options_.color_sigma);
    return options_.panelty * std::min<float>(abs(label0 - label1),
                                              options_.max_panelty);
  }

  void GetHWeights(std::vector<double> *weights) {
    CHECK(lab_image_.type() == CV_8UC3);
    int width = cost_volume_->GetWidth();
    int height = cost_volume_->GetHeight();
    weights->resize(width * height);
    int x, y, pix;
    for (y = 0; y < height; ++y) {
      for (x = 1; x < width; ++x) {
        pix = x + y * width;
        (*weights)[pix - 1] =
            options_.color_max_weight * exp(
                - Dist(*(lab_image_.ptr<cv::Vec3b>() + pix),
                       *(lab_image_.ptr<cv::Vec3b>() + pix - 1))
                / options_.color_sigma);
      }
    }
  }

  void GetVWeights(std::vector<double> *weights) {
    CHECK(lab_image_.type() == CV_8UC3);
    int width = cost_volume_->GetWidth();
    int height = cost_volume_->GetHeight();
    weights->resize(width * height);
    int x, y, pix;
    for (y = 1; y < height; ++y) {
      for (x = 0; x < width; ++x) {
        pix = x + y * width;
        (*weights)[pix - width] =
            options_.color_max_weight * exp(
                - Dist(*(lab_image_.ptr<cv::Vec3b>() + pix),
                       *(lab_image_.ptr<cv::Vec3b>() + pix - width))
                / options_.color_sigma);
      }
    }
  }

 private:

  template <typename T>
  float Dist(cv::Vec<T, 3> c0, cv::Vec<T, 3> c1) {
    return sqrt((c0[0] - c1[0]) * (c0[0] - c1[0]) +
                (c0[1] - c1[1]) * (c0[1] - c1[1]) +
                (c0[2] - c1[2]) * (c0[2] - c1[2]));
  }

 private:
  Options options_;
  cv::Mat lab_image_;
  const CostVolume *cost_volume_;
};

void SolveDiscreteDepth(const DiscreteMrfSettings &settings,
                        const CostVolume &cost_volume,
                        DepthMap *depthmap) {
  LOG(INFO) << "Solving by using Graph Cut";
  UnaryCost unary_cost(&cost_volume);
  BinaryCost::Options binary_options;
  binary_options.panelty = settings.fos().weight();
  binary_options.max_panelty = 0.15 *
      cost_volume.NumSamples();
  binary_options.color_max_weight = settings.fos().color_max_weight();
  binary_options.color_sigma = settings.fos().color_sigma();
  BinaryCost binary_cost(&cost_volume, binary_options);
  vector<double> v_weights, h_weights;
  binary_cost.GetVWeights(&v_weights);
  binary_cost.GetHWeights(&h_weights);
  DataCost data_term(unary_cost);
  // SmoothnessCost binary_term(binary_cost);
  SmoothnessCost binary_term(1,
                             binary_options.max_panelty,
                             binary_options.panelty,
                             h_weights.data(),
                             v_weights.data());
  // SmoothnessCost binary_term(1,
  //                            binary_options.max_panelty,
  //                            binary_options.panelty);
  EnergyFunction eng(&data_term, &binary_term);
  MRF* mrf = new Swap(cost_volume.GetWidth(), cost_volume.GetHeight(),
                           cost_volume.NumSamples(), &eng);
  // MRF* mrf = new MaxProdBP(cost_volume.GetWidth(), cost_volume.GetHeight(),
  //                          cost_volume.NumSamples(), &eng);
  mrf->initialize();
  mrf->clearAnswer();
  LOG(INFO) << "Finished initializing MRF";
  float time;
  for (int i = 0; i < settings.num_iterations(); ++i) {
    mrf->optimize(1, time);
    LOG(INFO) << i << " iteration took " << time << " sec";
    LOG(INFO) << "Total energy: " << mrf->totalEnergy();
  }

  LOG(INFO) << "Writing depth map";

  depthmap->Allocate(cost_volume.NumRows(), cost_volume.NumCols());
  for (int p = 0; p < cost_volume.NumRows() * cost_volume.NumCols(); ++p) {
    depthmap->SetDepthOfPixelIndex(p, mrf->getLabel(p));
  }

  delete mrf;
}

void SolveByDenseCrf(const DenseCrfSettings &settings,
                     const CostVolume &cost_volume,
                     DepthMap *depthmap) {
  LOG(INFO) << "Solving by using dense CRF";
  int width = cost_volume.GetWidth();
  int height = cost_volume.GetHeight();
  int num_labels = cost_volume.NumSamples();
  int max_panelty = num_labels * settings.max_panelty_ratio();
  // int max_panelty = 0;
  Mat image = cost_volume.GetImage();
  cvtColor(image, image, CV_BGR2Lab);
  MatrixXf unary;
  GetUnary(cost_volume, &unary);
  DenseCRF2D crf(width, height, num_labels);
  crf.setUnaryEnergy( unary );
  // add a color independent term (feature = pixel location 0..W-1, 0..H-1)
  // x_stddev = 3
  // y_stddev = 3
  // weight = 3
  if (settings.use_p()) {
    float x_stddev = settings.p_std();
    float y_stddev = settings.p_std();
    float weight = settings.p_weight();
    // crf.addPairwiseGaussian( x_stddev, y_stddev,
    //                          new PottsCompatibility( weight ) );
    crf.addPairwiseGaussian(
        x_stddev, y_stddev,
        // new PottsCompatibility( weight ) );
        new TruncatedLinearCompatibility(weight, max_panelty) );
  }

  if (settings.use_b()) {
    float x_stddev = settings.b_pstd();
    float y_stddev = settings.b_pstd();
    float r_stddev, g_stddev, b_stddev;
    r_stddev = g_stddev = b_stddev = settings.b_cstd();
    float weight = settings.b_weight();
    CHECK(image.isContinuous());
    CHECK(image.type() == CV_8UC3);
    crf.addPairwiseBilateral(
        x_stddev, y_stddev,
        r_stddev, g_stddev, b_stddev,
        image.ptr(),
        // new PottsCompatibility( weight ) );
        new TruncatedLinearCompatibility(weight, max_panelty));
  }

  if (settings.use_b2()) {
    float x_stddev = settings.b2_pstd();
    float y_stddev = settings.b2_pstd();
    float r_stddev, g_stddev, b_stddev;
    r_stddev = g_stddev = b_stddev = settings.b2_cstd();
    float weight = settings.b2_weight();
    CHECK(image.isContinuous());
    CHECK(image.type() == CV_8UC3);
    crf.addPairwiseBilateral(
        x_stddev, y_stddev,
        r_stddev, g_stddev, b_stddev,
        image.ptr(),
        // new PottsCompatibility( weight ) );
        new TruncatedLinearCompatibility(weight, max_panelty));
  }

  VectorXs map;
  VectorXf confidence;
  crf.map(settings.num_iterations(), &map, &confidence);
  // Store the result
  // unsigned char *res = colorize( map, W, H );
  depthmap->Allocate(height, width);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int index = i * width + j;
      if (confidence.coeff(index) >= settings.min_confidence() / num_labels)
        depthmap->SetDepth(i, j, map.coeff(index));
    }
  }
}

} // tiny
} // furry

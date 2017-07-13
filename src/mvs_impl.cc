#include "mvs.h"

#include <iostream>

#include <glog/logging.h>
#include <opencv2/opencv.hpp>

#include "furry/common/log.h"
#include "furry/common/parallel.h"
#include "furry/common/str.h"

#include "parallel.h"

using namespace Eigen;
using namespace std;
using namespace cv;

namespace {

bool IsInImage(const Mat& image, const Point2f& p) {
  if (p.x >= 0 && p.y >= 0 &&
      p.x < image.cols - 0.5 && p.y < image.rows - 0.5) {
    return true;
  } else {
    return false;
  }
}

template <typename T>
T* FirstIfNotNull(T* p0, T* p1) {
  if (p0 == nullptr) {
    return p1;
  } else {
    return p0;
  }
}

Point2f TransformPoint(const Matrix3d& H, const Point2f& p) {
  Vector3d pp = H * Vector3d(p.x, p.y, 1);
  return Point2f(pp.x() / pp.z(), pp.y() / pp.z());
}

double CalcNcc(const Mat& p0, const Mat& p1) {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  assert(p0.type() == CV_8UC1);
  double mu0, mu1, sigma0, sigma1, cov;
  mu0 = mu1 = sigma0 = sigma1 = cov = 0;
  int num_pixels = p0.cols * p0.rows;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      mu0 += p0.at<uint8_t>(r, c);
      mu1 += p1.at<uint8_t>(r, c);
    }
  }
  mu0 /= num_pixels;
  mu1 /= num_pixels;

  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      sigma0 += (p0.at<uint8_t>(r, c) - mu0) * (p0.at<uint8_t>(r, c) - mu0);
      sigma1 += (p1.at<uint8_t>(r, c) - mu1) * (p1.at<uint8_t>(r, c) - mu1);
      cov += (p0.at<uint8_t>(r, c) - mu0) * (p1.at<uint8_t>(r, c) - mu1);
    }
  }
  return cov / sqrt(sigma0) / sqrt(sigma1);
}

template <typename T>
void ExtractPatch(const cv::Mat& image,
                  const Matrix3d& Homography,
                  const cv::Point2f& p,
                  const cv::Size& size,
                  int stride,
                  Mat *patch_,
                  Mat *mask_ = nullptr,
                  uint8_t *completed_ = nullptr) {
  assert(image.type() == DataType<T>::type);
  Mat local_patch;
  Mat local_mask;
  uint8_t local_completed;
  Mat *patch, *mask;
  uint8_t *completed;
  patch = FirstIfNotNull(patch_, &local_patch);
  mask = FirstIfNotNull(mask_, &local_mask);
  completed = FirstIfNotNull(completed_, &local_completed);
  // Mat patch(size, image.type());
  patch->create(size, image.type());
  *patch = CV_RGB(0, 0, 0);
  mask->create(size, CV_8U);
  *completed = true;
  for (int i = 0; i < size.width; ++i) {
    for (int j = 0; j < size.height; ++j) {
      double x = (i - (size.width + 1) / 2) * stride + p.x;
      double y = (j - (size.height + 1) / 2) * stride + p.y;
      Point2f image_point = TransformPoint(Homography, Point2f(x, y));
      if (IsInImage(image, image_point)) {
        patch->at<T>(j, i) = image.at<T>(image_point);
        mask->at<uint8_t>(j, i) = 1;
      } else {
        mask->at<uint8_t>(j, i) = 0;
        *completed = false;
      }
    }
  }
}

}

namespace furry {
namespace tiny {

void MvsModel::SetCameras(const std::vector<Camera*>& all_cameras,
                           const std::vector<int>& camera_indexes,
                           int ref_camera_index) {
  ref_camera_index_ = ref_camera_index;
  cameras_.resize(0);
  for (int index : camera_indexes) {
    cameras_.push_back(all_cameras[index]);
  }
  for (auto& camera : cameras_) {
    // images_.push_back(cv::imread(StringPrintf("%s/%s", dirpath.c_str(),
    //                                          camera.filename().c_str())));
    // images_.push_back(camera->ReadUndistortedImage());
    images_.push_back(camera->ReadImage());
  }
  gray_images_.resize(images_.size());
  for (size_t i = 0; i < images_.size(); ++i) {
    cvtColor(images_[i], gray_images_[i], CV_BGR2GRAY);
  }
  ref_image_ = images_[ref_camera_index_];
}

int MvsModel::num_cameras() const {
  return cameras_.size();
}

const Camera* MvsModel::camera(int i) const {
  return cameras_[i];
}

const Camera* MvsModel::ref_camera() const {
  return cameras_[ref_camera_index_];
}

int MvsModel::ref_camera_index() const {
  return ref_camera_index_;
}

const cv::Mat& MvsModel::image(int i) const {
  return images_[i];
}

const std::vector<cv::Mat>& MvsModel::images() const {
  return images_;
}

vector<Eigen::Matrix3d> MvsModel::Homographies(double depth) const {
  vector<Matrix3d> homographies;
  for (int i = 0; i < num_cameras(); ++i) {
    homographies.push_back(camera(i)->HomographyFrom(*ref_camera(), depth));
  }
  return homographies;
}

void MvsModel::SetRefPoint(const cv::Point2f& p)  {
  ref_point_ = p;
}

void MvsModel::set_ref_depth(double d) {
  ref_depth_ = d;
  homographies_.resize(num_cameras());
  for (int i = 0; i < num_cameras(); ++i) {
    homographies_[i] = camera(i)->HomographyFrom(*ref_camera(), d);
    // cout << homographies_[i] << "\n\n";
  }
}

void MvsModel::SetDepthRange(float lower, float upper) {
  depth_range_[0] = lower;
  depth_range_[1] = upper;
}

void MvsModel::SetPatchSize(const cv::Size size) {
  patch_size_ = size;
}

void MvsModel::SetPatchRadius(int radius) {
  patch_size_ = Size((radius << 1) + 1, (radius << 1) + 1);
}

void MvsModel::set_patch_stride(double stride) {
  patch_stride_ = stride;
}

void MvsModel::SetDepthStep(float depth_step) {
  depth_step_ = depth_step;
}

const cv::Mat& MvsModel::GetRefImage() const {
  return images_[ref_camera_index_];
}

const cv::Point2f& MvsModel::GetRefPoint() const {
  return ref_point_;
}

double MvsModel::ref_depth() const {
  return ref_depth_;
}

Vector2f MvsModel::depth_range() const {
  return depth_range_;
}

float MvsModel::depth_step() const {
  return depth_step_;
}

std::vector<cv::Mat> MvsModel::MarkRefPointOnImages() const {
  vector<Mat> marked_images(num_cameras());
  auto points = ReprojectRefPointOnCameras();
  for (int i = 0; i < num_cameras(); ++i) {
    images_[i].copyTo(marked_images[i]);
    circle(marked_images[i], points[i], 20, CV_RGB(255, 0, 0), -1);
  }
  return marked_images;
}

std::vector<cv::Point2f> MvsModel::ReprojectRefPointOnCameras() const {
  vector<Point2f> points(num_cameras());
  for (int i = 0; i < num_cameras(); ++i) {
    // LOG(INFO) << "h: " << homographies_[i];
    points[i] = TransformPoint(homographies_[i], ref_point_);
  }
  return points;
}

std::vector<cv::Mat> MvsModel::ExtractPatches() const {
  std::vector<cv::Mat> patches(num_cameras());
  for (int i = 0; i < num_cameras(); ++i) {
    ExtractPatch<Vec3b>(images_[i], homographies_[i], ref_point_,
                        patch_size_, patch_stride_,
                        &patches[i]);
  }
  return patches;
}

std::vector<cv::Mat> MvsModel::ExtractGrayPatches() const {
  std::vector<cv::Mat> patches(num_cameras());
  for (int i = 0; i < num_cameras(); ++i) {
    ExtractPatch<uint8_t>(gray_images_[i], homographies_[i],
                          ref_point_, patch_size_, patch_stride_,
                          &patches[i]);
  }
  return patches;
}

std::vector<float> MvsModel::CalcPatchNcc() const {
  auto gray_patches = ExtractGrayPatches();
  vector<float> ncc_scores(num_cameras());
  for (int i = 0; i < num_cameras(); ++i) {
    ncc_scores[i] = CalcNcc(gray_patches[i], gray_patches[ref_camera_index()]);
  }
  return ncc_scores;
}

std::vector<std::vector<cv::Point2f>> MvsModel::CalcNccOnDepth() const {
  vector<vector<Point2f>> scores(num_cameras());
  int num_steps = (depth_range_[1] - depth_range_[0]) / depth_step_ + 1;
  // for (int i = 0; i < num_steps; ++i) {
  for (int i = 0; i < num_cameras(); ++i) {
    scores[i].resize(num_steps);
  }
  ParFor(0, num_steps, [&](int i) {
      float depth = i * depth_step_ + depth_range_[0];
      auto Hs = Homographies(depth);
      vector<Mat> patches(num_cameras());
      vector<uint8_t> completed(num_cameras());
      for (int j = 0; j < num_cameras(); ++j) {
        ExtractPatch<uint8_t>(
            gray_images_[j], Hs[j],
            ref_point_, patch_size_, patch_stride_,
            &patches[j], nullptr, &completed[j]);
      }
      if (!completed[ref_camera_index_]) {
        for (int j = 0; j < num_cameras(); ++j) {
          scores[j][i] = Point2f(depth, 0.0f);
        }
      } else {
        for (int j = 0; j < num_cameras(); ++j) {
          if (j == ref_camera_index_) {
            scores[j][i] = Point2f(depth, 1.0f);
          } else if (completed[j]) {
            scores[j][i] = Point2f(
                depth, CalcNcc(patches[j], patches[ref_camera_index_]));
          } else {
            scores[j][i] = Point2f(depth, 0.0f);
          }
        }
      }
    });
  return scores;
}

std::vector<std::vector<cv::Point2f>> MvsModel::TakeCostVolumeSamples() const {
  auto size = GetRefImage().size();
  Point2f ref_point = ref_point_;
  ref_point.x *= (float)cost_volume_.GetWidth() / size.width;
  ref_point.y *= (float)cost_volume_.GetHeight() / size.height;
  auto sample_ptr = cost_volume_.GetSamples(ref_point);
  vector<vector<cv::Point2f>> points(1);
  points[0].reserve(cost_volume_.NumSamples());
  // LOG(INFO) << "# Samples: " << cost_volume_.NumSamples() << ' ' << cost_volume_.GetSize(2);
  for (int i = 0; i < cost_volume_.NumSamples(); ++i) {
    points[0].push_back(Point2f(i, *sample_ptr++));
  }
  // LOG(INFO) << "Drawing "
  //     // << Join(points_to_draw.begin(), points_to_draw.end(), " ");
  //           << Join(points[0].begin(), points[0].end(), " ");
  return points;
}

std::vector<std::vector<cv::Point2f>> MvsModel::GetPhotoConsistency() const {
  if (use_cost_volume_) {
    return TakeCostVolumeSamples();
  } else {
    return CalcNccOnDepth();
  }
}

void MvsModel::SetCostVolume(CostVolume &&cost_volume) {
  cost_volume_ = move(cost_volume);
  SetDepthRange(cost_volume_.GetMinDepth(), cost_volume_.GetMaxDepth());
  SetDepthStep((depth_range_[1] - depth_range_[0]) / cost_volume_.NumSamples());
  // SetPatchRadius(cost_volume_.GetPatchRadius());
  use_cost_volume_ = true;
}

void MvsModel::GenerateCostVolume(CostVolume *cost_volume,
                                  DepthMap *depth_map) {
  int num_samples = (depth_range_[1] - depth_range_[0]) / depth_step_;
  cost_volume->Allocate(ref_image_.cols, ref_image_.rows, num_samples);
  (*cost_volume) = 0;
  cost_volume->SetMinDepth(depth_range_[0]);
  cost_volume->SetMaxDepth(depth_range_[1]);
  cost_volume->SetImage(ref_image_);
  depth_map->Allocate(ref_image_.rows, ref_image_.cols);
  int patch_radius = (patch_size_.width - 1) / 2;
  // for (int x = patch_radius; x < ref_image_.cols - patch_radius; ++x) {
  int num_cols = 1000;
  furry::ProgressBar generating_cost_volume("Generating Cost Volume",
                                            // ref_image_.cols - 2 * patch_radius);
                                            num_cols - patch_radius);
  // ParFor(patch_radius, ref_image_.cols - patch_radius, [&](int x) {
  // ParFor(patch_radius, 500, [&](int x) {
  for (int x = patch_radius; x < num_cols; ++x) {
      // LOG(INFO) << "X: " << x;
    for (int y = patch_radius; y < ref_image_.rows - patch_radius; ++y) {
      SetRefPoint(Point2f(x, y));
      auto points = GetPhotoConsistency();
      vector<Point2f> sum_points(points[0].size());
      for (size_t i = 0; i < points[0].size(); ++i) {
        Point2f p(points[0][i].x, 0);
        for (size_t j = 0; j < points.size(); ++j) {
          p.y += points[j][i].y;
        }
        p.y /= points.size();
        sum_points[i] = p;
      }
      float max_label = 0;
      double max_value = sum_points[0].y;
      for (int s = 0; s < num_samples; ++s) {
        cost_volume->Set(y, x, s, 1 - sum_points[s].y);
        if (sum_points[s].y > max_value) {
          max_label = s;
        }
      }
      depth_map->SetDepth(y, x, max_label);
    }
    generating_cost_volume.FinishOne();
    // });
  }
  generating_cost_volume.Done();
}

bool MvsModel::UseCostVolume() const {
  return use_cost_volume_;
}

} // tiny
} // furry

#include "cost_volume.h"

#include <algorithm>
#include <fstream>
#include <string>

#include <glog/logging.h>

using namespace std;

namespace {

float InversePerspectiveSampling(const float min, const float max,
                                 const float value) {
  return (max * min) / (max - (max - min) * value);
}

void GetMin(float *values, int size, float *first, float *second) {
  CHECK_GT(size, 1);
  float f, s;
  if (values[0] > values[1]) {
    f = values[1];
    s = values[0];
  } else {
    f = values[0];
    s = values[1];
  }

  float *last = values + size;
  values += 2;
  while (values < last) {
    if (*values < f) {
      f = *values;
    } else if (*values < s) {
      s = *values;
    }
    ++values;
  }

  if (first != nullptr) *first = f;
  if (second != nullptr) *second = s;
}

}

namespace furry {
namespace tiny {

CostVolume::CostVolume(CostVolume &&cost_volume) {
  *this = move(cost_volume);
}

CostVolume& CostVolume::operator = (CostVolume &&cost_volume) {
  if (this != &cost_volume) {
    Grid3f::operator = (move(cost_volume));
    psd_ = cost_volume.psd_;
    // min_depth_ = cost_volume.min_depth_;
    // max_depth_ = cost_volume.max_depth_;
    // image_ = cost_volume.image_;
    // cost_volume.min_depth_ = -1;
    // cost_volume.max_depth_ = -1;
  }
  return *this;
}

double CostVolume::GetDepth(int s) const {
  return InversePerspectiveSampling(
      psd_.min_depth, psd_.max_depth, s / (NumSamples() - 1.0f));
}

bool CostVolume::ReadFile(const std::string &filename) {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  if (file == nullptr) {
    LOG(ERROR) << "Can't open file " << filename << " to read";
    return false;
  }
  return ReadFile(file.Get());
}

bool CostVolume::WriteFile(const std::string &filename) const {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  if (!file) {
    LOG(WARNING) << "Can't open " << filename;
    return false;
  }
  return WriteFile(file.Get());
}

bool CostVolume::ReadFile(FILE *file) {
  if (!Grid3f::ReadBin(file)) return false;
  return psd_.ReadBin(file);
  // fread(&min_depth_, sizeof(double), 1, file);
  // fread(&max_depth_, sizeof(double), 1, file);
  // fread(intrinsic_matrix_.data(), sizeof(double), 9, file);
  // return ReadImage(file);
}

bool CostVolume::WriteFile(FILE *file) const {
  if (!Grid3f::WriteBin(file)) return false;
  return psd_.WriteBin(file);
  // fwrite(&min_depth_, sizeof(double), 1, file);
  // fwrite(&max_depth_, sizeof(double), 1, file);
  // fwrite(intrinsic_matrix_.data(), sizeof(double), 9, file);
  // return WriteImage(file);
}

void CostVolume::Modulate() {
  // PKRN method in
  // X. Hu and P. Mordohai
  // "A Quantitative Evaluation of Confidence Measures for Stereo Vision"
  float first, second, confidence;
  float *values;
  for (int r = 0; r < NumRows(); ++r) {
    for (int c = 0; c < NumCols(); ++c) {
      values = At(r, c);
      GetMin(values, NumSamples(), &first, &second);
      confidence = (second + 1) / (first + 1);
      // LOG_IF(INFO, first < 1e-10) << "first: " << first
      //                             << " second: " << second;
      for (int i = 0; i < NumSamples(); ++i) {
        values[i] = values[i] * confidence;
      }
    }
  }
}

void CostVolume::GetMinLabels(Eigen::MatrixXi *labels,
                              Eigen::MatrixXf *costs) const {
  labels->resize(NumRows(), NumCols());
  function<void (int, int, short, float)> record;
  if (costs == nullptr) {
    record = [labels](int r, int c, short label, float cost) {
      (*labels)(r, c) = label;
    };
  } else {
    record = [labels, costs](int r, int c, short label, float cost) {
      (*labels)(r, c) = label;
      (*costs)(r, c) = cost;
    };
    costs->resize(NumRows(), NumCols());
  }

  short best_label;
  float best_cost;
  const float *samples;
  for (int r = 0; r < NumRows(); ++r) {
    for (int c = 0; c < NumCols(); ++c) {
      samples = At(r, c);
      best_label = 0;
      best_cost = samples[0];
      for (int s = 0; s < NumSamples(); ++s) {
        if (samples[s] < best_cost) {
          best_label = s;
          best_cost = samples[s];
        }
      }
      record(r, c, best_label, best_cost);
    }
  }
}

void CostVolume::RemoveEndLabels() {
  Eigen::MatrixXi labels;
  GetMinLabels(&labels);
  int end0 = 0;
  int end1 = NumSamples() - 1;
  float *samples;
  for (int r = 0; r < labels.rows(); ++r) {
    for (int c = 0; c < labels.cols(); ++c) {
      if (labels(r, c) == end0 || labels(r, c) == end1) {
        samples = At(r, c);
        for (int s = 0; s < NumSamples(); ++s) {
          samples[s] = 0;
        }
      }
    }
  }
}

// void CostVolume::RemoveBorder() {
//   CostVolume new_volume;
//   int new_width = GetWidth() - 2 * patch_radius_;
//   int new_height = GetHeight() - 2 * patch_radius_;
//   new_volume.Allocate(new_width, new_height, NumSamples());

// }

// bool CostVolume::ReadImage(FILE *file) {
//   int type, rows, cols;
//   fread(&type, sizeof(type), 1, file);
//   fread(&rows, sizeof(rows), 1, file);
//   fread(&cols, sizeof(cols), 1, file);
//   LOG(INFO) << "type: " << type << " rows: " << rows << " cols: " << cols;
//   image_.create(rows, cols, type);
//   CHECK(image_.isContinuous());
//   fread(image_.ptr(), image_.elemSize(), image_.total(), file);
//   return true;
// }

// bool CostVolume::WriteImage(FILE *file) const {
//   int type = image_.type();
//   int rows = image_.rows;
//   int cols = image_.cols;
//   LOG(INFO) << "type: " << type << " rows: " << rows << " cols: " << cols;
//   fwrite(&type, sizeof(type), 1, file);
//   fwrite(&rows, sizeof(rows), 1, file);
//   fwrite(&cols, sizeof(cols), 1, file);
//   fwrite(image_.ptr(), image_.elemSize(), image_.total(), file);
//   return true;
// }

} // tiny
} // furry

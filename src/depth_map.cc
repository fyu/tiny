#include "depth_map.h"

#include <cstdio>
#include <string>
#include <fstream>
#include <limits>

#include <Eigen/Dense>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/3rdparty/ply/ply.h"
#include "furry/common/init.h"
#include "furry/common/path.h"
#include "furry/common/str.h"

#include "cost_volume.h"

using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

namespace {

float InversePerspectiveSampling(const float min, const float max,
                                 const float value) {
  return (max * min) / (max - (max - min) * value);
}

void ConvertLabelsToDepth(const cv::Mat &labels,
                          float min_depth,
                          float max_depth,
                          float num_samples,
                          cv::Mat *depth) {
  CHECK(labels.type() == cv::DataType<DepthMap::DepthType>::type);
  depth->create(labels.size(), CV_32F);
  for (int r = 0; r < labels.rows; ++r) {
    for (int c = 0; c < labels.cols; ++c) {
      if (labels.at<DepthMap::DepthType>(r, c) == DepthMap::kUnknownDepth) {
        depth->at<float>(r, c) = -1;
      } else {
        depth->at<float>(r, c) =
            InversePerspectiveSampling(
                min_depth, max_depth,
                labels.at<DepthMap::DepthType>(r, c) / (num_samples - 1.0f));
      }
    }
  }
}

void RenderDepth(const cv::Mat &depth_data, cv::Mat *depth_image) {
  double min_value;
  double max_value;
  // cv::minMaxLoc(depth_data, &min_value, &max_value, 0, 0, depth_data >= 0.0f);
  cv::minMaxLoc(depth_data, &min_value, &max_value, 0, 0, depth_data >= 0);
  LOG(INFO) << "Depthmap max: " << max_value << " min: " << min_value;
  double range = max_value - min_value;
  depth_image->create(depth_data.size(), CV_8UC3);
  for (int r = 0; r < depth_data.rows; ++r) {
    for (int c = 0; c < depth_data.cols; ++c) {
      if (depth_data.at<float>(r, c) >= 0) {
        float color = (depth_data.at<float>(r, c) - min_value) / range * 255;
        depth_image->at<cv::Vec3b>(r, c) = cv::Vec3b(color, color, color);
      } else {
        depth_image->at<cv::Vec3b>(r, c) = cv::Vec3b(0, 127, 255);
      }
    }
  }
}

struct Vertex {
  float x;
  float y;
  float z;
  uint8_t r;
  uint8_t g;
  uint8_t b;
};

bool WritePgm(const std::string &filename, const cv::Mat &image) {
  cv::Mat pgm_image = image;
  if (pgm_image.channels() == 3)
    cv::cvtColor(pgm_image, pgm_image, CV_BGR2GRAY);
  double min_value, max_value;
  minMaxLoc(pgm_image, &min_value, &max_value);
  CHECK(min_value >= 0 && max_value < 256);
  if (pgm_image.type() != CV_8U)
    pgm_image.convertTo(pgm_image, CV_8U);

  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  if (file == nullptr) return false;
  fprintf(file, "P5\n%d %d\n%d\n",
          pgm_image.cols, pgm_image.rows, int(max_value));
  CHECK(pgm_image.isContinuous());
  fwrite(pgm_image.ptr(), pgm_image.total(), 1, file);
  return true;
}

}

namespace furry {
namespace tiny {

int DepthMap::NumKnownDepths() const {
  int count = 0;
  for (int r = 0; r < NumRows(); ++r) {
    for (int c = 0; c < NumCols(); ++c) {
      if (GetDepth(r, c) != GetUnknownDepth()) ++count;
    }
  }
  return count;
}

bool DepthMap::WriteBin(const string &filename) const {
  CHECK(depth_.size() == psd_.image.size());
  CHECK(depth_.isContinuous());
  CHECK(psd_.image.isContinuous());
  FILE *file = fopen(filename.c_str(), "w");
  if (file == nullptr) return false;
  cv::Mat depth;
  ConvertLabelsToDepth(depth_, psd_.min_depth, psd_.max_depth, psd_.num_samples, &depth);
  // depth = depth_;
  int rows = depth.rows;
  int cols = depth.cols;
  fwrite(&rows, sizeof(rows), 1, file);
  fwrite(&cols, sizeof(rows), 1, file);
  int type = depth.type();
  fwrite(&type, sizeof(type), 1, file);
  fwrite(depth.ptr(), depth.elemSize(), depth.total(), file);
  type = psd_.image.type();
  fwrite(&type, sizeof(type), 1, file);
  cv::Mat rgb_image;
  cv::cvtColor(psd_.image, rgb_image, CV_BGR2RGB);
  // fwrite(psd_.image.ptr(), psd_.image.elemSize(), psd_.image.total(), file);
  fwrite(rgb_image.ptr(), psd_.image.elemSize(), psd_.image.total(), file);
  fclose(file);
  return true;
}

bool DepthMap::WriteDepthImage(const string &filename) const {
  cv::Mat depth_data, depth_image;
  // ConvertLabelsToDepth(depth_, psd_.min_depth, psd_.max_depth, psd_.num_samples,
  //                      &depth_data);
  string ext = GetExt(filename);
  if (ext == ".pgm") {
    return WritePgm(filename, depth_);
  } else {
    RenderDepth(depth_, &depth_image);
    return cv::imwrite(filename, depth_image);
  }
}

bool DepthMap::WritePly(const std::string &filename) const {
  PlyFile *ply;
  int num_elems;
  char **elist;
  int file_type;
  float version;
  // int num_vertexes = NumRows() * NumCols();
  int num_vertexes = NumKnownDepths();

  char *elem_names[] = {"vertex"};

  PlyProperty vertex_props[] = { /* list of property information for a vertex */
    {"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,x), 0, 0, 0, 0},
    {"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,y), 0, 0, 0, 0},
    {"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,z), 0, 0, 0, 0},
    {"red", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex,r), 0, 0, 0, 0},
    {"green", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex,g), 0, 0, 0, 0},
    {"blue", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex,b), 0, 0, 0, 0},
  };

  ply = ply_open_for_writing(
      (char*)filename.c_str(), 1, elem_names, PLY_ASCII, &version);

  ply_element_count(ply, "vertex", num_vertexes);
  for (int i = 0; i < 6; ++i) {
    ply_describe_property (ply, "vertex", &vertex_props[i]);
  }

  vector<Vertex> vertexes;
  Eigen::Matrix3d ik = psd_.intrinsic_matrix.inverse();
  // LOG(INFO) << "Inverse intrinsic matrix: " << ik;
  double actual_scale = pow(2, psd_.scale);
  for (int r = 0; r < NumRows(); ++r) {
    for (int c = 0; c < NumCols(); ++c) {
      if (GetDepth(r, c) == GetUnknownDepth()) continue;
      auto label = GetDepth(r, c) / (psd_.num_samples - 1);
      float depth = InversePerspectiveSampling(psd_.min_depth, psd_.max_depth, label);
      // float depth = label;
      // LOG(INFO) << "label: " << label
      //           << " max: " << psd_.max_depth
      //           << " min: " << psd_.min_depth
      //           << " depth: " << depth;
      Eigen::Vector3d dir = ik * Vector3d(c * actual_scale,
                                          r * actual_scale, 1);
      Eigen::Vector3d p = dir / dir(2) * depth;
      auto color = psd_.image.at<cv::Vec3b>(r, c);
      Vertex v;
      v.x = p[0];
      v.y = p[1];
      v.z = p[2];
      v.r = color[2];
      v.g = color[1];
      v.b = color[0];
      vertexes.push_back(v);
    }
  }

  ply_header_complete(ply);

  ply_put_element_setup (ply, "vertex");
  for (size_t i = 0; i < vertexes.size(); ++i) {
    ply_put_element(ply, (void*)&vertexes[i]);
  }
  ply_close (ply);
}

} // tiny
} // furry

#include <utility>
#include <iostream>

#include <Eigen/Dense>
#include <gflags/gflags.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "furry/common/str.h"
#include "furry/gui/window_QT.h"

#include "model.h"
#include "util_ui.h"

DEFINE_string(in_dir, "", "");
DEFINE_string(in_prefix, "", "");

using namespace cv;
using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

pair<Vector3d, Vector3d> GetEpipoles(const pair<Camera*, Camera*> &cameras) {
  return make_pair(cameras.first->ProjectPoint(cameras.second->C()),
                   cameras.second->ProjectPoint(cameras.first->C()));
}

Mat DrawEpipolarLines(const Mat &image, const Vector3d &epipole,
                      int num_lines = 10) {
  vector<Vector3d> start_points;
  Vector3d end_line;
  int step = image.rows / (num_lines + 1);
  double slop = (double)image.rows / image.cols;
  Vector3d diag_line1(slop, -1, 0);
  Vector3d diag_line2(slop, 1, -image.rows);
  double d1 = epipole.dot(diag_line1);
  double d2 = epipole.dot(diag_line2);
  if (d1 > 0) {
    if (d2 > 0) {
      for (int y = 0; y < image.rows; y += step) {
        start_points.push_back(Vector3d(0, y, 1));
      }
      end_line = Vector3d(1, 0, -image.cols);
    } else {
      for (int x = 0; x < image.cols; x += step) {
        start_points.push_back(Vector3d(x, image.rows - 1, 1));
      }
      end_line = Vector3d(0, 1, 0);
    }
  } else {
    if (d2 > 0) {
      for (int x = 0; x < image.cols; x += step) {
        start_points.push_back(Vector3d(x, 0, 1));
      }
      end_line = Vector3d(0, 1, -image.rows);
    } else {
      for (int y = 0; y < image.rows; y += step) {
        start_points.push_back(Vector3d(image.cols - 1, y, 1));
      }
      end_line = Vector3d(1, 0, 0);
    }
  }
  Mat show_image;
  image.copyTo(show_image);
  for (auto &start_point : start_points) {
    Vector3d line = start_point.cross(epipole);
    Vector3d end_point = line.cross(end_line);
    cv::line(show_image,
             Point2d(start_point[0], start_point[1]),
             Point2d(end_point[0] / end_point[2], end_point[1] / end_point[2]),
             kCvYellow,
             4,
             CV_AA);
  }
  return show_image;
}

int main(int argc, char *argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  vector<Camera*> all_cameras;
  ReadCameras(StringPrintf("%s/%s_cameras.txt",
                           FLAGS_in_dir.c_str(),
                           FLAGS_in_prefix.c_str()),
              &all_cameras);

  pair<Camera*, Camera*> cameras = make_pair(all_cameras[0], all_cameras[45]);
  pair<Mat, Mat> images;
  images.first = cameras.first->ReadImage();
  images.second = cameras.second->ReadImage();

  auto epipoles = GetEpipoles(cameras);
  Mat el_image1 = DrawEpipolarLines(images.first, epipoles.first);
  Mat el_image2 = DrawEpipolarLines(images.second, epipoles.second);
  // cout << "Epipoles:\n";
  // cout << epipoles.first.transpose() << ' '
  //      << epipoles.second.transpose() << '\n';

  QApplication app(argc, argv);
  furry::CvWindow image_window("epipolar");
  image_window.UpdateImage(VListImages({{el_image1, el_image2}}));
  app.exec();
  return 0;
}

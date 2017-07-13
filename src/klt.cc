#include <string>
#include <vector>
#include <cstdio>

#include "furry/gui/window_QT.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "furry/common/init.h"
#include "furry/common/str.h"
#include "furry/common/cv.h"
#include "furry/common/sequence.h"

#include "model.h"

using namespace cv;
using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_dir, "", "");
DEFINE_string(in_prefix, "", "");
DEFINE_string(images, "20,73", "");

void ChangeRotation(const Quaterniond &q, vector<Camera*> *cameras) {
  for (auto &camera : *cameras) {
    camera->SetQuaternion(q * camera->GetQuaternion());
  }
}

int main(int argc, char *argv[]) {
  furry::Init(&argc, &argv, "");
  vector<Camera*> all_cameras;
  ReadCameras(StringPrintf("%s/%s_cameras.txt",
                           FLAGS_in_dir.c_str(),
                           FLAGS_in_prefix.c_str()),
              &all_cameras);
  ChangeRotation(all_cameras[0]->GetQuaternion().conjugate(), &all_cameras);
  vector<string> image_indexes_str = furry::Split(FLAGS_images, ",");
  vector<int> image_indexes = To<vector<int>>(image_indexes_str);
  Camera *cameras[2];
  Mat images[2];
  Mat masks[2];
  for (int i = 0; i < 2; ++i) {
    cameras[i] = all_cameras[image_indexes[i] - 1];
    cout << "Reading " << cameras[i]->GetFilename() << '\n';
    images[i] = cameras[i]->RemoveRotation(
        cameras[i]->Undistort(cameras[i]->ReadImage()),
        &masks[i]);
    cvtColor(images[i], images[i], CV_BGR2GRAY);
  }

  vector<Point2f> corners;
  goodFeaturesToTrack(images[0], corners, 10000, 0.01, 11, masks[0], 11);
  cornerSubPix(images[0], corners, Size(5, 5), Size(-1, -1),
               TermCriteria(TermCriteria::EPS, 0, 0.001));
  cout << corners.size() << " corners detected\n";
  vector<Point2f> tracked_corners;
  vector<uint8_t> status;
  vector<float> err;
  calcOpticalFlowPyrLK(images[0], images[1], corners, tracked_corners,
                      status, err);
  corners = Filter(corners, status);
  tracked_corners = Filter(tracked_corners, status);
  vector<Point2f> canonical_corners;
  vector<Point2f> canonical_tracked_corners;
  perspectiveTransform(corners, canonical_corners,
                       To<Mat>(Matrix3d(cameras[0]->K().inverse())));
  perspectiveTransform(tracked_corners, canonical_tracked_corners,
                       To<Mat>(Matrix3d(cameras[1]->K().inverse())));
  Mat show_image = DrawPoints(images[0], corners,
                              vector<Scalar>(1, CV_RGB(0, 0, 0)), 5);
  imwrite("image1.png", show_image);
  show_image = DrawPoints(images[1], tracked_corners,
                          vector<Scalar>(1, CV_RGB(0, 0, 0)), 5);
  imwrite("image2.png", show_image);

  vector<Point2f> flow(corners.size());
  for (size_t i = 0; i < corners.size(); ++i) {
    flow[i].x = canonical_tracked_corners[i].x - canonical_corners[i].x;
    flow[i].y = canonical_tracked_corners[i].y - canonical_corners[i].y;
    // cout << flow[i] << '\n';
  }

  vector<Point2f> flow_inv(flow.size());
  for (size_t i = 0; i < flow.size(); ++i) {
    flow_inv[i].x = 1.0 / flow[i].x;
    flow_inv[i].y = 1.0 / flow[i].y;
  }

  FILE *file = fopen("depth_x.txt", "w");
  for (size_t i = 0; i < corners.size(); ++i) {
    fprintf(file, "%f %f %f\n", corners[i].x, corners[i].y, flow_inv[i].x);
  }
  fclose(file);

  file = fopen("depth_y.txt", "w");
  for (size_t i = 0; i < corners.size(); ++i) {
    fprintf(file, "%f %f %f\n", corners[i].x, corners[i].y, flow_inv[i].y);
  }
  fclose(file);


  return 0;
}

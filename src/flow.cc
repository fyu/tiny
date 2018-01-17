#include <string>
#include <iostream>
#include <vector>

#include <gflags/gflags.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/common/str.h"
#include "furry/common/cv.h"
#include "furry/gui/window_QT.h"

#include "model.h"
#include "mvs.h"
#include "mvs_ui.h"

using namespace cv;
using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_dir, "", "");
DEFINE_string(in_prefix, "", "");
DEFINE_string(images, "20,73", "");

void ChangeRotation(const Quaterniond &q, vector<Camera*> *cameras) {
  for (auto camera : *cameras) {
    camera->SetQuaternion(q * camera->GetQuaternion());
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  vector<Camera*> all_cameras;
  ReadCameras(StringPrintf("%s/%s_cameras.txt",
                           FLAGS_in_dir.c_str(),
                           FLAGS_in_prefix.c_str()),
              &all_cameras);
  ChangeRotation(all_cameras[0]->GetQuaternion().conjugate(), &all_cameras);
  vector<string> image_indexes_str = furry::Split(FLAGS_images, ",");
  vector<int> image_indexes = To<vector<int>>(image_indexes_str);
  {
  auto camera = all_cameras[image_indexes[0] - 1];
  auto image = camera->Undistort(camera->ReadImage());
  auto rect_image = camera->RemoveRotation(image);
  imwrite("rect_image_" + image_indexes_str[0] + ".png", rect_image);
  }

  {
  auto camera = all_cameras[image_indexes[1] - 1];
  auto image = camera->Undistort(camera->ReadImage());
  auto rect_image = camera->RemoveRotation(image);
  imwrite("rect_image_" + image_indexes_str[1] + ".png", rect_image);
  }

  return 0;
}

#include <string>
#include <algorithm>

#include "furry/gui/window_QT.h"

#include <gflags/gflags.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "furry/common/init.h"

#include "model.h"
#include "mvs.h"
#include "mvs_ui.h"

using namespace cv;
using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_model, "", "");
DEFINE_string(in_cv, "", "");
DEFINE_double(depth_step, 0.1, "");

Point2f TransformPoint(const Matrix3d& H, const Point2f& p) {
  Vector3d pp = H * Vector3d(p.x, p.y, 1);
  return Point2f(pp.x() / pp.z(), pp.y() / pp.z());
}

Mat GridImages(const vector<cv::Mat>& images) {
  int width = images[0].cols;
  int height = images[0].rows;
  int num_images = images.size();
  if (num_images < 4) {
    if (width > height) {
      cv::Mat show_image(height * num_images, width, images[0].type());
      for (int i = 0; i < num_images; ++i) {
        images[i].copyTo(show_image(cv::Range(i * height, (i + 1) * height),
                                    cv::Range::all()));
      }
      return show_image;
    } else {
      cv::Mat show_image(height, width * num_images, images[0].type());
      for (int i = 0; i < num_images; ++i) {
        images[i].copyTo(show_image(cv::Range::all(),
                                    cv::Range(i * width, (i + 1) * width)));
      }
      return show_image;
    }
  }
  return cv::Mat();
}

int main(int argc, char* argv[]) {
  Init(&argc, &argv, "");
  CHECK(!FLAGS_in_model.empty());
  vector<Camera*> all_cameras;
  ReadCameras(FLAGS_in_model + "_cameras.txt", &all_cameras);
  MvsModel model;
  // model.SetCameras(all_cameras, {{150, 310}}, dirpath, 0);
  vector<int> camera_indexes;
  for (int i = 1; i < 10; i += 1) {
    camera_indexes.push_back(i * 5);
  }
  model.SetCameras(all_cameras, camera_indexes, 0);
  model.set_ref_depth(1.5);
  model.SetRefPoint(cv::Point2f(model.camera(0)->internals().p[0],
                                model.camera(0)->internals().p[1]));
  model.SetDepthRange(0.1, 5.1);
  model.SetDepthStep(FLAGS_depth_step);

  if (!FLAGS_in_cv.empty()) {
    CostVolume cost_volume;
    // ReadImsCostVolume(FLAGS_in_cost_volume, &cost_volume);
    CHECK(cost_volume.ReadFile(FLAGS_in_cv));
    model.SetCostVolume(move(cost_volume));
  }

  // DepthMap depth_map;
  // CostVolume cost_volume;
  // model.GenerateCostVolume(&cost_volume, &depth_map);
  // cost_volume.WriteFile("mvs.cv");
  // depth_map.WriteDepthImage("mvs.jpg");

  // return 0;

  QApplication app(argc, argv);
  // furry::CvWindow patch_window("patch");
  // furry::CvWindow ref_window("Reference");
  // furry::CvWindow image_window("Images");
  // ref_window.updateImage(ref_image);
  // // patch_window.updateImage(image2);
  // image_window.updateImage(GridImages(images));
  // ref_window.setMouseCallback(
  //     [&](int event, int x, int y, int, void*) {
  //       if( event != EVENT_LBUTTONDOWN )
  //         return;
  //       vector<cv::Mat> patches(images.size());
  //       vector<cv::Mat> labeled_images(images.size());
  //       for (size_t i = 0; i < num_cameras; ++i) {
  //         images[i].copyTo(labeled_images[i]);
  //         cv::circle(labeled_images[i], Point(x, y), 10, CV_RGB(255, 0, 0), -1);
  //         patches[i] = ExtractPatch<Vec3b>(images[i], homographies[i],
  //                                   Point2f(x, y), Size(61, 61), 1);
  //       }
  //       image_window.updateImage(GridImages(labeled_images));
  //       patch_window.updateImage(GridImages(patches));
  //     });

  // // cv_window.displayStatusBar(QString("hello"), 0);
  MvsController controller(&model);
  controller.show();
  app.exec();

  // const char* win_name = "test";
  // cv::namedWindow(win_name, CV_GUI_EXPANDED | CV_WINDOW_NORMAL);
  // cv::imshow(win_name, image);
  // cv::waitKey();

  return 0;
}

const char *kUsageMessage = "Read bundle adjustment result and calculate the parallax information";

#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include "furry/common/cv.h"
#include "furry/common/init.h"

#include "model.h"

using namespace Eigen;
using namespace furry;
using namespace furry::tiny;
using namespace std;

DEFINE_string(in_model, "", "");
DEFINE_string(out_file, "", "");

int main (int argc, char *argv[]) {
  furry::Init(&argc, &argv, kUsageMessage);

  Model model;
  CHECK(model.ReadFile(FLAGS_in_model));

  vector<Eigen::Matrix3d> rotation_homographies(model.NumCameras());
  for (int i = 0; i < model.NumCameras(); ++i) {
    rotation_homographies[i] =
        model.GetCamera(i)->RotationHomography().inverse();
  }

  vector<vector<double>> parallax(model.NumTracks());
  for (int t = 0; t < model.NumTracks(); ++t) {
    parallax[t].resize(model.NumCameras() - 1);
    Vector2d p = To<Vector2d>(model.GetDetection(model.GetTrackId(t),
                                                 model.GetCameraId(0)));
    Vector3d p3 = model.GetTrack(t)->GetPoint();
    Vector3d center = model.GetCamera(0)->GetCenter();
    Vector3d d = (center - p3).normalized();
    for (int c = 1; c < model.NumCameras(); ++c) {
      // Vector2d pp = To<Vector2d>(model.GetDetection(model.GetTrackId(t),
      //                                               model.GetCameraId(c)));
      // pp = Vector3d(rotation_homographies[c] * pp.homogeneous()).hnormalized();
      center = model.GetCamera(c)->GetCenter();
      Vector3d dp = (center - p3).normalized();
      // parallax[t][c-1] = (p - pp).norm();
      parallax[t][c-1] = acos(d.dot(dp));
    }
  }

  ofstream ofs(FLAGS_out_file);
  // ofs << model.NumCameras() << ' ' << model.NumTracks() << '\n';
  for (int t = 0; t < model.NumTracks(); ++t) {
    for (int c = 0; c < model.NumCameras() - 1; ++c) {
      ofs << parallax[t][c] << ' ';
    }
    ofs << '\n';
  }
  return 0;
}

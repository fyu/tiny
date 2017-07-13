#include <iostream>
#include <unordered_map>
#include <random>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include "furry/common/init.h"

#include "numeric.h"
#include "model.h"
#include "tiny_bundle.h"
#include "tiny_bundle_one.h"
#include "tiny_bundle_util.h"

using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_model, "", "");
DEFINE_string(settings, "", "");
DEFINE_double(p_x, 0,
              "The projection of the 3d point in the reference "
              "view in retina coordinate");
DEFINE_double(p_y, 0, "");
DEFINE_double(min_depth, 1, "");
DEFINE_double(max_depth, 10, "");
DEFINE_int32(num_depth_samples, 100, "");
DEFINE_int32(min_cameras, 2, "");
DEFINE_int32(max_cameras, -1, "");
DEFINE_int32(camera_step, 1, "");
DEFINE_bool(change_coord, true, "");

bool ReadFileToProto(const string& filename, google::protobuf::Message* proto) {
  int file = open(filename.c_str(), O_RDONLY);
  google::protobuf::io::FileInputStream file_stream(file);
  return google::protobuf::TextFormat::Parse(&file_stream, proto);
}

void SampleDepth(const Point2 &projection,
                 double min_depth,
                 double max_depth,
                 int num_samples,
                 Model *model) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(0, 0.5);

  double step = (max_depth - min_depth) / num_samples;
  for (int i = 0; i < num_samples; ++i) {
    double depth = i * step + min_depth;
    auto track = new Track;
    track->SetPoint(projection.x * depth, projection.y * depth, depth);
    track->SetColor(255, 0, 0);
    vector<Point2> detections;
    for (int j = 0; j < model->NumCameras(); ++j) {
      auto camera = model->GetCamera(j);
      track->AddCameraId(camera->GetId());
      Vector2d p2 = camera->ProjectPoint(track->GetPoint()).hnormalized();
      detections.push_back(Point2(p2.x() + d(gen), p2.y() + d(gen)));
      // 2detections.push_back(Point2(p2.x(), p2.y()));
    }
    model->AddTrack(track, detections);
  }
}

void SynCameras(int num_cameras, Model *model) {
  auto camera = new Camera;
  double f[2] = {1781, 1781};
  double k[2] = {0, 0};
  camera->SetInternals(f, k);
  camera->SetCenter(0, 0, 0);
  camera->SetEulerAngles(0, 0, 0);
  model->AddCamera(camera);
  for (int i = 1; i < num_cameras; ++i) {
    camera = new Camera;
    camera->SetCenter(0.005, 0, 0);
    camera->SetEulerAngles(0, 0.01, 0);
    model->AddCamera(camera);
  }
}

int main (int argc, char *argv[]) {
  furry::Init(&argc, &argv, "");
  Model model;

  CHECK(!FLAGS_in_model.empty());
  CHECK(!FLAGS_settings.empty());

  model.ReadFile(FLAGS_in_model);
  model.RemoveTracks();
  // SynCameras(100, &model);

  BundleSettings settings;

  CHECK(ReadFileToProto(FLAGS_settings, &settings))
      << "Failed to read BA settings file " << FLAGS_settings;

  CHECK(settings.fix_internals());
  CHECK(settings.fix_externals());
  CHECK(!settings.fix_tracks());

  if (FLAGS_max_cameras < 0) FLAGS_max_cameras = model.NumCameras();

  if (FLAGS_change_coord) {
    ChangeCoord(model.GetCamera(settings.ref_camera_index())->GetExternals(),
                &model);
  }

  SampleDepth(Point2(FLAGS_p_x, FLAGS_p_y),
              FLAGS_min_depth,
              FLAGS_max_depth,
              FLAGS_num_depth_samples,
              &model);

  for (int num_cameras = FLAGS_max_cameras;
       num_cameras >= FLAGS_min_cameras;
       num_cameras -= FLAGS_camera_step) {
    model.RemoveCameras(num_cameras);

    BundleOne(settings, nullptr, &model);

    // for (int i = 0; i < model.NumTracks(); ++i) {
    for (int i = 0; i < 1; ++i) {
      auto track = model.GetTrack(i);
      printf("%d %lf %lf;\n",
             num_cameras,
             track->GetPoint().z(),
             sqrt(track->GetCovariance()[0]));
    }
  }

  return 0;
}

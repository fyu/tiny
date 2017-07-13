#include <iostream>
#include <limits>
#include <unordered_map>
#include <fstream>

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
DEFINE_string(in_image_dir, "", "");
DEFINE_string(in_image_pattern, "", "");
DEFINE_string(out_model, "", "");
DEFINE_string(out_uncertainty, "", "");
DEFINE_string(settings, "", "");
DEFINE_bool(set_c_zero, false, "");
DEFINE_bool(set_p_random, false, "");
DEFINE_bool(set_p_same, false, "");
DEFINE_double(projection_error_threshold, -1, "");
DEFINE_bool(one, true, "");
DEFINE_int32(num_cameras, -1, "");
DEFINE_bool(change_coord, true, "");
DEFINE_bool(triangulate, false, "");
DEFINE_bool(adjust, true, "");
DEFINE_bool(use_trim, true, "");
DEFINE_bool(use_separate, false, "");
DEFINE_string(camera_order, "", "");

bool ReadFileToProto(const string& filename, google::protobuf::Message* proto) {
  int file = open(filename.c_str(), O_RDONLY);
  google::protobuf::io::FileInputStream file_stream(file);
  return google::protobuf::TextFormat::Parse(&file_stream, proto);
}

void RandomInitPointsInViewDirection(Model *model) {
  for (int i = 0; i < model->NumTracks(); ++i) {
    Eigen::Vector2d mean_detection(0, 0);
    auto track = model->GetTrack(i);
    for (int j = 0; j < track->NumCameras(); ++j) {
      auto camera = model->GetCameraById(track->GetCameraId(j));
      auto K = camera->K();
      auto d = camera->GetDetectionOfTrack(track->GetId());
      Eigen::Vector3d detection(d.x, d.y, 1);
      Eigen::Vector3d norm_coord = K.colPivHouseholderQr().solve(detection);
      mean_detection += norm_coord.hnormalized();
    }
    mean_detection /= track->NumCameras();
    double rand_depth = uniform_rand() * 2 + 2;
    // double rand_depth = 1;
    double p[3];
    p[0] = mean_detection[0] * rand_depth;
    p[1] = mean_detection[1] * rand_depth;
    p[2] = rand_depth;
    // y = -1;
    // point.SetFrom3DPoint(point.type(), x, y, z);
    track->SetPoint(p);
  }
}

void RandomInitPointsInReferenceViewDirection(Model *model,
                                              int ref_camera_index) {
  auto ref_camera = model->GetCamera(ref_camera_index);
  for (int i = 0; i < model->NumTracks(); ++i) {
    // Eigen::Vector2d mean_detection(0, 0);
    auto track = model->GetTrack(i);
    // auto d = ref_camera->GetDetectionOfTrack(i);
    // Eigen::Vector3d detection(d.x, d.y, 1);
    // Eigen::Vector3d norm_coord = K.colPivHouseholderQr().solve(detection);
    // mean_detection = norm_coord.hnormalized();
    auto detection = ref_camera->GetDetectionOfTrackInRetinaCoord(track->GetId());
    // cout << mean_detection.transpose() << '\n';
    double rand_depth = uniform_rand() * 2 + 2;
    // double rand_depth = 1;
    double p[3];
    p[0] = detection.x * rand_depth;
    p[1] = detection.y * rand_depth;
    p[2] = rand_depth;
    // y = -1;
    // point.SetFrom3DPoint(point.type(), x, y, z);
    track->SetPoint(p);
  }
}

void SetPointsSameDepth(int ref_camera_index, double depth, Model *model) {
  auto ref_camera = model->GetCamera(ref_camera_index);
  for (int i = 0; i < model->NumTracks(); ++i) {
    auto track = model->GetTrack(i);
    auto detection = ref_camera->GetDetectionOfTrackInRetinaCoord(track->GetId());
    double p[3];
    p[0] = detection.x * depth;
    p[1] = detection.y * depth;
    p[2] = depth;
    track->SetPoint(p);
  }
}

void SetCameraZero(Model *model) {
  for (int i = 0; i < model->NumCameras(); ++i) {
    auto camera = model->GetCamera(i);
    camera->SetCenter(0, 0, 0);
    camera->SetQuaternion(1, 0, 0, 0);
  }
}

bool RemoveCenterTracks(int ref_camera_index, double threshold, Model *model) {
  vector<track_id_t> tracks_to_remove;
  auto camera = model->GetCamera(ref_camera_index);
  for (int i = 0 ; i < camera->NumDetections(); ++i) {
    auto d = camera->GetDetectionInRetinaCoord(i);
    if (fabs(d.x) < threshold || fabs(d.y) < threshold) {
      tracks_to_remove.push_back(camera->GetTrackId(i));
    }
  }
  model->RemoveTracksById(tracks_to_remove);
  return tracks_to_remove.size() > 0;
}

bool RemoveTracksInRange(double min, double max, Model *model) {
  vector<track_id_t> track_to_remove;
  for (int i = 0; i < model->NumTracks(); ++i) {
    auto track = model->GetTrack(i);
    if (track->GetPoint().z() >= min &&
        track->GetPoint().z() <= max) {
      track_to_remove.push_back(track->GetId());
    }
  }
  for (auto id : track_to_remove) {
    model->RemoveTrackById(id);
  }
  LOG(INFO) << "Removed " << track_to_remove.size() << " tracks";
  return track_to_remove.size() > 0;
}

bool RemoveTracksInPercentageRange(double pmin, double pmax, Model *model) {
  vector<double> depths(model->NumTracks());
  for (int i = 0; i < model->NumTracks(); ++i) {
    depths[i] = model->GetTrack(i)->GetPoint().z();
  }
  sort(depths.begin(), depths.end());
  double min, max;
  min = depths[pmin * depths.size()];
  max = depths[pmax * depths.size() - 1];
  return RemoveTracksInRange(min, max, model);
}

bool RemoveNegativeTracks(Model *model) {
  return RemoveTracksInRange(numeric_limits<double>::lowest(), 0, model);
}

bool RemoveOutlierTracks(double threshold, Model *model) {
  vector<track_id_t> tracks_to_remove;
  for (int i_track = 0; i_track < model->NumTracks(); ++i_track) {
    auto track = model->GetTrack(i_track);
    for (int i_camera = 0; i_camera < track->NumCameras(); ++i_camera) {
      auto camera_id = track->GetCameraId(i_camera);
      if (model->GetProjectionError(track->GetId(), camera_id) > threshold) {
        tracks_to_remove.push_back(track->GetId());
        break;
      }
    }
  }
  for (auto id : tracks_to_remove) {
    model->RemoveTrackById(id);
  }
  LOG(INFO) << "Removed " << tracks_to_remove.size() << " tracks";
  return tracks_to_remove.size() > 0;
}

void BundleSeparate(BundleSettings settings,
                    Model *model) {
  settings.set_zero_translation(true);
  settings.set_use_robustifier(true);
  settings.set_fix_tracks(true);
  BundleOne(settings, nullptr, model);
  settings.set_zero_translation(false);
  settings.set_zero_rotation(true);
  settings.set_fix_tracks(false);
  // RandomInitPointsInReferenceViewDirection(
  //     model, settings.ref_camera_index());
  settings.set_use_robustifier(false);
  BundleOne(settings, nullptr, model);
  settings.set_zero_rotation(false);
  BundleOne(settings, nullptr, model);
}

void ReorderCameras(const string &order_file, Model *model) {
  vector<int> order(model->NumCameras());
  ifstream ifs(order_file);
  for (size_t i = 0; i < order.size(); ++i) {
    ifs >> order[i];
  }
  model->ReorderCameras(order);
}

int main (int argc, char *argv[]) {
  furry::Init(&argc, &argv, "");
  CHECK(!FLAGS_in_model.empty());
  CHECK(!(FLAGS_set_p_random && FLAGS_set_p_same));
  Model model;

  model.ReadFile(FLAGS_in_model);

  if (!FLAGS_camera_order.empty()) {
    ReorderCameras(FLAGS_camera_order, &model);
  }
  // For 0008
  // model.RemoveCameras(90, 97);
  // model.RemoveCameras(78, 81);
  // RemoveTracksInRange(2, numeric_limits<double>::max(), &model);
  // RemoveTracksInRange(0, 2, &model);
  // RemoveTracksInPercentageRange(0.5, 1, &model);

  BundleSettings settings;
  CHECK(ReadFileToProto(FLAGS_settings, &settings))
      << "Failed to read BA settings file" << FLAGS_settings;

  // RemoveCenterTracks(settings.ref_camera_index(), 0.1, &model);

  if (FLAGS_num_cameras > 0 && FLAGS_num_cameras < model.NumCameras()) {
    model.RemoveCameras(FLAGS_num_cameras, model.NumCameras());
  }

  if (FLAGS_change_coord) {
    ChangeCoord(model.GetCamera(settings.ref_camera_index())->GetExternals(),
                &model);
  }

  if (FLAGS_triangulate) {
    for (int i = 0; i < model.NumTracks(); ++i) {
      model.TriangulateTrack(i);
    }
  }

  if (FLAGS_set_c_zero) {
    SetCameraZero(&model);
  }
  if (FLAGS_set_p_random) {
    if (FLAGS_one) {
      RandomInitPointsInReferenceViewDirection(&model,
                                               settings.ref_camera_index());
    } else {
      RandomInitPointsInViewDirection(&model);
    }
  } else if (FLAGS_set_p_same) {
    SetPointsSameDepth(settings.ref_camera_index(), 3.0f, &model);
  } else {
    if (FLAGS_projection_error_threshold > 0) {
      // ProgressBar removing_outliers("Removing Outliers", model.NumTracks());
      for (int i = 0; i < model.NumTracks(); ++i) {
        auto track = model.GetTrack(i);
        auto track_id = track->GetId();
        vector<camera_id_t> camera_ids_to_remove;
        for (int j = 0; j < track->NumCameras(); ++j) {
          auto camera_id = track->GetCameraId(j);
          if (model.GetProjectionError(track_id, camera_id) >
              FLAGS_projection_error_threshold) {
            // model.RemoveDetection(track_id, camera_id);
            camera_ids_to_remove.push_back(camera_id);
          }
        }
        for (auto id : camera_ids_to_remove) {
          model.RemoveDetection(track_id, id);
        }
      }
      model.RemoveEmptyTracks();
    }
  }

  vector<TrackPair> track_neighbors;
  if (settings.with_smooth_depth()) {
    if (FLAGS_one) {
      vector<Detection> detections;
      GetDetections(*model.GetCamera(settings.ref_camera_index()), &detections);
      ConnectTracksBy2dDistance(detections,
                                settings.smooth_depth().num_neighbors(),
                                settings.smooth_depth().range_squared(),
                                &track_neighbors);
    } else {
      unordered_map<camera_id_t, vector<Detection> > image_detections;
      GetDetections(model, &image_detections);
      ConnectTracksBy2dDistance(image_detections,
                                settings.smooth_depth().num_neighbors(),
                                settings.smooth_depth().range_squared(),
                                &track_neighbors);
    }
  }

  if (FLAGS_adjust) {
    if (FLAGS_one) {
      if (FLAGS_use_separate) {
        BundleSeparate(settings, &model);
      } else {
        BundleOne(settings, &track_neighbors, &model);
      }
    } else {
      Bundle(settings, &track_neighbors, &model);
    }
    int num_iterations = 1;
    while (FLAGS_use_trim &&
           (RemoveNegativeTracks(&model) ||
            (FLAGS_projection_error_threshold > 0 &&
             RemoveOutlierTracks(FLAGS_projection_error_threshold * pow(0.9, num_iterations), &model)))) {
      if (FLAGS_one) {
        if (FLAGS_use_separate) {
          BundleSeparate(settings, &model);
        } else {
          BundleOne(settings, &track_neighbors, &model);
        }
      } else {
        Bundle(settings, &track_neighbors, &model);
      }
      ++num_iterations;
      // if(num_iterations > 6) break;
    }
  }

  LOG(INFO) << "# Tracks: " << model.NumTracks();

  model.RelabelTracks();
  model.WriteFile(FLAGS_out_model, settings.calc_covariance());
  return 0;
}

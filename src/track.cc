#include "track.h"

using namespace std;

namespace furry {
namespace tiny {

Detection::Detection() {}
Detection::Detection(float x, float y) : Vector2f(x, y) {
}

void Detection::set_track_id(track_id_t track_id) {
  track_id_ = track_id;
}

track_id_t Detection::track_id() const {
  return track_id_;
}

float Distance(const Detection &d0, const Detection &d1) {
  return (d0 - d1).norm();
}

void GetDetections(
    const Model &model,
    unordered_map<camera_id_t, vector<Detection> > *image_detections) {
  for (int i = 0; i < model.NumCameras(); ++i) {
    vector<Detection> detections;
    auto camera = model.GetCamera(i);
    for (int j = 0; j < camera->NumDetections(); ++j) {
      Detection detection(camera->GetDetection(j).x, camera->GetDetection(j).y);
      detection.set_track_id(camera->GetTrackId(j));
      detections.push_back(detection);
    }
    image_detections->insert(make_pair(i, move(detections)));
  }
}

void GetDetections(const Camera &camera, vector<Detection> *detections) {
  for (int j = 0; j < camera.NumDetections(); ++j) {
    Detection detection(camera.GetDetection(j).x, camera.GetDetection(j).y);
    detection.set_track_id(camera.GetTrackId(j));
    detections->push_back(detection);
  }
}

} // tiny
} // furry

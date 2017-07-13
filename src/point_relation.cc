#include "point_relation.h"

#include "furry/3rdparty/kdtree/kdtree.h"

using namespace std;

namespace furry {
namespace tiny {

TrackPair::TrackPair(track_id_t track0, track_id_t track1, double score) {
  if (track0 > track1) {
    swap(track0, track1);
  }
  tracks_[0] = track0;
  tracks_[1] = track1;
  score_ = score;
}

track_id_t TrackPair::track(int i) const {
  return tracks_[i];
}

void TrackPair::set_score(double score) {
  score_ = score;
}

double TrackPair::score() const {
  return score_;
}

bool TrackPair::operator < (const TrackPair &p) const {
  return tracks_[0] < p.tracks_[0]
      || (tracks_[0] == p.tracks_[0] && tracks_[1] < p.tracks_[1]);
}

void ConnectTracksBy2dDistance(const vector<Detection> &detections,
                               int num_neighbors, float range_sq,
                               vector<TrackPair> *track_pairs) {
  int num_points = detections.size();
  vector<float> points(num_points * 2);
  for (size_t i = 0; i < num_points; ++i) {
    points[i << 1] = detections[i][0];
    points[(i << 1) + 1] = detections[i][1];
  }
  kdtree::kdtree<2, float> tree;
  tree.initialize(points.data(), num_points, false);
  vector<int> heap(num_neighbors + 1);
  Detection d0, d1;
  float dist;
  for (size_t i = 0; i < num_points; ++i) {
    d0 = detections[i];
    int num_found = tree.nearestNeighbors(d0.data(), heap.data(),
                                          num_neighbors + 1, range_sq);
    for (int j = 0; j < num_found; ++j) {
      d1 = detections[heap[j]];
      dist = Distance(d0, d1);
      if (d0.track_id() == d1.track_id() || dist < 1e-20) continue;
      track_pairs->push_back(TrackPair(d0.track_id(), d1.track_id(),
                                       // 1.0 / dist));
                                       dist));
    }
  }
}

void ConnectTracksBy2dDistance(
    const unordered_map<camera_id_t, vector<Detection>> &image_detections,
    int num_neighbors, float range_sq,
    vector<TrackPair> *track_pairs) {
  for (auto &detections : image_detections) {
    ConnectTracksBy2dDistance(detections.second, num_neighbors, range_sq,
                              track_pairs);
  }
}

} // tiny
} // furry

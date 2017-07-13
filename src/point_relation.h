#ifndef POINT_RELATION_H
#define POINT_RELATION_H

#include <vector>
#include <cstdint>
#include <unordered_map>

#include "track.h"

namespace furry {
namespace tiny {

class TrackPair {
 public:
  TrackPair(track_id_t track0, track_id_t track1, double score);
  track_id_t track(int i) const;
  void set_score(double score);
  double score() const;

  bool operator < (const TrackPair &p) const;

 private:
  track_id_t tracks_[2];
  double score_;
};


void ConnectTracksBy2dDistance(const std::vector<Detection> &detections,
                               int num_neighbors, float range_sq,
                               std::vector<TrackPair> *track_pairs);

void ConnectTracksBy2dDistance(
    const std::unordered_map<camera_id_t, std::vector<Detection>>
    &image_detections,
    int num_neighbors, float range_sq,
    std::vector<TrackPair> *track_pairs);

} // tiny
} // furry

static inline void mix(uint64_t& a, uint64_t& b, uint64_t& c) {     // 64bit version
  a -= b; a -= c; a ^= (c>>43);
  b -= c; b -= a; b ^= (a<<9);
  c -= a; c -= b; c ^= (b>>8);
  a -= b; a -= c; a ^= (c>>38);
  b -= c; b -= a; b ^= (a<<23);
  c -= a; c -= b; c ^= (b>>5);
  a -= b; a -= c; a ^= (c>>35);
  b -= c; b -= a; b ^= (a<<49);
  c -= a; c -= b; c ^= (b>>11);
  a -= b; a -= c; a ^= (c>>12);
  b -= c; b -= a; b ^= (a<<18);
  c -= a; c -= b; c ^= (b>>22);
}

inline uint64_t Hash64NumWithSeed(uint64_t num, uint64_t c) {
  uint64_t b = 0xe08c1d668b756f82ULL;   // more of the golden ratio
  mix(num, b, c);
  return c;
}
namespace std {
template <>
struct hash<furry::tiny::TrackPair> {
  size_t operator () (const furry::tiny::TrackPair &p) const {
    return Hash64NumWithSeed(p.track(0), p.track(1));
  }
};
} // std

#endif // POINT_RELATION_H

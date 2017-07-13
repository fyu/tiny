#ifndef TRACK_H
#define TRACK_H

#include <unordered_map>

#include <Eigen/Dense>

#include "model.h"

using Eigen::Vector2f;

namespace furry {
namespace tiny {

class Detection : public Eigen::Vector2f {
 public:
  Detection();
  Detection(float x, float y);
  void set_track_id(track_id_t track_id);
  track_id_t track_id() const;
 private:
  track_id_t track_id_;
};

float Distance(const Detection &d0, const Detection &d1);

void GetDetections(const Camera &camera, std::vector<Detection> *detections);

void GetDetections(
    const Model &model,
    std::unordered_map<camera_id_t, std::vector<Detection> > *image_detections);

} // tiny
} // furry

#endif

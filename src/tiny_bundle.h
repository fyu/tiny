#ifndef TINY_BUNDLE_H
#define TINY_BUNDLE_H

#include <vector>

#include "model.h"
#include "ba_settings.pb.h"
#include "point_relation.h"

namespace furry {
namespace tiny {

// void ChangeCoord(const lightfield_sfm::CameraExternals& externals,
//                  lightfield_sfm::Model* model);

bool Bundle(const tiny::BundleSettings& settings,
            const std::vector<TrackPair> *track_neighbors,
            Model* model);

} // tiny
} // furry

#endif // TINY_BUNDLE_H

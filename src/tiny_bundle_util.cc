#include "tiny_bundle_util.h"

namespace furry {
namespace tiny {

void ChangeCoord(const Camera::Externals &externals, Model *model) {
  auto R = externals.R();
  auto C = externals.C();
  Eigen::Matrix3d RT = R.transpose();

  for (int i = 0; i < model->NumCameras(); ++i) {
    auto camera = model->GetCamera(i);
    camera->SetRotationMatrix(camera->R() * RT);
    camera->SetCenter(R * (camera->C() - C));
  }

  for (int i = 0; i < model->NumTracks(); ++i) {
    auto track = model->GetTrack(i);
    track->SetPoint(R * (track->GetPoint() - C));
  }
}

} // tiny
} // furry

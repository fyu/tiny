#ifndef TINY_UTIL_H
#define TINY_UTIL_H

#include <string>

#include "geo/lightfield/sfm/model.pb.h"
#include "geo/lightfield/sfm/model.h"

namespace tiny {

struct CameraTracks {
  lightfield_sfm::docid_hash_t docid;
  lightfield_sfm::Camera camera;
  vector<lightfield_sfm::TrackMapValue> tracks;
};

struct CameraTracksMap : hash_map<lightfield_sfm::docid_hash_t, CameraTracks> {};

void BuildCameraTracksMap(const lightfield_sfm::Model& model,
                          CameraTracksMap* ct_map);


string GetFilename(const string& path);

string GetDir(const string& path);

void ExportCameras(const lightfield_sfm::Model& model, const string& filename);

void ExportTracks(const lightfield_sfm::Model& model, const string& filename);

bool ExportModel(const lightfield_sfm::Model &model, const string &prefix);

void InitGlobalDocidToFilepath(const lightfield_sfm::Model& model);

void RemoveCameras(const string &camera_pattern, int track_size_threshold,
                   lightfield_sfm::Model* model);

}

#endif

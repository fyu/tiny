#include "util.h"

#include <vector>
#include <utility>
#include <cstdlib>

#include "util/hash/hash.h"
#include "util/regexp/re2/re2.h"

using namespace lightfield_sfm;

namespace tiny {

void BuildCameraTracksMap(const Model& model,
                         CameraTracksMap* ct_map) {
  auto& tracks = model.GetTracks();
  auto docids = model.GetAllDocIDHashesSet();
  for (auto docid : docids) {
    (*ct_map)[docid].camera = model.GetCamera(docid);
    (*ct_map)[docid].docid = docid;
    vector<TrackMapValue> visible_tracks;
    for (auto track : tracks) {
      if (track.second.track.HasDocID(docid)) {
        // (*ct_map)[docid].tracks.push_back(track.second);
        visible_tracks.push_back(track.second);
      }
    }
    (*ct_map)[docid].tracks = move(visible_tracks);
    LOG(INFO) << docid << " have " << (*ct_map)[docid].tracks.size()
              << " tracks\n";
  }
}

double uniform_rand() {
  return (double)std::rand() / RAND_MAX;
}

string GetFilename(const string& path) {
  auto n = path.rfind("/");
  if (n == string::npos) {
    return path;
  } else {
    return path.substr(n + 1, string::npos);
  }
}

string GetDir(const string& path) {
  auto n = path.rfind("/");
  if (n == string::npos) {
    return "";
  } else {
    return path.substr(0, n);
  }
}

bool CompNamedCamera(const pair<string, Camera>& c1,
                     const pair<string, Camera>& c2) {
  return c1.first < c2.first;
}

void ExportCameras(const Model& model, const string& filename) {
  vector<pair<string,Camera>> cameras;
  auto doc_ids = model.GetAllDocIDHashesSet();
  string dir = GetDir(model.GetCamera(*doc_ids.begin()).GetDebugUrl());
  for (auto docid : doc_ids) {
    auto camera = model.GetCamera(docid);
    auto path = camera.GetDebugUrl();
    cameras.push_back(make_pair(GetFilename(path), camera));
  }
  sort(cameras.begin(), cameras.end(), CompNamedCamera);
  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "%d %s\n", model.NumCameras(), dir.c_str());
  for (auto& cam_pair : cameras) {
    auto& cam = cam_pair.second;
    auto i = cam.internals();
    auto e = cam.externals();
    // Write internal parameters
    fprintf(file, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f ",
            i.k[0], i.k[1], i.f[0], i.f[1], i.p[0], i.p[1], i.skew);
    // Write external parameters
    fprintf(file, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f ",
            e.c[0], e.c[1], e.c[2], e.q[0], e.q[1], e.q[2], e.q[3]);
    fprintf(file, "%s\n", cam_pair.first.c_str());
  }
  fclose(file);
}

void ExportTracks(const Model& model, const string& filename) {
  auto tracks = model.GetTracks();
  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "%lu\n", tracks.size());
  // export the track 3d positions
  for (auto track : tracks) {
    Point& p = track.second.point;
    fprintf(file, "%f %f %f\n", p[0], p[1], p[2]);
  }
  // export track colors (only RGB)
  int r, g, b;
  for (auto track : tracks) {
    auto argb = track.second.track.Color();
    r = (argb >> 16) & 0xFF;
    g = (argb >>  8) & 0xFF;
    b = (argb >>  0) & 0xFF;
    fprintf(file, "%d %d %d\n", r, g, b);
  }
  fclose(file);
}

void InitGlobalDocidToFilepath(const lightfield_sfm::Model& model) {
  auto& cameras = model.GetCameras();
  for (auto& camera : cameras) {
    global_docid_hash_to_filepath[camera.first] =
        camera.second.camera.GetDebugUrl();
  }
}

void RemoveCameras(const string &camera_pattern, int track_size_threshold,
                   Model* model) {
  auto doc_ids = model->GetAllDocIDHashesSet();
  // RE2 pattern("long_.*.png");
  // RE2 pattern("stable_00[0-2][0-9].png");
  RE2 pattern(camera_pattern);
  for (auto& id : doc_ids) {
    auto &camera = model->GetCamera(id);
    auto &url = camera.GetDebugUrl();
    if (!RE2::PartialMatch(url, pattern)) {
      model->EraseCamera(id);
    }
  }
  auto& track_map = model->GetTracks();
  vector<docid_hash_t> track_ids_to_delete;
  for (auto it_track = track_map.begin();
       it_track != track_map.end();
       ++it_track) {
    if (it_track->second.track.size() < track_size_threshold) {
      track_ids_to_delete.push_back(it_track->first);
    }
  }
  for (auto id : track_ids_to_delete) {
    model->EraseTrack(id);
  }
}

struct RichCamera {
  string filename;
  docid_hash_t docid;
  Camera camera;
};

bool CompRichCamera(const RichCamera &c0, const RichCamera &c1) {
  return c0.filename < c1.filename;
}

bool ExportModel(const Model &model, const string &prefix) {
  vector<RichCamera> cameras;
  hash_map<docid_hash_t, int> camera_index;
  hash_map<track_id_t, int> track_index;
  auto doc_ids = model.GetAllDocIDHashesSet();
  auto& tracks = model.GetTracks();
  string dir = GetDir(model.GetCamera(*doc_ids.begin()).GetDebugUrl());
  for (auto docid : doc_ids) {
    auto camera = model.GetCamera(docid);
    auto path = camera.GetDebugUrl();
    cameras.push_back({GetFilename(path), docid, camera});
  }
  sort(cameras.begin(), cameras.end(), CompRichCamera);
  for (size_t i = 0; i < cameras.size(); ++i) {
    camera_index[cameras[i].docid] = i;
  }
  int track_count = 0;
  for (auto track : tracks) {
    track_index[track.first] = track_count++;
  }
  string filename = prefix + "_cameras.txt";
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL) {
    LOG(ERROR) << "Can't open " << filename;
    return false;
  }
  fprintf(file, "%d %s\n", model.NumCameras(), dir.c_str());
  for (auto& camera : cameras) {
    auto& cam = camera.camera;
    auto i = cam.internals();
    auto e = cam.externals();
    // Write internal parameters
    fprintf(file, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f ",
            i.k[0], i.k[1], i.f[0], i.f[1], i.p[0], i.p[1], i.skew);
    // Write external parameters
    fprintf(file, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f ",
            e.c[0], e.c[1], e.c[2], e.q[0], e.q[1], e.q[2], e.q[3]);
    fprintf(file, "%s\n", camera.filename.c_str());
  }

  for (auto& camera : cameras) {
    vector<Detection> detections;
    vector<track_id_t> track_ids;
    for (auto track : tracks) {
      if (track.second.track.HasDocID(camera.docid)) {
        detections.push_back(
            track.second.track.GetItem(camera.docid).detection);
        track_ids.push_back(track_index[track.first]);
      }
    }
    fprintf(file, "%lu ", detections.size());
    for (size_t i = 0; i < detections.size(); ++i) {
      fprintf(file, "%f %f ", detections[i].x, detections[i].y);
    }
    for (size_t i = 0; i < track_ids.size(); ++i) {
      fprintf(file, "%llu ", track_ids[i]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  // Write Tracks
  filename = prefix + "_tracks.txt";
  file = fopen(filename.c_str(), "w");
  if (file == NULL) {
    LOG(ERROR) << "Can't open " << filename;
    return false;
  }
  fprintf(file, "%lu\n", tracks.size());
  // export the track 3d positions
  for (auto track : tracks) {
    Point& p = track.second.point;
    fprintf(file, "%f %f %f\n", p[0], p[1], p[2]);
  }
  // export track colors (only RGB)
  int r, g, b;
  for (auto track : tracks) {
    auto argb = track.second.track.Color();
    r = (argb >> 16) & 0xFF;
    g = (argb >>  8) & 0xFF;
    b = (argb >>  0) & 0xFF;
    fprintf(file, "%d %d %d\n", r, g, b);
  }

  for (auto track : tracks) {
    fprintf(file, "%d ", track.second.track.size());
    for (auto it = track.second.track.begin();
         it != track.second.track.end(); ++it) {
      fprintf(file, "%d ", camera_index[it->docid_hash]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  return true;
}

} // tiny

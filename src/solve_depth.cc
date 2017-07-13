#include <cstdio>
#include <string>
#include <fstream>
#include <limits>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <Eigen/Dense>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "furry/3rdparty/mrf/GCoptimization.h"
#include "furry/3rdparty/mrf/TRW-S.h"
#include "furry/3rdparty/mrf/MaxProdBP.h"
#include "furry/3rdparty/mrf/ICM.h"
#include "furry/common/init.h"
#include "furry/common/path.h"
#include "furry/common/protobuf.h"
#include "furry/common/str.h"

#include "cost_volume.h"
#include "depth_map.h"
#include "mvs_opt.h"
#include "smooth.h"

using namespace Eigen;
using namespace std;
using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_cv, "", "");
DEFINE_string(settings, "", "");
DEFINE_bool(use_cont, false, "");
DEFINE_bool(use_disc, false, "use graph cut");
DEFINE_bool(use_densecrf, false, "use dense crf");
DEFINE_string(out_dir, "", "");
DEFINE_string(prefix, "", "");
DEFINE_string(smooth_settings, "", "");

int main(int argc, char *argv[]) {
  furry::Init(&argc, &argv, "");
  CostVolume cost_volume;
  DepthMap depthmap;
  CHECK(cost_volume.ReadFile(FLAGS_in_cv));
  CHECK(!FLAGS_out_dir.empty());
  CHECK(!FLAGS_prefix.empty());

  string prefix =
      StringPrintf("%s/%s/", FLAGS_out_dir.c_str(), FLAGS_prefix.c_str());
  MakeDir(prefix);
  prefix += FLAGS_prefix;

  auto image = cost_volume.GetImage();
  if (!FLAGS_smooth_settings.empty()) {
    cv::Mat smooth_image;
    SmoothSettings settings;
    ReadTxtFileToProto(FLAGS_smooth_settings, &settings);
    SmoothImage(settings, cost_volume.GetImage(), &smooth_image);
    imwrite(prefix + "_smooth.png", smooth_image);
    WriteProtoToTxtFile(prefix + "_smooth_settings.txt", settings);
    cost_volume.SetImage(smooth_image);
  }

  // CHECK(!cost_volume.HasNaN());
  // cost_volume.RemoveEndLabels();

  if (FLAGS_use_disc) {
    LOG(INFO) << "Discrete optimization";
    DiscreteMrfSettings settings;
    CHECK(ReadTxtFileToProto(FLAGS_settings, &settings));
    SolveDiscreteDepth(settings, cost_volume, &depthmap);
    prefix += "_disc";
    WriteProtoToTxtFile(prefix + "_settings.txt", settings);
  } else if (FLAGS_use_cont) {
    LOG(INFO) << "Continuous optimization";
    ContinuousMrfSettings settings;
    CHECK(ReadTxtFileToProto(FLAGS_settings, &settings));
    SolveContinuousDepth(settings, cost_volume, &depthmap);
    prefix += "_cont";
    WriteProtoToTxtFile(prefix + "_settings.txt", settings);
  } else if (FLAGS_use_densecrf) {
    DenseCrfSettings settings;
    CHECK(ReadTxtFileToProto(FLAGS_settings, &settings));
    SolveByDenseCrf(settings, cost_volume, &depthmap);
    prefix += "_crf";
    WriteProtoToTxtFile(prefix + "_settings.txt", settings);
  } else {
    WtaDepth(cost_volume, &depthmap);
    prefix += "_wta";
  }
  cost_volume.SetImage(image);
  PlaneSweepingData psd = cost_volume.GetPlaneSweepingData();
  depthmap.SetPlaneSweepingData(psd);
  psd.WriteTxt(prefix + "_info.txt");
  LOG(INFO) << "min depth: " << psd.min_depth;
  depthmap.WriteDepthImage(prefix + ".png");
  depthmap.WriteDepthImage(prefix + ".pgm");
  depthmap.WritePly(prefix + ".ply");

  return 0;
}

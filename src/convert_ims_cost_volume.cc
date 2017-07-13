const char *kUsageMessage = "Convert IMS (refocus) cost volume to furry " \
    "cost volume";

#include <fstream>
#include <string>

#include <opencv2/opencv.hpp>

#include "furry/common/init.h"

#include "cost_volume.h"

using namespace cv;
using namespace furry;
using namespace furry::tiny;
using namespace std;

DEFINE_string(in_cv, "", "");
DEFINE_double(min_depth, -1, "");
DEFINE_double(max_depth, -1, "");
DEFINE_double(depth_step, -1, "");
DEFINE_string(image_path, "", "");
DEFINE_string(out_cv, "", "");

void ReadImsCostVolume(const string &filename, CostVolume *cost_volume) {
  int width, height, num_samples;
  ifstream dim_file(filename + ".dim");
  dim_file >> num_samples >> width >> height;
  dim_file.close();
  LOG(INFO) << "Cost volume dimension: "
            << width << ' ' << height << ' ' << num_samples;
  cost_volume->Allocate(width, height, num_samples);
  ifstream ima_file(filename + ".ima");
  ima_file.read((char*)cost_volume->GetData(),
                width * height * num_samples * sizeof(float));
  ima_file.close();
  cost_volume->SetMinDepth(1);
  cost_volume->SetMaxDepth(num_samples + 1);
}

int main(int argc, char *argv[]) {
  Init(&argc, &argv, kUsageMessage);

  CostVolume cost_volume;
  ReadImsCostVolume(FLAGS_in_cv, &cost_volume);
  auto image = imread(FLAGS_image_path);
  resize(image, image, Size(cost_volume.GetWidth(), cost_volume.GetHeight()));
  cost_volume.SetImage(image);
  if (FLAGS_min_depth > 0) cost_volume.SetMinDepth(FLAGS_min_depth);
  if (FLAGS_max_depth > 0) cost_volume.SetMaxDepth(FLAGS_max_depth);
  cost_volume.WriteFile(FLAGS_out_cv);

  return 0;
}

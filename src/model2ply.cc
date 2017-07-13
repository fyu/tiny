const char *kUsageMessage = "Convert Model tracks to ply";

#include "furry/common/init.h"

#include "model.h"

using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_model, "", "");
DEFINE_string(out_path, "", "");

int main(int argc, char *argv[]) {
  Init(&argc, &argv, kUsageMessage);
  Model model;
  model.ReadFile(FLAGS_in_model);
  Model2Ply(FLAGS_out_path, model);
  return 0;
}
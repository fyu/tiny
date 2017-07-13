const char *kUsageMessage = "Convert Model to PMVS input";

#include "furry/common/init.h"

#include "model.h"

using namespace furry;
using namespace furry::tiny;

DEFINE_string(in_model, "", "");
DEFINE_string(out_dir, "", "");

int main(int argc, char *argv[]) {
  Init(&argc, &argv, kUsageMessage);
  Model model;
  model.ReadFile(FLAGS_in_model);
  Model2Pmvs(FLAGS_out_dir, model);
  return 0;
}

const kUsageMessage = "Optimize pixel depth with MRF";

#include "furry/common/init.h"

#include "mvs_opt_cont.h"

using namespace furry;
using namespace furry::tiny;

int main(int argc, char *argv[]) {
  Init(&argc, &argv, kUsageMessage);
  return 0;
}

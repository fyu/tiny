#include <time.h>

#include <glog/logging.h>
#include <Eigen/Dense>

#include "furry/common/parallel.h"

namespace furry {

inline int Init(int *argc, char ***argv, const char *usage_message) {
  tbb::task_scheduler_init tbb_init(FLAGS_p);
  google::InitGoogleLogging(**argv);
  gflags::SetUsageMessage(usage_message);
  gflags::ParseCommandLineFlags(argc, argv, true);
  // Eigen::initParallel();
  LOG(INFO) << **argv << " built time: " << __DATE__ << ' ' << __TIME__;
  return 0;
}

inline int init(int *argc, char ***argv, const char *usage_message) {
  tbb::task_scheduler_init tbb_init(FLAGS_p);
  google::InitGoogleLogging(**argv);
  gflags::SetUsageMessage(usage_message);
  gflags::ParseCommandLineFlags(argc, argv, true);
  // Eigen::initParallel();
  LOG(INFO) << **argv << " built time: " << __DATE__ << ' ' << __TIME__;
  return 0;
}

} // furry

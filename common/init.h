#ifndef FURRY_COMMON_INIT_H
#define FURRY_COMMON_INIT_H

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "furry/common/parallel.h"

DECLARE_int32(p);

namespace furry {
int Init(int *argc, char ***argv, const char *usage_message = nullptr);
int init(int *argc, char ***argv, const char *usage_message = nullptr);

#define INIT(argc, argv, usage_message)                                 \
  tbb::task_scheduler_init tbb_init(FLAGS_p);                           \
  do {                                                                  \
    google::InitGoogleLogging(*argv);                                   \
    google::SetUsageMessage(usage_message);                             \
    google::ParseCommandLineFlags(&argc, &argv, true);                    \
    LOG(INFO) << *argv << " built time: " << __DATE__ << ' ' << __TIME__; \
  } while (0)

} // furry

#include "init-inl.h"

#endif

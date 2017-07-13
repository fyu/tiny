#ifndef FURRY_COMMON_LOG_H
#define FURRY_COMMON_LOG_H

#include <iostream>
#include <string>

#include "gflags/gflags.h"

#include "furry/common/parallel.h"
#include "furry/common/clock.h"

DECLARE_int32(verbose);

namespace furry {

namespace internal {

template <typename T>
std::ostream& log(std::ostream& os, const T &arg) {
  os << arg << '\n';
  return os;
}

inline std::ostream& log(std::ostream& os) {
  return os;
}

template <typename T, typename... Args>
std::ostream& log(std::ostream& os, const T &arg, Args... args) {
  os << arg << ' ';
  return log(os, args...);
}

} // internal

#define F_LOG(os, ...) furry::internal::log(os, "[", __func__, "]", __VA_ARGS__)
#define F_VLOG(v, ...) do {if (FLAGS_verbose >= v) F_LOG(std::cout, __VA_ARGS__);} while(0)

#define F_ELOG(...) do {F_LOG(std::cerr, "ERROR:", __VA_ARGS__); exit(1);} while(0)
#define F_WLOG(...) F_LOG(std::cerr, "WARNING:", __VA_ARGS__)

class ProgressBar {
 public:
  ProgressBar(const std::string& title, int num_jobs);
  void FinishOne();
  void Done();
 private:
  std::string title_;
  int num_jobs_;
  int num_finished_jobs_;
  furry::Mutex mutex_;
  timeval clock_;
};

} // furry

#endif

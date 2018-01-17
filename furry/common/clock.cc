#include "clock.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

#include <sys/resource.h>

#include <glog/logging.h>

#include "furry/common/log.h"

namespace furry
{

namespace
{
timeval tictv;
timeval toctv;
}

timeval
tic()
{
  gettimeofday(&tictv, NULL);
  return tictv;
}

double
toc(timeval tv)
{
  if (tv.tv_sec == 0 && tv.tv_usec == 0) {
    tv = tictv;
  }
  gettimeofday(&toctv, NULL);
  return toctv.tv_sec - tv.tv_sec + (toctv.tv_usec - tv.tv_usec) / 1e6;
}

void PrintEventTime(const std::string &event_name, timeval tv) {
  if (event_name == "")
    LOG(INFO) << toc(tv) << " seconds elapsed\n";
  else
    LOG(INFO) << event_name << " took " << toc(tv) << " seconds\n";
}

void PrintEventTime(int verbose, const std::string &event_name, timeval tv) {
  if (FLAGS_verbose >= verbose)
    PrintEventTime(event_name, tv);
}

Timer::Timer() {
  gettimeofday(&begin_, NULL);
}

double Timer::Lap() {
  timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec - begin_.tv_sec + (time.tv_usec - begin_.tv_usec) / 1e6;
}

std::string print_memory_usage() {
  rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  double maxrss;
#ifdef __APPLE__
  maxrss = r_usage.ru_maxrss / 1024.0 / 1024.0;
#else
  maxrss = r_usage.ru_maxrss / 1024.0;
#endif
  std::stringstream ss;
  ss << maxrss << " MB physical memory";
  return ss.str();
}

} // furry

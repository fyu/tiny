#ifndef _FURRY_CORE_TIME_H_
#define _FURRY_CORE_TIME_H_

#include <sys/time.h>

#include <string>

#include <glog/logging.h>

namespace furry
{

// return time in second
timeval tic();
double toc(timeval tv = {0, 0});

std::string print_memory_usage();

void PrintEventTime(const std::string &event_name = "", timeval tv = {0, 0});
void PrintEventTime(int verbose, const std::string &event_name = "",
                    timeval tv = {0, 0});

#define print_event_time(event_name, timer) \
  LOG(INFO) << (event_name) << " took " << toc(timer) << " seconds and " \
  << print_memory_usage()

#define PRINT_EVENT_TIME print_event_time

class Timer {
 public:
  Timer();
  double Lap();
 private:
  timeval begin_;
};

} // furry

#endif // _FURRY_CORE_TIME_H_

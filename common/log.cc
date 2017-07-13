#include "furry/common/log.h"

#include <iostream>
#include <cstdio>

#include "furry/common/str.h"

using namespace std;

DEFINE_int32(verbose, 0, "");

namespace furry {

ProgressBar::ProgressBar(const std::string& title, int num_jobs)
    : title_(title), num_jobs_(num_jobs), num_finished_jobs_(0) {
  clock_ = tic();
  printf("%s: %2.2lf%% Elapsed: 0 sec ETA: N/A sec",
         title_.c_str(), (double) num_finished_jobs_ / num_jobs_ * 100);
  fflush(stdout);
}
void ProgressBar::FinishOne() {
  furry::Mutex::scoped_lock lock(mutex_);
  ++num_finished_jobs_;
  double time = toc(clock_);
  printf("\r%s: %2.2lf%% Elapsed: %0.2lf sec ETA: %0.2lf sec",
         title_.c_str(), (double) num_finished_jobs_ / num_jobs_ * 100,
         time, time / num_finished_jobs_ * num_jobs_ - time);
  fflush(stdout);
}

void ProgressBar::Done() {
  // printf("\n");
  printf(" Done\n");
  fflush(stdout);
}

} // furry

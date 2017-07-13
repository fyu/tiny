#include "furry/common/tbb.h"

#include <tbb/tbb.h>

namespace furry {

void ParFor(int begin, int end, std::function<void (int)> func) {
  tbb::parallel_for(
      tbb::blocked_range<int>(begin, end),
      [&](const tbb::blocked_range<int> &range) {
        for (int i = range.begin(); i < range.end(); ++i) {
          func(i);
        }
      });
}

void parfor(int begin, int end, std::function<void (int)> func) {
  ParFor(begin, end, func);
}

} // furry

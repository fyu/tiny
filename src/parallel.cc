#include "parallel.h"

#include <tbb/tbb.h>

namespace furry {
namespace tiny {

void ParFor(int begin, int end, std::function<void (int)> func) {
  tbb::parallel_for(
      tbb::blocked_range<int>(begin, end),
      [&](const tbb::blocked_range<int> &range) {
        for (int i = range.begin(); i < range.end(); ++i) {
          func(i);
        }
      });
}

} // tiny
} // furry

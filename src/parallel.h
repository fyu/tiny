#ifndef FURRY_TINY_PARALLEL_H
#define FURRY_TINY_PARALLEL_H

#include <functional>

namespace furry {
namespace tiny {

void ParFor(int begin, int end, std::function<void (int)> func);

} // tiny
} // furry

#endif

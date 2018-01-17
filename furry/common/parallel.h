#ifndef FURRY_COMMON_PARALLEL_H
#define FURRY_COMMON_PARALLEL_H

#include <tbb/tbb.h>

#include "furry/common/tbb.h"

namespace furry {

typedef tbb::mutex Mutex;
typedef tbb::mutex::scoped_lock ScopedLock;

} // furry

#endif

#ifndef FURRY_ALGORITHM_UTIL_H
#define FURRY_ALGORITHM_UTIL_H

#include <vector>
#include <functional>
#include <algorithm>

namespace furry {

template <typename T, typename CompFunc = std::less<T>>
T Max(const std::vector<T> &v, CompFunc comp_func = CompFunc());

template <typename T>
T Mean(const std::vector<T> &v);

template <typename T>
T Clamp(T v, const T &min, const T &max);

} // furry

#include "furry/algorithm/util-inl.h"

#endif // FURRY_ALGORITHM_UTIL_H

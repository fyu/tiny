#ifndef FURRY_COMMON_STL_H_
#define FURRY_COMMON_STL_H_

#include <utility>
#include <iostream>

#include "furry/common/hash.h"

namespace furry {

} // furry

namespace std {

template <typename T0, typename T1>
struct hash<std::pair<T0, T1>> {
  size_t operator () (const std::pair<T0, T1> &p) const {
    size_t seed = 0;
    seed = furry::CombineHash(seed, hash<T0>()(p.first));
    seed = furry::CombineHash(seed, hash<T1>()(p.second));
    return seed;
  }
}; // hash

template <typename T0, typename T1>
struct less<std::pair<T0, T1>> {
  size_t operator () (const std::pair<T0, T1> &p0,
                      const std::pair<T0, T1> &p1) const {
    less<T0> comp_t0;
    if (comp_t0(p0.first, p1.first))
      return true;
    else if (comp_t0(p1.first, p0.first))
      return false;
    if (less<T1>()(p0.second, p1.second))
      return true;
    return false;
  }
};

template <typename T0, typename T1>
std::ostream& operator << (std::ostream &os, const std::pair<T0, T1> &p) {
  os << p.first << ' ' << p.second;
  return os;
}

// template <typename T0, typename T1>
// void swap(pair<T0, T1> &a, pair<T0, T1> &b) {
//   swap(a.first, b.first);
//   swap(a.second, b.second);
// }

} // std

#endif

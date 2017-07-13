#ifndef FURRY_COMMON_HASH_H
#define FURRY_COMMON_HASH_H

#include <functional>
#include <utility>

namespace furry {

template <typename T>
size_t CombineHash(size_t seed, const T& v) {
  // From boost::hash_combine
  seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

template <typename T>
size_t combine_hash(size_t seed, const T& v) {
  // From boost::hash_combine
  seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

} // furry

#endif

#ifndef FURRY_COMMON_SEQUENCE_H
#define FURRY_COMMON_SEQUENCE_H

#include <vector>
#include <functional>

namespace furry {

template <typename ValueType, typename PredicateType>
std::vector<ValueType> Filter(const std::vector<ValueType> &values,
                              const std::vector<PredicateType> &predicates) {
  std::vector<ValueType> results;
  for (size_t i = 0; i < predicates.size(); ++i) {
    if (predicates[i])
      results.push_back(values[i]);
  }
  return results;
}

template <typename ValueType, typename PredicateFunc>
std::vector<ValueType> Filter(const std::vector<ValueType> &values,
                              PredicateFunc func) {
  std::vector<ValueType> results;
  for (size_t i = 0; i < values.size(); ++i) {
    if (func(values[i])) results.push_back(values[i]);
  }
  return results;
}

template <typename ValueType, typename IndexType>
std::vector<ValueType> Slice(const std::vector<ValueType> &values,
                             const std::vector<IndexType> &indexes) {
  std::vector<ValueType> results;
  results.reserve(indexes.size());
  for (auto index: indexes) {
    results.push_back(values[index]);
  }
  return results;
}

template <typename T, typename CompFunc>
std::vector<bool> Compare(const std::vector<T> &values, CompFunc comp_func) {
  std::vector<bool> result(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    result[i] = comp_func(values[i]);
  }
  return result;
}

} // furry

#endif // FURRY_COMMON_SEQUENCE_H

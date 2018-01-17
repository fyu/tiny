#ifndef FURRY_ALG_SELECT_H_
#define FURRY_ALG_SELECT_H_

#include <algorithm>
#include <functional>
#include <iostream>

namespace furry {

template <typename InputIterator,
          typename CompFunc = std::less<typename InputIterator::value_type>>
    InputIterator
SelectKth(InputIterator first,
          InputIterator last,
          int k,
          CompFunc comp = CompFunc()) {
  int last_pivot_rank, pivot_rank = 0;
  InputIterator pivot, big, small;
  while (first != last) {
    pivot = first;
    big = last;
    small = pivot + 1;
    last_pivot_rank = pivot_rank;
    while (big != small) {
      if (comp(*small, *pivot)) {
        std::swap(*small, *pivot);
        ++pivot;
        ++pivot_rank;
        ++small;
      } else {
        std::swap(*small, *(--big));
      }
    }
    if (pivot_rank == k)
      return pivot;
    if (pivot_rank < k) {
      first = small;
      ++pivot_rank;
    } else {
      last = pivot;
      pivot_rank = last_pivot_rank;
    }
  }
  return last;
} // SelectKth

template <typename RandomIt,
          typename CompFunc = std::less<typename RandomIt::value_type>>
    void
SelectTopK(RandomIt first,
           RandomIt last,
           int k,
           CompFunc comp = CompFunc()) {
  sort(first, last, comp);
  return;
  int last_pivot_rank, pivot_rank = 0;
  RandomIt pivot, big, small;
  ++k;
  RandomIt begin = first;
  while (first != last) {
    pivot = first;
    big = last;
    small = pivot + 1;
    last_pivot_rank = pivot_rank;
    while (big != small) {
      if (comp(*small, *pivot)) {
        std::swap(*small, *pivot);
        ++pivot;
        ++pivot_rank;
        ++small;
      } else {
        std::swap(*small, *(--big));
      }
    }
    if (pivot_rank == k) {
      last = pivot;
      break;
    }
    if (pivot_rank < k) {
      first = small;
      ++pivot_rank;
    } else {
      last = pivot;
      pivot_rank = last_pivot_rank;
    }
  }
  sort(begin, last, comp);
  return;
} // Select Top k

} // furry

#endif // FURRY_ALG_SELECT_H_

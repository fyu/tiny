#ifndef FURRY_GRID_INL_H
#define FURRY_GRID_INL_H

#include <cassert>
#include <cmath>
#include <cstdint>

#include <glog/logging.h>
#include <opencv2/core/core.hpp>

#include "furry/common/memory.h"

namespace furry {

template <typename T>
inline Grid3<T>::Grid3(Grid3<T> &&grid) {
  *this = grid;
}

template <typename T>
inline Grid3<T>::Grid3(const Grid3<T> &grid) {
  *this = grid;
}

template <typename T>
inline bool Grid3<T>::Allocate(int d0, int d1, int d2) {
  if (!(d0 > 0 && d1 > 0 && d2 > 0))
    return false;
  if (data_ != nullptr)
    delete data_;
  data_ = (T*)malloc(d0 * d1 * d2 * sizeof(T));
  if (data_ == nullptr)
    return false;
  sizes_[0] = d0;
  sizes_[1] = d1;
  sizes_[2] = d2;
  return true;
}

template <typename T>
inline void Grid3<T>::Deallocate() {
  if (data_ != nullptr) delete data_;
  data_ = nullptr;
  memset(sizes_, 0, sizeof(sizes_));
}

template <typename T>
inline bool Grid3<T>::IsEmpty() const {
  return data_ == nullptr;
}

template <typename T>
inline bool Grid3<T>::HasNaN() const {
  if (IsEmpty()) return false;
  int size = GetTotalSize();
  for (int i = 0; i < size; ++i) {
    if (isnan(data_[i])) return true;
  }
  return false;
}

template <typename T>
inline int Grid3<T>::GetSize(int d) const {
  assert(d >= 0 && d < 3);
  return sizes_[d];
}

template <typename T>
inline int Grid3<T>::GetTotalSize() const {
  return sizes_[0] * sizes_[1] * sizes_[2];
}

template <typename T>
inline const T& Grid3<T>::At(int i0, int i1, int i2) const {
  return const_cast<Grid3<T>*>(this)->At(i0, i1, i2);
}

template <typename T>
inline T& Grid3<T>::At(int i0, int i1, int i2) {
  assert(i0 < sizes_[0] && i1 < sizes_[1] && i2 < sizes_[2]);
  return data_[Grid2Array(i0, i1, i2)];
}

template <typename T>
inline const T* Grid3<T>::At(int i0, int i1) const {
  return &At(i0, i1, 0);
}

template <typename T>
inline T* Grid3<T>::At(int i0, int i1) {
  return &At(i0, i1, 0);
}

template <typename T>
inline T* Grid3<T>::GetData() {
  return data_;
}

template <typename T>
inline const T* Grid3<T>::GetData() const {
  return data_;
}

template <typename T>
inline int Grid3<T>::NumRows() const {
  return sizes_[0];
}

template <typename T>
inline int Grid3<T>::NumCols() const {
  return sizes_[1];
}

template <typename T>
inline int Grid3<T>::NumBytes() const {
  return GetTotalSize() * sizeof(T);
}


template <typename T>
inline void Grid3<T>::Set(int i0, int i1, int i2, T v) {
  data_[Grid2Array(i0, i1, i2)] = v;
}

template <typename T>
inline Grid3<T>::~Grid3() {
  if (data_)
    free(data_);
}

template <typename T>
inline Grid3<T>& Grid3<T>::operator = (Grid3<T> &&grid) {
  if (this == &grid)
    return *this;
  if (data_ != nullptr) delete data_;
  data_ = grid.data_;
  memcpy(sizes_, grid.sizes_, sizeof(sizes_));
  grid.data_ = nullptr;
  grid.sizes_[0] = 0;
  grid.sizes_[1] = 0;
  grid.sizes_[2] = 0;
  return *this;
}

template <typename T>
inline Grid3<T>& Grid3<T>::operator = (const Grid3<T> &grid) {
  if (this == &grid)
    return *this;
  memcpy(sizes_, grid.sizes_, sizeof(sizes_));
  Allocate(sizes_[0], sizes_[1], sizes_[2]);
  memcpy(data_, grid.data_, NumBytes());
  return *this;
}

template <typename T>
inline Grid3<T>& Grid3<T>::operator = (T value) {
  CHECK(!IsEmpty());
  int size = GetTotalSize();
  for (int i = 0; i < size; ++i) {
    data_[i] = value;
  }
  return *this;
}

template <typename T>
inline bool Grid3<T>::ReadBin(const std::string &filename) {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "r"));
  if (!file) return false;
  return ReadBin(file.Get());
}

template <typename T>
inline bool Grid3<T>::WriteBin(const std::string &filename) const {
  auto file = ScopedPtr<FILE>(fopen(filename.c_str(), "w"));
  if (!file) return false;
  return WriteBin(file.Get());
}

template <typename T>
inline bool Grid3<T>::ReadBin(FILE *file) {
  int type;
  fread(&type, sizeof(type), 1, file);
  if (type != cv::DataType<T>::type)
    return false;
  fread(sizes_, sizeof(sizes_[0]), 3, file);
  LOG(INFO) << "Sizes " << sizes_[0] << ' '
            << sizes_[1] << ' '
            << sizes_[2] << ' ';
  Allocate(sizes_[0], sizes_[1], sizes_[2]);
  fread(data_, sizeof(T), GetTotalSize(), file);
  return true;
}

template <typename T>
inline bool Grid3<T>::WriteBin(FILE *file) const {
  int type = cv::DataType<T>::type;
  fwrite(&type, sizeof(type), 1, file);
  fwrite(sizes_, sizeof(sizes_[0]), 3, file);
  fwrite(data_, sizeof(T), GetTotalSize(), file);
  return true;
}

template <typename T>
inline int Grid3<T>::Grid2Array(int i0, int i1, int i2) const {
  return (i0 * sizes_[1] + i1) * sizes_[2] + i2;
}

} // furry

#endif // FURRY_GRID_INL_H

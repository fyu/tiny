#ifndef FURRY_GRID_H
#define FURRY_GRID_H

#include <cstdio>

namespace furry {

// Grid is in row major if you treat it as an image
template <typename T>
class Grid3 {
 public:
  Grid3() {}
  Grid3(Grid3<T> &&grid);
  Grid3(const Grid3<T> &grid);
  ~Grid3();

  Grid3<T>& operator = (Grid3<T> &&grid);
  Grid3<T>& operator = (const Grid3<T> &grid);
  Grid3<T>& operator = (T value);

  bool Allocate(int d0, int d1, int d2);
  void Deallocate();

  bool IsEmpty() const;
  bool HasNaN() const;

  // Accessor
  T& At(int i0, int i1, int i2);
  const T& At(int i0, int i1, int i2) const;
  T* At(int i0, int i1);
  const T* At(int i0, int i1) const;
  T* GetData();
  const T* GetData() const;
  int GetSize(int d) const;
  int GetTotalSize() const;
  int NumRows() const;
  int NumCols() const;
  int NumBytes() const;

  // Mutator
  void Set(int i0, int i1, int i2, T v);

  bool ReadBin(const std::string &filename);
  bool WriteBin(const std::string &filename) const;

 protected:
  bool ReadBin(FILE *file);
  bool WriteBin(FILE *file) const;

  int Grid2Array(int i0, int i1, int i2) const;

 private:
  T *data_ = nullptr;
  int sizes_[3] = {0, 0, 0};
};

typedef Grid3<double> Grid3d;
typedef Grid3<float> Grid3f;

} // furry

#include "grid-inl.h"

#endif // FURRY_GRID_H

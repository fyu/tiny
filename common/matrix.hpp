#ifndef _FURRY_COMMON_MAT_H_
#define _FURRY_COMMON_MAT_H_

//#include "opencv2/core/core.hpp"
#include <iostream>
#include "defs.hpp"

#include <algorithm>
#include <cstring>
#include <vector>
#include <cassert>
#include <algorithm>
#include <limits>
#include <initializer_list>

namespace furry
{

/*td::ostream& operator << (std::ostream& os, const cv::Point2d& p);
std::ostream& operator << (std::ostream& os, const cv::Point3d& p);
template<typename T, int n>
std::ostream& operator << (std::ostream& os, const cv::Vec<T, n>& p);

////////////////////////////////////////////////////////////
// inline implementation
////////////////////////////////////////////////////////////

inline std::ostream& operator << (std::ostream& os, const cv::Point2d& p)
{
  os << p.x << ' ' << p.y;
  return os;
}

inline std::ostream& operator << (std::ostream& os, const cv::Point3d& p)
{
  os << p.x << ' ' << p.y << ' ' << p.z;
  return os;
}

template<typename DataType, int n>
inline std::ostream& operator << (std::ostream& os, const cv::Vec<DataType, n>& v)
{
  os << v[0];
  for (int i = 1; i < n; ++i)
  {
    os << ' ' << v[i];
  }
  return os;
  }*/

////////////////////////////////////////////////////////////
// Size
//////////////////////////////////////////////////////////

template <typename DataType, int Dims>
class Size_
{
public:
  Size_();
  Size_(std::initializer_list<DataType> d);

  DataType& operator [] (int index);
  const DataType& operator [] (int index) const;

private:
  DataType data_ = {0};
};

typedef Size_<int, 2> Size2i;
typedef Size2i Size;

// Data Traits
template <typename DataType>
struct data_traits
{
  typedef DataType value_type;
  static const int type_size = sizeof(DataType);
  static const int num_channels = 1;
};

class Range
{
public:

  static Range all();

  bool is_all();
  bool has_all_begin();
  bool has_all_end();


  Range(int begin);
  Range(int begin, int end);
  bool empty() const;
  int size() const;
  bool operator == (const Range& r) const;
  bool operator != (const Range& r) const;

  // Check whether the begin or end is all
  // if it is, give the concrete value
  void set_concrete(int begin, int end);

  int begin = Range::MIN;
  int end = Range::MAX;

private:

  static const int MAX = std::numeric_limits<int>::max();
  static const int MIN = std::numeric_limits<int>::min();
};

// The base of planar matrix
// Provide the basic functions common to 2 dimentional matrix
template <typename Derived>
class PlanarMatBase
{
public:
  typedef typename data_traits<Derived>::value_type value_type;

  value_type& operator () (int row, int col);
  const value_type& operator () (int row, int col) const;

  int num_rows() const;
  int num_cols() const;
};

// Matrix
// Stored in row major
template <typename DataType, int... dims> class Mat;

template <typename DataType, int nrows, int ncols>
class Mat<DataType, nrows, ncols> : PlanarMatBase<Mat<DataType, nrows, ncols>>
{
public:
  typedef DataType value_type;

  Mat();
  Mat(DataType* data);
  Mat(std::initializer_list<DataType> data);
  Mat(const DataType& v);
  Mat(const std::vector<DataType>& v);
  Mat(const Mat<DataType, nrows, ncols>& m);
  Mat(const Size& s);

  DataType& operator () (int index);
  DataType& operator () (int row, int col);
  DataType& operator [] (int index);
  const DataType& operator () (int index) const;
  const DataType& operator () (int row, int col) const;
  const DataType& operator [] (int index) const;

  int num_rows() const;
  int num_cols() const;

protected:
  static const int kTypeSize = sizeof(DataType);
  DataType data_[nrows * ncols] = {0};
};

template <typename DataType, int nrows, int ncols> Mat<DataType, nrows, ncols>&
operator *= (Mat<DataType, nrows, ncols>& m, double s);

template <typename DataType, int nrows, int ncols> Mat<DataType, nrows, ncols>
operator * (Mat<DataType, nrows, ncols>& m, double s);

template <typename DataType, int nrows, int ncols> Mat<DataType, nrows, ncols>
operator * (double s, Mat<DataType, nrows, ncols>& m);

template <typename DataType, int nrows, int ncols> Mat<DataType, nrows, ncols>&
operator += (Mat<DataType, nrows, ncols>& m1, Mat<DataType, nrows, ncols>& m2);

template <typename DataType, int nrows, int ncols> Mat<DataType, nrows, ncols>
operator + (Mat<DataType, nrows, ncols>& m1, Mat<DataType, nrows, ncols>& m2);

template <typename DataType, int n>
using Vec = Mat<DataType, n, 1>;

// dynamic N dimentional array
// optimized for 2 dimentional operation
template <typename DataType>
class Mat<DataType> : PlanarMatBase<Mat<DataType>>
{
public:
  typedef DataType value_type;
  typedef typename DataType::value_type channel_type;

  static const int kTypeSize = sizeof(DataType);
  static const int kAutoStep = -1;

  Mat();
  Mat(int ndims,
      const int* sizes,
      DataType* data,
      int row_step = kAutoStep,
      bool flip = false,
      bool copy = false);

  Mat(int rows,
      int cols,
      DataType* data,
      int row_step = kAutoStep,
      bool flip = false,
      bool copy = false);

  Mat(const Size& size);

  Mat(const Mat<DataType>& m,
      Range row_range = Range::all(),
      Range col_range = Range::all());

  ~Mat();

  // Turn this matrix into null matrix
  // delete the data if reference count becomes 0
  void release();

  Mat<DataType>& operator = (const Mat<DataType>& m);

  DataType& operator () (int index);
  DataType& operator () (int row, int col);
  DataType& operator [] (int index);
  const DataType& operator () (int index) const;
  const DataType& operator () (int row, int col) const;
  const DataType& operator [] (int index) const;

  int get_size(int dim) const;
  int num_dims() const;
  const DataType* get_row(int row) const;
  Mat<DataType> submat(Range row_range = Range::all(),
                       Range col_range = Range::all());

  int num_rows() const;
  int num_cols() const;

protected:
  void _reset();
  void _init_by_data(DataType* data,
                     int row_step,
                     bool flip,
                     bool copy);

protected:

  struct MatSize
  {
    MatSize(int* p);
    const int& operator [] (int i) const;
    int& operator [] (int i);
    int* p = NULL;
  };

  struct MatStep
  {
    MatStep();
    MatStep(int row_step);
    MatStep(int row_step, int col_step);
    const int& operator [] (int i) const;
    int& operator [] (int i);
    int buf[2] = {1, 1};
  };

protected:
  // number of dimentions of the matrix
  int num_dims_ = 0;
  // Size of each matrix dimension
  int num_rows_ = 0;
  int num_cols_ = 0;
  MatSize _size = MatSize(&num_rows_);

  // steps between rows and cols
  MatStep _step;
  // number of total elements
  //int _data_size = 0;
  // reference count of the data
  // if reference count is -1, the data is not owned by the matrix
  int* _refcount = NULL;
  // Data array
  DataType* data_ = NULL;
  DataType* data_limit_ = NULL;
  DataType* data_begin_ = NULL;
  DataType* data_end_ = NULL;
  // Whether the data is fliped on y-axis
  bool _flip = false;
};



template <typename DataType, int nchannels>
Mat<Vec<DataType, nchannels>> data_to_mat(int rows,
                                          int cols,
                                          const DataType* data,
                                          bool flip = false);

template <typename Derived> typename PlanarMatBase<Derived>::value_type
bilinear_at(const PlanarMatBase<Derived>& m, double row, double col);

template <typename InputMat, typename SampleFunction, typename OutputMat = InputMat>
OutputMat sample(const InputMat& m, Size outsize, SampleFunction func = bilinear_at);

#include "mat-inl.hpp"

} // furry

#endif // _FURRY_COMMON_MAT_H_


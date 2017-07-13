////////////////////////////////////////////////////////////
// Data Traits
////////////////////////////////////////////////////////////

template <typename DataType, int n>
struct data_traits<Vec<DataType, n>>
{
  typedef DataType value_type;
  static const int type_size = sizeof(Vec<DataType, n>);
  static const int num_channels = n;
};

////////////////////////////////////////////////////////////
// Size_
////////////////////////////////////////////////////////////

template <typename DataType, int Dims> inline
Size_<DataType, Dims>::Size_()
{
}

template <typename DataType, int Dims> inline
Size_<DataType, Dims>::Size_(std::initializer_list<DataType> d)
{
  assert(d.size() == Dims);
  copy(d.begin(), d.end(), data_);
}

template <typename DataType, int Dims> inline DataType&
Size_<DataType, Dims>::operator [] (int index)
{
  return data_[index];
}

template <typename DataType, int Dims> inline const DataType&
Size_<DataType, Dims>::operator [] (int index) const
{
  return const_cast<Size_<DataType, Dims>*>(this)->operator [] (index);
}

////////////////////////////////////////////////////////////
// Range
////////////////////////////////////////////////////////////

Range::Range(int b, int e): begin(b), end(e)
{
}

Range::Range(int b) : begin(b)
{
}

inline Range
Range::all()
{
  return Range(MIN, MAX);
}

inline bool
Range::is_all()
{
  return has_all_begin() && has_all_end();
}

inline bool
Range::has_all_begin()
{
  return begin == MIN;
}

inline bool
Range::has_all_end()
{
  return end == MAX;
}

inline bool
Range::empty() const
{
  return begin >= end;
}

inline int
Range::size() const
{
  return end - begin;
}

inline void
Range::set_concrete(int b, int e)
{
  begin = (begin == MIN)? b : begin;
  end = (end == MAX)? e : end;
}

inline bool
Range::operator == (const Range& r) const
{
  return begin == r.begin && end == r.end;
}

inline bool
Range::operator != (const Range& r) const
{
  return !(*this == r);
}

////////////////////////////////////////////////////////////
// PlanarMatBase
////////////////////////////////////////////////////////////

template <typename Derived> inline typename PlanarMatBase<Derived>::value_type&
PlanarMatBase<Derived>::operator () (int row, int col)
{
  return static_cast<Derived*>(this)->operator () (row, col);
}

template <typename Derived> inline const typename PlanarMatBase<Derived>::value_type&
PlanarMatBase<Derived>::operator () (int row, int col) const
{
  return static_cast<const Derived*>(this)->operator () (row, col);
}

template <typename Derived> inline int
PlanarMatBase<Derived>::num_rows() const
{
  return static_cast<const Derived*>(this)->num_rows();
}

template <typename Derived> inline int
PlanarMatBase<Derived>::num_cols() const
{
  return static_cast<const Derived*>(this)->num_cols();
}

////////////////////////////////////////////////////////////
// two dim mat
////////////////////////////////////////////////////////////

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat()
{
}

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat(DataType* data)
{
  memcpy(data_, data, rows * cols * kTypeSize);
}

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat(std::initializer_list<DataType> data)
{
  assert(data.size() == rows * cols);
  std::copy(data.begin(), data.end(), data_);
}

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat(const DataType& v)
{
  std::for_each(data_, data_ + rows*cols,
                [&v](DataType& e) {
                  e = v;
                });
}

template <typename DataType, int nrows, int ncols> inline
Mat<DataType, nrows, ncols>::Mat(const Mat<DataType, nrows, ncols>& m)
{
  memcpy(data_, m.data_, nrows * ncols * data_traits<DataType>::type_size);
}

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat(const std::vector<DataType>& v)
{
  assert(v.size() > rows * cols);
  memcpy(data_, &v[0], rows * cols * kTypeSize);
}

template <typename DataType, int rows, int cols> inline
Mat<DataType, rows, cols>::Mat(const Size& s)
{
  assert(s[0] == rows && s[1] == cols);
}

template <typename DataType, int rows, int cols> inline DataType&
Mat<DataType, rows, cols>::operator () (int index)
{
  return data_[index];
}

template <typename DataType, int rows, int cols> inline DataType&
Mat<DataType, rows, cols>::operator () (int row, int col)
{
  return data_[row * cols + col];
}

template <typename DataType, int rows, int cols> inline DataType&
Mat<DataType, rows, cols>::operator [] (int index)
{
  return data_[index];
}

template <typename DataType, int rows, int cols> inline const DataType&
Mat<DataType, rows, cols>::operator () (int index) const
{
  return data_[index];
}

template <typename DataType, int rows, int cols> inline const DataType&
Mat<DataType, rows, cols>::operator () (int row, int col) const
{
  return data_[row * cols + col];
}

template <typename DataType, int rows, int cols> inline const DataType&
Mat<DataType, rows, cols>::operator [] (int index) const
{
  return data_[index];
}

template <typename DataType, int Rows, int Cols> inline int
Mat<DataType, Rows, Cols>::num_rows() const
{
  return Rows;
}

template <typename DataType, int Rows, int Cols> inline int
Mat<DataType, Rows, Cols>::num_cols() const
{
  return Cols;
}

////////////////////////////////////////////////////////////
// Two dim mat operation
////////////////////////////////////////////////////////////


template <typename DataType, int nrows, int ncols>
inline Mat<DataType, nrows, ncols>&
operator *= (Mat<DataType, nrows, ncols>& m, double s)
{
  for (int i = 0; i < nrows * ncols; ++i)
    m[i] *= s;
  return m;
}

template <typename DataType, int nrows, int ncols>
inline Mat<DataType, nrows, ncols>
operator * (Mat<DataType, nrows, ncols>& m, double s)
{
  Mat<DataType, nrows, ncols> tmp(m);
  return tmp *= s;
}

template <typename DataType, int nrows, int ncols>
inline Mat<DataType, nrows, ncols>
operator * (double s, Mat<DataType, nrows, ncols>& m)
{
  return m * s;
}

template <typename DataType, int nrows, int ncols>
inline Mat<DataType, nrows, ncols>&
operator += (Mat<DataType, nrows, ncols>& m1, Mat<DataType, nrows, ncols>& m2)
{
  for (int i = 0; i < nrows * ncols; ++i)
    m1[i] += m2[i];
  return m1;
}

template <typename DataType, int nrows, int ncols>
inline Mat<DataType, nrows, ncols>
operator + (Mat<DataType, nrows, ncols>& m1, Mat<DataType, nrows, ncols>& m2)
{
  Mat<DataType, nrows, ncols> tmp(m1);
  return tmp += m2;
}

////////////////////////////////////////////////////////////
// Dynamic matrix inline
////////////////////////////////////////////////////////////

template <typename DataType> inline
Mat<DataType>::Mat()
{
  _reset();
}

template <typename DataType> inline
Mat<DataType>::Mat(int ndims,
                   const int* sizes,
                   DataType* data,
                   int row_step,
                   bool flip,
                   bool copy): num_dims_(ndims)
{
  if (num_dims_ >= 1)
    num_rows_ = sizes[0];
  if (num_dims_ >= 2)
    num_cols_ = sizes[1];
  if (num_dims_ > 2)
  {
    _size.p = new int[num_dims_];
    memcpy(_size, sizes, num_dims_ * sizeof(int));
  }
  _init_by_data(data, row_step, copy, flip);
}

template <typename DataType> inline
Mat<DataType>::Mat(int rows,
                   int cols,
                   DataType* data,
                   int row_step,
                   bool flip,
                   bool copy): num_dims_(2), num_rows_(rows), num_cols_(cols)
{
  _init_by_data(data, row_step, copy, flip);
}

template <typename DataType> inline
Mat<DataType>::Mat(const Size& size) : num_dims_(2),
                                       num_rows_(size[0]),
                                       num_cols_(size[1]),
                                       _step(num_cols, 1),
                                       _refcount(new int),
                                       data_(new DataType[_step[0] * num_rows]),
                                       data_limit_(data_ + _step[0] * num_rows),
                                       data_begin_(data_),
                                       data_end_(data_limit_)
{
  *_refcount = 1;
}

template <typename DataType> inline
Mat<DataType>::Mat(const Mat<DataType>& m,
                   Range row_range,
                   Range col_range)
  : num_dims_(m.num_dims_), _step(m._step) ,data_(m.data_), data_limit_(m.data_limit_)
{
  row_range.set_concrete(0, num_rows_);
  col_range.set_concrete(0, num_cols_);
  num_rows_ = row_range.size();
  num_cols_ = col_range.size();
  if (num_dims_ > 2)
  {
    _size.p = new int[num_dims_];
    _size[0] = num_rows_;
    _size[1] = num_cols_;
    memcpy(_size.p+2, m._size.p+2, (num_dims_-1) * sizeof(int));
  }
  _refcount = m._refcount;
  if (_refcount != NULL)
    *_refcount += 1;
  data_begin_ = data_ + row_range.begin * _step[0] + col_range.begin * _step[1];
  data_end_ = data_ + row_range.end * _step[0] + col_range.begin * _step[1];
}

template <typename DataType> inline
Mat<DataType>::~Mat()
{
  release();
}

template <typename DataType> inline void
Mat<DataType>::_init_by_data(DataType* data,
                             int row_step,
                             bool flip,
                             bool copy)
{
  _flip = flip;
  int data_size = 1;
  for (int i = 0; i < num_dims_; ++i)
  {
    data_size *= _size[i];
  }

  _step[0] = data_size / _size[0];
  if (row_step == kAutoStep)
    row_step = _step[0];

  if (copy)
  {
    _refcount = new int;
    *_refcount = 1;
    data_ = new DataType[data_size];
    if (row_step == _step[0])
      memcpy(data_, data, data_size * kTypeSize);
    else
    {
      for (int i = 0; i < _size[0]; ++i)
        memcpy(data_ + _step[0] * i, data + row_step * i, _step[0] * kTypeSize);
    }
  }
  else
  {
    _step[0] = row_step;
    data_ = data;
  }
  data_begin_ = data_;
  data_end_ = data_ + data_size;
  data_limit_ = data_end_;
}

template <typename DataType> inline void
Mat<DataType>::_reset()
{
  num_dims_ = 0;
  _size.p = &num_rows_;
  _step = MatStep(1, 1);
  _refcount = NULL;
  data_ = NULL;
  data_begin_ = NULL;
  data_end_ = NULL;
  data_limit_ = NULL;
  _flip = false;
}

template <typename DataType> inline void
Mat<DataType>::release()
{
  if (_refcount != NULL)
  {
    *_refcount -= 1;
    if (*_refcount == 0)
    {
      delete data_;
      delete _refcount;
    }
  }
  if(_size.p != &num_rows_)
    delete _size.p;
  _reset();
}

template <typename DataType> inline Mat<DataType>&
Mat<DataType>::operator = (const Mat<DataType>& m)
{
  if (this != &m)
  {
    // Release the current matrix
    release();

    num_dims_ = m.num_dims_;
    num_rows_ = m.num_rows_;
    num_cols_ = m.num_cols_;
    //memcpy(_sizes, m._sizes, num_dims_ * sizeof(int));
    if (num_dims_ > 2)
    {
      _size.p = new int[num_dims_];
      memcpy(_size.p, m._size.p, num_dims_ * sizeof(int));
    }
    _step = m._step;
    _refcount = m._refcount;
    if (_refcount != NULL)
      *_refcount += 1;

    data_ = m.data_;
    data_begin_ = m.data_begin_;
    data_end_ = m.data_end_;
    data_limit_ = m.data_limit_;
    _flip = m._flip;
  }
  return *this;
}

template <typename DataType> inline DataType&
Mat<DataType>::operator () (int index)
{
  return data_begin_[index];
}

template <typename DataType> inline const DataType&
Mat<DataType>::operator () (int index) const
{
  return const_cast<Mat<DataType>*>(this)->operator() (index);
}

template <typename DataType> inline DataType&
Mat<DataType>::operator [] (int index)
{
  return *this(index);
}

template <typename DataType> inline const DataType&
Mat<DataType>::operator [] (int index) const
{
  return *this(index);
}

template <typename DataType> inline DataType&
Mat<DataType>::operator () (int row, int col)
{
  return data_begin_[row * _step[0] + col * _step[1]];
}

template <typename DataType> inline const DataType&
Mat<DataType>::operator () (int row, int col) const
{
  return const_cast<Mat<DataType>*>(this)->operator() (row, col);
}

template <typename DataType> inline int
Mat<DataType>::get_size(int dim) const
{
  assert(dim > 0 && dim < num_dims_);
  return _size[dim];
}

template <typename DataType> inline int
Mat<DataType>::num_dims() const
{
  return num_dims_;
}

template <typename DataType> inline const DataType*
Mat<DataType>::get_row(int row) const
{
  return data_begin_ + row * _step[0];
}

template <typename DataType> inline Mat<DataType>
Mat<DataType>::submat(Range row_range, Range col_range)
{
  return Mat<DataType>(*this, row_range, col_range);
}

// Mat::MatSize
template <typename DataType>
Mat<DataType>::MatSize::MatSize(int* p)
{
  this->p = p;
}

template <typename DataType> const int&
Mat<DataType>::MatSize::operator [] (int i) const
{
  return p[i];
}

template <typename DataType> int&
Mat<DataType>::MatSize::operator [] (int i)
{
  return p[i];
}

template <typename DataType>
Mat<DataType>::MatStep::MatStep()
{
}

// Mat::MatStep
template <typename DataType>
Mat<DataType>::MatStep::MatStep(int row_step)
{
  buf[0] = row_step;
}

template <typename DataType>
Mat<DataType>::MatStep::MatStep(int row_step, int col_step)
{
  buf[1] = col_step;
}

template <typename DataType> const int&
Mat<DataType>::MatStep::operator [] (int i) const
{
  return buf[i];
}

template <typename DataType> int&
Mat<DataType>::MatStep::operator [] (int i)
{
  return buf[i];
}

template <typename DataType> inline int
Mat<DataType>::num_rows() const
{
  return num_rows_;
}

template <typename DataType> inline int
Mat<DataType>::num_cols() const
{
  return num_cols_;
}

////////////////////////////////////////////////////////////
// Implementation of functions
////////////////////////////////////////////////////////////

template <typename DataType, int nchannels> inline Mat<Vec<DataType, nchannels>>
data_to_mat(int rows, int cols, const DataType* data, bool flip)
{
  return Mat<Vec<DataType, nchannels>>(rows,
                                       cols,
                                       (Vec<DataType, nchannels>*)data,
                                       Mat<DataType>::kAutoStep,
                                       flip,
                                       true);
}

template <typename Derived> inline typename PlanarMatBase<Derived>::value_type
bilinear_at(const PlanarMatBase<Derived>& m, double row, double col)
{
  typedef typename PlanarMatBase<Derived>::value_type value_type;

  if (row < 0 || row > m.num_rows() - 1 || col < 0 || col > m.num_cols() - 1)
    return value_type();

  int x0 = row;
  int x1 = x0 + 1;
  int y0 = col;
  int y1 = y0 + 1;

  // in case the sampled position is on the border of the matrix
  x1 = (x1 == m.num_rows()) ? x0 : x1;
  y1 = (y1 == m.num_cols()) ? y0 : y1;

  value_type v0 = m(x0, y0) * (row - x0) + m(x1, y0) * (x1 - row);
  value_type v1 = m(x0, y1) * (row - x0) + m(x1, y1) * (x1 - row);
  return v0 * (col - y0) + v1 * (y1 - col);
}

template <typename InputMat, typename SampleFunction, typename OutputMat = InputMat>
OutputMat sample(const InputMat& m, Size outsize, SampleFunction func)
{
  OutputMat output_mat(outsize);
  Size insize = {m.num_rows, m.num_rows};
  double row_step = ((double)insize[0]) / outsize[0];
  double col_step = ((double)insize[1]) / outsize[1];
  for (int i = 0; i < outsize[0]; ++i)
  {
    for (int j = 0; j < outsize[1]; ++j)
    {
      output_mat(i, j) = func(m, i * row_step, j * col_step);
    }
  }
  return output_mat;
}

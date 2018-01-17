/* ********************************************************** matrix.cpp *** *
 * 行列クラス
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 * 
 *                                    Time-stamp: <2008-03-07 13:51:54 sugaya>
 * ************************************************************************* */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <cmath>
#include "macros.h"
#include "matrix.h"

/* ************************************************************************* */
static bool
IsZero (double	val,
	double	epsilon = -1) {
  if (epsilon == -1) epsilon = MATH_DBL_EPSILON;
  if (fabs (val) < epsilon) {
    return true;
  } else {
    return false;
  }
}

/* ************************************************************************* *
 * コンストラクタ & デストラクタ
 * ************************************************************************* */

/* デフォルトコンストラクタ ************************************************ */
Matrix::Matrix (void) {
  row = 0;
  col = 0;
  nelements = 0;
  data = 0;
  data_ptr = 0;
  reference = false;
}

/* 次元を指定したコンストラクタ ******************************************** */
Matrix::Matrix (int 	_row,
		int	_col) {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  row = _row;
  col = _col;
  nelements = row * col;
  reference = false;
  
  data = new double [row * col];
  data_ptr = new double_ptr [row];
  
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0; n < nelements; n++) data[n] = 0.0;
}

/* 次元と値を指定したコンストラクタ **************************************** */
Matrix::Matrix (int	_row,
		int	_col,
		double	a, ...) {
  va_list	ptr;
  
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  row = _row;
  col = _col;
  nelements = row * col;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];
  
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  
  *data = a;
  va_start (ptr, a);
  for (int n = 1; n < nelements; n++) {
    data[n] = (double) va_arg (ptr, double);
  }
  va_end (ptr);
}

/* 次元と値を指定したコンストラクタ **************************************** */
Matrix::Matrix (int	_row,
		int	_col,
		double	*a) {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  row = _row;
  col = _col;
  nelements = row * col;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];
  
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0; n < nelements; n++) data[n] = a[n];
}

/* 次元と値を指定したコンストラクタ **************************************** */
Matrix::Matrix (int	_row,
		int	_col,
		double	**a) {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  row = _row;
  col = _col;
  nelements = row * col;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];
  
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0; n < row; n++) {
    for (int m = 0; m < col; m++) data_ptr[n][m] = a[n][m];
  }
}

/* 行列を指定したコンストラクタ ******************************************** */
Matrix::Matrix (const Matrix&	a) {
  if (a.row < 1 || a.col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  row  = a.row;
  col  = a.col;
  nelements = a.nelements;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];

  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0; n < nelements; n++) data[n] = a.data[n];
}

/* ベクトルを指定したコンストラクタ **************************************** */
Matrix::Matrix (const Vector&	v,
		int    		_row,
		int		_col) {
  if (v.Dimension() != _row * _col) {
    fprintf (stderr, "Specified dimension is invalid.\n");
    return;
  }
  row  = _row;
  col  = _col;
  nelements = row * col;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];

  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0; n < nelements; n++) data[n] = v[n];
}

/* 行列を指定したコンストラクタ ******************************************** */
Matrix::Matrix (const Matrix&	a,
		int		start_row,
		int		start_col,
		int		_row,
		int		_col) {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  int end_row = start_row + _row - 1;
  int end_col = start_col + _col - 1;

  if (start_row < 0 || end_row >= a.row ||
      start_col < 0 || end_col >= a.col) {      
    fprintf (stderr, "Specified index is out of range.\n");
    return;
  }
  row = _row;
  col = _col;
  nelements = row * col;
  reference = false;

  data = new double[row * col];
  data_ptr = new double_ptr [row];
  
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
  for (int n = 0, i = start_row; n < row; n++, i++) {
    for (int m = 0, j = start_col; m < col; m++, j++) {
      data_ptr[n][m] = a.data_ptr[i][j];
    }
  }
}

/* ベクトルを指定したコンストラクタ **************************************** */
Matrix::Matrix (const Vector&	a,
		const Vector&	b) {
  if (a.Dimension() < 1 || b.Dimension() < 1) {
    fprintf (stderr, "Vector dimension must be greater than zero.\n");
    return;
  } else if (a.Dimension() != b.Dimension()) {
    fprintf (stderr, "Dimensions of two vectors do not match.\n");
    return;
  } else {
    row  = a.Dimension();
    col  = a.Dimension();
    nelements = row * col;
    reference = false;
    
    data = new double[row * col];
    data_ptr = new double_ptr [row];
    
    for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
    for (int n = 0; n < row; n++) {
      for (int m = 0; m < col; m++) data_ptr[n][m] = a[n] * b[m];
    }
  }
}

/* ベクトルを指定したコンストラクタ **************************************** */
Matrix::Matrix (const Vector&	v,
		double		angle) {
  if (v.Dimension() != 3) {
    fprintf (stderr, "Vector dimension must be three.\n");
    return;
  }
  Vector v_(v.Dimension());
  double c, s;
  
  v_ = Normalize (v);
  angle = angle * M_PI / 180.0;
  c = cos (angle);
  s = sin (angle);
  
  row = 3;
  col = 3;
  nelements = 9;
  reference = false;
  
  data = new double[row * col];
  data_ptr = new double_ptr [row];
    
  for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;

  *(data    ) = c + v_[0] * v_[0] * (1.0 - c);
  *(data + 1) = v_[0] * v_[1] * (1.0 - c) - v_[2] * s;
  *(data + 2) = v_[0] * v_[2] * (1.0 - c) + v_[1] * s;
  *(data + 3) = v_[1] * v_[0] * (1.0 - c) + v_[2] * s;
  *(data + 4) = c + v_[1] * v_[1] * (1.0 - c);
  *(data + 5) = v_[1] * v_[2] * (1.0 - c) - v_[0] * s;
  *(data + 6) = v_[2] * v_[0] * (1.0 - c) - v_[1] * s;
  *(data + 7) = v_[2] * v_[1] * (1.0 - c) + v_[0] * s;
  *(data + 8) = c + v_[2] * v_[2] * (1.0 - c);
}	

/* デストラクタ ************************************************************ */
Matrix::~Matrix (void) {
  if (!reference && data) delete [] data;
  if (!reference && data_ptr) delete [] data_ptr;
  data = 0;
  data_ptr = 0;
  row = 0;
  col = 0;
  nelements = 0;
  reference = false;
}

/* ************************************************************************* *
 * メンバアクセス関数
 * ************************************************************************* */

/* 行列の行数 ************************************************************** */
int
Matrix::Row (void) const {
  return row;
}

/* 列列の行数 ************************************************************** */
int
Matrix::Column (void) const {
  return col;
}

/* n列のベクトルを取り出す ************************************************* */
Vector
Matrix::operator() (int	index) const {
  Vector	v(row);

  if (index >= 0 && index < col) {
    for (int n = 0; n < row; n++) v[n] = data_ptr[n][index];
    return v;    
  } else {
    fprintf (stderr, "Specified index is out of range.\n");
    return 0;
  }
}

/* n行のベクトルを取り出す ************************************************* */
double*
Matrix::operator[] (int	index) {
  if (index >= 0 && index < row) {
    return data_ptr[index];
  } else {
    fprintf (stderr, "Specified index is out of range.\n");    
    return 0;
  }
}

/* n行のベクトルを取り出す ************************************************* */
const double*
Matrix::operator[] (int	index) const {
  if (index >= 0 && index < row) {
    return data_ptr[index];
  } else {
    fprintf (stderr, "Specified index is out of range.\n");    
    return 0;
  }
}

/* n行のベクトルを取り出す ************************************************* */
Vector
Matrix::RowVector (int	index) const {
  Vector	v(col);
  
  if (index >= 0 && index < row) {
    for (int n = 0; n < col; n++) v[n] = data_ptr[index][n];
    return v;
  } else {
    fprintf (stderr, "Specified index is out of range.\n");    
    return 0;
  }
}

/* n列のベクトルを取り出す ************************************************* */
Vector
Matrix::ColumnVector (int	index) const {
  Vector	v(row);

  if (index >= 0 && index < col) {
    for (int n = 0; n < row; n++) v[n] = data_ptr[n][index];
    return v;    
  } else {
    fprintf (stderr, "Specified index is out of range.\n");
    return 0;
  }
}

/* ************************************************************************* *
 *　公開メソッド
 * ************************************************************************* */

/* コピーメソッド ********************************************************** */
void
Matrix::Copy (Matrix&	dst,
	      int	start_row,
	      int	start_col,
	      int	_row,
	      int	_col) const {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  int end_row = start_row + _row - 1;
  int end_col = start_col + _col - 1;
  
  if (start_row < 0 || end_row >= row || start_col < 0 || end_col >= col) {
    fprintf (stderr, "Specified index is out of range.\n");
    return;
  }
  for (int n = 0, i = start_row; n < _row; n++, i++) {
    for (int m = 0, j = start_col; m < _col; m++, j++) {
      dst.data_ptr[n][m] = data_ptr[i][j];
    }
  }
}

/* コピーメソッド ********************************************************** */
void
Matrix::ReferenceMatrix (Matrix&	ref,
			 int		start_row,
			 int		start_col,
			 int		_row,
			 int		_col) const {
  if (_row < 1 || _col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return;
  }
  int end_row = start_row + _row - 1;
  int end_col = start_col + _col - 1;
  
  if (start_row < 0 || end_row >= row || start_col < 0 || end_col >= col) {
    fprintf (stderr, "Specified index is out of range.\n");
    return;
  }
  if (ref.data) delete[] ref.data;
  if (ref.data_ptr) delete[] ref.data_ptr;
  
  ref.row  = _row;
  ref.col  = _col;
  ref.data = 0;
  ref.nelements = _row * _col;
  ref.reference = true;  
  ref.data_ptr = new double_ptr [row];
  
  for (int n = 0, i = start_row; n < _row; n++, i++) {
    ref.data_ptr[n] = data_ptr[i] + start_col;
  }
}

/* 次元の再設定 ************************************************************ */
void
Matrix::Initialize (int	_row,
		    int	_col) {
  if (!reference) {
    if (data) delete[] data;
    if (data_ptr) delete[] data_ptr;

    row = _row;
    col = _col;
    nelements = row * col;

    data = new double [_row * _col];
    data_ptr = new double_ptr [_row];
    for (int n = 0; n < row; n++) data_ptr[n] = data + n * col;
    for (int n = 0; n < nelements; n++) data[n] = 0.0;
  } else {
    fprintf (stderr, "Can not apply Initialize for reference matrix.\n");
  }
}

/* 全ての要素を指定した値でうめるメソッド ********************************** */
void
Matrix::Fill (double val) {
  if (!reference) {
    for (int n = 0; n < nelements; n++) *(data + n) = val;
  } else {
    for (int n = 0; n < nelements; n++) data_ptr[n / col][n % col] = val;
  }
}

/* 全ての要素を0でうめるメソッド ******************************************* */
void
Matrix::Clear (void) {
  if (!reference) {
    for (int n = 0; n < nelements; n++) *(data + n) = 0.0;
  } else {
    for (int n = 0; n < nelements; n++) data_ptr[n / col][n % col] = 0.0;
  }
}

/* 単位行列 **************************************************************** */
void
Matrix::Unit (void) {
  for (int n = 0; n < row; n++) {
    for (int m = 0; m < col; m++) {
      data_ptr[n][m] = (n == m) ? 1.0 : 0.0;
    }
  }
}

/* 対角行列 **************************************************************** */
void
Matrix::Diagonal (double	a, ...) {
  va_list	ptr;

  for (int n = 0; n < nelements; n++) *(data + n) = 0.0;    
  data_ptr[0][0] = a;
  va_start (ptr, a);
  for (int n = row + 1; n < nelements; n+= (row + 1)) {
    *(data + n) = (double) va_arg (ptr, double);
  }
  va_end (ptr);
}

/* ************************************************************************* */
void
Matrix::Set (double	a, ...) {
  va_list	ptr;

  data_ptr[0][0] = a;
  va_start (ptr, a);
  
  if (!reference) {
    for (int n = 1; n < nelements; n++) {
      *(data + n) = (double) va_arg (ptr, double);
    }
  } else {
    for (int n = 1; n < nelements; n++) {
      data_ptr[n / col][n % col] = (double) va_arg (ptr, double);
    }
  } 
  va_end (ptr);
}

/* ************************************************************************* */
void
Matrix::Set (double	*a) {
  if (!reference) {
    for (int n = 0; n < nelements; n++) data[n] = a[n];
  } else {
    for (int n = 0; n < nelements; n++) data_ptr[n / col][n % col] = a[n];
  }
}

/* ************************************************************************* */
void
Matrix::Set (double	**a) {
  for (int n = 0; n < row; n++) {
    for (int m = 0; m < col; m++) data_ptr[n][m] = a[n][m];
  }
}

/* ************************************************************************* */
void
Matrix::Set (int	index_row,
	     int	index_col,
	     double	val) {
  if (CheckRange (index_row, index_col)) {
    data_ptr[index_row][index_col] = val;
  } else {
    fprintf (stderr, "Specified index is out of range.\n");
  }
}

/* ************************************************************************* */
void
Matrix::SetColumn (int			index_col,
		   int			start_row,
		   const Vector&	v) {
  if (index_col < 0 || index_col >= col ||
      start_row < 0 || start_row >= row ||
      start_row + v.Dimension() > row) {
    fprintf (stderr, "Specified index is out of range.\n");
  } else {
    for (int n = 0, i = start_row; n < v.Dimension(); n++, i++) {
      data_ptr[i][index_col] = v[n];
    }
  }
}

/* ************************************************************************* */
void
Matrix::SetRow (int		index_row,
		int		start_col,
		const Vector&	v) {
  if (index_row < 0 || index_row >= row ||
      start_col < 0 || start_col >= col ||
      start_col + v.Dimension() > col) {
    fprintf (stderr, "Specified index is out of range.\n");
  } else {
    for (int n = 0, i = start_col; n < v.Dimension(); n++, i++) {
      data_ptr[index_row][i] = v[n];
    }
  }
}

/* 列の入れ換え ************************************************************ */
void
Matrix::SwapColumn (int	index_i,
		    int	index_j) {
  if (index_i < 0 || index_i >= col || index_j < 0 || index_j >= col) {
    fprintf (stderr, "Specified index is out of range.\n");
  } else {
    Vector temp(row);  

    for (int n = 0; n < row; n++) temp[n] = data_ptr[n][index_i];
    for (int n = 0; n < row; n++) data_ptr[n][index_i] = data_ptr[n][index_j];
    for (int n = 0; n < row; n++) data_ptr[n][index_j] = temp[n];
  }
}

/* 行の入れ換え ************************************************************ */
void
Matrix::SwapRow (int	index_i,
		 int	index_j) {
  if (index_i < 0 || index_i >= row || index_j < 0 || index_j >= row) {
    fprintf (stderr, "Specified index is out of range.\n");
  } else {
    Vector temp(row);  

    for (int n = 0; n < col; n++) temp[n] = data_ptr[index_i][n];
    for (int n = 0; n < col; n++) data_ptr[index_i][n] = data_ptr[index_j][n];
    for (int n = 0; n < col; n++) data_ptr[index_j][n] = temp[n];
  }
}

/* ************************************************************************* */
inline bool
Matrix::IsSquare (void) {
  return (row == col);
}

/* ************************************************************************* */
bool
Matrix::IsSymmetric (void) {
  if (row != col) return false;
  for (int n = 0; n < row; n++) {
    for (int m = 0; m < col; m++) {
      if (!IsZero (fabs (data_ptr[n][m] - data_ptr[m][n]))) return false;
    }
  }
  return true;
}

/* 出力メソッド ************************************************************ */
void
Matrix::Print (FILE *fp, const char *format, const char *message) const {
  int	msg_length = 0;

  if (message) msg_length = strlen (message);

  for (int n = 0; n < row; n++) {
    if (message) {
      if (n == ((row -1 ) / 2)) {
	fprintf (fp, "%s| ", message);
      } else {
	for (int l = 0; l < msg_length; l++) fprintf (fp, " ");
	fprintf (fp, "| ");
      }
    }
    for (int m = 0; m < col; m++) fprintf (fp, format, data_ptr[n][m]);
    if (message) fprintf (fp, " |");
    fprintf (fp, "\n");
  }
}

/* ************************************************************************* */
double
Matrix::CoFactor (int		index_i,
		  int		index_j) const {
  if (CheckRange (index_j, index_i)) {
    Matrix B (row-1, col-1);
    for (int i = 0, j = 0; i < nelements; i++) {
      int c = i % col;
      int r = i / row;
      if (c == index_i || r == index_j) continue;
      B[j/2][j%2] = data[i];
      j++;
    }
    double det = Determinant (B);
    return ((index_i + index_j) % 2 == 0) ? det : -det;
  } else {
    fprintf (stderr, "Specified index is out of range.\n");
    return -1;
  }
}

/* ************************************************************************* */
static void
_jacobi_rotate (Matrix&	A,
		int	i,
		int 	j,
		int 	k,
		int 	l,
		double	*g,
		double 	*h,
		double 	s,
		double	tau) {
  *g = A[i][j];
  *h = A[k][l];
  A[i][j] = *g - s * (*h + *g * tau);
  A[k][l] = *h + s * (*g - *h * tau);
}

/* ************************************************************************* */
static bool
_eigen_sort (Vector&	eval,
	     Matrix&	evec,
	     int	sort_type) {
  int	N = eval.Dimension();

  for (int i = 0; i < N - 1; i++) {
    int	k = i;
    double ek = eval[i];

    for (int j = i + 1; j < N; j++) {
      bool	test;
      double	ej = eval[j];
	  
      switch (sort_type) {
      case 0: // Ascending order
	test = (ej < ek);
	break;
      case 2: // Ascending order (use absolute value)
	test = (fabs (ej) < fabs (ek));
	break;
      case 3: // Descending order (use absolute value)
	test = (fabs (ej) > fabs (ek));
	break;
      case 1: // Descending order
      default:
	test = (ej > ek);
	break;
      }
      if (test) {
	k = j;
	ek = ej;
      }
    }
    if (k != i) {
      eval.Swap (i, k);
      evec.SwapColumn (i, k);
    }
  }
  return true;
}

/* LU分解 ****************************************************************** */
bool
Matrix::LUDecomp (Matrix&	LU,
		  int		*permutation,
		  double	*determinant) const {
  if (row != col) {
    fprintf (stderr, "Matrix size must be square.\n");
    return false;
  }
#if 0
  int		N = row;
  int		i, j, k, ii, ik;
  double	t, u;
  double	*weight = new double[N];
	
  LU = *this;
  *determinant = 0;

  for (i = 0; i < N; i++) permutation[i] = i;

  for (k = 0; k < N; k++) {
    u = 0;
    permutation[k] = k;
    for (j = 0; j < N; j++) {
      t = fabs (LU[k][j]);
      if (t > u) u = t;
    }
    if (u == 0) {
      delete[] weight;
      return false;
    }
    weight[k] = 1.0 / u;
  }
  *determinant = 1;

  for (k = 0; k < N; k++) {
    u = -1;
    for (i = k; i < N; i++) {
      ii = permutation[i];
      t = fabs (LU[ii][k]) * weight[ii];
      if (t > u) {
	u = t;
	j = i;
      }
    }
    ik = permutation[j];
    if (j != k) {
      permutation[j] = permutation[k];
      permutation[k] = ik;
      *determinant = -(*determinant);
    }
    u = LU[ik][k];
    *determinant *= u;
    if (u == 0) {
      delete[] weight;
      return false;
    }
    for (i = k + 1; i < N; i++) {
      ii = permutation[i];
      t = (LU[ii][k] /= u);
      for (j = k + 1; j < N; j++) LU[ii][j] -= t * LU[ik][j];
    }
  }
  delete[] weight;
  
  return true;
#else
  int		N = row;
  int		i, j, k;

  LU = *this;
  *determinant = 1;

  for (i = 0; i < N; i++) permutation[i] = i;

  for (j = 0; j < N - 1; j++) {
    double	ajj;
    double	max = fabs (LU[j][j]);
    int		i_pivot = j;
    
    for (i = j + 1; i < N; i++) {
      double	aij = fabs (LU[i][j]);
      if (aij > max) {
	max = aij;
	i_pivot = i;
      }
    }
    if (i_pivot != j) {
      int tmp;
      LU.SwapRow (j, i_pivot);
      tmp = permutation[j];
      permutation[j] =  permutation[i_pivot];
      permutation[i_pivot] = tmp;
      *determinant = -(*determinant);
    }
    ajj = LU[j][j];

    if (ajj != 0.0) {
      for (i = j + 1; i < N; i++) {
	double aij = LU[i][j] / ajj;
	LU[i][j] = aij;
	for (k = j + 1; k < N; k++) {
	  double	aik = LU[i][k];
	  double	ajk = LU[j][k];
	  LU[i][k] = aik - aij * ajk;
	}
      }
    }
  }
  return true;
#endif
}

/* ************************************************************************* *
 * オペレータ
 * ************************************************************************* */

/* 行列の代入 ************************************************************** */
Matrix&
Matrix::operator= (const Matrix& a) {
  if (this != &a) {
    if (row == a.row && col == a.col) {
      if (!reference) {
	for (int n = 0; n < nelements; n++) data[n] = a.data[n];
      } else {
	for (int n = 0; n < nelements; n++) {
	  data_ptr[n / col][n % col] = a.data_ptr[n / col][n % col];
	}
      }
    } else {
      fprintf (stderr, "Invarid dimension.\n");
    }
  }
  return *this;
}

/* 行列の加算 ************************************************************** */
Matrix&
Matrix::operator+= (const Matrix&	a) {
  if (row == a.row && col == a.col) {
    if (!reference) {
      for (int n = 0; n < nelements; n++) data[n] += a.data[n];
    } else {
      for (int n = 0; n < nelements; n++) {
	data_ptr[n / col][n % col] += a.data_ptr[n / col][n % col];
      }
    }
  } else {
    fprintf (stderr, "Dimensions of two matrices do not match.\n");
  }	
  return *this;
}

/* 行列の減算 ************************************************************** */
Matrix&
Matrix::operator-= (const Matrix&	a) {
  if (row == a.row && col == a.col) {
    if (!reference) {
      for (int n = 0; n < nelements; n++) data[n] -= a.data[n];
    } else {
      for (int n = 0; n < nelements; n++) {
	data_ptr[n / col][n % col] -= a.data_ptr[n / col][n % col];
      }
    }
  } else {
    fprintf (stderr, "Dimensions of two matrices do not match.\n");    
  }
  return *this;
}

/* 行列要素の乗算 ********************************************************** */
Matrix&
Matrix::operator*= (double	val) {
  if (!reference) {
    for (int n = 0; n < nelements; n++) data[n] *= val;
  } else {
    for (int n = 0; n < nelements; n++) {
      data_ptr[n / col][n % col] *= val;
    }
  }
  return *this;
}

/* 行列要素の除算 ********************************************************** */
Matrix&
Matrix::operator/= (double	val) {
  if (!IsZero (val)) {
    if (!reference) {
      for (int n = 0; n < nelements; n++) data[n] /= val;
    } else {
      for (int n = 0; n < nelements; n++) {
	data_ptr[n / col][n % col] /= val;
      }
    }
  } else {
    fprintf (stderr, "Specified value is nearly zero.\n");    
  }	
  return *this;
}

/* 行列の加算 ************************************************************** */
Matrix
operator+ (const Matrix&	a,
	   const Matrix&	b) {
  Matrix	c(a.row, a.col);  

  if (a.row == b.row && a.col == b.col) {
    if (!a.reference && !b.reference) {
      for (int n = 0; n < a.nelements; n++) c.data[n] = a.data[n] + b.data[n];
    } else {
      for (int n = 0; n < a.nelements; n++) {
	c.data_ptr[n / c.col][n % c.col] =
	  a.data_ptr[n / a.col][n % a.col] + b.data_ptr[n / b.col][n % b.col];
      }
    }
  } else {
    fprintf (stderr, "Dimensions of two matrices do not match.\n");
  }
  return c;  
}

/* 行列の減算 ************************************************************** */
Matrix
operator- (const Matrix&	a) {
  Matrix	b(a.row, a.col);  

  if (!a.reference) {
    for (int n = 0; n < a.nelements; n++) b.data[n] = -a.data[n];
  } else {
    for (int n = 0; n < a.nelements; n++) {
      b.data_ptr[n / b.col][n % b.col] = -a.data_ptr[n / a.col][n % a.col];
    }
  }
  return b;
}

/* 行列の減算 ************************************************************** */
Matrix
operator- (const Matrix&	a,
	   const Matrix&	b) {
  Matrix	c(a.row, a.col);

  if (a.row == b.row && a.col == b.col) {
    if (!a.reference && !b.reference) {    
      for (int n = 0; n < a.nelements; n++) c.data[n] = a.data[n] - b.data[n];
    } else {
      for (int n = 0; n < a.nelements; n++) {
	c.data_ptr[n / c.col][n % c.col] =
	  a.data_ptr[n / a.col][n % a.col] - b.data_ptr[n / b.col][n % b.col];
      }
    }
  } else {
    fprintf (stderr, "Invalid dimension.\n");
  }
  return c;
}

/* 行列の乗算 ************************************************************** */
Matrix
operator* (const Matrix&	a,
	   double		val) {
  Matrix	b(a.row, a.col);

  if (!a.reference) {
    for (int n = 0; n < a.nelements; n++) b.data[n] = a.data[n] * val;
  } else {
    for (int n = 0; n < a.nelements; n++) {
      b.data_ptr[n / b.col][n % b.col] = a.data_ptr[n / a.col][n % a.col]*val;
    }
  }
  return b;
}

/* 行列の乗算 ************************************************************** */
Matrix
operator* (double		val,
	   const Matrix&	a) {
  Matrix	b(a.row, a.col);

  if (!a.reference) {  
    for (int n = 0; n < a.nelements; n++) b.data[n] = a.data[n] * val;
  } else {
    for (int n = 0; n < a.nelements; n++) {
      b.data_ptr[n / b.col][n % b.col] = a.data_ptr[n / b.col][n % b.col]*val;
    }
  }
  return b;
}

/* 行列の乗算 ************************************************************** */
Vector
operator* (const Matrix&	a,
	   const Vector&	v) {
  Vector	b(a.row);

  if (a.col == v.Dimension()) {
    for (int n = 0; n < a.row; n++) {
      for (int m = 0; m < a.col; m++) {
	b[n] += a.data_ptr[n][m] * v[m];
      }
    }
  } else {
    fprintf (stderr, "Dimensions of matrix and vector do not match.\n");
  }
  return b;
}

/* 行列の乗算 ************************************************************** */
Matrix
operator* (const Matrix&	a,
	   const Matrix&	b) {
  Matrix	c(a.row, b.col);

  if (a.col == b.row) {
    for (int n = 0; n < c.row; n++) {
      for (int m = 0; m < c.col; m++) {
	for (int l = 0; l < a.col; l++) { 
	  c.data_ptr[n][m] += a.data_ptr[n][l] * b.data_ptr[l][m];
	}
      }
    }
  } else {
    fprintf (stderr, "Dimensions of two matrices do not match.\n");
  }
  return c;
}

/* 行列の除算 ************************************************************** */
Matrix
operator/ (const Matrix&	a,
	   double		val) {
  Matrix	b(a.row, a.col);

  if (!IsZero (val)) {
    if (!a.reference) {
      for (int n = 0; n < a.nelements; n++) b.data[n] = a.data[n] / val;
    } else {
      for (int n = 0; n < a.nelements; n++) {
	b.data_ptr[n / b.col][n % b.col]
	  = a.data_ptr[n / a.col][n % a.col] / val;
      }
    }
  } else {
    fprintf (stderr, "Specified value is nearly zero.\n");    
  }
  return b;
}

/* 行列の外積 ************************************************************** */
Matrix
operator% (const Matrix&	a,
	   const Vector&	v) {
  return Transpose (v % Transpose (a));  
}

/* 行列の外積 ************************************************************** */
Matrix
operator% (const Vector&	v,
	   const Matrix&	a) {
  Matrix	b(a.row, a.col);

  if (v.Dimension() < 3 || a.row < 3) {
    fprintf (stderr,
	     "Vector dimension and row of matrix must be greater "
	     "than three.\n");
  } else if (v.Dimension() != a.row) {
    fprintf (stderr, "Vector dimension and row of matrix do not match.\n");
  } else {
    for (int n = 0; n < a.col; n++) b.SetColumn (n, 0, v % a(n));
  }
  return b;
}

/* 行列の内積 ************************************************************** */
double
operator, (const Matrix&	a,
	   const Matrix&	b) {
  double	val = 0.0;
  
  if (a.row == b.row && a.col || b.col) {

    if (!a.reference && !b.reference) {
      for (int n = 0; n < a.nelements; n++) val += a.data[n] * b.data[n];
    } else {
      for (int n = 0; n < a.nelements; n++) {
	val += a.data_ptr[n / a.col][n % a.col]
	  * b.data_ptr[n / b.col][n % b.col];
      }
    }
  } else {
    fprintf (stderr, "Dimensions of two matrices do not match.\n");
  }
  return val;
}

/* 行列の比較 ************************************************************** */
bool
operator== (const Matrix&	a,
	    const Matrix&	b) {
  if (a.row == b.row && a.col == b.col && IsZero (fabs (Norm (a - b)))) {
    return true;
  } else {
    return false;
  }
}

/* 行列の比較 ************************************************************** */
bool
operator!= (const Matrix&	a,
	    const Matrix&	b) {
  if (a.row == b.row && a.col == b.col && IsZero (fabs (Norm (a - b)))) {
    return false;
  } else {
    return true;
  }
}

/* 行列のノルム ************************************************************ */
double
Norm (const Matrix&	a) {
  double	norm = 0.0;

  if (!a.reference) {
    for (int n = 0; n < a.nelements; n++) norm += a.data[n] * a.data[n];
  } else {
    for (int n = 0; n < a.nelements; n++) {
      double val = a.data_ptr[n / a.col][n % a.col];
      norm += val * val;
    }
  }
  return sqrt (norm);
}

/* ************************************************************************* */
Matrix
Normalize (const Matrix&	a) {
  Matrix	b(a.row, a.col);
  double	norm = Norm (a);

  if (!a.reference) {
    for (int n = 0; n < a.nelements; n++) b.data[n] = a.data[n] / norm;
  } else {
    for (int n = 0; n < a.nelements; n++) {
      b.data_ptr[n / b.col][n % b.col]
	= a.data_ptr[n / a.col][n % a.col] / norm;
    }
  }
  return b;
}

/* ************************************************************************* */
Matrix
Transpose (const Matrix&	a) {
  Matrix	b (a.col, a.row);

  for (int n = 0; n < a.row; n++) {
    for (int m = 0; m < a.col; m++) b.data_ptr[m][n] = a.data_ptr[n][m];
  }
  return b;
}

/* ************************************************************************* */
Matrix
ProjectionMatrix (const Vector&	v) {
  Matrix	I(v.Dimension(), v.Dimension());

  I.Unit();
  return I - Matrix (Normalize (v) , Normalize (v));
}

/* ************************************************************************* */
double
Trace (const Matrix&	a) {
  double	trace = 0.0;

  if (a.row == a.col) {
    for (int n = 0; n < a.row; n++) trace += a.data_ptr[n][n];
  } else {
    fprintf (stderr, "Matrix size must be square.\n");
  }
  return trace;
}


/* ************************************************************************* */
double
Determinant (const Matrix&	a) {
  if (a.row != a.col) {
    fprintf (stderr, "Matrix size must be square.\n");
    return 0.0;
  }
  int		i, j, k, l;
  double	det = 0;
  Matrix	temp;

  if (a.Row() != a.Column()) return 0;
  if (a.Row() == 0) return 0;
  if (a.Row() == 1) return a[0][0];

  temp.Initialize (a.Row() - 1, a.Column() - 1);
  
  /*  第1列に関する展開  */
  for(i = 0; i < a.Row(); i++) {
    if (a[i][0] != 0) {    
      l = 0;
      for (j = 0; j < a.Row(); j++) {
	if (i != j) {
	  for (k = 1; k < a.Column(); k++) {
	    temp[l / temp.Column()][l % temp.Column()] = a[j][k];
	    l++;
	  }
	}
      }
      /*  iが偶数か奇数か．(-1)^{i-1}  */
      if ((i + 1) % 2 == 1) {
	det += a[i][0] * Determinant (temp);
      } else {
	det -= a[i][0] * Determinant(temp);
      }
    }
  }
  return det;
}

/* ************************************************************************* */
int
Rank (const Matrix&	a) {
  Matrix	evec;
  Vector	eval;
  int		rank = 0;
  
  (a * Transpose (a)).Eigen (eval, evec);
  for (int n = 0; n < eval.Dimension(); n++) {
    if (!IsZero (eval[n])) rank++;
  }
  return rank;
}

/* 逆行列 ****************************************************************** */
Matrix
Inverse (const Matrix&	a) {
  Matrix	Inv (a);
  int		i, j, k, N = a.row;
  double	t, u, determinant;
  
  if (a.row != a.col) {
    fprintf (stderr, "Matrix size must be square.\n");
    return Inv;
  }

  determinant = 1;
  for (k = 0; k < N; k++) {
    t = Inv[k][k];
    determinant *= t;
    for (i = 0; i < N; i++) Inv[k][i] /= t;
    Inv[k][k] = 1 / t;
    for (j = 0; j < N; j++)
      if (j != k) {
	u = Inv[j][k];
	for (i = 0; i < N; i++)
	  if (i != k) {
	    Inv[j][i] -= Inv[k][i] * u;
	  } else {
	    Inv[j][i] = -u / t;
	  }
      }
  }
  return Inv;
}

/* ************************************************************************* */
Matrix
GeneralInverse (const Matrix&	a,
		int		rank) {
  if (a.row != a.col) {
    fprintf (stderr, "Matrix size must be square.\n");
    return Matrix (a);
  }
  Matrix	evec, Inv (a.row, a.col);
  Vector	eval;

  if (rank == 0) rank = Rank (a);

  a.Eigen (eval, evec);
  for (int n = 0; n < rank; n++) {
    if (!IsZero (eval[n]))  Inv += 1.0 / eval[n] * Matrix (evec(n), evec(n));
  }
  return Inv;
}

/* 対角要素を取り出したベクトルを生成 ************************************** */
Vector
Diagonal (const Matrix&	a) {
  int		N = (a.row < a.col) ? a.row : a.col;
  Vector	v(N);

  for (int n = 0; n < a.nelements; n++) v[n] = a[n][n];

  return v;
}

/* 与えたベクトル要素を対角成分とする対角行列を生成 ************************ */
Matrix
Diagonal (const Vector&	v) {
  Matrix	a(v.Dimension(), v.Dimension());
  
  for (int n = 0; n < v.Dimension(); n++) a[n][n] = v[n];
  return a;
}

/* ************************************************************************* */
Vector
ToVector (const Matrix&	a) {
  if (a.row < 1 || a.col < 1) {
    fprintf (stderr, "Matrix dimension must be greater than zero.\n");
    return Vector (1);
  }
  Vector	v(a.row * a.col);

  for (int n = 0, l = 0; n < a.row; n++) {
    for (int m = 0; m < a.col; m++) {
      v[l++] = a[n][m];
    }
  }
  return v;
}

/* ************************************************************************* */
bool
Matrix::CheckRange (int	index_row,
		    int	index_col) const {
  if (index_row >= 0 && index_row < row &&
      index_col >= 0 && index_col < col) {
    return true;
  } else {
    return false;
  }
}

/* *************************************************** End of matrix.cpp *** */

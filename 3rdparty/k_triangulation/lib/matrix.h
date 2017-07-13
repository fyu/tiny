/* ************************************************************ matrix.h *** *
 * 行列クラス ヘーダファイル
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <2008-03-07 13:50:30 sugaya>
 * ************************************************************************* */
#ifndef	__MATRIX_H__
#define	__MATRIX_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "vector.h"

typedef double	*double_ptr;

/* 行列クラス ************************************************************** */
class Matrix {
 private:
  int		row;
  int		col;
  int		nelements;
  bool		reference;
  double	*data;
  double	**data_ptr;

 public:

  /* コンストラクタ & デストラクタ */
  Matrix	(void);
  Matrix	(int _row, int _col);
  Matrix	(int _row, int _col, double a, ...);
  Matrix	(int _row, int _col, double *a);
  Matrix	(int _row, int _col, double **a);
  Matrix	(const Matrix& a);
  Matrix	(const Vector& v, int _row, int _col);  
  Matrix	(const Matrix& a,
		 int start_row, int start_col, int _row, int _col);
  Matrix	(const Vector& a, const Vector& b);
  Matrix	(const Vector& a, double angle);

  ~Matrix	(void);  

  /* メンバアクセスメソッド */
  int		Row		(void) const;
  int		Column		(void) const;
  double*	operator[]	(int 		index);
  const double*	operator[]	(int 		index) const;
  Vector	operator() 	(int		index) const;
  Vector	RowVector	(int		index) const;
  Vector	ColumnVector	(int		index) const;  
    
  /* 公開メソッド */
  void		Copy		(Matrix&	dst,
				 int 		start_row,
				 int		start_col,
				 int		_row,
				 int		_col) const;
  void		ReferenceMatrix	(Matrix&	ref,
				 int 		start_row,
				 int		start_col,
				 int		_row,
				 int		_col) const;
  void		Initialize	(int 		_row,
				 int		_col);
  void		Fill		(double 	val);
  void		Clear		(void);
  void		Unit		(void);
  void		Diagonal	(double		a, ...);
  void		Set		(double		a, ...);
  void		Set		(double		*a);
  void		Set		(double		**a);
  void		Set		(int		index_row,
				 int		index_col,
				 double		val);
  void		SetColumn	(int		index_col,
				 int		start_row,
				 const Vector&	v);
  void		SetRow		(int		index_row,
				 int		start_col,
				 const Vector&	v);
  void		SwapColumn	(int		index_i,
				 int		index_j);
  void		SwapRow		(int		index_i,
				 int		index_j);
  bool		IsSquare	(void);
  bool		IsSymmetric	(void);
  void		Print		(FILE		*fp,
				 const char 	*format,
				 const char 	*message = 0) const;
  double	CoFactor 	(int		index_i,
				 int		index_j) const ;
 
  bool		LUDecomp	(Matrix&	LU,
				 int		*permutation,
				 double		*determinant) const;
  bool		Eigen	 	(Vector		&eval,
				 Matrix		&evec,
				 double		eps = 1.0e-15) const;
  bool		SvDecomp 	(Matrix&	U,
				 Matrix&	W,
				 Matrix&	V) const;
  
  /* オペレータ */
  Matrix&	operator=	(const Matrix&	a);
  Matrix&	operator+= 	(const Matrix&	a);
  Matrix&	operator-= 	(const Matrix&	a);
  Matrix&	operator*= 	(double		val);
  Matrix&	operator/= 	(double		val);
  friend Matrix	operator+	(const Matrix& a, const Matrix& b);
  friend Matrix	operator-	(const Matrix& a);
  friend Matrix	operator-	(const Matrix& a, const Matrix& b);
  friend Matrix	operator*	(const Matrix& a, double val);
  friend Matrix	operator*	(double val, const Matrix& a);  
  friend Vector	operator* 	(const Matrix& a, const Vector&	v);
  friend Matrix	operator* 	(const Matrix& a, const Matrix&	b);
  friend Matrix	operator/	(const Matrix& a, double val);
  friend Matrix	operator% 	(const Matrix& a, const Vector&	v);
  friend Matrix	operator% 	(const Vector& v, const Matrix&	a);
  //friend matrix	operator% 	(const Matrix& a, const Matrix&	b);
  friend double operator,	(const Matrix& a, const Matrix& b);
  friend bool	operator== 	(const Matrix& a, const Matrix& b);
  friend bool	operator!= 	(const Matrix& a, const Matrix& b);
    
  /* フレンド関数 */    
  friend double	Norm		(const Matrix&	a);
  friend Matrix	Normalize	(const Matrix&	a);
  friend Matrix	Transpose 	(const Matrix&	a);
  friend Matrix ProjectionMatrix(const Vector&	v);
  friend double	Trace 		(const Matrix&	a);
  friend double	Determinant	(const Matrix&	a);
  friend int	Rank 		(const Matrix&	a);
  friend Matrix	Inverse		(const Matrix&	a);
  friend Matrix	GeneralInverse	(const Matrix&	a,
				 int		rank = 0);
  friend Vector Diagonal	(const Matrix&	a);
  friend Matrix Diagonal	(const Vector&	v);
  friend Vector ToVector	(const Matrix&	a);
  bool		CheckRange	(int		index_row,
				 int		index_col) const;
  friend int	HouseholderSolve(const Matrix&	a,
				 const Vector&	b,
				 Vector&	x);		

};

#endif /* __MATRIX_H__ */

/* ***************************************************** End of matrix.h *** */

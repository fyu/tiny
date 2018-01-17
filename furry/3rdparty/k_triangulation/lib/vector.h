/* ************************************************************ vector.h *** *
 * ベクトルクラス ヘーダファイル
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <06/12/27 22:26:35 sugaya>
 * ************************************************************************* */
#ifndef	__VECTOR_H__
#define	__VECTOR_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>

/* ベクトルクラス ********************************************************** */
class Vector {
 private:
  int		dimension;
  bool		reference;
  double	*data;
 public:

  /* コンストラクタ & デストラクタ */
  Vector	(void);
  Vector	(int _dimension);
  Vector	(int _dimension, double a, ...);
  Vector	(int _dimension, double *a);
  Vector	(const Vector&	v);
  Vector 	(const Vector& v, int start_index, int _dimension);
  
  ~Vector	(void);

  /* メンバアクセスメソッド */
  int		Dimension	(void) const;
  double&	operator[]	(int 		index);
  const double&	operator[]	(int 		index) const;  
  
  /* 公開メソッド */
  void		Copy		(Vector&	dst,
				 int 		start_index,
				 int 		_dimension) const;
  void		ReferenceVector	(Vector&	ref,
				 int 		start_index,
				 int 		_dimension) const;
  void		Initialize	(int 		_dimension);
  void		Fill		(double 	val);
  void		Clear		(void);
  void		Set		(double		a, ...);
  void		Set		(double		*a);
  void		Set		(int 		index,
				 double 	val);
  void		Swap		(int 		index_i,
				 int 		index_j);
  void		Print		(FILE		*fp,
				 const char 	*format,
				 const char 	*message = 0) const;

  /* オペレータ */
  Vector&	operator=	(const Vector&	v);
  Vector&	operator+=	(const Vector&	v);
  Vector&	operator-=	(const Vector&	v);
  Vector&	operator*=	(double	val);
  Vector&	operator/=	(double	val);
  friend Vector operator+	(const Vector& a, const Vector& b);
  friend Vector operator-	(const Vector& v);
  friend Vector operator-	(const Vector& a, const Vector& b);
  friend Vector operator*	(const Vector& v, double x);
  friend Vector operator*	(double x, const Vector& v);
  friend Vector operator/	(const Vector& v, double x);
  friend double operator,	(const Vector& 	a, const Vector& b);
  friend Vector	operator%	(const Vector&	a, const Vector& b);
  friend bool	operator== 	(const Vector&	a, const Vector& b);
  friend bool	operator!= 	(const Vector&	a, const Vector& b);

  /* フレンド関数 */  
  friend double	Norm		(const Vector&	v);
  friend Vector	Normalize	(const Vector&	v);
  friend Vector	ZOperator	(const Vector&	v);
  friend double	TripleProduct	(const Vector&	a,
				 const Vector&	b,
				 const Vector&	c);

  bool		CheckRange	(int		index) const;
};

#endif /* __VECTOR_H__ */

/* ***************************************************** End of vector.h *** */

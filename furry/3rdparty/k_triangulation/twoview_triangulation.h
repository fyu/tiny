/* ********************************************* twoview_triangulation.h *** *
 * 2画像からの最適三角測量に関する関数
 *
 * Copyright (C) 2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/09 17:32:35 sugaya>
 * ************************************************************************* */
#ifndef	__TWOVIEW_TRIANGULATION_H__
#define	__TWOVIEW_TRIANGULATION_H__

#include "lib/libmatrix.h"

void calc_fundamental_matrix (const Matrix&	P1,
			      const Matrix&	P2,
			      Matrix&		F);

double two_view_triangulation (const Vector&	x1,
			       const Vector&	x2,
			       const Matrix&	F,
			       Vector&		_x1,
			       Vector&		_x2,
			       double		f0,
			       int		iter_max,
			       double		convergence);

double two_view_triangulation_fast (const Vector&	x1,
				    const Vector&	x2,
				    const Matrix&	F,
				    Vector&		_x1,
				    Vector&		_x2,
				    double		f0,
				    int			iter_max,
				    double		convergence);

#endif	/* __TWOVIEW_TRIANGULATION_H__ */

/* ************************************** End of twoview_triangulation.h *** */

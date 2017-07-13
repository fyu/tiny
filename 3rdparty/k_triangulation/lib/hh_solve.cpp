/* ********************************************************** hh_solve.c *** *
 * ハウスホルダー変換により線形問題を解く関数
 *
 * the original source code comes from gsl-1.6/linalg/hh.c
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <06/11/13 14:39:03 sugaya>
 * ************************************************************************* */

/* ************************************************************************* *
 * linalg/hh.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, #Include <., 675 Mass Ave, Cambridge, MA 02139, USA.
 * ************************************************************************* */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "macros.h"
#include "matrix.h"

/* ************************************************************************* */
#define _SIGN(a) (((a) == 0) ? 0 : (((a) > 0) ? 1 : -1))

/* ************************************************************************* */
static int
_gsl_linalg_HH_svx (Matrix		&A,
		    Vector		&x) {
  const int	N = A.Row();
  const int	M = A.Column();
  int		i, j, k;
  double	*d = new double [N];

  for (i = 0; i < N; i++) {
    const double aii = A[i][i];
    double	 alpha, f, ak;
    double	 max_norm = 0.0;
    double	 r = 0.0;

    for (k = i; k < M; k++) {
      double aki = A[k][i];
      r += aki * aki;
    }
    if (r == 0) {
      delete [] d;
      fprintf (stderr, "Matrix is rank deficient.\n");
      return 0;
    }
    alpha = sqrt (r) * _SIGN(aii);
    ak = 1.0 / (r + alpha * aii);
    A[i][i] = aii + alpha;
    d[i] = -alpha;

    for (k = i + 1; k < N; k++) {
      double	norm = 0.0;
      f = 0.0;
      for (j = i; j < M; j++) {
	double	ajk = A[j][k];
	double	aji = A[j][i];
	norm += ajk * ajk;
	f += ajk * aji;
      }
      max_norm = (max_norm < norm) ? norm : max_norm;
      f *= ak;

      for (j = i; j < M; j++) {
	double	ajk = A[j][k];
	double	aji = A[j][i];
	A[j][k] = ajk - f * aji;
      }
    }
    if (fabs (alpha) < 2.0 * MATH_DBL_EPSILON * sqrt (max_norm)) {
      delete [] d;
      fprintf (stderr, "Apparent sigularity detected.\n");
      return 0;
    }
    f = 0.0;
    for (j = i; j < M; j++) f += x[j] * A[j][i];
    f *= ak;
    for (j = i; j < M; j++) {
      double	xj = x[j];
      double	aji = A[j][i];
      x[j] = xj - f * aji;
    }
  }
  for (i = N; i > 0 && i--;) {
    double	xi = x[i];
    double	sum = 0.0;
    for (k = i + 1; k < N; k++) sum += A[i][k] * x[k];
    x[i] = (xi - sum) / d[i];
  }
  delete [] d;
  return 1;
}

/* ************************************************************************* */
int
HouseholderSolve (const Matrix&	a,
		  const Vector&	b,
		  Vector&	x) {
  if (a.row > a.col) {
    fprintf (stderr, "System is undertermined.\n");
    return 0;
  } else if (a.col != x.Dimension()) {
    fprintf (stderr, "Matrix and vector sizes must be equal.\n");
    return 0;
  }
  Matrix	A(a);
  x = b;
  return _gsl_linalg_HH_svx (A, x);
}

/* ************************************************** End of hh_solve.c *** */

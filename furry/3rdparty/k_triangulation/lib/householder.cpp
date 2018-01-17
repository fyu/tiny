/* ***************************************************** householder.cpp *** *
 *
 *                                    Time-stamp: <06/11/13 11:11:14 sugaya>
 * ************************************************************************* */

/* ************************************************************************* *
 * linalg/householder.
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
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
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 * ************************************************************************* */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "matrix.h"

/* ************************************************************************* */
double
gsl_linalg_householder_transform (Vector&	v) {
  /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */

  const int n = v.Dimension();

  if (n == 1) {
    return 0.0; /* tau = 0 */
  } else {
    double alpha, beta, tau;
    Vector x;
    v.ReferenceVector (x, 1, n - 1);
    double xnorm = Norm (x);
    
    if (xnorm == 0) return 0.0; /* tau = 0 */
      
    alpha = v[0];
    beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
    tau = (beta - alpha) / beta;

    x *= (1.0 / (alpha - beta));
    v[0] = beta;
    
    return tau;
  }
}

/* ************************************************************************* */
bool
gsl_linalg_householder_hm (double		tau,
			   const Vector&	v,
			   Matrix&		A) {
  /* applies a householder transformation v,tau to matrix m */

  int i, j;

  if (tau == 0.0) return true;

  for (j = 0; j < A.Column(); j++) {
    double wj = A[0][j];
    for (i = 1; i < A.Row(); i++) wj += A[i][j] * v[i];
    double A0j = A[0][j];
    A[0][j] = A0j - tau * wj;
    for (i = 1; i < A.Row(); i++) {
      double Aij = A[i][j];
      double vi = v[i];
      A[i][j] = Aij - tau * vi * wj;
    }
  }
  return true;
}

/* ************************************************************************* */
bool
gsl_linalg_householder_mh (double		tau,
			   const Vector&	v,
			   Matrix&		A) {
  /* applies a householder transformation v,tau to matrix m from the
     right hand side in order to zero out rows */

  int i, j;

  if (tau == 0) return true;

  for (i = 0; i < A.Row(); i++) {
    double wi = A[i][0];
    for (j = 1; j < A.Column(); j++) wi += A[i][j] * v[j];
    double Ai0 = A[i][0];
    A[i][0] = Ai0 - tau * wi;
    for (j = 1; j < A.Column(); j++) {
      double vj = v[j];
      double Aij = A[i][j];
      A[i][j] = Aij - tau * wi * vj;
    }
  }
  return true;
}

/* ************************************************************************* */
bool
gsl_linalg_householder_hv (double		tau,
			   const Vector&	v,
			   Vector&		w) {
  /* applies a householder transformation v to vector w */
  const int N = v.Dimension();
 
  if (tau == 0) return true;

  /* compute d = v'w */

  double d0 = w[0];
  double d1, d;
  
  Vector v1(v, 1, N - 1);
  Vector w1;
  w.ReferenceVector (w1, 1, N - 1);

  d1 = (v1, w1);
  d = d0 + d1;

  /* compute w = w - tau (v) (v'w) */
  
  double w0 = w[0];
  w[0] = w0 - tau * d;
  w1 += v1 * (-tau * d);
       
  return true;
}

/* ************************************************************************* */
bool
gsl_linalg_householder_hm1 (double	tau,
			    Matrix&	A) {
  /* applies a householder transformation v,tau to a matrix being
     build up from the identity matrix, using the first column of A as
     a householder vector */

  int i, j;

  if (tau == 0) {
    A[0][0] = 1.0;
      
    for (j = 1; j < A.Column(); j++) A[0][j] = 0.0;
    for (i = 1; i < A.Row(); i++) A[i][0] = 0.0;

    return true;
  }
  /* w = A' v */
  Matrix A1;
  A.ReferenceMatrix (A1, 1, 0, A.Row() - 1, A.Column());
  Vector v1(A1(0));

  for (j = 1; j < A.Column(); j++) {
    double wj = 0.0;   /* A0j * v0 */
    Vector A1j(A1(j));
    
    wj = (A1j, v1);

    /* A = A - tau v w' */
    A[0][j] = -tau * wj;

    A1j += v1 * (-tau * wj);    
    A1.SetColumn (j, 0, A1j);
  }
  v1 *= -tau;
  A1.SetColumn (0, 0, v1);
  A[0][0] = 1.0 - tau;
  
  return true;
}

/* ********************************************** End of householder.cpp *** */

/* ********************************************************** bidiag.cpp *** *
 *
 *                                    Time-stamp: <07/09/12 09:47:10 sugaya>
 * ************************************************************************* */

/* ************************************************************************* *
 * linalg/bidiag.c
 * 
 * Copyright (C) 2001 Brian Gough
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
 */

/* Factorise a matrix A into
 *
 * A = U B V^T
 *
 * where U and V are orthogonal and B is upper bidiagonal. 
 *
 * On exit, B is stored in the diagonal and first superdiagonal of A.
 *
 * U is stored as a packed set of Householder transformations in the
 * lower triangular part of the input matrix below the diagonal.
 *
 * V is stored as a packed set of Householder transformations in the
 * upper triangular part of the input matrix above the first
 * superdiagonal.
 *
 * The full matrix for U can be obtained as the product
 *
 *       U = U_1 U_2 .. U_N
 *
 * where 
 *
 *       U_i = (I - tau_i * u_i * u_i')
 *
 * and where u_i is a Householder vector
 *
 *       u_i = [0, .. , 0, 1, A(i+1,i), A(i+3,i), .. , A(M,i)]
 *
 * The full matrix for V can be obtained as the product
 *
 *       V = V_1 V_2 .. V_(N-2)
 *
 * where 
 *
 *       V_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [0, .. , 0, 1, A(i,i+2), A(i,i+3), .. , A(i,N)]
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Algorithm 5.4.2 
 *
 * Note: this description uses 1-based indices. The code below uses
 * 0-based indices 
 * ************************************************************************* */
#include <cstdlib>
#include <cstring>
#include "matrix.h"
#include "householder.h"

/* ************************************************************************* */
bool
gsl_linalg_bidiag_decomp (Matrix&	A,
			  Vector&	tau_U,
			  Vector&	tau_V) {
  if (A.Row() < A.Column()) {
    fprintf (stderr, "Bidiagonal decomposition requires M >= N.\n");
  } else if (tau_U.Dimension() != A.Column()) {
    fprintf (stderr, "Size of tau_U must be N.\n");
  } else if (tau_V.Dimension() + 1 != A.Column()) {
    fprintf (stderr, "Size of tau_V must be (N - 1).\n");
  } else {
    const int M = A.Row();
    const int N = A.Column();
    int	 i;

    for (i = 0 ; i < N; i++) {
      /* Apply Householder transformation to current column */
      {
	Vector c(A(i));
	Vector v;
	c.ReferenceVector(v, i, M - i);
	double tau_i = gsl_linalg_householder_transform (v);
	A.SetColumn (i, 0, c);

	/* Apply the transformation to the remaining columns */
            
	if (i + 1 < N) {
	  Matrix m;
	  A.ReferenceMatrix(m, i, i + 1, M - i, N - (i + 1));
	  gsl_linalg_householder_hm (tau_i, v, m);
	  A.SetColumn (i, 0, c);
	}
	tau_U[i] = tau_i;            
      }
      
      /* Apply Householder transformation to current row */
          
      if (i + 1 < N) {
	Vector r(A.Column(), A[i]);
	Vector v;
	r.ReferenceVector (v, i + 1, N - (i + 1));
	double tau_i = gsl_linalg_householder_transform (v);
	A.SetRow (i, 0, r);
	
	/* Apply the transformation to the remaining rows */
	if (i + 1 < M) {
	  Matrix m;
	  A.ReferenceMatrix(m, i + 1, i + 1, M - (i + 1), N - (i + 1));
	  gsl_linalg_householder_mh (tau_i, v, m);
	  A.SetRow (i, 0, r);
	}
	tau_V[i] = tau_i;
      }
    }
  }
  return true;
}

/* ************************************************************************* */
bool
gsl_linalg_bidiag_unpack2 (Matrix&	A, 
                           Vector&	tau_U, 
                           Vector&	tau_V,
                           Matrix&	V) {
  const int M = A.Row();
  const int N = A.Column();
  const int K = (M < N) ? M : N;

  if (M < N) {
    fprintf (stderr, "Matrix A must have M >= N.\n");
  } else if (tau_U.Dimension() != K) {
    fprintf (stderr, "Size of tau must be MIN(M, N).\n");
  } else if (tau_V.Dimension() + 1 != K) {
    fprintf (stderr, "Size of tau must be MIN(M, N) - 1.\n");
  } else if (V.Row() != N || V.Column() != N) {
    fprintf (stderr, "Size of V must be N x N.\n");
  } else {
    int i, j;

    /* Initialize V to the identity */

    V.Unit();

    for (i = N - 1; i > 0 && i--;) {

      /* Householder row transformation to accumulate V */
      Vector r(A.Row(), A[i]);
      Vector h (r, i + 1, N - (i + 1));

      double ti = tau_V[i];

      Matrix m;
      V.ReferenceMatrix (m, i + 1, i + 1, N - (i + 1), N - (i + 1));
      
      gsl_linalg_householder_hm (ti, h, m);
    }

    /* Copy superdiagonal into tau_v */

    for (i = 0; i < N - 1; i++) {
      double Aij = A[i][i + 1];
      tau_V[i] = Aij;
    }

    /* Allow U to be unpacked into the same memory as A, copy
       diagonal into tau_U */

    for (j = N; j > 0 && j--;) {

      /* Householder column transformation to accumulate U */
	double tj = tau_U[j];
	double Ajj = A[j][j];
	Matrix m;
	A.ReferenceMatrix(m, j, j, M - j, N - j);

	tau_U[j] = Ajj;
	gsl_linalg_householder_hm1 (tj, m);
    }
    return true;
  }
  return false;
}

/* *************************************************** End of bidiag.cpp *** */

/* ************************************************************* svd.cpp *** *
 *
 *                                    Time-stamp: <06/10/09 19:55:37 sugaya>
 * ************************************************************************* */

/* ************************************************************************* *
 * linalg/svd.c
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
#include <cstdlib>
#include <cstring>
#include "matrix.h"
#include "bidiag.h"
#include "svdstep.h"
#include "givens.h"

/* Factorise a general M x N matrix A into,
 *
 *   A = U D V^T
 *
 * where U is a column-orthogonal M x N matrix (U^T U = I), 
 * D is a diagonal N x N matrix, 
 * and V is an N x N orthogonal matrix (V^T V = V V^T = I)
 *
 * U is stored in the original matrix A, which has the same size
 *
 * V is stored as a separate matrix (not V^T). You must take the
 * transpose to form the product above.
 *
 * The diagonal matrix D is stored in the vector S,  D_ii = S_i
 */
bool
Matrix::SvDecomp (Matrix&	U,
		  Matrix&	W,
		  Matrix&	V) const {
  U.Initialize (row, col);
  W.Initialize (col, col);
  V.Initialize (col, col);
  U.Clear();
  W.Clear();
  V.Clear();
  
  int		a, b, i, j;
  const int	M = U.Row();
  const int	N = U.Column();
  const int	K = (M < N) ? M : N;
  Vector	S(W.Row());
  Vector	work(N);

  
  if (M < N) {
    fprintf (stderr, "SvDecomp of M x N matrix, M < N, is not implemented.\n");
  } else if (V.Row() != N) {
    fprintf (stderr,
	     "Square matrix V must match second dimension of matrix U.\n");
  } else if (V.Row() != V.Column()) {
    fprintf (stderr, "Matrix V must be square.\n");
  } else if (S.Dimension() != N) {
    fprintf (stderr,
	     "Length of vector S must match second dimension of matrix U.\n");
  }
  /* Handle the case of N = 1 (SVD of a column vector) */

  U = *this;
  if (N == 1) {
    Vector column(U(0));
    double norm = Norm (column);

    S[0] = norm; 
    V[0][0] = 1.0;
      
    if (norm != 0.0) column *= (1.0 / norm);
    U.SetColumn(0, 0, column);

    return true;
  }
  Vector f;
  work.ReferenceVector (f, 0, K - 1);
    
  /* bidiagonalize matrix A, unpack A into U S V */
  gsl_linalg_bidiag_decomp (U, S, f);

  gsl_linalg_bidiag_unpack2 (U, S, f, V);
    
  /* apply reduction steps to B=(S,Sd) */
  chop_small_elements (S, f);
    
  /* Progressively reduce the matrix until it is diagonal */
  b = N - 1;
    
  while (b > 0) {

    double fbm1 = f[b - 1];

    if (fbm1 == 0.0 || isnan (fbm1)) {
      b--;
      continue;
    }
    /* Find the largest unreduced block (a,b) starting from b
       and working backwards */
        
    a = b - 1;
        
    while (a > 0) {
      double fam1 = f[a - 1];

      if (fam1 == 0.0 || isnan (fam1)) break;
      a--;
    }
        
    const int n_block = b - a + 1;
    Vector S_block;
    Vector f_block;
    S.ReferenceVector (S_block, a, n_block);
    f.ReferenceVector (f_block, a, n_block - 1);

    Matrix U_block;
    Matrix V_block;
    U.ReferenceMatrix (U_block, 0, a, U.Row(), n_block);
    V.ReferenceMatrix (V_block, 0, a, V.Row(), n_block);
          
    qrstep (S_block, f_block, U_block, V_block);
#if 0
    {
    U.Print (stderr, "%.16e ");
    fprintf (stderr, "\n");
    S.Print (stderr, "%.16e ");
    fprintf (stderr, "\n");    
    V.Print (stderr, "%.16e ");
    fprintf (stderr, "\n");    
  }
#endif
    /* remove any small off-diagonal elements */
          
    chop_small_elements (S_block, f_block);
  }

  /* Make singular values positive by reflections if necessary */
  
  for (j = 0; j < K; j++) {

    double Sj = S[j];
      
    if (Sj < 0.0) {
      for (i = 0; i < N; i++) {
	double Vij = V[i][j];
	V[i][j] = -Vij;
      }
      S[j] = -Sj;
    }
  }
  
  /* Sort singular values into decreasing order */
  
  for (i = 0; i < K; i++) {
    
    double S_max = S[i];
    int i_max = i;
      
    for (j = i + 1; j < K; j++) {

      double Sj = S[j];
          
      if (Sj > S_max) {
	S_max = Sj;
	i_max = j;
      }
    }
      
    if (i_max != i) {

      /* swap eigenvalues */
      S.Swap (i, i_max);
          
      /* swap eigenvectors */
      U.SwapColumn (i, i_max);
      V.SwapColumn (i, i_max);
    }
  }
  for (i = 0; i < S.Dimension(); i++) W[i][i] = S[i];
  
  return true;
}

/* ****************************************************** End of svd.cpp *** */

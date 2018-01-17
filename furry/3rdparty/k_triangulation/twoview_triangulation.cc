/* ********************************************* twoview_triangulation.c *** *
 * 2画像からの三角測量プログラム
 *
 * Copyright (C) 2007-2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/09 18:03:12 sugaya>
 * ************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include "lib/libmatrix.h"

/* 基礎行列の計算 ********************************************************** */
void
calc_fundamental_matrix (const Matrix& P1,
			 const Matrix& P2,
			 Matrix&       F) {
  Matrix W(4, 4);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[1][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[1][i];
    W[3][i] = P2[2][i];
  }
  F[0][0] = Determinant (W);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[1][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[2][i];
  }
  F[0][1] = -Determinant (W);
  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[1][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[1][i];
  }
  F[0][2] = Determinant (W);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[1][i];
    W[3][i] = P2[2][i];
  }
  F[1][0] = -Determinant (W);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[2][i];
  }
  F[1][1] = Determinant (W);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[2][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[1][i];
  }
  F[1][2] = -Determinant (W);
  
  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[1][i];
    W[2][i] = P2[1][i];
    W[3][i] = P2[2][i];
  }
  F[2][0] = Determinant (W);

  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[1][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[2][i];
  }
  F[2][1] = -Determinant (W);
  
  for (int i = 0; i < 4; i++) {
    W[0][i] = P1[0][i];
    W[1][i] = P1[1][i];
    W[2][i] = P2[0][i];
    W[3][i] = P2[1][i];
  }
  F[2][2] = Determinant (W);

  F = Normalize (F);
}

double
two_view_triangulation (const Vector&	x1,
			const Vector&	x2,
			const Matrix&	F,
			Vector&		_x1,
			Vector&		_x2,
			double		f0,
			int		iter_max,
			double		convergence) {
  Matrix Pk(3, 3), FT(3, 3);  
  Pk[0][0] = Pk[1][1] = 1.0;
  FT = Transpose (F);

  _x1 = x1;
  _x2 = x2;

  Vector dx1(3), dx2(3);  
  double J0 = 1.0e+10, J = 0;
  int    iters = 0;  
  while (1) {

    if (++iters == iter_max) break;
    
    double param1;
    param1 = (_x1, F * _x2) + (F * _x2, dx1) + (FT * _x1, dx2);
    double param2;
    param2 = (F * _x2, Pk * F * _x2) + (FT * _x1, Pk * FT * _x1);

    dx1 = (param1 / param2) * Pk * F * _x2;
    dx2 = (param1 / param2) * Pk * FT * _x1;

    J = (dx1[0] * dx1[0] + dx1[1] * dx1[1] +
	 dx2[0] * dx2[0] + dx2[1] * dx2[1]) * f0 * f0;
    
    if (fabs (J - J0) > convergence) {
      J0 = J;
      _x1 = x1 - dx1;
      _x2 = x2 - dx2;
    } else {
      break;
    }
  }
  return J / 2.0;
}

double
two_view_triangulation_fast(const Vector&	x1,
			    const Vector&	x2,
			    const Matrix&	F,
			    Vector&		_x1,
			    Vector&		_x2,
			    double		f0,
			    int			iter_max,	
			    double		convergence) {
  Vector u(9);
  u[0] = F[0][0];
  u[1] = F[0][1];
  u[2] = F[0][2];  
  u[3] = F[1][0];
  u[4] = F[1][1];
  u[5] = F[1][2];  
  u[6] = F[2][0];
  u[7] = F[2][1];
  u[8] = F[2][2];  

  /* Step 1 */
  double f0f0 = f0 * f0;
  double Q00  = f0f0 * (u[2] * u[2] + u[5] * u[5] + u[6] * u[6] + u[7] * u[7]);
  double p9   = f0f0 * u[8];
  
  /* Step 2 */
  double F1[2][2], F2[2][2];
  
  F1[0][0] = F2[0][0] = u[0];
  F1[0][1] = F2[1][0] = u[1];
  F1[1][0] = F2[0][1] = u[3];
  F1[1][1] = F2[1][1] = u[4];

  double f1[2], f2[2];
  
  f1[0] = u[2] * f0;
  f1[1] = u[5] * f0;
  f2[0] = u[6] * f0;
  f2[1] = u[7] * f0;

  /* Step 3 */
  double Ox1 = x1[0] * f0;
  double Oy1 = x1[1] * f0;  
  double Ox2 = x2[0] * f0;
  double Oy2 = x2[1] * f0;
  double Tx1 = Ox1;
  double Ty1 = Oy1;
  double Tx2 = Ox2;
  double Ty2 = Oy2;
  double Dx1 = 0.0;
  double Dy1 = 0.0;
  double Dx2 = 0.0;
  double Dy2 = 0.0;

  double J0 = 1.0e+10, J = 0;
  int    iters = 0;  

  while (1) {

    if (++iters >= iter_max) break;
    
    /* Step 4*/
    double p1 = (Tx1 * Tx2 + Tx2 * Dx1 + Tx1 * Dx2) * u[0];
    double p2 = (Tx1 * Ty2 + Ty2 * Dx1 + Tx1 * Dy2) * u[1];
    double p3 = f0 * (Tx1 + Dx1) * u[2];
    double p4 = (Ty1 * Tx2 + Tx2 * Dy1 + Ty1 * Dx2) * u[3];
    double p5 = (Ty1 * Ty2 + Ty2 * Dy1 + Ty1 * Dy2) * u[4];
    double p6 = f0 * (Ty1 + Dy1) * u[5];
    double p7 = f0 * (Tx2 + Dx2) * u[6];
    double p8 = f0 * (Ty2 + Dy2) * u[7];

    double xx1 = Tx1 * Tx1;
    double xx2 = Tx2 * Tx2;
    double yy1 = Ty1 * Ty1;	
    double yy2 = Ty2 * Ty2;	
    double xy1 = Tx1 * Ty1;
    double xy2 = Tx2 * Ty2;   
    double f0x1 = f0 * Tx1;
    double f0x2 = f0 * Tx2;    
    double f0y1 = f0 * Ty1;
    double f0y2 = f0 * Ty2;

    double Q11 = (xx1 + xx2) * u[0] * u[0];
    double Q22 = (xx1 + yy2) * u[1] * u[1];
    double Q44 = (yy1 + xx2) * u[3] * u[3];
    double Q55 = (yy1 + yy2) * u[4] * u[4];
    double Q12 = xy2 * u[0] * u[1];
    double Q13 = f0x2 * u[0] * u[2];
    double Q14 = xy1 * u[0] * u[3];
    double Q17 = f0x1 * u[0] * u[6];
    double Q23 = f0y2 * u[1] * u[2];
    double Q25 = xy1 * u[1] * u[4];
    double Q28 = f0x1 * u[1] * u[7];
    double Q45 = xy2 * u[3] * u[4];
    double Q46 = f0x2 * u[3] * u[5];
    double Q47 = f0y1 * u[3] * u[6];
    double Q56 = f0y2 * u[4] * u[5];
    double Q58 = f0y1 * u[4] * u[7];

    double uxi = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;    
    double uxiu = (Q00 + Q11 + Q22 + Q44 + Q55 +
		   2.0 * (Q12 + Q13 + Q14 + Q17 + Q23 + Q25 + Q28 +
			  Q45 + Q46 + Q47 + Q56 + Q58));
    double C = uxi / uxiu;

    /* Step 5 */
    Dx1 = C * (F1[0][0] * Tx2 + F1[0][1] * Ty2 + f1[0]);
    Dy1 = C * (F1[1][0] * Tx2 + F1[1][1] * Ty2 + f1[1]);    
    Dx2 = C * (F2[0][0] * Tx1 + F2[0][1] * Ty1 + f2[0]);
    Dy2 = C * (F2[1][0] * Tx1 + F2[1][1] * Ty1 + f2[1]);    

    /* Step 6 */
    J = Dx1 * Dx1 + Dy1 * Dy1 + Dx2 * Dx2 + Dy2 * Dy2;

    if (fabs (J - J0) > convergence) {
      J0= J;
      Tx1 = Ox1 - Dx1;
      Ty1 = Oy1 - Dy1;
      Tx2 = Ox2 - Dx2;
      Ty2 = Oy2 - Dy2;
    } else {
      break;
    }
  }
  _x1[0] = Tx1 / f0;
  _x1[1] = Ty1 / f0;
  _x2[0] = Tx2 / f0;
  _x2[1] = Ty2 / f0;  
  
  return J / 2.0;
}

/* ************************************** End of twoview_triangulation.c *** */

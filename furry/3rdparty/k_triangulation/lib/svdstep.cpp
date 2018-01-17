/* ********************************************************* svdstep.cpp *** *
 *
 *                                    Time-stamp: <06/02/23 11:18:52 sugaya>
 * ************************************************************************* */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "macros.h"
#include "matrix.h"
#include "givens.h"

/* ************************************************************************* */
void
chop_small_elements (const Vector&	d,
		     Vector&		f) {
  const int	N = d.Dimension();
  double 	d_i = d[0];

  for (int i = 0; i < N - 1; i++) {
    double f_i   = f[i];
    double d_ip1 = d[i + 1];

    if (fabs (f_i) < MATH_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1))) {
      f[i] = 0.0;
    }
    d_i = d_ip1;
  }
}

/* ************************************************************************* */
double
trailing_eigenvalue (const Vector&	d,
		     const Vector&	f) {
  const int	n = d.Dimension();

  double da  = d[n - 2];
  double db  = d[n - 1];
  double fa  = (n > 2) ? f[n - 3] : 0.0;
  double fb  = f[n - 2];
  double ta  = da * da + fa * fa;
  double tb  = db * db + fb * fb;
  double tab = da * fb;
  double dt  = (ta - tb) / 2.0;
  double mu;

  if (dt >= 0) {
    mu = tb - (tab * tab) / (dt + hypot (dt, tab));
  } else {
    mu = tb + (tab * tab) / ((-dt) + hypot (dt, tab));
  }
  return mu;
}

/* ************************************************************************* */
void
create_schur (double	d0,
	      double	f0,
	      double	d1,
	      double	*c,
	      double 	*s) {
  double apq = 2.0 * d0 * f0;
  
  if (apq != 0.0) {
    double t;
    double tau = (f0 * f0 + (d1 + d0) * (d1 - d0)) / apq;
    
    if (tau >= 0.0) {
      t = 1.0 / (tau + hypot (1.0, tau));
    } else {
      t = -1.0 / (-tau + hypot (1.0, tau));
    }
    *c = 1.0 / hypot (1.0, t);
    *s = t * (*c);
  } else {
    *c = 1.0;
    *s = 0.0;
  }
}

/* ************************************************************************* */
void
svd2 (Vector&	d,
      Vector&	f,
      Matrix&	U,
      Matrix&	V) {
  int		i;
  double 	c, s, a11, a12, a21, a22;
  const int	M = U.Row();
  const int	N = V.Row();
  double 	d0 = d[0];
  double 	d1 = d[1];
  double 	f0 = f[0];

  if (d0 == 0.0) {
    /* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */
    create_givens (f0, d1, &c, &s);

    /* compute B <= G^T B X,  where X = [0,1;1,0] */
    d[0] = c * f0 - s * d1;
    f[0] = s * f0 + c * d1;
    d[1] = 0.0;
      
    /* Compute U <= U G */
    for (i = 0; i < M; i++) {
      double Uip = U[i][0];
      double Uiq = U[i][1];
      U[i][0] = c * Uip - s * Uiq;
      U[i][1] = s * Uip + c * Uiq;
    }
    /* Compute V <= V X */
    V.SwapColumn (0, 1);

    return;
  } else if (d1 == 0.0) {
    /* Eliminate off-diagonal element in [d0,f0;0,0] */
    create_givens (d0, f0, &c, &s);

    /* compute B <= B G */
    d[0] = d0 * c - f0 * s;
    f[0] = 0.0;

    /* Compute V <= V G */
    for (i = 0; i < N; i++) {
      double Vip = V[i][0];
      double Viq = V[i][1];
      V[i][0] = c * Vip - s * Viq;
      V[i][1] = s * Vip + c * Viq;
    }
    return;
  } else {
    /* Make columns orthogonal, A = [d0, f0; 0, d1] * G */
    create_schur (d0, f0, d1, &c, &s);
      
    /* compute B <= B G */
    a11 = c * d0 - s * f0;
    a21 = - s * d1;
    a12 = s * d0 + c * f0;
    a22 = c * d1;
      
    /* Compute V <= V G */
    for (i = 0; i < N; i++) {
      double Vip = V[i][0];
      double Viq = V[i][1];
      V[i][0] = c * Vip - s * Viq;
      V[i][1] = s * Vip + c * Viq;
    }
    /* Eliminate off-diagonal elements, bring column with largest
       norm to first column */
    if (hypot(a11, a21) < hypot(a12,a22)) {
      double t1, t2;

      /* B <= B X */
      t1  = a11;
      a11 = a12;
      a12 = t1;
      t2  = a21;
      a21 = a22;
      a22 = t2;

      /* V <= V X */
      V.SwapColumn (0, 1);
    } 
    create_givens (a11, a21, &c, &s);
      
    /* compute B <= G^T B */
    d[0] = c * a11 - s * a21;
    d[1] = s * a12 + c * a22;
    f[0] = c * a12 - s * a22;
      
    /* Compute U <= U G */
    for (i = 0; i < M; i++) {
      double Uip = U[i][0];
      double Uiq = U[i][1];
      U[i][0] = c * Uip - s * Uiq;
      U[i][1] = s * Uip + c * Uiq;
    }
    return;
  }
}

/* ************************************************************************* */
void
chase_out_intermediate_zero (Vector&	d,
			     Vector&	f,
			     Matrix&	U,
			     int	k0) {
  const int	M = U.Row();
  const int	n = d.Dimension();
  double 	c, s;
  double 	x, y;
  int	 	i, k;
  
  x = f[k0];
  y = d[k0 + 1];

  for (k = k0; k < n - 1; k++) {
    create_givens (y, -x, &c, &s);
      
    /* Compute U <= U G */
    for (i = 0; i < M; i++) {
      double Uip = U[i][k0];
      double Uiq = U[i][k + 1];
      U[i][k0]    = c * Uip - s * Uiq;
      U[i][k + 1] = s * Uip + c * Uiq;
    }
    /* compute B <= G^T B */
    d[k + 1] = s * x + c * y;

    if (k == k0) f[k] = c * x - s * y;
    if (k < n - 2) {
      double z = f[k + 1];
      f[k + 1] =  c * z; 
      x = -s * z ;
      y = d[k + 2]; 
    }
  }
}

/* ************************************************************************* */
void
chase_out_trailing_zero (Vector&	d,
			 Vector&	f,
			 Matrix&	V) {
  const int	N = V.Row();
  const int	n = d.Dimension();
  double 	c, s;
  double 	x, y;
  int	 	i, k;

  x = d[n - 2];
  y = f[n - 2];

  for (k = n - 1; k > 0 && k--;) {
    create_givens (x, y, &c, &s);

    /* Compute V <= V G where G = [c, s ; -s, c] */
    for (i = 0; i < N; i++) {
      double Vip = V[i][k];
      double Viq = V[i][n - 1];
      V[i][k]     = c * Vip - s * Viq;
      V[i][n - 1] = s * Vip + c * Viq;
    }
    /* compute B <= B G */
    d[k] = c * x - s * y;
    if (k == n - 2) f[k] = s * x + c * y;
    if (k > 0) {
      double z = f[k - 1];
      f[k - 1] = c * z; 
      x = d[k - 1]; 
      y = s * z;
    }
  }
}

/* ************************************************************************* */
void
qrstep (Vector&	d,
	Vector& f,
	Matrix& U,
	Matrix& V) {
  const int	M = U.Row();
  const int	N = V.Row();
  const int	n = d.Dimension();
  double 	y, z;
  double 	ak, bk, zk, ap, bp, aq, bq;
  int		i, k;

  if (n == 1) return;  /* shouldn't happen */

  /* Compute 2x2 svd directly */
  if (n == 2) {
    svd2 (d, f, U, V);
    return;
  }
  /* Chase out any zeroes on the diagonal */
  for (i = 0; i < n - 1; i++) {
    double d_i = d[i];
    if (d_i == 0.0) {
      chase_out_intermediate_zero (d, f, U, i);
      return;
    }
  }
  /* Chase out any zero at the end of the diagonal */
  {
    double d_nm1 = d[n - 1];
    if (d_nm1 == 0.0) {
      chase_out_trailing_zero (d, f, V);
      return;
    }
  }
  /* Apply QR reduction steps to the diagonal and offdiagonal */
  {
    double d0 = d[0];
    double d1 = d[1];
    double f0 = f[0];
    double f1 = f[1];
    {
      double mu = trailing_eigenvalue (d, f);
      y = d0 * d0 - mu;
      z = d0 * f0;
    }
    /* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    ak = 0;
    bk = 0;
    ap = d0;
    bp = f0;
    aq = d1;
    bq = f1;
  }
  for (k = 0; k < n - 1; k++) {
    double c, s;
    create_givens (y, z, &c, &s);

    /* Compute V <= V G */
    for (i = 0; i < N; i++) {
      double Vip = V[i][k];
      double Viq = V[i][k + 1];
      V[i][k]     = c * Vip - s * Viq;
      V[i][k + 1] = s * Vip + c * Viq;
    }
    /* compute B <= B G */
    {
      double bk1 = c * bk - s * z;
      double ap1 = c * ap - s * bp;
      double bp1 = s * ap + c * bp;
      double zp1 = -s * aq;
      double aq1 = c * aq;

      if (k > 0) f[k - 1] = bk1;

      ak = ap1;
      bk = bp1;
      zk = zp1;
      ap = aq1;

      if (k < n - 2) {
	bp = f[k + 1];
      } else {
	bp = 0.0;
      }
      y = ak;
      z = zk;
    }
    create_givens (y, z, &c, &s);

    /* Compute U <= U G */
    for (i = 0; i < M; i++) {
      double Uip = U[i][k];
      double Uiq = U[i][k + 1];
      U[i][k]     = c * Uip - s * Uiq;
      U[i][k + 1] = s * Uip + c * Uiq;
    }
    /* compute B <= G^T B */
    {
      double ak1 = c * ak - s * zk;
      double bk1 = c * bk - s * ap;
      double zk1 = -s * bp;
      double ap1 = s * bk + c * ap;
      double bp1 = c * bp;

      d[k] = ak1;

      ak = ak1;
      bk = bk1;
      zk = zk1;
      ap = ap1;
      bp = bp1;

      if (k < n - 2) {
	aq = d[k + 2];
      } else {
	aq = 0.0;
      }
      y = bk;
      z = zk;
    }
  }
  f[n - 2] = bk;
  d[n - 1] = ap;
}

/* ************************************************** End of svdstep.cpp *** */

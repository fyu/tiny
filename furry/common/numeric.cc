#include "furry/common/numeric.h"
#include <cmath>

namespace furry
{

/* The code is from Ken Turkowski
 * interface changed to meet the c++ stype
 * The original copyright notice is below
 */
/* Copyright (C) 1997-2001 Ken Turkowski. <turk@computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.
 */

/*******************************************************************************
 * cubic_roots
 *
 *  Solve:
 *      coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *  returns:
 *      3 - 3 real roots
 *      1 - 1 real root (2 complex conjugate)
 *******************************************************************************/

std::vector<double>
cubic_roots(const double coeff[4])
{
  std::vector<double> x;
  double a1 = coeff[2] / coeff[3];
  double a2 = coeff[1] / coeff[3];
  double a3 = coeff[0] / coeff[3];

  double Q = (a1 * a1 - 3 * a2) / 9;
  double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
  double Qcubed = Q * Q * Q;
  double d = Qcubed - R * R;

  /* Three real roots */
  if (d >= 0) {
    x.resize(3);
    double theta = acos(R / sqrt(Qcubed));
    double sqrtQ = sqrt(Q);
    x[0] = -2 * sqrtQ * cos( theta             / 3) - a1 / 3;
    x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
    x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
    return x;
  }

  /* One real root */
  else {
    x.resize(1);
    double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
    if (R > 0)
      e = -e;
    x[0] = (e + Q / e) - a1 / 3.;
    return x;
  }
}

/*  Fast cube root by Ken Turkowski */
float
cube_root(float x)
{
  float fr, r;
  int ex, shx;

  /* Argument reduction */
  fr = frexp(x, &ex);     /* separate into mantissa and exponent */
  shx = ex % 3;
  if (shx > 0)
    shx -= 3;           /* compute shx such that (ex - shx) is divisible by 3 */
  ex = (ex - shx) / 3;    /* exponent of cube root */
  fr = ldexp(fr, shx);                                        /* 0.125 <= fr < 1.0 */

  /* Use quadric rational polynomial with error < 2^(-24) */
  fr = ((((45.2548339756803022511987494 * fr +
           192.2798368355061050458134625) * fr +
          119.1654824285581628956914143) * fr +
         13.43250139086239872172837314) * fr +
        0.1636161226585754240958355063) /
    ((((14.80884093219134573786480845 * fr +
        151.9714051044435648658557668) * fr +
       168.5254414101568283957668343) * fr +
      33.9905941350215598754191872) * fr + 1.0);
  r = ldexp(fr, ex);                                          /* 24 bits of precision */

  return(r);
}

double
cube_root(double x)
{
  double fr, r;
  int ex, shx;

  /* Argument reduction */
  fr = frexp(x, &ex);     /* separate into mantissa and exponent */
  shx = ex % 3;
  if (shx > 0)
    shx -= 3;           /* compute shx such that (ex - shx) is divisible by 3 */
  ex = (ex - shx) / 3;    /* exponent of cube root */
  fr = ldexp(fr, shx);                                        /* 0.125 <= fr < 1.0 */

  /* Compute seed with a quadratic qpproximation */
  fr = (-0.46946116F * fr + 1.072302F) * fr + 0.3812513F;     /* 0.5   <= fr < 1.0 */
  r = ldexp(fr, ex);                                          /* 6 bits of precision */

  /* Newton-Raphson iterations */
  r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);  /* 12 bits of precision */
  r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);  /* 24 bits of precision */
  r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);  /* 48 bits of precision */
  r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);  /* 96 bits of precision */

  return(r);
}

} // furry

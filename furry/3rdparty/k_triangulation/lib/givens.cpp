/* ************************************************************ givens.c ***
 *
 *                                    Time-stamp: <06/02/23 11:16:36 sugaya>
 * ************************************************************************* */

/* ************************************************************************* *
 * linalg/givens.c
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
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 * Generate a Givens rotation (cos,sin) which takes v=(x,y) to (|v|,0) 
 *
 * From Golub and Van Loan, "Matrix Computations", Section 5.1.8
 * ************************************************************************* */
#include <cstdio>
#include <cstdlib>
#include <cmath>

/* ************************************************************************* */
void
create_givens (const double	a,
	       const double	b,
	       double		*c,
	       double 		*s) {
  if (b == 0) {
    *c = 1;
    *s = 0;
  } else if (fabs (b) > fabs (a)) {
    double t = -a / b;
    double s1 = 1.0 / sqrt (1 + t * t);
    *s = s1;
    *c = s1 * t;
  } else {
    double t = -b / a;
    double c1 = 1.0 / sqrt (1 + t * t);
    *c = c1;
    *s = c1 * t;
  }
}

/* *************************************************** End of givens.cpp *** */

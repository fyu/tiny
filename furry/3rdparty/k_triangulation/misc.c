/* ************************************************************** misc.c *** *
 * その他の補助関数
 *
 * Copyright (C) 2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/09 18:06:11 sugaya>
 * ************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include "libmatrix.h"

void
read_projection_matrix (char*	filename,
			double	F0,
			Matrix&	P) {
  Matrix Pk(3, 3);
  Pk[0][0] = Pk[1][1] = 1.0;
  Pk[2][2] = F0;

  FILE* fp = fopen (filename, "r");
  fscanf (fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	  &P[0][0], &P[0][1], &P[0][2], &P[0][3],
	  &P[1][0], &P[1][1], &P[1][2], &P[1][3],
	  &P[2][0], &P[2][1], &P[2][2], &P[2][3]);
  fclose (fp);
  P = Pk * P;
}

void
read_image_data (char		*filename,
		 double		F0,
		 Matrix&	X) {
  

  FILE* fp = fopen (filename, "r");
  int npoints;
  fscanf (fp, "%d", &npoints);

  X.Initialize (3, npoints);
  for (int n = 0; n < npoints; n++) {
    fscanf (fp, "%lf %lf", &X[0][n], &X[1][n]);
    X[0][n] /= F0;
    X[1][n] /= F0;
    X[2][n] = 1.0;
  }
  fclose (fp);
}

void
write_image_data (char		*filename,
		  double	F0,
		  const Matrix&	X) {
  
  FILE* fp = fopen (filename, "w");
  fprintf (fp, "%d\n\n", X.Column());
  for (int n = 0; n < X.Column(); n++) {
    fprintf (fp, "%.16e\t%.16e\n", X[0][n] * F0, X[1][n] * F0);
  }
  fclose (fp);
}
		  
/* ******************************************************* End of misc.c *** */

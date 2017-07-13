/* ************************************************************ matutil.h *** *
 * ヘーダファイル
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <06/12/27 13:51:28 sugaya>
 * ************************************************************************* */
#ifndef	__MATUTIL_H__
#define	__MATUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#ifndef SCALAR
    #define SCALAR double
#endif
typedef SCALAR *vector, **matrix;

void error(char *message);
vector newvec(int n);
matrix newmat(int nrow, int ncol);
vector new_vector(int n);
matrix new_matrix(int nrow, int ncol);
void free_vector(vector v);
void free_matrix(matrix a);
double innerproduct(int n, vector u, vector v);
void vecprint(vector v, int n, int perline, char *format);
void matprint(matrix a, int ncol, int perline, char *format);

#endif	/* __MATUTIL_H__ */

/* **************************************************** End of matutil.h *** */

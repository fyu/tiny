/* ******************************************************** extensions.h *** *
 * お助け関数 ヘーダファイル
 *
 * Copyright (C) 2004-2006 Yasuyuki SUGAYA <sugaya@suri.it.okayama-u.ac.jp>
 *
 *                                    Time-stamp: <06/02/23 11:31:49 sugaya>
 * ************************************************************************* */
#ifndef	__MATRIX_EXTENSIONS_H__
#define	__MATRIX_EXTENSIONS_H__

#define	I_TO_X(i, i0)	((i0) - (i))
#define	J_TO_Y(j, j0)	((j) - (j0))
#define	X_TO_I(x, i0)	((i0) - (x))
#define	Y_TO_J(y, j0)	((j0) + (y))

double**	convert_from_matrix 		(const Matrix&	a);
void		free_array_data 		(double		**a,
						 int		row);
void		convert_coordinate_IJ_to_XY	(Matrix&	X,
						 double		i0 = 240.0,
						 double		j0 = 320.0,
						 double		F0 = 1.0);

void		convert_coordinate_XY_to_IJ 	(Matrix&	X,
						 double		i0 = 240.0,
						 double		j0 = 320.0,
						 double		F0 = 1.0);

#endif	/* __MATRIX_EXTENSIONS_H__ */

/* ************************************************* End of extensions.h *** */

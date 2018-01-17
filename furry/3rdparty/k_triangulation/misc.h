/* ************************************************************** misc.h *** *
 * その他の補助関数 ヘッダファイル
 *
 * Copyright (C) 2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/09 17:39:55 sugaya>
 * ************************************************************************* */
#ifndef	__MISC_H__
#define	__MISC_H__

#include "libmatrix.h"

void    read_projection_matrix	(char	 *filename,
				 double	 F0,
				 Matrix& P);

void  	read_image_data 	(char	 *filename,
				 double	 F0,
				 Matrix& X);

void    write_image_data 	(char	       *filename,
				 double	       F0,
				 const Matrix& X);

#endif	/* __MISC_H__ */

/* ******************************************************* End of misc.h *** */

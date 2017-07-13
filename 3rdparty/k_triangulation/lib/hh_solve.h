/* *********************************************************** hh_solve.h *** *
 * ハウスホルダー変換により線形問題を解く関数 ヘーダファイル
 *
 * Copyright (C) 2006 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <06/11/13 11:39:26 sugaya>
 * ************************************************************************* */
#ifndef	__HH_SOLVE_H__
#define	__HH_SOLVE_H__

int	_gsl_linalg_HH_svx (Matrix		&A,
			    Vector		&x);

#endif	/* __HH_SOLVE_H__ */

/* *************************************************** End of hh_solve.h *** */

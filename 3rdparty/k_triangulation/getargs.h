/* *********************************************************** getargs.h *** *
 * 引数解析関数 ヘッダファイル
 *
 * Copyright (C) 2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/09 17:51:06 sugaya>
 * ************************************************************************* */
#ifndef	__COMMON_FUNC_GETARGS_H__
#define	__COMMON_FUNC_GETARGS_H__

/* ************************************************************************* *
 * 引数用の構造体
 * 	- プログラムに必要な引数に応じてメンバを変更してください -
 * ************************************************************************* */
typedef struct _Argument {
  char		*P_filename[2];
  char		*data_filename[2];
  char		*out_filename[2];
  int		max_iters;
  double	F0;
  double	convergence;
} Argument;

/* ************************************************************************* *
 * 関数
 * ************************************************************************* */

void		usage		(Argument	*arg,
				 char		*cmd);         /* 使い方表示 */
void		version 	(Argument	*arg); /* バージョン情報表示 */
void		argument_free	(Argument	*arg);         /* メモリ解放 */
Argument*	argument_new 	(void);                  /* 引数構造体の生成 */
int		getargs 	(int		argc,        /* 引数解析関数 */
				 char		**argv,
				 Argument	**arg);

#endif	/* __COMMON_FUNC_GETARGS_H__ */

/* **************************************************** End of getargs.h *** */

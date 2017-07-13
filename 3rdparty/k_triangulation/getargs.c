/* *********************************************************** getargs.c *** *
 * 引数解析関数
 *
 * Copyright (C) 2008 Yasuyuki SUGAYA <sugaya@iim.ics.tut.ac.jp>
 *
 *                                    Time-stamp: <08/12/10 12:30:51 sugaya>
 * ************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "getargs.h"

const char*  PROGRAM = "triangulation2";
const char*  VERSION = "1.0.0";
const int    DEFAULT_ITER_MAX = 1000;
const double DEFAULT_F0 = 600.0;
const double DEFAULT_CONVERGENCE = 1.0e-10;

/* 仕様表示 **************************************************************** */
void
usage (Argument	*arg,
       char	*cmd) {
  fprintf (stderr,
	   "\n"
	   "* ***************************************************************"
	   "********** *\n");
  fprintf (stderr,
	   "Usage:%s P1.dat P2.dat data1.dat data2.dat out1.dat out2.dat\n",
	   cmd);
  fprintf (stderr,
	   "\n"
	   "\t --convergence #: Threshold for convergence (Default 1.0e-10).\n"
	   "\t -F0 #\t\t: Default focal length (Default 600.0).\n"
	   "\t -i #\t\t: maximum iterations (Default 1000).\n"
	   "\t -h\t\t: Show Usage.\n"
	   "\t -v\t\t: Show Version.\n"
	   "* ***************************************************************"
	   "********** *\n");
  
  argument_free (arg);
  exit (1);
}

/* バージョン表示 ********************************************************** */
void
version (Argument	*arg) {
  fprintf (stderr, "%s Ver.%s\n", PROGRAM, VERSION);
  argument_free (arg);
  exit (1);
}

/* 引数用構造体の生成 ****************************************************** */
Argument*
argument_new (void) {
  Argument	*arg;

  arg = (Argument *) malloc (sizeof (Argument));

  arg->P_filename[0]    = NULL;
  arg->P_filename[1]    = NULL;  
  arg->data_filename[0] = NULL;
  arg->data_filename[1] = NULL;  
  arg->out_filename[0]  = NULL;
  arg->out_filename[1]  = NULL;    
  arg->max_iters        = DEFAULT_ITER_MAX;
  arg->convergence      = DEFAULT_CONVERGENCE;
  arg->F0	        = DEFAULT_F0;
  
  return arg;
}

/* 引数用構造体のメモリ解放 ************************************************ */
void
argument_free (Argument	*arg) {
  free (arg);
}

/* 引数のチェック ********************************************************** */
int
getargs (int		argc,
 	 char		**argv,
	 Argument	**arg) {
  int 	args;

  *arg = argument_new ();
  
  for (args = 1; args < argc; ++args) {
    if ((*argv[args] == '-') && *(argv[args]+1)) {
      if ((strcmp ("h", ++argv[args])) == 0) {
	usage (*arg, argv[0]);
      } else if ((strcmp ("v", argv[args])) == 0) {
	version (*arg);
      } else if ((strcmp ("-convergence", argv[args])) == 0) {
	if (++args >= argc) return 0;
	(*arg)->convergence = atof (argv[args]);
      } else if ((strcmp ("F0", argv[args])) == 0) {
	if (++args >= argc) return 0;	
	(*arg)->F0 = atof (argv[args]);
      } else if ((strcmp ("i", argv[args])) == 0) {
	if (++args >= argc) return 0;	
	(*arg)->max_iters = atoi (argv[args]);
      }
    } else {
      if (!((*arg)->P_filename[0])) {
	(*arg)->P_filename[0] = argv[args];
      } else if (!((*arg)->P_filename[1])) {
	(*arg)->P_filename[1] = argv[args];
      } else if (!((*arg)->data_filename[0])) {
	(*arg)->data_filename[0] = argv[args];
      } else if (!((*arg)->data_filename[1])) {
	(*arg)->data_filename[1] = argv[args];
      } else if (!((*arg)->out_filename[0])) {
	(*arg)->out_filename[0] = argv[args];
      } else if (!((*arg)->out_filename[1])) {
	(*arg)->out_filename[1] = argv[args];
      } else {
	return 0;
      }
    }
  }
  if (!((*arg)->P_filename[0]) || !((*arg)->P_filename[1]) ||
      !((*arg)->data_filename[0]) || !((*arg)->data_filename[1]) ||
      !((*arg)->out_filename[0]) || !(*arg)->out_filename[1]) {
    return 0;
  }
  return 1;
}

/* **************************************************** End of getargs.c *** */

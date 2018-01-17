/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to UMFPACK sparse linear solver
//////  Copyright (C) 2008  Manolis Lourakis (lourakis at ics.forth.gr)
//////  Institute of Computer Science, Foundation for Research & Technology - Hellas
//////  Heraklion, Crete, Greece.
//////
//////  This program is free software; you can redistribute it and/or modify
//////  it under the terms of the GNU General Public License as published by
//////  the Free Software Foundation; either version 2 of the License, or
//////  (at your option) any later version.
//////
//////  This program is distributed in the hope that it will be useful,
//////  but WITHOUT ANY WARRANTY; without even the implied warranty of
//////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//////  GNU General Public License for more details.
//////
/////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <umfpack.h>


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


struct solvstat{ // solver state
  int n;
  double Control[UMFPACK_CONTROL];
  void *Symbolic;
};

/* solves Ax=B using UMFPACK sparse solver.
 * The function returns 0 in case of error, 1 if successfull
 * Note also that the contents of A and B are not modified on exit.
 *
 * Since this function will be called repeatedly for matrices with
 * identical zero patterns, it will pay off to save some computations
 * by retaining information between calls. Argument "state" is a handle
 * to previously computed such information. Argument "what" specifies
 * the action to be taken by the solver:
 *    0: perform all steps without relying upon any previous computations;
 *       this ignores "state" and is the slowest option overall
 *    1: initialize parameters, allocate memory & perform symbolic factorization
 *       & solve system; store necessary info in "state"
 *    2: numeric factorization and system solution using previous info in "state"
 *    3: free internal state accessed through "state"
 */
int splm_Axb_UMFPACK(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
int i, ret;
struct solvstat *stat;
void *Numeric;
double Info[UMFPACK_INFO];

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_UMFPACK()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;

      /* set default parameters */
      umfpack_di_defaults(stat->Control); // print level is at Control[UMFPACK_PRL]
      //umfpack_di_report_control(Control); // print the control parameters

      /* perform symbolic manipulation of sparsity pattern */
      i=umfpack_di_symbolic(stat->n, stat->n, A->colptr, A->rowidx, A->val, &stat->Symbolic, stat->Control, Info);
      if(i<0){
        umfpack_di_report_info(stat->Control, Info);
        umfpack_di_report_status(stat->Control, i);
        fprintf(stderr, "umfpack_di_symbolic() failed\n");
        exit(1);
      }
      //printf("\nSymbolic factorization: ");
      //umfpack_di_report_symbolic(stat->Symbolic, stat->Control); // print the symbolic factorization

      /* prepare for solving linear system */
      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      /* factorize matrix */
      i=umfpack_di_numeric(A->colptr, A->rowidx, A->val, stat->Symbolic, &Numeric, stat->Control, Info);
      if(i<0){
        umfpack_di_report_info(stat->Control, Info);
        umfpack_di_report_status(stat->Control, i);
        fprintf(stderr, "umfpack_di_numeric() failed\n");

        //exit(1);
        ret=0;
        break;
      }
      //printf("\nNumeric factorization: ");
      //umfpack_di_report_numeric(Numeric, stat->Control); // print the numeric factorization

      /* solve equations */
      i=umfpack_di_solve(UMFPACK_A, A->colptr, A->rowidx, A->val, x, B, Numeric, stat->Control, Info);
      if(i<0){
        umfpack_di_report_info(stat->Control, Info);
        umfpack_di_report_status(stat->Control, i);
        fprintf(stderr, "umfpack_di_solve() failed\n");

        ret=0;
        break;
      }
      umfpack_di_free_numeric(&Numeric);

      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_UMFPACK()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    umfpack_di_free_symbolic(&stat->Symbolic);

    free(stat);
    *state=NULL;
  }

  return ret;
}

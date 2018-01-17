/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to CHOLMOD sparse linear solver
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

#include <cholmod.h>


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


struct solvstat{ // solver state
  int n;
  cholmod_common c;
  cholmod_sparse cmA;
  cholmod_dense cmb;
  cholmod_factor *L;
};

/* solves Ax=B (A symmetric positive definite) using CHOLMOD sparse solver.
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
int splm_Axb_CHOLMOD(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i;
int ret;
cholmod_dense *cmx;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_CHOLMOD()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;

      /* set default parameters */
      cholmod_start(&stat->c);

      stat->cmA.nrow=stat->cmA.ncol=A->nr;
      stat->cmA.nzmax=A->nnz;
      stat->cmA.p=A->colptr;
      stat->cmA.i=A->rowidx;
      stat->cmA.x=A->val;
      stat->cmA.nz=stat->cmA.z=NULL;
      stat->cmA.stype=-1; // square symmetric matrix, use lower triangular part - ignore upper
      stat->cmA.itype=CHOLMOD_INT;
      stat->cmA.xtype=CHOLMOD_REAL;
      stat->cmA.dtype=CHOLMOD_DOUBLE;
      stat->cmA.sorted=stat->cmA.packed=1;

      /* select ordering */
      /* when factorizing many matrices with the same nonzero pattern, the following
       * line instructs cholmod to spend a great deal of time finding a good permutation
       */
      //stat->c.nmethods=CHOLMOD_MAXMETHODS;

#if 0
      /* use AMD or METIS only */
      stat->c.nmethods=1;
      stat->c.method[0].ordering=CHOLMOD_AMD;
      //stat->c.method[0].ordering=CHOLMOD_METIS;
#endif

      /* perform symbolic manipulation of sparsity pattern */
      stat->L=cholmod_analyze(&stat->cmA, &stat->c);

#if 0
      printf("Ordering method %d selected by cholmod_analyze()\n", stat->c.method[stat->c.selected].ordering);
#endif

      if(!stat->L){
        fprintf(stderr, "Unsuccessful termination of cholmod_analyze()\n");
        exit(1);
      }

      /* prepare for solving linear system */
      stat->cmb.nrow=stat->cmb.nzmax=stat->cmb.d=stat->n;
      stat->cmb.ncol=1;
      stat->cmb.z=NULL;
      stat->cmb.xtype=CHOLMOD_REAL;
      stat->cmb.dtype=CHOLMOD_DOUBLE;
      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      stat->cmA.p=A->colptr;
      stat->cmA.i=A->rowidx;
      stat->cmA.x=A->val;
      /* factorize matrix */
      i=cholmod_factorize(&stat->cmA, stat->L, &stat->c);
      if(!i){
        fprintf(stderr, "Unsuccessful termination of cholmod_factorize()\n");

        //exit(1);
        ret=0;
        break;
      }

      /* check for nonzero diagonal */
      if(stat->L->minor<stat->L->n){
        //fprintf(stderr, "Indefinite or singular matrix in splm_Axb_CHOLMOD() [ L(%d, %d)=0 ]\n", stat->L->minor, stat->L->minor);

        //exit(1);
        ret=0;
        break;
      }

      stat->cmb.x=B;
      /* solve equations */
      cmx=cholmod_solve(CHOLMOD_A, stat->L, &stat->cmb, &stat->c);
#if 0
      {
      double *ptr;
      ptr=(double *)cmx->x;
      for(i=stat->n; i-->0; ) //for(i=0; i<stat->n; ++i)
        x[i]=ptr[i];
      }
#endif
      memcpy(x, cmx->x, stat->n*sizeof(double));
      cholmod_free_dense(&cmx, &stat->c);
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_CHOLMOD()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    cholmod_free_factor(&stat->L, &stat->c);
    cholmod_finish(&stat->c);

    free(stat);
    *state=NULL;
  }

  return ret;
}

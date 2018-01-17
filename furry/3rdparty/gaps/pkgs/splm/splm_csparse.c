/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to CSparse sparse linear solver
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

#include "CSparse/cs.h"


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


struct solvstat{ // solver state
  css *S;
  double *x;
  int dim;
};

/* solves Ax=B (A symmetric positive definite) using Davis' sparse Cholesky factorization.
 * The function returns 0 in case of error, 1 if successfull
 * Note also that the contents of A and B are destroyed on exit. // CHECKME
 *
 * Since this function will be called repeatedly for matrices with
 * identical zero patterns, it will pay off to save some computations
 * by retaining information between calls. Argument "state" is a handle
 * to previously computed such information. Argument "what" specifies
 * the action to be taken by the solver:
 *    0: perform all steps without relying upon any previous computations;
 *       this ignores "state" and is the slowest option overall
 *    1: initialize parameters, allocate memory, perform symbolic factorization
 *       & solve system; store necessary info in "state"
 *    2: numeric factorization and system solution using previous info in "state"
 *    3: free internal state accessed through "state"
 */
int splm_Axb_CSPARSE(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
int ret;
struct solvstat *stat;
cs csA;
csn *N;
int order=1; /* 0: natural, 1: amd(A+A') */

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_CSPARSE()!\n");
    exit(1);
  }

  if(what<=2){ /* initialize CSparse matrix */
    csA.nzmax=A->nnz;
    csA.m=csA.n=A->nr;
    csA.p=A->colptr;
    csA.i=A->rowidx;
    csA.x=A->val;
    csA.nz=-1;
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->dim=A->nr;

      stat->S=cs_schol(order, &csA); /* ordering and symbolic analysis */
      if(!stat->S){
        fprintf(stderr, "cs_schol() failed in splm_Axb_CSPARSE()!\n");
        exit(1);
      }

      /* prepare for solving linear system */
      stat->x=(double *)emalloc(stat->dim*sizeof(double));
      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      N=cs_chol(&csA, stat->S); /* numeric Cholesky factorization */
      if(N){
        cs_ipvec(stat->S->pinv, B, stat->x, stat->dim) ; /* x = P*b */
        cs_lsolve(N->L, stat->x) ;   /* x = L\x */
        cs_ltsolve(N->L, stat->x) ;    /* x = L'\x */
        cs_pvec(stat->S->pinv, stat->x, x, stat->dim) ;  /* x = P'*x */

        cs_nfree(N);
        ret=1;
      } else{
        fprintf(stderr, "cs_chol() failed in splm_Axb_CSPARSE()!\n");
        ret=0;
      }
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_CSPARSE()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->x);
    cs_sfree(stat->S);
    free(stat);
    *state=NULL;
  }

  return ret;
}

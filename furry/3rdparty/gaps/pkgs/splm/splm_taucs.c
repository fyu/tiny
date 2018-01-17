/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to TAUCS sparse linear solver
//////  Copyright (C) 2005  Manolis Lourakis (lourakis at ics.forth.gr)
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
#include <math.h>

/* TAUCS prototypes */
#include <taucs.h>


#include "splm.h"
#include "splm_priv.h"

struct solvstat{ // solver state
  taucs_ccs_matrix *A;
  void *F;
};

#undef TAUCS_COPY_MATRIX

/* solves Ax=B (A symmetric) using TAUCS LL^T direct solver.
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
 *    1: initialize parameters, allocate memory, perform symbolic factorization
 *       & solve system; store necessary info in "state"
 *    2: numeric factorization and system solution using previous info in "state"
 *    3: free internal state accessed through "state"
 */
int splm_Axb_TAUCS(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
int ret, tmp;
#ifdef TAUCS_COPY_MATRIX
int trinnz;
register int i, j, k;
#endif

/* Auxiliary variables */
// request a Cholesky multifrontal factorization
/* ordering below can be one of metis, amd,md,mmd, genmmd; metis seems faster */
char *options[6]={"taucs.factor.LLT=true", "taucs.factor.ordering=metis", "taucs.factor.mf=true", NULL, NULL, NULL}; // trailing NULLs used as placeholders here

/* options for out-of-core processing */
//char *ooc[]={"taucs.factor.LLT=true", "taucs.ooc=true", "taucs.ooc.basename=taucs-test", "taucs.ooc.memory=512", NULL};

struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_TAUCS()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)malloc(sizeof(struct solvstat));
      if(!stat){
        fprintf(stderr, "Memory allocation request for \"stat\" failed in splm_Axb_TAUCS()!\n");
        exit(1);
      }

#if 0
      taucs_logfile("./taucs_solver.log"); // logfile name
#endif

      stat->F=NULL;
#ifdef TAUCS_COPY_MATRIX
      /* check if the whole A or its lower triangular part has been supplied */
      if(A->colptr[A->nc]-A->colptr[A->nc-1]==1) // only one element in last column, thus A is triangular
        trinnz=A->nnz;
      else{ // full A
        /* Since A is spd, all its diagonal elements are nonzero.
         * Thus the strictly upper and lower triangles contain nnz-n nonzeros in total
         */
        trinnz=((A->nnz-A->nr)>>1) + A->nr; // (A->nnz-A->nr)/2 + A->nr
      }
      stat->A=taucs_ccs_create(A->nr, A->nr, trinnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);

      /* copy lower part of input matrix */
      for(j=k=0; j<A->nc; ++j){
        stat->A->colptr[j]=k;
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(j<=A->rowidx[i]){
            stat->A->rowind[k]=A->rowidx[i];
            stat->A->values.d[k]=A->val[i];
            ++k;
          }
        }
      }
      stat->A->colptr[A->nc]=k; // k==trinnz here
#else
      /* avoid copying matrix, make TAUCS's matrix reuse A's components */
      stat->A=(taucs_ccs_matrix *)malloc(sizeof(taucs_ccs_matrix));
      if(!stat->A){
        fprintf(stderr, "Memory allocation request for \"stat.A\" failed in splm_Axb_TAUCS()!\n");
        exit(1);
      }
      stat->A->m=stat->A->n=A->nr;
      stat->A->flags=TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
      stat->A->colptr=A->colptr;
      stat->A->rowind=A->rowidx;
      stat->A->values.d=A->val;
#endif /* TAUCS_COPY_MATRIX */

      /* reordering and symbolic factorization. This step also allocates
       * all memory that is necessary for the factorization
       */
      options[3]="taucs.factor.numeric=false"; // no numeric factorization
      options[4]="taucs.factor.symbolic=true"; // perform symbolic factorization
      ret=taucs_linsolve(stat->A, &(stat->F), 0, NULL, NULL, options, NULL);
      if(ret!=TAUCS_SUCCESS){
        fprintf(stderr, "TAUCS error during symbolic factorization: %d\n", ret);
        //exit(1);
        ret=0;
        break;
      }

      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); // retrieve solver state

#ifdef TAUCS_COPY_MATRIX
      /* copy matrix */
      for(j=k=0; j<A->nc; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(j<=A->rowidx[i]){
            stat->A->values.d[k++]=A->val[i];
          }
        }
#else
      stat->A->values.d=A->val;
#endif /* TAUCS_COPY_MATRIX */

      options[3]="taucs.factor.numeric=true"; // perform numeric factorization
      options[4]="taucs.factor.symbolic=false"; // reuse symbolic factorization
      ret=taucs_linsolve(stat->A, &(stat->F), 1, x, B, options, NULL);
      if(ret!=TAUCS_SUCCESS){
        fprintf(stderr, "TAUCS error during numerical factorization: %d\n", ret);
        //exit(2);
        ret=0;
        break;
      }

      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); // retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_TAUCS()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    /* termination and release of memory */
    tmp=taucs_linsolve(NULL, &(stat->F), 0, NULL, NULL, NULL, NULL); // free the factorization
    if(tmp!=TAUCS_SUCCESS)
      fprintf(stderr, "TAUCS error while releasing memory: %d\n", tmp);
#ifdef TAUCS_COPY_MATRIX
    taucs_ccs_free(stat->A);
#else
    free(stat->A);
#endif
    free(stat);
    *state=NULL;
  }

  return ret;
}

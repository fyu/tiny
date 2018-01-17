/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to SPOOLES sparse linear solver
//////  Copyright (C) 2009  Manolis Lourakis (lourakis at ics.forth.gr)
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


/* solve A*X=Y using LU factorization from SPOOLES, code based on LinSol/drivers/testWrapper.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* SPOOLES includes */
#include <LinSol/Bridge.h>


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


struct solvstat{ // solver state
  int n;
  DenseMtx *mtxY, *mtxX;
  Bridge *bridge;
};

/* solves Ax=B (A symmetric) using SPOOLES sparse solver.
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
int splm_Axb_SPOOLES(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j;
int ret;
InpMtx *mtxA=NULL;
FILE   *msgFile;
struct solvstat *stat;
int    error, msglvl, permuteflag, rc;
register double *ptr;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_SPOOLES()!\n");
    exit(1);
  }

  /* set parameters */
  msglvl=0;
  msgFile=stdout; // NULL
  permuteflag=1;

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;

      /* create the InpMtx object */
      mtxA=InpMtx_new();
      InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, A->nnz, A->nr);
      for(j=0; j<A->nc; ++j)
        InpMtx_inputRealColumn(mtxA, j, A->colptr[j+1]-A->colptr[j], A->rowidx+A->colptr[j], A->val+A->colptr[j]);
      InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);

      /* create the rhs object Y */
      stat->mtxY=DenseMtx_new();
      DenseMtx_init(stat->mtxY, SPOOLES_REAL, 0, 0, A->nr, 1, 1, A->nr);

      /* create a bridge object */
      stat->bridge=Bridge_new();
      Bridge_setMatrixParams(stat->bridge, A->nr, SPOOLES_REAL, SPOOLES_SYMMETRIC);
      Bridge_setMessageInfo(stat->bridge, msglvl, msgFile);
      rc=Bridge_setup(stat->bridge, mtxA);
      if(rc!=1){
        fprintf(stderr, "error return %d from Bridge_setup() in splm_Axb_SPOOLES()", rc);
        ret=0;
        break;
      }

      /* prepare for solving linear system */
      stat->mtxX=DenseMtx_new();
      DenseMtx_init(stat->mtxX, SPOOLES_REAL, 0, 0, stat->n, 1, 1, stat->n);
      //DenseMtx_zero(stat->mtxX);

      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      if(!mtxA){ // not allocated in case 1 above
        /* allocate matrix & copy values */
        mtxA=InpMtx_new();
        InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, A->nnz, A->nr);
        for(j=0; j<stat->n; ++j)
          InpMtx_inputRealColumn(mtxA, j, A->colptr[j+1]-A->colptr[j], A->rowidx+A->colptr[j], A->val+A->colptr[j]);
        InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
      }

      /* copy rhs values */
#if 1
      ptr=DenseMtx_entries(stat->mtxY);
      for(i=stat->n; i-->0;  )
        ptr[i]=B[i];
#else
      for(i=0; i<stat->n; ++i)
        DenseMtx_setRealEntry(stat->mtxY, i, 0, B[i]);
#endif /* 1 */

      /* factor the matrix */
      /* Bridge_factor() permutes mtxA into new ordering, thus it is afterwards
       * free'd and reallocated during the next iteration
       */
      rc=Bridge_factor(stat->bridge, mtxA, permuteflag, &error);
      InpMtx_free(mtxA);
      if(rc!=1){
        fprintf(stderr, "error return %d from Bridge_factor() in splm_Axb_SPOOLES()", rc);
        ret=0;
        break;
      }

      /* solve the linear system */
      rc=Bridge_solve(stat->bridge, permuteflag, stat->mtxX, stat->mtxY);
      if(rc!=1){
        fprintf(stderr, "error return %d from Bridge_solve() in splm_Axb_SPOOLES()", rc);
        ret=0;
        break;
      }

#if 1
      ptr=DenseMtx_entries(stat->mtxX);
      /*
      for(i=stat->n; i-->0;  )
        x[i]=ptr[i];
       */
      memcpy(x, ptr, stat->n*sizeof(double));
#else
      for(i=0; i<stat->n; ++i)
        DenseMtx_realEntry(stat->mtxX, i, 0, x+i);
#endif /* 1 */
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_SPOOLES()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    DenseMtx_free(stat->mtxX);
    DenseMtx_free(stat->mtxY);
    Bridge_free(stat->bridge);

    free(stat);
    *state=NULL;
  }

  return ret;
}

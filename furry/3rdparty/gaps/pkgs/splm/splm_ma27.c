/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MA27 sparse linear solver
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
#include <string.h>
#include <math.h>

#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"

extern int ma27id_(int *icntl, double *cntl);

extern int ma27ad_(int *n, int *nz, int *irn, int *icn, int *iw, int *liw, int *ikeep, int *iw1,
                   int *nsteps, int *iflag, int *icntl, double *cntl, int *info, double *ops);

extern int ma27bd_(int *n, int *nz, int *irn, int *icn, double *a, int *la, int *iw, int *liw,
                   int *ikeep, int *nsteps, int *maxfrt, int *iw1, int *icntl, double *cntl, int *info);

extern int ma27cd_(int *n, double *a, int *la, int *iw, int *liw, double *w, int *maxfrt, double *rhs,  
                   int *iw1, int *nsteps, int *icntl, int *info);

#define CNTLSZ  5
#define ICNTLSZ 30
#define INFOSZ  20


struct solvstat{ // solver state
  int n, *irn, *icn, *iw, liw, *ikeep, la, *iw1, nz, icntl[ICNTLSZ], nsteps;
  double *a, *w, cntl[CNTLSZ];
};


/* solves Ax=B (A symmetric) using HSL' MA27 sparse solver.
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
int splm_Axb_MA27(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j, k;
int info[INFOSZ], ret, iflag, maxfrt;
double ops;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_MA27()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;
      
      /* set default parameters */
      ma27id_(stat->icntl, stat->cntl);

#if 0
      /* count nonzero elements in the lower triangle */
      for(j=k=0; j<stat->n; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i)
          if(A->rowidx[i]>=j) ++k;
      stat->nz=k;
#endif

      if(A->colptr[A->nc]-A->colptr[A->nc-1]==1) // A is triangular
        stat->nz=A->nnz;
      else
        stat->nz=((A->nnz-A->nr)>>1) + A->nr; // (A->nnz-A->nr)/2 + A->nr

      /* set up indices */
      stat->irn=(int *)emalloc(stat->nz*sizeof(int));
      stat->icn=(int *)emalloc(stat->nz*sizeof(int));
      for(j=k=0; j<stat->n; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j){ // only copy lower triangle
            stat->irn[k]=A->rowidx[i]+1; // note: fortran indices start from 1!
            stat->icn[k++]=j+1;
          }
        }

      stat->liw=(int)((2*stat->nz+3*stat->n+1)*1.25);
      stat->iw=(int *)emalloc(stat->liw*sizeof(int));
      stat->ikeep=(int *)emalloc((3*stat->n)*sizeof(int));
      //memset(stat->ikeep, 0, (3*stat->n)*sizeof(int)); // set to zero
      stat->iw1=(int *)emalloc((2*stat->n+2)*sizeof(int)); // max(2*n, n, 2*n+2)
      iflag=0; // automatic pivot order
      /* perform symbolic manipulation of sparsity pattern */
      ma27ad_(&(stat->n), &(stat->nz), stat->irn, stat->icn, stat->iw, &(stat->liw), stat->ikeep, stat->iw1, &(stat->nsteps), &iflag,
                  stat->icntl, stat->cntl, info, &ops);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA27AD(), return code %d\n", info[0]);

        if(info[0]==-1) fprintf(stderr, "Value of N out of range (N=%d)\n", stat->n);
        else if(info[0]==-2) fprintf(stderr, "Value NZ out of range (NZ=%d)\n", stat->nz);
        else if(info[0]==-3) fprintf(stderr, "Insufficient space allocated for array IW (need at least %d)\n", info[1]);
        //else if(info[0]==-5) fprintf(stderr, "Indices in KEEP do not constitute a permutation\n");
        //else if(info[0]==-6) fprintf(stderr, "A negative index in KEEP is not followed by another negative one\n");
        else if(info[0]==1) fprintf(stderr, "An index in IRN or ICN is out of range\n");
        //else if(info[0]==2) fprintf(stderr, "Duplicate entries found\n");
        //else if(info[0]>=4) fprintf(stderr, "Matrix is rank deficient\n");
        exit(1);
      }

      if(info[11]>10) fprintf(stderr, "%d compresses performed by MA27AD(), performance can be improved by increasing "
                            "the length LIW of IW; LIW should be at least %d\n", info[11], info[6]);
      if(info[5]>stat->liw){
        fprintf(stderr, "Size LIW of input matrix IW in MA27BD() should be at least %d (now %d)\n", info[5], stat->liw);
        exit(1);
      }

      /* prepare for solving linear system */
      stat->la=(int)(info[4]*1.5); // 50% more than suggested by MA27AD
      stat->a=(double *)emalloc(stat->la*sizeof(double));
      stat->w=(double *)emalloc(stat->n*sizeof(double)); // required size is maxfrt<=n
      /* fall through */

    case 2: /* numerical factorization & system solution */  
      stat=(struct solvstat *)(*state); // retrieve solver state
      /* copy matrix and rhs vector */
      for(j=k=0; j<stat->n; ++j){
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j) // only copy lower triangle
            stat->a[k++]=A->val[i];
        }

        x[j]=B[j];
      }

      /* factorize matrix */
      ma27bd_(&stat->n, &stat->nz, stat->irn, stat->icn, stat->a, &stat->la, stat->iw, &stat->liw, stat->ikeep,
                  &(stat->nsteps), &maxfrt, stat->iw1, stat->icntl, stat->cntl, info);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA27BD(), return code %d\n", info[0]);

        if(info[0]==-1) fprintf(stderr, "Value of N out of range (N=%d)\n", stat->n);
        else if(info[0]==-2) fprintf(stderr, "Too few non zero elements (NZ=%d)\n", stat->nz);
        else if(info[0]==-3) fprintf(stderr, "Insufficient space allocated for array IW (need at least %d)\n", info[1]);
        else if(info[0]==-4) fprintf(stderr, "Insufficient space allocated for array A (need at least %d)\n", info[1]);
        else if(info[0]==-5) fprintf(stderr, "Matrix is singular\n");
        else if(info[0]==-6) fprintf(stderr, "Change of pivots sign detected when U was negative\n");
        else if(info[0]==-7) fprintf(stderr, "Value of NSTEPS out of range (NSTEPS=%d)\n", stat->nsteps);
        else if(info[0]==2) fprintf(stderr, "Pivots of different signs whan factorizing supposedly definite matrix\n");
        else if(info[0]>=3) fprintf(stderr, "Matrix is rank deficient\n");
        //exit(1);
        ret=0;
        break;
      }

      /* solve equations */
      ma27cd_(&stat->n, stat->a, &stat->la, stat->iw, &stat->liw, stat->w, &maxfrt, x, stat->iw1, &(stat->nsteps), stat->icntl, info);
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); // retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_MA27()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->w);
    free(stat->iw1);
    free(stat->ikeep);
    free(stat->iw);
    free(stat->a);
    free(stat->icn);
    free(stat->irn);
    free(stat);
    *state=NULL;
  }

  return ret;
}

/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MA47 sparse linear solver
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

extern int ma47id_(double *cntl, int *icntl);

extern int ma47ad_(int *n, int *ne, int *irn, int *jcn, int *iw, int *liw, int *keep,
                   int *icntl, double *rinfo, int *info);

extern int ma47bd_(int *n, int *ne, int *jcn, double *a, int *la, int *iw, int *liw, int *keep,
                   double *cntl, int *icntl, int *iw1, double *rinfo, int *info);

extern int ma47cd_(int *n, double *a, int *la, int *iw, int *liw, double *w,
                   double *rhs, int *iw1, int *icntl);

#define CNTLSZ  2
#define RINFOSZ 4
#define ICNTLSZ 7
#define INFOSZ  24


struct solvstat{ // solver state
  int n, *irn, *jcn, *iw, liw, *keep, la, *iw1, ne, icntl[ICNTLSZ];
  double *a, *w, cntl[CNTLSZ];
};

/* solves Ax=B (A symmetric) using HSL' MA47 sparse indefinite solver.
 * The function returns 0 in case of error, 1 if successfull
 *  Note also that the contents of A and B are not modified on exit.
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
int splm_Axb_MA47(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j, k;
double rinfo[RINFOSZ];
int info[INFOSZ], ret;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_MA47()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;

      /* set default parameters */
      ma47id_(stat->cntl, stat->icntl);

#if 0
      /* count nonzero elements in the lower triangle */
      for(j=k=0; j<stat->n; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i)
          if(A->rowidx[i]>=j) ++k;
      stat->ne=k;
#endif

      if(A->colptr[A->nc]-A->colptr[A->nc-1]==1) // A is triangular
        stat->ne=A->nnz;
      else
        stat->ne=((A->nnz-A->nr)>>1) + A->nr; // (A->nnz-A->nr)/2 + A->nr

      /* set up indices */
      stat->irn=(int *)emalloc(stat->ne*sizeof(int));
      stat->jcn=(int *)emalloc(stat->ne*sizeof(int));
      for(j=k=0; j<stat->n; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j){ // only copy lower triangle
            stat->irn[k]=A->rowidx[i]+1; // note: fortran indices start from 1!
            stat->jcn[k++]=j+1;
          }
        }

      stat->liw=(int)((stat->ne*2+5*stat->n+4)*1.25);
      stat->iw=(int *)emalloc(stat->liw*sizeof(int));
      stat->keep=(int *)emalloc((stat->ne+5*stat->n+2)*sizeof(int));
      memset(stat->keep, 0, (stat->ne+5*stat->n+2)*sizeof(int)); // set to zero

      /* perform symbolic manipulation of sparsity pattern */
      ma47ad_(&(stat->n), &(stat->ne), stat->irn, stat->jcn, stat->iw, &(stat->liw), stat->keep, stat->icntl, rinfo, info);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA47AD(), return code %d\n", info[0]);

        if(info[0]==-1) fprintf(stderr, "Value of N out of range (N=%d)\n", stat->n);
        else if(info[0]==-2) fprintf(stderr, "Too few non zero elements (NE=%d)\n", stat->ne);
        else if(info[0]==-3) fprintf(stderr, "Insufficient space allocated for array IW (need at least %d)\n", info[1]);
        else if(info[0]==-5) fprintf(stderr, "Indices in KEEP do not constitute a permutation\n");
        else if(info[0]==-6) fprintf(stderr, "A negative index in KEEP is not followed by another negative one\n");
        else if(info[0]==1) fprintf(stderr, "An index in IRN or JCN is out of range\n");
        else if(info[0]==2) fprintf(stderr, "Duplicate entries found\n");
        else if(info[0]>=4) fprintf(stderr, "Matrix is rank deficient\n");
        exit(1);
      }

      if(info[11]>10) fprintf(stderr, "%d compresses performed by MA47AD(), performance can be improved by increasing "
                            "the length LIW of IW; LIW should be at least %d\n", info[11], info[6]);
      if(info[6]>stat->liw){
        fprintf(stderr, "Size LIW of input matrix IW in MA47BD() should be at least %d (now %d)\n", info[6], stat->liw);
        exit(1);
      }

      /* prepare for solving linear system */
      stat->la=(int)(info[5]*1.5); // 50% more than suggested by MA47AD
      stat->a=(double *)emalloc(stat->la*sizeof(double));
      stat->iw1=(int *)emalloc((2*stat->n+2)*sizeof(int));
      stat->w=(double *)emalloc(stat->n*sizeof(double));
      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state
      /* copy matrix and rhs vector */
      for(j=k=0; j<stat->n; ++j){
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j) // only copy lower triangle
            stat->a[k++]=A->val[i];
        }

        x[j]=B[j];
      }

      /* factorize matrix */
      ma47bd_(&(stat->n), &(stat->ne), stat->jcn, stat->a, &(stat->la), stat->iw, &(stat->liw), stat->keep,
                  stat->cntl, stat->icntl, stat->iw1, rinfo, info);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA47BD(), return code %d\n", info[0]);

        if(info[0]==-2) fprintf(stderr, "Too few non zero elements (NE=%d)\n", stat->ne);
        else if(info[0]==-3) fprintf(stderr, "Insufficient space allocated for array IW (need at least %d)\n", info[1]);
        else if(info[0]==-4) fprintf(stderr, "Insufficient space allocated for array A (need at least %d)\n", info[1]);
        else if(info[0]>=4) fprintf(stderr, "Matrix is rank deficient\n");
        //exit(1);
        ret=0;
        break;
      }

      if(info[17]>10) fprintf(stderr, "%d compresses on real data performed by MA47BD(), factorization performance "
                        "can be improved by increasing the length LA of A\n", info[17]);
      if(info[18]>10) fprintf(stderr, "%d compresses on integer data performed by MA47BD(), factorization performance "
                        "can be improved by increasing the length LIW of IW\n", info[18]);

      /* solve equations */
      ma47cd_(&(stat->n), stat->a, &(stat->la), stat->iw, &(stat->liw), stat->w, x, stat->iw1, stat->icntl);
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_MA47()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->w);
    free(stat->iw1);
    free(stat->keep);
    free(stat->iw);
    free(stat->a);
    free(stat->jcn);
    free(stat->irn);
    free(stat);
    *state=NULL;
  }

  return ret;
}

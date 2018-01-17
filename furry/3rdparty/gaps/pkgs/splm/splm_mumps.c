/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MUMPS sparse linear solver
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


/* solve A*X=Y using LU factorization from MUMPS, code based on LinSol/drivers/testWrapper.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* MUMPS includes */
#include <mpi.h>
#include <dmumps_c.h>

#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"

#define USE_COMM_WORLD -987654


struct solvstat{ // solver state
  DMUMPS_STRUC_C id;
};

/* solves Ax=B (A symmetric) using MUMPS sparse solver.
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
int splm_Axb_MUMPS(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j, k;
int ret;
int myid, ierr;
struct solvstat *stat;
register double *ptr;

int argc=1;
char **argv, *name="splm_Axb_MUMPS";

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_MUMPS()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));

      argv=&name;
      ierr=MPI_Init(&argc, &argv);
      ierr=MPI_Comm_rank(MPI_COMM_WORLD, &myid);

      /* init a MUMPS instance */
      stat->id.par=1; // host processor is working
      stat->id.sym=1; // A is SPD
      stat->id.comm_fortran=USE_COMM_WORLD;
      stat->id.job=-1;
      dmumps_c(&stat->id);
      if(stat->id.info[0]){
        if(stat->id.info[0]<0){
          fprintf(stderr, "error (code %d) during initialization in dmumps_c()!\n", stat->id.info[1]);
          ret=0;
          break;
        }
        fprintf(stderr, "warning (code %d) during initialization in dmumps_c()!\n", stat->id.info[1]);
      }

      /* specify parameters */
      stat->id.icntl[0]=-1; stat->id.icntl[1]=-1; stat->id.icntl[2]=-1; stat->id.icntl[3]=0; // no outputs
      stat->id.icntl[4]=0; stat->id.icntl[17]=0; // centralized assembled matrix input
      stat->id.icntl[5]=0; // don't compute column permutation
      //stat->id.icntl[6]=0; // AMD ordering
      //stat->id.icntl[6]=2; // AMF ordering
      //stat->id.icntl[6]=4; // PORD ordering
      //stat->id.icntl[6]=5; // METIS ordering
      stat->id.icntl[6]=7; // automatic choice for ordering

      stat->id.n=A->nr;
      if(A->colptr[A->nc]-A->colptr[A->nc-1]==1) // A is triangular
        stat->id.nz=A->nnz;
      else
        stat->id.nz=((A->nnz-A->nr)>>1) + A->nr; // (A->nnz-A->nr)/2 + A->nr

      /* set up indices */
      stat->id.irn=(int *)emalloc(stat->id.nz*sizeof(int));
      stat->id.jcn=(int *)emalloc(stat->id.nz*sizeof(int));
      stat->id.a=NULL;
      for(j=k=0; j<A->nr; ++j)
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j){ // only copy lower triangle
            stat->id.irn[k]=A->rowidx[i]+1; // note: fortran indices start from 1!
            stat->id.jcn[k++]=j+1;
          }
      }

      /* symbolic analysis */
      stat->id.job=1;
      dmumps_c(&stat->id);
      if(stat->id.info[0]){
        if(stat->id.info[0]<0){
          fprintf(stderr, "error (code %d) during analysis in dmumps_c()!\n", stat->id.info[1]);
          ret=0;
          break;
        }
        fprintf(stderr, "warning (code %d) during analysis in dmumps_c()!\n", stat->id.info[1]);
      }

      /* prepare for solving linear system */
      stat->id.a=(double *)emalloc(stat->id.nz*sizeof(double));
      stat->id.rhs=(double *)emalloc(A->nr*sizeof(double));

      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      /* copy matrix and rhs vector */
      ptr=stat->id.rhs;
      for(j=k=0; j<stat->id.n; ++j){
        for(i=A->colptr[j]; i<A->colptr[j+1]; ++i){
          if(A->rowidx[i]>=j) // only copy lower triangle
            stat->id.a[k++]=A->val[i];
        }

        ptr[j]=B[j];
      }

      /* factor the matrix & solve system */
      stat->id.job=5; // combines actions of job=2 and job=3
      dmumps_c(&stat->id);
      if(stat->id.info[0]){
        if(stat->id.info[0]<0){
          fprintf(stderr, "error (code %d) during factor & solve in dmumps_c()!\n", stat->id.info[1]);
          ret=0;
          break;
        }
        fprintf(stderr, "warning (code %d) during factor & solve in dmumps_c()!\n", stat->id.info[1]);
      }

#if 0
      ptr=stat->id.rhs;
      for(i=stat->id.n; i-->0;  )
        x[i]=ptr[i];
#endif
      memcpy(x, stat->id.rhs, stat->id.n*sizeof(double));
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_MUMPS()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->id.irn);
    free(stat->id.jcn);
    free(stat->id.a);
    free(stat->id.rhs);

    stat->id.job=-2; // kill instance
    dmumps_c(&stat->id);
    if(stat->id.info[0]){
      if(stat->id.info[0]<0){
        fprintf(stderr, "error (code %d) during finalize in dmumps_c()!\n", stat->id.info[1]);
        ret=0;
      }
      fprintf(stderr, "warning (code %d) during finalize in dmumps_c()!\n", stat->id.info[1]);
    }
    ierr=MPI_Finalize();

    free(stat);
    *state=NULL;
  }

  return ret;
}

/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MKL's DSS sparse linear solver simplified PARDISO interface.
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

#include "splm.h"
#include "splm_priv.h"

#include <mkl_dss.h>

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif


struct solvstat{ // solver state
  _MKL_DSS_HANDLE_t handle;
  MKL_INT opt, sym, typ;
  int n, nrhs, *ia, *ja;
  double *trivals;
};

/* solves Ax=B (A symmetric) using MKL's DSS solver PARDISO interface.
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
int splm_Axb_DSS(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
int trinnz, ret;

/* Auxiliary variables */
register int i, j, k;
int *rowptr, *colidx;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_DSS()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      if(!stat){
        fprintf(stderr, "Memory allocation request for \"stat\" failed in splm_Axb_DSS()!\n");
        exit(1);
      }

      /* initialize solver */
      stat->opt=MKL_DSS_DEFAULTS;
      /* symmetric spd matrix */
      stat->sym=MKL_DSS_SYMMETRIC;
      stat->typ=MKL_DSS_POSITIVE_DEFINITE;
      stat->nrhs=1; // single rhs

      if(dss_create(stat->handle, stat->opt)!=MKL_DSS_SUCCESS){
        fprintf(stderr, "DSS initialization failed in splm_Axb_DSS()!\n");
        exit(1);
      }
      stat->n=A->nr;

      /* DSS PARDISO interface accepts CRS matrices. Furthermore, for symmetric matrices,
       * only the upper triangular part + diagonal should be provided. Input matrix is in CCS,
       * thus it can be regarded as A^T in CRS, which is identical to A due to symmetry.
       * Note that 1 has to be added to the elements of colptr & rowidx since DSS
       * uses Fortran's index convention
       */

      /* Build the upper triangle in CRS */
      /* check if the whole A or its lower triangular part has been supplied */
      if(A->colptr[A->nc]-A->colptr[A->nc-1]==1) // only one element in last column, thus A is triangular
        trinnz=A->nnz;
      else{ // full A

        /* Since A is spd, all its diagonal elements are nonzero.
         * Thus the strictly upper and lower triangles contain nnz-n nonzeros in total
         */
        trinnz=((A->nnz-stat->n)>>1) + stat->n; // (A->nnz-stat->n)/2 + stat->n
      }

      /* set up indices */
      stat->ia=(int *)emalloc((stat->n+1)*sizeof(int));
      stat->ja=(int *)emalloc(trinnz*sizeof(int));
      stat->trivals=(double *)emalloc(trinnz*sizeof(double));
      if(!stat->ia || !stat->ja || !stat->trivals){
        fprintf(stderr, "Memory allocation request failed in splm_Axb_DSS()!\n");
        exit(1);
      }

      rowptr=A->colptr; colidx=A->rowidx;
      for(i=k=0; i<stat->n; ++i){
        stat->ia[i]=k+1; // fortran!
        for(j=rowptr[i]; j<rowptr[i+1]; ++j){
          if(i<=colidx[j]){
            stat->ja[k]=colidx[j]+1; // fortran!
            //stat->trivals[k]=A->val[j]; // NOTE: matrix values are copied below!
            ++k;
          }
        }
      }
      stat->ia[stat->n]=k+1; // fortran! // k==trinnz here

      /* define the non-zero structure of the matrix */
      if(dss_define_structure(stat->handle, stat->sym, stat->ia, stat->n, stat->n, stat->ja, trinnz)!=MKL_DSS_SUCCESS){
        fprintf(stderr, "DSS matrix structure definition failed in splm_Axb_DSS()!\n");
        exit(1);
      }

      /* reorder matrix */
      if(dss_reorder(stat->handle, stat->opt, 0)!=MKL_DSS_SUCCESS){
        fprintf(stderr, "DSS matrix reordering failed in splm_Axb_DSS()!\n");
        exit(1);
      }

      //printf("\nReordering completed ... ");

      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); // retrieve solver state

      /* copy matrix */
      rowptr=A->colptr; colidx=A->rowidx;
      for(i=k=0; i<stat->n; ++i)
        for(j=rowptr[i]; j<rowptr[i+1]; ++j){
          if(i<=colidx[j]){
            stat->trivals[k++]=A->val[j];
          }
        }

      /* factor matrix */
      if(dss_factor_real(stat->handle, stat->typ, stat->trivals)!=MKL_DSS_SUCCESS){
        fprintf(stderr, "DSS matrix factoring failed in splm_Axb_DSS()!\n");
        //exit(2);
        ret=0;
        break;
      }
      //printf("\nFactorization completed ...\n ");

      /* retrieve solution vector */
      if(dss_solve_real(stat->handle, stat->opt, B, stat->nrhs, x)!=MKL_DSS_SUCCESS){
        fprintf(stderr, "DSS error during solution in splm_Axb_DSS()!\n");
        //exit(3);
        ret=0;
      }
      else ret=1;
      //printf("\nSolve completed ... ");
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); // retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_DSS()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    /* termination and release of memory */
    dss_delete(stat->handle, stat->opt);
    free(stat->ia);
    free(stat->ja);
    free(stat->trivals);
    free(stat);
    *state=NULL;
  }

  return ret;
}

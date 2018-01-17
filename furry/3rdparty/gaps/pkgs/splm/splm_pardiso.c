/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to PARDISO sparse linear solver
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

/*#include <mkl_pardiso.h>*/

/* PARDISO prototypes */
extern int pardisoinit_(void *pt[], int *mtype, int *iparm);

extern int pardiso_(void *pt[], int *maxfct, int *mnum, int *mtype, int *phase, int *n, double *a, int *ia, int *ja,
                   int *perm, int *nrhs, int *iparm, int *msglvl, double *b, double *x, int *error);

#define PSIZE  64


struct solvstat{ // solver state
  /* internal solver memory pointer pt,                     */
  /* 32-bit: int pt[PSIZE]; 64-bit: long int pt[PSIZE]      */
  /* or void *pt[PSIZE] should be OK on both architectures  */ 
  void *pt[PSIZE]; 
  /* Pardiso control parameters */
  int iparm[PSIZE];
  int n, mtype, nrhs, *ia, *ja, maxfct, mnum, msglvl;
  double *trivals;
};

/* solves Ax=B (A symmetric) using PARDISO solvers. Adapted from pardiso_sym.c,
 * see http://www.computational.unibas.ch/cs/scicomp
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
int splm_Axb_PARDISO(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
int phase, trinnz, ret;
static int num_procs=-1; // number of processors, determined on 1st call only

/* Auxiliary variables */
register int i, j, k;
int *rowptr, *colidx, error=0;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_PARDISO()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      if(!stat){
        fprintf(stderr, "Memory allocation request for \"stat\" failed in splm_Axb_PARDISO()!\n");
        exit(1);
      }

      stat->n=A->nr;
      stat->mtype=2; // real symmetric spd matrix; use 11 for nonsymmetric matrices
      stat->nrhs=1; // single rhs
      stat->maxfct=1;	// max. numerical factorizations
      stat->mnum=1; // factorization to use
    
      stat->msglvl=0; // do not print statistical info


      /* Numbers of processors, value of OMP_NUM_THREADS */
      if(num_procs<0){
        char *var=getenv("OMP_NUM_THREADS");

        if(var!=NULL)
          sscanf(var, "%d", &num_procs);
        else{
          printf("Environment variable OMP_NUM_THREADS not set, assuming 1\n");
          num_procs=1;
        }
      }

      /* setup Pardiso control parameters */
      pardisoinit_(stat->pt, &(stat->mtype), stat->iparm); 
      stat->iparm[2]=num_procs; //omp_get_max_threads(); /* Number of processors */
      stat->iparm[7]=1; // max. iterative refinement steps

#if 0
      for(i=0; i<PSIZE; ++i){
        stat->iparm[i]=0;
        stat->pt[i]=0;
      }
      stat->iparm[0]=1;     /* No solver default */
      stat->iparm[1]=2;     /* Fill-in reordering from METIS */
      stat->iparm[2]=num_procs; //omp_get_max_threads(); /* Number of processors */
      stat->iparm[3]=0;     /* No iterative-direct algorithm */
      stat->iparm[4]=0;     /* No user fill-in reducing permutation */
      stat->iparm[5]=0;     /* Write solution into x */
      stat->iparm[6]=16;    /* Default logical fortran unit number for output */
      stat->iparm[7]=2;     /* Max numbers of iterative refinement steps */
      stat->iparm[8]=0;     /* Not in use */
      stat->iparm[9]=13;    /* Perturb the pivot elements with 1E-13 */
      stat->iparm[10]=1;    /* Use nonsymmetric permutation and scaling MPS */
      stat->iparm[11]=0;    /* Not in use */
      stat->iparm[12]=0;    /* Not in use */
      stat->iparm[13]=0;    /* Output: Number of perturbed pivots */
      stat->iparm[14]=0;    /* Not in use */
      stat->iparm[15]=0;    /* Not in use */
      stat->iparm[16]=0;    /* Not in use */
      stat->iparm[17]=-1;   /* Output: Number of nonzeros in the factor LU */
      stat->iparm[18]=-1;   /* Output: Mflops for LU factorization */
      stat->iparm[19]=0;    /* Output: Numbers of CG Iterations */
#endif /* 0 */

      /* PARDISO accepts CRS matrices. Furthermore, for symmetric matrices, only the
       * upper triangular part + diagonal should be provided. Input matrix is in CCS,
       * thus it can be regarded as A^T in CRS, which is identical to A due to symmetry.
       * Note that 1 has to be added to the elements of colptr & rowidx since PARDISO
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
        fprintf(stderr, "Memory allocation request failed in splm_Axb_PARDISO()!\n");
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

      /* reordering and symbolic factorization. This step also allocates
       * all memory that is necessary for the factorization
       */
      phase=11; 
      pardiso_(stat->pt, &(stat->maxfct), &(stat->mnum), &(stat->mtype), &phase, &(stat->n), NULL /*trivals*/, stat->ia, stat->ja, NULL,
                    &(stat->nrhs), stat->iparm, &(stat->msglvl), NULL, NULL, &error);
      if(error!=0) {
        fprintf(stderr, "Pardiso error during symbolic factorization: %d\n", error);
        //exit(1);
        ret=0;
        break;
      }

      //printf("\nReordering completed ... ");
      //printf("\nNumber of nonzeros in factors  = %d", stat->iparm[17]);
      //printf("\nNumber of factorization MFLOPS = %d\n", stat->iparm[18]);

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

      phase=22;
      pardiso_(stat->pt, &(stat->maxfct), &(stat->mnum), &(stat->mtype), &phase, &(stat->n), stat->trivals, stat->ia, stat->ja, NULL,
                      &(stat->nrhs), stat->iparm, &(stat->msglvl), NULL, NULL, &error);
      if(error!=0) {
        fprintf(stderr, "Pardiso error during numerical factorization: %d\n", error);
        //exit(2);
        ret=0;
        break;
      }
      //printf("\nFactorization completed ...\n ");

      /* back substitution and iterative refinement */
      phase=33;
      pardiso_(stat->pt, &(stat->maxfct), &(stat->mnum), &(stat->mtype), &phase, &(stat->n), stat->trivals, stat->ia, stat->ja, NULL,
                      &(stat->nrhs), stat->iparm, &(stat->msglvl), B, x, &error);
      if(error!=0) {
        fprintf(stderr, "Pardiso error during solution: %d\n", error);
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
      fprintf(stderr, "Unknown mode %d in splm_Axb_PARDISO()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    /* termination and release of memory */
    phase=-1; // release internal memory
    pardiso_(stat->pt, &(stat->maxfct), &(stat->mnum), &(stat->mtype), &phase, &(stat->n), NULL, stat->ia, stat->ja, NULL,
                    &(stat->nrhs), stat->iparm, &(stat->msglvl), NULL, NULL, &error);
    free(stat->ia);
    free(stat->ja);
    free(stat->trivals);
    free(stat);
    *state=NULL;
  }

  return ret;
}

/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to SuperLU sparse linear solver
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

/* SuperLU prototypes */
#include <slu_ddefs.h>


#include "splm.h"
#include "splm_priv.h"

/* The following functions are based on dCreate_CompCol_Matrix() & dCreate_Dense_Matrix() in SuperLU's dutil.c
 * They copy the supplied values to their SuperMatrix arguments, without allocating any memory.
 */
static void
_dSet_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, 
		       double *nzval, int *rowind, int *colptr,
		       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    Astore = A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

static void
_dSet_Dense_Matrix(SuperMatrix *X, int m, int n, double *x, int ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (double *) x;
}


struct solvstat{ // solver state
  superlu_options_t options;
  SuperMatrix A, B, X;
  int *perm_c; /* column permutation vector */
  int *perm_r; /* row permutations from partial pivoting */
  int *etree, lwork;
  double *R, *C, *work;
  /* double *rhs; */
};

/* solves Ax=B (A symmetric) using SuperLU solvers. Adapted from dlinsol1.c/dlinsolx2.c
 * in SuperLU v. 3.0, see http://crd.lbl.gov/~xiaoye/SuperLU/
 * The function returns 0 in case of error, 1 if successfull
 * Note also that the contents of A and B are destroyed on exit.
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
 *
 */
int splm_Axb_SuperLU(struct splm_ccsm *inA, double *inB, void **state, int what, double *x)
{
SuperMatrix L, U; /* factors L & U */
SCformat *Lstore;
NCformat *Ustore;
int nrhs=1, info, ret, verbose=0;
double *aux, ferr, berr, rpg, rcond;
mem_usage_t mem_usage;
SuperLUStat_t statistics;
register int i;
char equed;
struct solvstat *stat;
    
  if(inA && inA->nr!=inA->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_SuperLU()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)malloc(sizeof(struct solvstat));
      if(!stat){
        fprintf(stderr, "Memory allocation request for \"stat\" failed in splm_Axb_SuperLU()!\n");
        exit(1);
      }

      /* Set the default input options */
      stat->options.Fact=DOFACT;
      stat->options.Equil=YES;
      stat->options.ColPerm=COLAMD;
      stat->options.DiagPivotThresh=1.0;
      stat->options.Trans=NOTRANS;
      stat->options.IterRefine=NOREFINE;
      stat->options.SymmetricMode=NO;
      stat->options.PivotGrowth=NO;
      stat->options.ConditionNumber=NO;
      stat->options.PrintStat=YES;
      //set_default_options(&options);

      /* Now modify the default options to use the symmetric mode. */
      stat->options.SymmetricMode=YES; // A is expected to be diagonally dominant
      stat->options.ColPerm=MMD_AT_PLUS_A;
      stat->options.DiagPivotThresh=0.001;
      stat->options.PrintStat=NO;
      stat->options.IterRefine=SINGLE;

      stat->work=NULL; stat->lwork=0;

      /* Create a matrix in Harwell-Boeing format. */
      //m=inA->nr; n=inA->nc;

      dCreate_CompCol_Matrix(&(stat->A), inA->nr, inA->nc, inA->nnz, inA->val, inA->rowidx, inA->colptr, SLU_NC, SLU_D, SLU_GE);
    
      /* if(!(stat->rhs=doubleMalloc(inA->nr * nrhs))) ABORT("Malloc fails for rhs[] in splm_Axb_SuperLU()!"); */
      dCreate_Dense_Matrix(&(stat->B), inA->nr, nrhs, inB /*stat->rhs*/, inA->nr, SLU_DN, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&(stat->X), inA->nr, nrhs, x, inA->nr, SLU_DN, SLU_D, SLU_GE);

      if(!(stat->perm_c=intMalloc(inA->nc))) ABORT("Malloc fails for perm_c[] in splm_Axb_SuperLU()!");
      if(!(stat->perm_r=intMalloc(inA->nr))) ABORT("Malloc fails for perm_r[] in splm_Axb_SuperLU()!");
      if(!(stat->etree=intMalloc(inA->nc))) ABORT("Malloc fails for etree[] in splm_Axb_SuperLU()!");
      if(!(stat->R=(double *)SUPERLU_MALLOC(inA->nr*sizeof(double)))) ABORT("SUPERLU_MALLOC fails for R[] in splm_Axb_SuperLU()!");
      if(!(stat->C=(double *)SUPERLU_MALLOC(inA->nc*sizeof(double)))) ABORT("SUPERLU_MALLOC fails for C[] in splm_Axb_SuperLU()!");

      /* Initialize the statistics variables. */
      StatInit(&statistics);

#if 0
      /* copy rhs vector */
      aux=(double *)((DNformat *)(stat->B).Store)->nzval;
      for(i=0; i<inA->nr; ++i)
        aux[i]=inB[i];
#endif

      /* first time solution of the linear system */
      dgssvx(&(stat->options), &(stat->A), stat->perm_c, stat->perm_r, stat->etree, &equed, stat->R, stat->C,
                     &L, &U, stat->work, stat->lwork, &(stat->B), &(stat->X), &rpg, &rcond, &ferr, &berr, &mem_usage, &statistics, &info);

      if(info==0 || info==inA->nc+1){
        /* This is how you could access the solution matrix. */
        aux=(double *)((DNformat *)(stat->X).Store)->nzval;
        for(i=0; i<inA->nr; ++i)
          x[i]=aux[i];
        ret=1;

        if(verbose){
          if(stat->options.PivotGrowth) printf("Recip. pivot growth = %e\n", rpg);
          if(stat->options.ConditionNumber) printf("Recip. condition number = %e\n", rcond);

          Lstore=(SCformat *) L.Store;
          Ustore=(NCformat *) U.Store;
          printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
          printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
          printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - inA->nc);
          printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6, mem_usage.expansions);
          if(stat->options.IterRefine){
            printf("Iterative Refinement:\n%8s%16s%16s\n", "Steps", "FERR", "BERR");
            printf("%8d%16e%16e\n", statistics.RefineSteps, ferr, berr);
          }
          fflush(stdout);
        }
      }else{
        printf("dgssvx() error returns INFO= %d\n", info);
        if(info>0 && stat->lwork==-1)
          printf("** Estimated memory: %d bytes\n", info - inA->nc);

        ret=0;
      }

      if(stat->options.PrintStat) StatPrint(&statistics);

      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
      StatFree(&statistics);
      break;
      /* NOTE: we do NOT fall through in this case! */

    case 2:
      stat=(struct solvstat *)(*state); //retrieve solver state

      /* the following ensures correctness even if the address of inA, inB & x changes among calls */
      _dSet_CompCol_Matrix(&(stat->A), inA->nr, inA->nc, inA->nnz, inA->val, inA->rowidx, inA->colptr, SLU_NC, SLU_D, SLU_GE);
      _dSet_Dense_Matrix(&(stat->B), inA->nr, nrhs, inB /*stat->rhs*/, inA->nr, SLU_DN, SLU_D, SLU_GE);
      _dSet_Dense_Matrix(&(stat->X), inA->nr, nrhs, x, inA->nr, SLU_DN, SLU_D, SLU_GE);
                                    
      /* solve the linear system for the N-th time, N>1 */
      stat->options.Fact=SamePattern; // sparsity pattern identical to previous one 
      StatInit(&statistics); /* Initialize the statistics variables. */

#if 0
      /* copy rhs vector */
      aux=(double *)((DNformat *)(stat->B).Store)->nzval;
      for(i=0; i<inA->nr; ++i)
        aux[i]=inB[i];
#endif
    
      dgssvx(&(stat->options), &(stat->A), stat->perm_c, stat->perm_r, stat->etree, &equed, stat->R, stat->C,
                     &L, &U, stat->work, stat->lwork, &(stat->B), &(stat->X), &rpg, &rcond, &ferr, &berr, &mem_usage, &statistics, &info);
     
      if(info==0 || info==stat->A.ncol+1){
        /* This is how you could access the solution matrix. */
        aux=(double *)((DNformat *)(stat->X).Store)->nzval;
#if 0
        for(i=0; i<inA->nr; ++i)
          x[i]=aux[i];
#endif
        memcpy(x, aux, inA->nr*sizeof(double));
        ret=1;

        if(verbose){
          if(stat->options.PivotGrowth) printf("Recip. pivot growth = %e\n", rpg);
          if(stat->options.ConditionNumber) printf("Recip. condition number = %e\n", rcond);

          Lstore = (SCformat *) L.Store;
          Ustore = (NCformat *) U.Store;
          printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
          printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
          printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - inA->nc);
          printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6, mem_usage.expansions);
          if(stat->options.IterRefine){
            printf("Iterative Refinement:\n%8s%16s%16s\n", "Steps", "FERR", "BERR");
            printf("%8d%16e%16e\n", statistics.RefineSteps, ferr, berr);
          }
          fflush(stdout);
        }
      }else{
        printf("dgssvx() error returns INFO= %d\n", info);
        if(info>0 && stat->lwork==-1)
          printf("** Estimated memory: %d bytes\n", info - inA->nc);

        ret=0;
      }

      if(stat->options.PrintStat) StatPrint(&statistics);

      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
      StatFree(&statistics);
      break;

      case 3: /* release memory */
        stat=(struct solvstat *)(*state); //retrieve solver state
        ret=1;
        break;

      default:
        fprintf(stderr, "Unknown mode %d in splm_Axb_SuperLU()!\n", what);
        exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    /* SUPERLU_FREE(stat->rhs); */
    SUPERLU_FREE(stat->perm_r);
    SUPERLU_FREE(stat->perm_c);
    SUPERLU_FREE(stat->etree);
    SUPERLU_FREE(stat->R);
    SUPERLU_FREE(stat->C);
    //Destroy_CompCol_Matrix(&(stat->A));
    SUPERLU_FREE((stat->A).Store);
    /* Destroy_SuperMatrix_Store(&(stat->B)); */
    SUPERLU_FREE((stat->B).Store);
    SUPERLU_FREE((stat->X).Store);
    free(stat);
    *state=NULL;
  }

  return ret;
}

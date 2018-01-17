/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to LDL sparse linear solver
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

#define _LDL_USE_AMD

#include <ldl.h>
#ifdef _LDL_USE_AMD
#include <amd.h>
#endif


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


struct solvstat{ // solver state
  int *Li, *Lp, *Parent, *Lnz, *Flag, *Pattern, *P, *Pinv, dim;
  double *Lx, *D, *Y;
};

/* solves Ax=B (A symmetric) using Davis' sparse LDL factorization.
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
 *    1: initialize parameters, allocate memory, perform symbolic factorization
 *       & solve system; store necessary info in "state"
 *    2: numeric factorization and system solution using previous info in "state"
 *    3: free internal state accessed through "state"
 */
int splm_Axb_LDLP(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j;
int ret;
struct solvstat *stat;

#ifdef _LDL_USE_AMD
double Info[AMD_INFO], *Control=NULL;
#endif

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_LDLP()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->dim=A->nr;

#ifdef _LDL_USE_AMD
      stat->P=(int *)emalloc(stat->dim*sizeof(int));
      stat->Pinv=(int *)emalloc(stat->dim*sizeof(int));
      if((i=amd_order(stat->dim, A->colptr, A->rowidx, stat->P, Control, Info))!=AMD_OK){
        fprintf(stderr, "Call to AMD failed [%d]!\n", i);
        exit(1);
      }

      //amd_control(Control);
      //amd_info(Info);
#else
      stat->P=stat->Pinv=NULL;
#endif /* _LDL_USE_AMD */

      stat->D=(double *)emalloc(stat->dim*sizeof(double));
      stat->Y=(double *)emalloc(stat->dim*sizeof(double));
      stat->Lp=(int *)emalloc((stat->dim+1)*sizeof(int));
      stat->Parent=(int *)emalloc(stat->dim*sizeof(int));
      stat->Lnz=(int *)emalloc(stat->dim*sizeof(int));
      stat->Flag=(int *)emalloc(stat->dim*sizeof(int));
      stat->Pattern=(int *)emalloc(stat->dim*sizeof(int));

      /* factorize A into LDL' */
      LDL_symbolic(stat->dim, A->colptr, A->rowidx, stat->Lp, stat->Parent, stat->Lnz, stat->Flag, stat->P, stat->Pinv);

      /* prepare for solving linear system */
      stat->Lx=(double *)emalloc(stat->Lp[stat->dim]*sizeof(double));
      stat->Li=(int *)emalloc(stat->Lp[stat->dim]*sizeof(int));
      /* fall through */

    case 2: /* numerical factorization & system solution */
      stat=(struct solvstat *)(*state); //retrieve solver state

      j=LDL_numeric(stat->dim, A->colptr, A->rowidx, A->val, stat->Lp, stat->Parent, stat->Lnz, stat->Li, stat->Lx,
                               stat->D, stat->Y, stat->Pattern, stat->Flag, stat->P, stat->Pinv);
      if(j==stat->dim){
        /* solve Ax=B, storing the solution in x */
#ifdef _LDL_USE_AMD
        LDL_perm(stat->dim, stat->Y, B, stat->P);                       /* Y = PB */
        LDL_lsolve(stat->dim, stat->Y, stat->Lp, stat->Li, stat->Lx);   /* Y = L\Y */
        LDL_dsolve(stat->dim, stat->Y, stat->D);                        /* Y = D\Y */
        LDL_ltsolve(stat->dim, stat->Y, stat->Lp, stat->Li, stat->Lx);  /* Y = L'\Y */
        LDL_permt(stat->dim, x, stat->Y, stat->P);                      /* x = P'Y */
#else
        for(i=0; i<stat->dim; ++i) x[i]=B[i]; // copy RHS to x to avoid destroying it

        LDL_lsolve(stat->dim, x, stat->Lp, stat->Li, stat->Lx);   /* x = L\x */
        LDL_dsolve(stat->dim, x, stat->D);                        /* x = D\x */
        LDL_ltsolve(stat->dim, x, stat->Lp, stat->Li, stat->Lx);        /* x = L'\x */
#endif /* _LDL_USE_AMD */
        ret=1;
      } else{
        fprintf(stderr, "LDL_numeric() failed, D(%d,%d) is zero\n", j, j);
        ret=0;
      }
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_LDLP()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->Lx);
    free(stat->D);
    free(stat->Y);
    free(stat->Li);
    free(stat->Lp);
    free(stat->Parent);
    free(stat->Lnz);
    free(stat->Flag);
    free(stat->Pattern);
#ifdef _LDL_USE_AMD
    free(stat->P);
    free(stat->Pinv);
#endif
    free(stat);
    *state=NULL;
  }

  return ret;
}

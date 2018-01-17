/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MA57 sparse linear solver
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

extern int ma57id_(double *cntl, int *icntl);

extern int ma57ad_(int *n, int *ne, int *irn, int *jcn, int *lkeep, int *keep,
        int *iwork, int *icntl, int *info, double *rinfo);

extern int ma57bd_(int *n, int *ne, double *a, double *fact, int *lfact, int *ifact, int *lifact,
        int *lkeep, int *keep, int *ppos, int *icntl, double *cntl, int *info, double *rinfo);

extern int ma57cd_(int *job, int *n, double *fact, int *lfact, int *ifact, int *lifact,
        int *nrhs, double *rhs, int *lrhs, double *w, int *lw, int *iw1, int *icntl, int *info);


static void ma57_errmsg(int *info, double *cntl)
{
  if(info[0]>0) fprintf(stderr, "MA57 Warning: ");

  switch(info[0]){
    case -1: fprintf(stderr, "Value of N is out of range (%d)\n", info[1]); 
             break;

    case -2: fprintf(stderr, "Value of NE is out of range (%d)\n", info[1]);
             break;

    case -3: fprintf(stderr, "Failure due to insufficient REAL space on call to MA57BD. %d may suffice for LFACT (currently %d)\n", info[16], info[1]);
             fprintf(stderr, "User can allocate a larger array and copy the contents of FACT into it using MA57ED, and recall MA57BD\n");
             break;

    case -4: fprintf(stderr, "Failure due to insufficient INTEGER space on call to MA57BD. %d may suffice for for LIFACT (currently %d)\n", info[17], info[1]);
             fprintf(stderr, "User allocate a larger array and copy the contents of IFACT into it using MA57ED, and recall MA57BD\n");
             break;

    case -5: fprintf(stderr, "A pivot with magnitude less than or equal to %g was found at pivot step %d when calling "
                    "MA57BD with ICNTL(7) = 2 or 3, or the correction obtained when using matrix modification does not give a "
                    "pivot greater than %g when ICNTL(7) = 4\n", cntl[1], info[1], cntl[1]);
             break;

    case -6: fprintf(stderr, "A change of sign of pivots has been detected when ICNTL(7) = 2. Change was detected at pivot step %d on call to MA57BD\n", info[1]);
             break;

    case -7: fprintf(stderr, "Either LNEW < LFACT or LINEW < LIFACT on call to MA57ED. LNEW/LINEW is %d\n", info[1]);
             break;

    case -8: fprintf(stderr, "Iterative refinement fails to converge in specified number of iterations on call to MA57DD\n");
             break;

    case -9: fprintf(stderr, "Error in permutation array when ICNTL(6)=1 on call to MA57A/AD. %d is first component at which error was detected\n", info[1]);
             break;

    case -10: fprintf(stderr, "Value of ICNTL(7) out of range on call to MA57BD. Value given is %d\n", info[1]);
              break;

    case -11: fprintf(stderr, "LRHS < N on call to MA57C/CD. LRHS is %d\n", info[1]);
              break;

    case -12: fprintf(stderr, "Invalid value for JOB on call to MA57DD (%d)\n", info[1]);
              break;

    case -13: fprintf(stderr, "Invalid value of ICNTL(9) on call to MA57DD (%d)\n", info[1]);
              break;

    case -14: fprintf(stderr, "Failure of MC71AD on call to MA57DD with ICNTL(10)> 0.\n");
              break;

    case -15: fprintf(stderr, "LKEEP less than 5*N+NE+MAX(N,NE)+42 on call to MA57AD or MA57BD (%d)\n", info[1]);
              break;

    case -16: fprintf(stderr, "NRHS less than 1 on call to MA57C/CD (%d)\n", info[1]);
              break;


    case -17: fprintf(stderr, "LWORK too small on entry to MA57C/CD, minimum length required is %d\n", info[1]);
              break;

    case -18: fprintf(stderr, "An ordering using MeTiS was requested in the call to MA57AD but the MeTiS package was not linked in\n");
              break;

    case +1: fprintf(stderr, "Index (in IRN or JCN) out of range on call to MA57AD or MA57DD. Action taken by subroutine is to ignore "
                              "any such entries and continue. %d out-of-range entries", info[2]);
             break;

    case +2: fprintf(stderr, "Duplicate indices on call to MA57AD or MA57DD. Action taken by subroutine is to keep the duplicates and "
                             "then to sum corresponding reals when MA57BD is called. %d duplicate entries\n", info[3]);
             break;

    case +3: fprintf(stderr, "Both out-of-range indices and duplicates exist.\n");
             break;

    case +4: fprintf(stderr, "Matrix is rank deficient on exit from MA57BD. In this case, a decomposition will still have been produced that "
                              "will enable the subsequent solution of consistent equations. Rank of the factorized matrix is %d\n", info[24]);
             break;

    case +5: fprintf(stderr, "Pivots have different signs when factorizing a supposedly definite matrix (ICNTL(7) = 3) on call to MA57BD. %d sign changes.\n", info[25]);
             break;

    case +8: fprintf(stderr, "During error analysis the infinity norm of the computed solution was found to be zero\n");
             break;

    case +10: fprintf(stderr, "Insufficient real space to complete factorization when MA57BD called with ICNTL(8) != 0. User can copy real "
                              "values to a longer array using MA57ED and recall MA57BD using this longer array to continue the factorization\n");
              break;

    case +11: fprintf(stderr, "Insufficient integer space to complete factorization when MA57BD called with ICNTL(8) != 0. User can copy "
                              "integer values to a longer array using MA57ED and recall MA57BD using this longer array to continue the factorization\n");
              break;

    default: fprintf(stderr, "Unknown error code %d in ma57_errmsg()!\n", info[0]);
  }
}

#define CNTLSZ  5
#define RINFOSZ 20
#define ICNTLSZ 20
#define INFOSZ  40


struct solvstat{ // solver state
  int n, *irn, *jcn, *keep, lkeep, *iwork, ne, icntl[ICNTLSZ], lfact, lifact, *ifact, lwork;
  double *a, cntl[CNTLSZ], *fact, *work;
};

/* solves Ax=B (A symmetric) using HSL' MA57 sparse indefinite solver.
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
int splm_Axb_MA57(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i, j, k;
double rinfo[RINFOSZ];
int info[INFOSZ], ret;
int job, nrhs;
struct solvstat *stat;

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_MA57()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      *state=stat=(struct solvstat *)emalloc(sizeof(struct solvstat));
      stat->n=A->nr;

      /* set default parameters */
      ma57id_(stat->cntl, stat->icntl);
      //stat->icntl[5]=0; // AMD ordering without dense row detection
      //stat->icntl[5]=2; // AMD ordering
      //stat->icntl[5]=3; // minimum degree ordering as in MA27
      stat->icntl[5]=4; // METIS ordering
      stat->icntl[5]=5; // automatic choice between METIS and AMD ordering (default)

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

      stat->lkeep=6*stat->n + stat->ne + ((stat->n > stat->ne)? stat->n : stat->ne) + 42; // N more than minimum
      stat->keep=(int *)emalloc(stat->lkeep*sizeof(int));
      memset(stat->keep, 0, (stat->lkeep)*sizeof(int)); // set to zero
      stat->iwork=(int *)emalloc(5*stat->n*sizeof(int)); // max(N, 5*N)

      /* perform symbolic manipulation of sparsity pattern */
      ma57ad_(&(stat->n), &(stat->ne), stat->irn, stat->jcn, &(stat->lkeep), stat->keep, stat->iwork, stat->icntl, info, rinfo);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA57AD(), return code %d\n", info[0]);

        ma57_errmsg(info, stat->cntl);
        if(info[0]<0) exit(1);
      }

      if(info[12]>10) fprintf(stderr, "%d data compresses performed by MA57AD(), performance may be improved by increasing "
                                      "the length LKEEP of KEEP (currently %d)\n", info[12], stat->lkeep);

      /* prepare for solving linear system */
      stat->a=(double *)emalloc(stat->ne*sizeof(double));
      stat->lfact=(int)(info[8]*1.5); // 50% more than suggested by MA57AD
      stat->fact=(double *)emalloc(stat->lfact*sizeof(double));
      stat->lifact=(int)(info[9]*1.2); // 20% more than suggested by MA57AD
      stat->ifact=(int *)emalloc(stat->lifact*sizeof(int));
      stat->lwork=(int)(stat->n*1.20); // 20% more than needed
      stat->work=(double *)emalloc(stat->lwork*sizeof(double));
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
      ma57bd_(&(stat->n), &(stat->ne), stat->a, stat->fact, &(stat->lfact), stat->ifact, &(stat->lifact), &(stat->lkeep),
                  stat->keep, stat->iwork, stat->icntl, stat->cntl, info, rinfo);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA57BD(), return code %d\n", info[0]);

        ma57_errmsg(info, stat->cntl);
        //exit(1);
        ret=0;
        break;
      }

      /* solve equations */
      job=1;
      nrhs=1;
      ma57cd_(&job, &(stat->n), stat->fact, &(stat->lfact), stat->ifact, &(stat->lifact), &nrhs, x,
              &(stat->n), stat->work, &(stat->lwork), stat->iwork, stat->icntl, info);
      if(info[0]!=0){
        fprintf(stderr, "Unsuccessful termination of MA57CD(), return code %d\n", info[0]);

        ma57_errmsg(info, stat->cntl);
        ret=0;
        break;
      }
      ret=1;
      break;

    case 3: /* release memory */
      stat=(struct solvstat *)(*state); //retrieve solver state
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_MA57()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
    free(stat->work);
    free(stat->ifact);
    free(stat->fact);
    free(stat->a);
    free(stat->iwork);
    free(stat->keep);
    free(stat->jcn);
    free(stat->irn);
    free(stat);
    *state=NULL;
  }

  return ret;
}

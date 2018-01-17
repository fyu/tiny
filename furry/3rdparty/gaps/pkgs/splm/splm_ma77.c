/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Interface to MA77 sparse linear solver
//////  Copyright (C) 2010  Manolis Lourakis (lourakis at ics.forth.gr)
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
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <cholmod.h>


#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"


/* converts a vector of C strings into FORTRAN strings.
 * Based on a similar function from cfortran.h
 */
static char *c2fstrv(char* cstr, char *fstr, int eltlen, int sizeofcstr)
{
register int i,j;

  /* eltlen includes a trailing zero for C strings, Fortran strings lack it. */
  for(i=0; i<sizeofcstr/eltlen; i++){
    for(j=1; j<eltlen && *cstr; j++) *fstr++=*cstr++;
    cstr+=1+eltlen-j;
    for( ; j<eltlen; j++) *fstr++=' ';
  }

  return fstr-sizeofcstr+sizeofcstr/eltlen; 
}

/* solves Ax=B (A symmetric positive definite) using MA77 sparse solver.
 * The function returns 0 in case of error, 1 if successfull
 * Note also that the contents of A and B are not modified on exit.
 *
 * Since this function will be called repeatedly for matrices with
 * identical zero patterns, it will pay off to save some computations
 * by retaining information between calls. Argument "state" is unused here
 * since previously computed information is kept internally by ma77_wrapper().
 * Argument "what" specifies the action to be taken by the solver:
 *    0: perform all steps without relying upon any previous computations;
 *       this is the slowest option overall
 *    1: initialize parameters, allocate memory & perform symbolic factorization
 *       & solve system
 *    2: numeric factorization and system solution using previous info
 *    3: free internal data structures
 */
int splm_Axb_MA77(struct splm_ccsm *A, double *B, void **state, int what, double *x)
{
register int i;
int wwhat, ret;
static char fnames[4][128];
char *fstr=NULL;
unsigned int len=0;
double dummy;
int zero=0;

extern int
ma77_wrapper_(int *what, int *n, char *filenames,
              int *nnz, int *colptr, int *rowidx, double *val, double *x,
              unsigned int filenameslen);

  if(A && A->nr!=A->nc){
    fprintf(stderr, "Non-square matrix passed to splm_Axb_MA77()!\n");
    exit(1);
  }

  switch(what){
    case 0:
    case 1: /* initialize parameters, allocate memory & perform symbolic manipulation */
      tmpnam(fnames[0]); // integers
      tmpnam(fnames[1]); // reals
      tmpnam(fnames[2]); // reals
      tmpnam(fnames[3]); // integers

      len=((sizeof(fnames[0])==1 ? sizeof(fnames) : sizeof(fnames[0])))-1;
      fstr=c2fstrv((char *)fnames, (char *)fnames, len+1, sizeof(fnames));

      //printf("filenames: %s %s %s %s\n", fnames[0], fnames[1], fnames[2], fnames[3]);

      wwhat=1;
      /* FORTRAN array indices start from 1 */
      for(i=A->nnz; i-->0; ) A->rowidx[i]++;
      for(i=A->nr+1; i-->0; ) A->colptr[i]++;

      ret=ma77_wrapper_(&wwhat, &A->nr,(char *)fstr, &A->nnz,A->colptr,A->rowidx,A->val, x, len);

      for(i=A->nnz; i-->0; ) A->rowidx[i]--;
      for(i=A->nr+1; i-->0; ) A->colptr[i]--;

      if(ret) break;

      ret=1;
      /* fall through */

    case 2: /* loading of array elements, numerical factorization & system solution */

      wwhat=2;
      /* FORTRAN array indices start from 1 */
      for(i=A->nnz; i-->0; ) A->rowidx[i]++;
      for(i=A->nr+1; i-->0; ) A->colptr[i]++;

      for(i=A->nr; i-->0; ) x[i]=B[i];

      ret=ma77_wrapper_(&wwhat, &A->nr,(char *)fstr, &A->nnz,A->colptr,A->rowidx,A->val, x, len);

      for(i=A->nnz; i-->0; ) A->rowidx[i]--;
      for(i=A->nr+1; i-->0; ) A->colptr[i]--;

      if(ret) break;

      ret=1;
      break;

    case 3: /* release memory */
      ret=1;
      break;

    default:
      fprintf(stderr, "Unknown mode %d in splm_Axb_MA77()!\n", what);
      exit(1);
  }

  /* cleanup */
  if(what==0 || what==3){
      wwhat=3;
      ret=ma77_wrapper_(&wwhat, &zero,(char *)fstr, &zero,&zero,&zero,&dummy, &dummy, len);

      unlink(fnames[0]);
      unlink(fnames[1]);
      unlink(fnames[2]);
      unlink(fnames[3]);
  }

  return ret;
}

/////////////////////////////////////////////////////////////////////////////////
//// 
////  CRS sparse matrices manipulation routines
////  Copyright (C) 2004-2011  Manolis Lourakis (lourakis at ics.forth.gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

/* NOTE: most functions assume that CRS matrices are "sorted", i.e. for each row,
 * colidx contains indices in ascending order
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "splm.h"

/* allocate a sparse CRS matrix */
void splm_crsm_alloc(struct splm_crsm *sm, int nr, int nc, int nnz)
{
  sm->nr=nr;
  sm->nc=nc;
  sm->nnz=nnz;

  sm->val=(double *)malloc(nnz*sizeof(double));

  sm->colidx=(int *)malloc(nnz*sizeof(int));
  sm->rowptr=(int *)malloc((nr+1)*sizeof(int));
  if(!sm->val || !sm->colidx || !sm->rowptr){
    fprintf(stderr, "memory allocation request failed in splm_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }
}

/* following functions can be used to allocate a CRS matrix in steps;
 * they can be employed when its actual nonzero pattern is unknown beforehand
 */

/* allocate all fields except values */
void splm_crsm_alloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz)
{
  sm->nr=nr;
  sm->nc=nc;
  sm->nnz=nnz;

  sm->val=NULL;
  sm->colidx=(int *)malloc(nnz*sizeof(int));
  sm->rowptr=(int *)malloc((nr+1)*sizeof(int));
  if(!sm->colidx || !sm->rowptr){
    fprintf(stderr, "memory allocation request failed in splm_crsm_alloc_novalues() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }
}

/* allocate values, remaining fields assumed already allocated. #nonzeros assumed contained in sm */
void splm_crsm_alloc_values(struct splm_crsm *sm)
{
  sm->val=(double *)malloc(sm->nnz*sizeof(double));
  if(!sm->val){
    fprintf(stderr, "memory allocation request failed in splm_crsm_alloc_values() [nr=%d, nc=%d, nnz=%d]\n", sm->nr, sm->nc, sm->nnz);
    exit(1);
  }
}

/* new dimensions & nonzero size, does not allocate values */
void splm_crsm_realloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz)
{
  sm->nc=nc;
  sm->nnz=nnz;

  sm->colidx=(int *)realloc(sm->colidx, nnz*sizeof(int));
  if(!sm->colidx){
    fprintf(stderr, "memory allocation request failed in splm_crsm_realloc_novalues() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }

  if(nr!=sm->nr){
    sm->nr=nr;
    sm->rowptr=(int *)realloc(sm->rowptr, (nr+1)*sizeof(int));
    if(!sm->rowptr){
      fprintf(stderr, "memory allocation request failed in splm_crsm_realloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
      exit(1);
    }
  }
}

/* free a sparse CRS matrix */
void splm_crsm_free(struct splm_crsm *sm)
{
  if(sm->val) free(sm->val);
  sm->val=NULL;

  sm->nr=sm->nc=sm->nnz=-1;
  free(sm->colidx);
  free(sm->rowptr);
  sm->colidx=sm->rowptr=NULL;
}

#if 0
static void splm_crsm_print(struct splm_crsm *sm, FILE *fp);
static void splm_crsm_build(struct splm_crsm *sm, int *m, int nr, int nc);

static void splm_crsm_print(struct splm_crsm *sm, FILE *fp)
{
register int i;

  fprintf(fp, "matrix is %dx%d, %d non-zeros\nval: ", sm->nr, sm->nc, sm->nnz);
  for(i=0; i<sm->nnz; ++i)
    fprintf(fp, "%g ", sm->val[i]);
  fprintf(fp, "\ncolidx: ");
  for(i=0; i<sm->nnz; ++i)
    fprintf(fp, "%d ", sm->colidx[i]);
  fprintf(fp, "\nrowptr: ");
  for(i=0; i<=sm->nr; ++i)
    fprintf(fp, "%d ", sm->rowptr[i]);
  fprintf(fp, "\n");
}

/* build an int sparse CRS matrix from a dense one. intended to serve as an example for sm creation */
static void splm_crsm_build(struct splm_crsm *sm, int *m, int nr, int nc)
{
int nnz;
register int i, j, k;

  /* count nonzeros */
  for(i=nnz=0; i<nr; ++i)
    for(j=0; j<nc; ++j)
      if(m[i*nc+j]!=0) ++nnz;

  splm_crsm_alloc(sm, nr, nc, nnz);

  /* fill up the sm structure */
  for(i=k=0; i<nr; ++i){
    sm->rowptr[i]=k;
    for(j=0; j<nc; ++j)
      if(m[i*nc+j]!=0){
        sm->val[k]=m[i*nc+j];
        sm->colidx[k++]=j;
      }
  }
  sm->rowptr[nr]=nnz;
}
#endif

/* returns the index of the (i, j) element. No bounds checking! */
int splm_crsm_elmidx(struct splm_crsm *sm, int i, int j)
{
register int low, high, mid, diff, *Colidx;

  low=sm->rowptr[i];
  high=sm->rowptr[i+1]-1;
  Colidx=sm->colidx;

  /* binary search for finding the element at column j */
  while(low<=high){
    if(j<Colidx[low] || j>Colidx[high]) return -1; /* not found */
    
    mid=(low + high) >> 1; //(low+high)/2;
    //mid=low+((high-low)>>1); /* ensures no index overflows */
    diff=j-Colidx[mid];
    if(diff<0)
        high=mid-1;
    else if(diff>0)
        low=mid+1;
    else
      return mid;
  }

  return -1; /* not found */
}

/* returns the row of element with index idx; inverse of splm_crsm_elmidx().
 * Note that the element's column is simply sm->colidx[idx]
 */
int splm_crsm_elmrow(struct splm_crsm *sm, int idx)
{
register int low, high, mid, diff, *Rowptr;
int leq=-1;

  low=0;
  high=sm->nr;
  Rowptr=sm->rowptr;

  /* binary search for finding the row whose first element is
   * at an index that is the largest value smaller than idx
   */
  while(low<=high){
    mid=(low + high) >> 1; //(low+high)/2;
    //mid=low+((high-low)>>1); /* ensures no index overflows */
    diff=idx-Rowptr[mid];
    if(diff<0)
      high=mid-1;
    else if(diff>0){
      low=mid+1;
      leq=mid;
    }
    else
      return mid; /* this is the 1st element at row mid */
  }

  return leq; /* element is between Rowptr[leq], Rowptr[leq+1] */
}

/* returns the number of nonzero elements in row i and
 * fills up the vidxs and jidxs arrays with the val and column
 * indexes of the elements found, respectively.
 * vidxs and jidxs are assumed preallocated and of max. size sm->nc
 */
int splm_crsm_row_elmidxs(struct splm_crsm *sm, int i, int *vidxs, int *jidxs)
{
register int j, k;

  for(j=sm->rowptr[i], k=0; j<sm->rowptr[i+1]; ++j, ++k){
    vidxs[k]=j;
    jidxs[k]=sm->colidx[j];
  }

  return k;
}

/* returns the maximum number of nonzero elements across all rows */
int splm_crsm_row_maxnelms(struct splm_crsm *sm)
{
register int i, n, max;

  for(i=sm->nr, max=-1; i-->0;  )
    if((n=sm->rowptr[i+1]-sm->rowptr[i])>max) max=n;

  return max;
}

/* returns the number of nonzero elements in col j and
 * fills up the vidxs and iidxs arrays with the val and row
 * indexes of the elements found, respectively.
 * vidxs and iidxs are assumed preallocated and of max. size sm->nr
 */
int splm_crsm_col_elmidxs(struct splm_crsm *sm, int j, int *vidxs, int *iidxs)
{
register int *nextrowptr=sm->rowptr+1;
register int i, l;
register int low, high, mid, diff;

  for(i=l=0; i<sm->nr; ++i){
    low=sm->rowptr[i];
    high=nextrowptr[i]-1;

    /* binary search attempting to find an element at column j */
    while(low<=high){
      if(j<sm->colidx[low] || j>sm->colidx[high]) break; /* not found */

      mid=(low + high) >> 1; //(low+high)/2;
      //mid=low+((high-low)>>1); /* ensures no index overflows */
      diff=j-sm->colidx[mid];
      if(diff<0)
          high=mid-1;
      else if(diff>0)
          low=mid+1;
      else{ /* found */
        vidxs[l]=mid;
        iidxs[l++]=i;
        break;
      }
    }
  }

  return l;
}

#if 0
/* compute y=A*x for a matrix A in CRS format and a vector x */
void splm_crsm_vec_multd(struct splm_crsm *A, double *x, double *y)
{
register int i, j, k1, k2;
register double sum;

  for(i=0; i<A->nr; ++i){
    k1=A->rowptr[i];
    k2=A->rowptr[i+1]-1;
    for(j=k1, sum=0.0; j<=k2; ++j)
      sum+=A->val[j]*x[A->colidx[j]]; // dotprod(val(k1:k2), x(colidx(k1:k2)))
    y[i]=sum;
  }
}
#endif

/* convert a matrix from the CRS format to CCS. If crs->val is NULL, only the nonzero pattern is converted */
void splm_crsm2ccsm(struct splm_crsm *crs, struct splm_ccsm *ccs)
{
register int i, j, k, l;
int nr, nc, nnz, jmax;
int *colidx, *rowptr, *rowidx, *colptr;
int *colcounts; // counters for the number of nonzeros in each column

  nr=crs->nr; nc=crs->nc;
  nnz=crs->nnz;

  /* allocate only if there are values in crs and ccs->val is NULL */
  if(crs->val && !ccs->val) splm_ccsm_alloc(ccs, nr, nc, nnz);

  ccs->nr=nr; ccs->nc=nc; // ensure that ccs has the correct dimensions

  if((colcounts=(int *)calloc(nc, sizeof(int)))==NULL){ /* initialize to zero */
    fprintf(stderr, "memory allocation request failed in splm_crsm2ccsm() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }

  colidx=crs->colidx; rowptr=crs->rowptr;
  rowidx=ccs->rowidx; colptr=ccs->colptr;

  /* 1st pass: count #nonzeros in each column */
  /*
  for(i=0; i<nr; ++i)
    for(j=rowptr[i]; j<rowptr[i+1]; ++j)
  */
  for(j=rowptr[nr]; j-->0;  )
    ++(colcounts[colidx[j]]);

  /* 2nd pass: copy every nonzero to its right position into the CCS structure */
  for(j=k=0; j<nc; ++j){
    colptr[j]=k;
    k+=colcounts[j];
    colcounts[j]=0; // clear to avoid separate loop below
  }
  colptr[nc]=nnz;

  /* colcounts[j] will count the #nonzeros in col. j seen before the current row */

  if(crs->val){ // are there any values to copy?
    register double *crsv, *ccsv;

    crsv=crs->val; ccsv=ccs->val;
    for(i=0; i<nr; ++i){
      jmax=rowptr[i+1];
      for(j=rowptr[i]; j<jmax; ++j){
        l=colidx[j];
        k=colptr[l];
        k+=colcounts[l]++;

        rowidx[k]=i;

        /* copy values */
        ccsv[k]=crsv[j];
      }
    }
  }
  else{ // no values, copy just structure
    for(i=0; i<nr; ++i){
      jmax=rowptr[i+1];
      for(j=rowptr[i]; j<jmax; ++j){
        l=colidx[j];
        k=colptr[l];
        k+=colcounts[l]++;

        rowidx[k]=i;
      }
    }
  }

  free(colcounts);
}

/* sorts the column indices in every row using radix sort */
void splm_crsm_row_sort(struct splm_crsm *sm)
{
register int i, j, k;
int *rowptr, *colidx, *tmpi, nc;
double *val, *tmpv;

  nc=sm->nc;
  rowptr=sm->rowptr; colidx=sm->colidx;
  val=sm->val;

  tmpv=(double *)malloc(nc*sizeof(double));
  tmpi=(int *)malloc(nc*sizeof(int));
  if(!tmpv || !tmpi){
    fprintf(stderr, "memory allocation request failed in splm_crsm_row_sort()\n");
    exit(1);
  }

  for(j=nc; j-->0;  ) tmpi[j]=0;

  for(i=sm->nr; i-->0;  ){
    for(j=rowptr[i]; j<rowptr[i+1]; ++j){
      k=colidx[j];
      tmpv[k]=val[j];
      tmpi[k]=1;
    }

    for(j=0, k=rowptr[i]; j<nc; ++j){
      if(tmpi[j]==0) continue;

      tmpi[j]=0; // reset
      val[k]=tmpv[j];
      colidx[k++]=j;
    }
  }

  free(tmpi);
  free(tmpv);
}

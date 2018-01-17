/////////////////////////////////////////////////////////////////////////////////
////// 
//////  CCS sparse matrices manipulation routines
//////  Copyright (C) 2005-2011  Manolis Lourakis (lourakis at ics.forth.gr)
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

/* NOTE: most functions assume that CCS matrices are "sorted", i.e. for each column,
 * rowidx contains indices in ascending order
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "splm.h"
#include "compiler.h"

/* allocate a sparse CCS matrix */
void splm_ccsm_alloc(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
  sm->nr=nr;
  sm->nc=nc;
  sm->nnz=nnz;

  sm->val=(double *)malloc(nnz*sizeof(double));
  sm->rowidx=(int *)malloc(nnz*sizeof(int));
  sm->colptr=(int *)malloc((nc+1)*sizeof(int));
  if(!sm->val || !sm->rowidx || !sm->colptr){
    fprintf(stderr, "memory allocation request failed in splm_ccsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }
}

/* following functions can be used to allocate a CCS matrix in steps;
 * they can be employed when its actual nonzero pattern is unknown beforehand
 */

/* allocate all fields except values */
void splm_ccsm_alloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
  sm->nr=nr;
  sm->nc=nc;
  sm->nnz=nnz;

  sm->val=NULL;
  sm->rowidx=(int *)malloc(nnz*sizeof(int));
  sm->colptr=(int *)malloc((nc+1)*sizeof(int));
  if(!sm->rowidx || !sm->colptr){
    fprintf(stderr, "memory allocation request failed in splm_ccsm_alloc_novalues() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }
}

/* allocate values, remaining fields assumed already allocated. #nonzeros assumed contained in sm */
void splm_ccsm_alloc_values(struct splm_ccsm *sm)
{
  sm->val=(double *)malloc(sm->nnz*sizeof(double));
  if(!sm->val){
    fprintf(stderr, "memory allocation request failed in splm_ccsm_alloc_values() [nr=%d, nc=%d, nnz=%d]\n", sm->nr, sm->nc, sm->nnz);
    exit(1);
  }
}

/* new dimensions & nonzero size, does not allocate values */
void splm_ccsm_realloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
  sm->nr=nr;
  sm->nnz=nnz;

  sm->rowidx=(int *)realloc(sm->rowidx, nnz*sizeof(int));
  if(!sm->rowidx){
    fprintf(stderr, "memory allocation request failed in splm_ccsm_realloc_novalues() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }

  if(nc!=sm->nc){
    sm->nc=nc;
    sm->colptr=(int *)realloc(sm->colptr, (nc+1)*sizeof(int));
    if(!sm->colptr){
      fprintf(stderr, "memory allocation request failed in splm_ccsm_realloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
      exit(1);
    }
  }
}

/* free a sparse CCS matrix */
void splm_ccsm_free(struct splm_ccsm *sm)
{
  if(sm->val) free(sm->val);
  sm->val=NULL;

  sm->nr=sm->nc=sm->nnz=-1;
  free(sm->rowidx);
  free(sm->colptr);
  sm->rowidx=sm->colptr=NULL;
}

/* returns the index of the (i, j) element. No bounds checking! */
int splm_ccsm_elmidx(struct splm_ccsm *sm, int i, int j)
{
register int low, high, mid, *Rowidx, diff;

  low=sm->colptr[j];
  high=sm->colptr[j+1]-1;
  Rowidx=sm->rowidx;

  /* binary search for finding the element at row i */
  while(low<=high){
    //if(i<Rowidx[low] || i>Rowidx[high]) return -1; /* not found */
    
    mid=(low + high) >> 1; //(low+high)/2;
    //mid=low+((high-low)>>1); /* ensures no index overflows */
    diff=i-Rowidx[mid];
    if(diff<0)
        high=mid-1;
    else if(diff>0)
        low=mid+1;
    else
      return mid;
  }

  return -1; /* not found */
}

/* returns the column of element with index idx; inverse of splm_ccsm_elmidx().
 * Note that the element's row is simply sm->rowidx[idx]
 */
int splm_ccsm_elmcol(struct splm_ccsm *sm, int idx)
{
register int low, high, mid, diff, *Colptr;
int leq=-1;

  low=0;
  high=sm->nc;
  Colptr=sm->colptr;

  /* binary search for finding the column whose first element is
   * at an index that is the largest value smaller than idx
   */
  while(low<=high){
    mid=(low + high) >> 1; //(low+high)/2;
    //mid=low+((high-low)>>1); /* ensures no index overflows */
    diff=idx-Colptr[mid];
    if(diff<0)
      high=mid-1;
    else if(diff>0){
      low=mid+1;
      leq=mid;
    }
    else
      return mid; /* this is the 1st element at column mid */
  }

  return leq; /* element is between Colptr[leq], Colptr[leq+1] */
}

/* returns the number of nonzero elements in row i and
 * fills up the vidxs and jidxs arrays with the val and column
 * indexes of the elements found, respectively.
 * vidxs and jidxs are assumed preallocated and of max. size sm->nc
 */
int splm_ccsm_row_elmidxs(struct splm_ccsm *sm, int i, int *vidxs, int *jidxs)
{
register int *nextcolptr=sm->colptr+1;
register int j, l;
register int low, high, mid, diff;

  for(j=l=0; j<sm->nc; ++j){
    low=sm->colptr[j];
    high=nextcolptr[j]-1;

    /* binary search attempting to find an element at row i */
    while(low<=high){
      //if(i<sm->rowidx[low] || i>sm->rowidx[high]) break; /* not found */

      mid=(low + high) >> 1; //(low+high)/2;
      //mid=low+((high-low)>>1); /* ensures no index overflows */
      diff=i-sm->rowidx[mid];
      if(diff<0)
          high=mid-1;
      else if(diff>0)
          low=mid+1;
      else{ /* found */
        vidxs[l]=mid;
        jidxs[l++]=j;
        break;
      }
    }
  }

  return l;
}

/* returns the number of nonzero elements in col j and
 * fills up the vidxs and iidxs arrays with the val and row
 * indexes of the elements found, respectively.
 * vidxs and iidxs are assumed preallocated and of max. size sm->nr
 */
int splm_ccsm_col_elmidxs(struct splm_ccsm *sm, int j, int *vidxs, int *iidxs)
{
register int i, k;
int low, high, *rowidx=sm->rowidx;

  low=sm->colptr[j];
  high=sm->colptr[j+1];
  for(i=low, k=0; i<high; ++i, ++k){
    vidxs[k]=i;
    iidxs[k]=rowidx[i];
  }

  return k;
}

/* returns the maximum number of nonzero elements across all columns */
int splm_ccsm_col_maxnelms(struct splm_ccsm *sm)
{
register int j, n, max;

  for(j=sm->nc, max=-1; j-->0;  )
    if((n=sm->colptr[j+1]-sm->colptr[j])>max) max=n;

  return max;
}

#if 0
/* returns the minimum number of nonzero elements across all columns */
int splm_ccsm_col_minnelms(struct splm_ccsm *sm)
{
register int j, n, min;

  j=sm->nc-1;
  min=sm->colptr[sm->nc]-sm->colptr[j];
  while(j-->0)
    if((n=sm->colptr[j+1]-sm->colptr[j])<min) min=n;

  return min;
}

static void splm_ccsm_print(struct splm_ccsm *sm, FILE *fp);
static void splm_ccsm_build(struct splm_ccsm *sm, int *m, int nr, int nc);

static void splm_ccsm_print(struct splm_ccsm *sm, FILE *fp)
{
register int i;

  fprintf(fp, "matrix is %dx%d, %d non-zeros\nval: ", sm->nr, sm->nc, sm->nnz);

  for(i=0; i<sm->nnz; ++i)
    fprintf(fp, "%g ", sm->val[i]);
  fprintf(fp, "\nrowidx: ");
  for(i=0; i<sm->nnz; ++i)
    fprintf(fp, "%d ", sm->rowidx[i]);
  fprintf(fp, "\ncolptr: ");
  for(i=0; i<=sm->nc; ++i)
    fprintf(fp, "%d ", sm->colptr[i]);
  fprintf(fp, "\n");
}

/* build an int sparse CCS matrix from a dense one. intended to serve as an example for sm creation */
static void splm_ccsm_build(struct splm_ccsm *sm, int *m, int nr, int nc)
{
int nnz;
register int i, j, k;

  /* count nonzeros */
  for(i=nnz=0; i<nr; ++i)
    for(j=0; j<nc; ++j)
      if(m[i*nc+j]!=0.0) ++nnz;

  splm_ccsm_alloc(sm, nr, nc, nnz);

  /* fill up the sm structure */
  for(j=k=0; j<nc; ++j){
    sm->colptr[j]=k;
    for(i=0; i<nr; ++i)
      if(m[i*nc+j]!=0){
        sm->val[k]=m[i*nc+j];
        sm->rowidx[k++]=i;
      }
  }
  sm->colptr[nc]=nnz;
}
#endif

/* drop the first ncols columns of a CCS array. This is achieved by adjusting
 * the colptr and shifting the val & rowidx arrays. Note since ncols is expected
 * to be much smaller compared to A->nc, no memory reallocation takes place. Thus,
 * the new array occupies the same amount of storage as the original one, without
 * actually using it all. Note also that the the val & rowidx pointers are moved
 * forward so that they point to the first nonzero of column ncols, whereas colptr
 * is advanced by "ncols" columns. Additionally, the contents of colptr are recuced
 * to account for the fact that column "ncols" becomes column 0.
 *
 * Returns the number of nonzeros found in the dropped columns, which is
 * necessary for undoing the drop action later.
 */
int splm_ccsm_drop_cols(struct splm_ccsm *A, int ncols)
{
register int i;
int ncnnz;

  if(ncols==0) return 0; // nothing to do

  /* number of nonzeros in the first "ncols" columns */
  ncnnz=A->colptr[ncols];

  /* adjust the contents of the column pointers */
  for(i=ncols; i<A->nc+1; ++i)
    A->colptr[i]-=ncnnz;

  A->nnz-=ncnnz;
  A->nc-=ncols;

  /* advance pointers */
  A->colptr+=ncols;
  A->rowidx+=ncnnz;
  A->val+=ncnnz;

  return ncnnz;
}

/* undo the action of dropping "ncols" columns (which contained "ncnnz" nonzero elements */
void splm_ccsm_restore_cols(struct splm_ccsm *A, int ncols, int ncnnz)
{
register int i;

  if(ncols==0) return; // nothing to do

  /* decrement pointers */
  A->colptr-=ncols;
  A->rowidx-=ncnnz;
  A->val-=ncnnz;

  A->nnz+=ncnnz;
  A->nc+=ncols;

  /* adjust the contents of the column pointers */
  for(i=ncols; i<A->nc+1; ++i)
    A->colptr[i]+=ncnnz;
}

/* convert a matrix from the CCS format to CRS. If ccs->val is NULL, only the nonzero pattern is converted */
void splm_ccsm2crsm(struct splm_ccsm *ccs, struct splm_crsm *crs)
{
register int i, j, k, l;
int nr, nc, nnz;
int *colidx, *rowptr, *rowidx, *colptr;
int *rowcounts; // counters for the number of nonzeros in each row

  nr=ccs->nr; nc=ccs->nc;
  nnz=ccs->nnz;

  /* allocate only if there are values in ccs and crs->val is NULL */
  if(ccs->val && !crs->val) splm_crsm_alloc(crs, nr, nc, nnz);

  crs->nr=nr; crs->nc=nc; // ensure that crs has the correct dimensions

  if((rowcounts=(int *)calloc(nr, sizeof(int)))==NULL){ /* initialize to zero */
    fprintf(stderr, "memory allocation request failed in splm_ccsm2crsm() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }

  colidx=crs->colidx; rowptr=crs->rowptr;
  rowidx=ccs->rowidx; colptr=ccs->colptr;

  /* 1st pass: count #nonzeros in each row */
/***
  for(j=0; j<nc; ++j)
    for(i=colptr[j]; i<colptr[j+1]; ++i)
      ++(rowcounts[rowidx[i]]);
***/
  for(i=colptr[nc]; i-->0;  )
    ++(rowcounts[rowidx[i]]);

  /* 2nd pass: copy every nonzero to its right position into the CRS structure */
  for(i=k=0; i<nr; ++i){
    rowptr[i]=k;
    k+=rowcounts[i];
    rowcounts[i]=0; // clear to avoid separate loop below
  }
  rowptr[nr]=nnz;

  /* rowcounts[i] will count the #nonzeros in row i seen before the current column */
/****
  for(i=nr; i-->0;  )
    rowcounts[i]=0;
****/

  if(ccs->val){ // are there any values to copy?
    register double *crsv, *ccsv;

    crsv=crs->val; ccsv=ccs->val;
    for(j=0; j<nc; ++j){
      for(i=colptr[j]; i<colptr[j+1]; ++i){
        l=rowidx[i];
        k=rowptr[l];
        k+=rowcounts[l]++;

        colidx[k]=j;

        /* copy values */
        crsv[k]=ccsv[i];
      }
    }
  }
  else{ // no values, copy just structure
    for(j=0; j<nc; ++j){
      for(i=colptr[j]; i<colptr[j+1]; ++i){
        l=rowidx[i];
        k=rowptr[l];
        k+=rowcounts[l]++;

        colidx[k]=j;
      }
    }
  }

  free(rowcounts);
}

/* sorts the row indices in every column using radix sort */
void splm_ccsm_col_sort(struct splm_ccsm *sm)
{
register int i, j, k;
int *colptr, *rowidx, *tmpi, nr;
double *val, *tmpv;

  nr=sm->nr;
  colptr=sm->colptr; rowidx=sm->rowidx;
  val=sm->val;

  tmpv=(double *)malloc(nr*sizeof(double));
  tmpi=(int *)malloc(nr*sizeof(int));
  if(!tmpv || !tmpi){
    fprintf(stderr, "memory allocation request failed in splm_ccsm_col_sort()\n");
    exit(1);
  }

  for(i=nr; i-->0;  ) tmpi[i]=0;

  for(j=sm->nc; j-->0;  ){
    for(i=colptr[j]; i<colptr[j+1]; ++i){
      k=rowidx[i];
      tmpv[k]=val[i];
      tmpi[k]=1;
    }

    for(i=0, k=colptr[j]; i<nr; ++i){
      if(tmpi[i]==0) continue;

      tmpi[i]=0; // reset
      val[k]=tmpv[i];
      rowidx[k++]=i;
    }
  }

  free(tmpi);
  free(tmpv);
}

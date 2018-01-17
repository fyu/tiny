/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Set up and compute CCS/CRS approximate sparse Hessians
//////  Copyright (C) 2005-2010  Manolis Lourakis (lourakis at ics.forth.gr)
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

#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"

#define USE_BTA_MAP // faster splm_setup_AtA_ccsA()/splm_setup_AtA_crsA at the cost of some extra memory

#define USE_LINSRCH // use linear search when binsearching short arrays?
#undef USE_MKL_ATX // use MKL for computing A*x for a CCS array A?

#ifdef USE_ATA_QUADS
//#define USE_VBYTE_COMPRESSION
#define USE_S9_COMPRESSION
#endif

#ifdef  USE_LINSRCH
#define LINSRCH_MAXLEN  64 // max. array size for linear vs. binary search
#endif  /* USE_LINSRCH */
static int binsearch_geq(const int *mat, int n, int val)
{
register int low=0, high=n, mid;

  if(val>mat[high-1]) return -1; /* no larger element */

#ifdef USE_LINSRCH
  /* binary search while the list is big enough */
  while(low+LINSRCH_MAXLEN<high){
#else 
  while(low<high){
#endif /* USE_LINSRCH */
    mid=(low+high)>>1; //(low+high)/2;
    //mid=low+((high-low)>>1); /* ensure no index overflows */
    if(val>mat[mid])
      low=mid+1;
    else
      high=mid;
  }

#ifndef USE_LINSRCH
  /* at this point low==high */
  return high;
#else /* USE_LINSRCH */
  /* if present, val must be in the short slice mat[low]...mat[high-1],
   * which is searched linearly using an unrolled loop. Note that val
   * is known to be <= mat[high-1], so there is no need to check whether
   * mid is at the end of the array in the loop below
   */
  for(mid=low; ; mid+=8){
    if(mat[mid  ]>=val) return mid;
    if(mat[mid+1]>=val) return mid+1;
    if(mat[mid+2]>=val) return mid+2;
    if(mat[mid+3]>=val) return mid+3;
    if(mat[mid+4]>=val) return mid+4;
    if(mat[mid+5]>=val) return mid+5;
    if(mat[mid+6]>=val) return mid+6;
    if(mat[mid+7]>=val) return mid+7;
  }

  /* this point is never reached */
#endif /* !USE_LINSRCH */
}
#undef LINSRCH_MAXLEN

/* bit array operations */
#define SPLM_BTA_SET(map, elt)         (map)[(elt)>>3]|=1<<((elt)&7)   // map[(elt)/8]|=1<<(elt%8)
#define SPLM_BTA_ISSET(map, elt)       (((map)[(elt)>>3]&(1<<((elt)&7)))!=0)
/* as SPLM_BTA_ISSET() but avoids some bit operations if whole byte is zero; useful for very sparse bit arrays */
#define SPLM_BTA_ISSET_SPRS(map, elt)  ((map)[(elt)>>3]!=0 && ((map)[(elt)>>3]&(1<<((elt)&7)))!=0) // short-circuiting
#define SPLM_BTA_ISSET_SPRS2(map, aux, elt)  ((aux)=(map)[(elt)>>3], ((aux)!=0 && ((aux)&(1<<((elt)&7)))!=0)) // short-circuiting
#define SPLM_BTA_CLEAR(map, elt)       (map)[(elt)>>3]&=~(1<<((elt)&7))
#define SPLM_BTA_CLEARALL(map, elt)    (map)[(elt)>>3]=0
#define SPLM_BTA_TOGGLE(map, elt)      (map)[(elt)>>3]^=(1<<((elt)&7))

/* as above but forcing evaluation of elt to the supplied auxiliary variable */
#define SPLM_BTA_SET_FE(map, aux, elt)         ((aux)=(elt), SPLM_BTA_SET((map), (aux)))
#define SPLM_BTA_ISSET_FE(map, aux, elt)       ((aux)=(elt), SPLM_BTA_ISSET((map), (aux)))
#define SPLM_BTA_ISSET_SPRS_FE(map, aux, elt)  ((aux)=(elt), SPLM_BTA_ISSET_SPRS((map), (aux)))
#define SPLM_BTA_CLEAR_FE(map, aux, elt)       ((aux)=(elt), SPLM_BTA_CLEAR((map), (aux)))
#define SPLM_BTA_TOGGLE_FE(map, aux, elt)      ((aux)=(elt), SPLM_BTA_TOGGLE((map), (aux)))


/************************** At*A and A^t*x for a CCS matrix A **************************/

/* compute the CCS structure of A^T*A for a matrix A in CCS format. if AtA->nnz>0, assumes that AtA
 * is large enough to hold the product. Otherwise, it computes its number of nonzeros and allocates
 * it accordingly.
 *
 * If job is 1, computes just the structure of the lower triangular & diagonal part of the product
 * If job is > 1, computes the structure of the full matrix product (i.e., all L, D & U parts)
 */
void splm_setup_AtA_ccsA(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job)
{
register int i, j, k;
register int *AColptr, *ARowidx, *AtARowidx;
int nc, *AtAColptr, *ARowidxM1, *AColptrP1;
int ilow, ihigh, jlow, jhigh;
int allocating=0;

#ifdef USE_BTA_MAP
unsigned char *rowmap; // bit array marking the nonzero rows of A's column j
register int *ridx;
int *lastp;
#endif

  nc=A->nc;
  AColptr=A->colptr; AColptrP1=AColptr+1; // for all i, AColptrP1[i]==AColptr[i+1] 
  ARowidx=A->rowidx; ARowidxM1=ARowidx-1; // for all i, ARowidxM1[i]==ARowidx[i-1]
  AtAColptr=AtA->colptr; AtARowidx=AtA->rowidx;

#ifdef USE_BTA_MAP
  i=(A->nr>>3)+1; // A->nr/8+1
  rowmap=(unsigned char *)emalloc(i*sizeof(unsigned char));
  memset(rowmap, 0, i*sizeof(unsigned char));
#endif

  /* compute the nonzero pattern of A^T*A. The code below is very similar to that
   * used for the actual computation of the matrix product in splm_calc_AtA_ccsA()
   */
  if(AtA->nnz<=0){ // AtA not preallocated, start by assuming it has as many nonzeros as A
    allocating=1;
    splm_ccsm_alloc_novalues(AtA, nc, nc, A->nnz);
    /* re-init pointers */
    AtAColptr=AtA->colptr; AtARowidx=AtA->rowidx;
  }

  for(j=k=0; j<nc; ++j){
    AtAColptr[j]=k;

    jlow=ARowidx[AColptr[j]];
    jhigh=ARowidxM1[AColptrP1[j]];

#ifdef USE_BTA_MAP
    for(ridx=ARowidx+AColptr[j], lastp=ARowidx+AColptrP1[j]; ridx!=lastp; ++ridx)
      SPLM_BTA_SET_FE(rowmap, i, *ridx);
#endif

    for(i=(job==1)? j : 0; i<nc; ++i){
      ilow=ARowidx[AColptr[i]];
      if(jhigh<ilow) continue; /* skip columns with non overlapping row index ranges */

      ihigh=ARowidxM1[AColptrP1[i]];
      if(jlow>ihigh) continue; /* skip columns with non overlapping row index ranges */

#ifdef USE_BTA_MAP
      for(ridx=ARowidx+AColptr[i], lastp=ARowidx+AColptrP1[i]; ridx!=lastp; ++ridx){
        register int aux;

        if(SPLM_BTA_ISSET_FE(rowmap, aux, *ridx)){
          /* assume no zero cancellations */
          if(allocating && k>=AtA->nnz){ // more memory needed, double current size
            splm_ccsm_realloc_novalues(AtA, nc, nc, AtA->nnz<<1); // 2*AtA->nnz
            AtARowidx=AtA->rowidx; // re-init
          }
          AtARowidx[k++]=i;
          break;
        }
      }
#else
      /* following code avoids the use of rowmap; it is based on searching and is hence slower */
      {
      register int ii, jj, diff;
      int idx;

      for(ii=AColptr[i], jj=AColptr[j]; 1 /*ii<AColptrP1[i] && jj<AColptrP1[j]*/;  ){ // ii & jj are modified within the loop
        diff=ARowidx[ii]-ARowidx[jj];
        if(diff<0){
          idx=binsearch_geq(ARowidx+ii, AColptrP1[i]-ii, ARowidx[jj]);
          if(idx<0) break;
          ii+=idx;
          continue;
        }
        else if(diff>0){
          idx=binsearch_geq(ARowidx+jj, AColptrP1[j]-jj, ARowidx[ii]);
          if(idx<0) break;
          jj+=idx;
          continue;
        }
        else{ // ARowidx[ii]==ARowidx[jj]
          /* assume no zero cancellations */
          if(allocating && k>=AtA->nnz){ // more memory needed, double current size
            splm_ccsm_realloc_novalues(AtA, nc, nc, AtA->nnz<<1); // 2*AtA->nnz
            AtARowidx=AtA->rowidx; // re-init
          }
          AtARowidx[k++]=i;
          break;
        }
      }
      }
#endif /* USE_BTA_MAP */
    }
#ifdef USE_BTA_MAP
    for(ridx=ARowidx+AColptr[j], lastp=ARowidx+AColptrP1[j]; ridx!=lastp; ++ridx) // restore to all zeros
      SPLM_BTA_CLEARALL(rowmap, *ridx);
#endif
  }
  AtAColptr[nc]=k;

#ifdef USE_BTA_MAP
  free(rowmap);
#endif

  if(allocating){
    splm_ccsm_realloc_novalues(AtA, nc, nc, k); // adjust to actual size...
    splm_ccsm_alloc_values(AtA); // ...and allocate values
    /* re-init pointers */
    AtARowidx=AtA->rowidx; // re-init
  }
  else{
    if(k!=AtA->nnz){
      int enough=AtA->nnz>k;

      fprintf(stderr, "%sreallocated size of array A^T*A does not match the number "
              "of detected nonzeros in splm_setup_AtA_ccsA()! [%d != %d]\n", enough? "Warning: p" : "P", AtA->nnz, k);
      if(!enough) exit(1);
    }
  }
}

/* compute A^T*A in CCS for a matrix A in CCS format. It is assumed that the structure
 * of AtA has been determined with a previous call to splm_setup_AtA_ccsA().
 *
 * If job is 1, computes just the lower triangular & diagonal part of the product
 * If job is > 1, computes the full matrix product (i.e., all L, D, & U parts)
 */
void splm_calc_AtA_ccsA(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job)
{
register int i, j, ii, jj, diff;
register int *AColptr, *AColptrP1, *ARowidx;
int *AtAColptr, *AtAColptrP1, *AtARowidx;
double *AVal, *AtAVal;
int nc, idx, k;
register double sum;

  if(AtA->nnz<=0){
    fprintf(stderr, "product matrix must have been preallocated in splm_calc_AtA_ccsA()!\n");
    exit(1);
  }
  AtAVal=AtA->val;
  AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1] 
  AtARowidx=AtA->rowidx;

  nc=A->nc;
  AColptr=A->colptr; AColptrP1=AColptr+1; // for all i, AColptrP1[i]==AColptr[i+1] 
  ARowidx=A->rowidx; //ARowidxM1=ARowidx-1; // for all i, ARowidxM1[i]==ARowidx[i-1]
  AVal=A->val;

  /* compute (column by column) AtA[i][j] as \sum_k A[k][i]*A[k][j] */
  /* Symmetry in A^T*A is exploited in order to reduce computations in half */
  for(j=0; j<nc; ++j){
    //jlow=ARowidx[AColptr[j]];
    //jhigh=ARowidxM1[AColptrP1[j]];

    if(job==1){
      /* assume that splm_setup_AtA_ccsA() has also been called with job==1,
       * so that the first element on column j is the one on the diagonal
       */
      //assert(AtARowidx[AtAColptr[j]]==j);
      k=AtAColptr[j];
    }
    else{
      /* ignore i<j, i.e. compute elements below or on the diagonal only */
      k=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
      k+=AtAColptr[j];
    }
    for(  ; k<AtAColptrP1[j]; ++k){ // loop over nonzero [A^T*A]_ij, i=0...nc-1 in the lower triangular part (i.e., i>=j)
      i=AtARowidx[k];

      //ilow=ARowidx[AColptr[i]];
      //ihigh=ARowidxM1[AColptrP1[i]];

      for(ii=AColptr[i], jj=AColptr[j], sum=0.0; ii<AColptrP1[i] && jj<AColptrP1[j];  ){ // ii & jj are modified within the loop
        diff=ARowidx[ii]-ARowidx[jj];
        if(diff==0){
          /* the possibility of an inner product summing up to zero is ignored */
          sum+=AVal[ii]*AVal[jj];
          ++ii; ++jj;
          continue;
        } else if(diff<0){
          idx=binsearch_geq(ARowidx+ii, AColptrP1[i]-ii, ARowidx[jj]);
          if(idx<0) break;
          ii+=idx;
          continue;
        }
        else{ // diff>0
          idx=binsearch_geq(ARowidx+jj, AColptrP1[j]-jj, ARowidx[ii]);
          if(idx<0) break;
          jj+=idx;
          continue;
        }
      }

      AtAVal[k]=sum;
    }
  }

  if(job==1) return; /* no upper part */

  /* copy the elements below the diagonal to those above it */
  for(j=0; j<nc; ++j){
    int ilow, ihigh;

    /* find nonzero [A^T*A]_ij, i=0...nc-1 above the diagonal (i.e., i<j) */
    /* note that element _jj on the diagonal is guaranteed to be nonzero */
    jj=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
    jj+=AtAColptr[j];
    for(k=AtAColptr[j]; k<jj; ++k){
      i=AtARowidx[k];

      /* for all elements (i, j) copy the corresponding elements (j, i) to them */

      /* a binary search over all rows of column i is performed next, some speed
       * might be gained by searching below the diagonal only by setting ilow to
       *
       * ilow=binsearch_geq(AtARowidx+AtAColptr[i], AtAColptrP1[i]-AtAColptr[i], i) + AtAColptr[i];
       */
      ilow=AtAColptr[i];
      ihigh=AtAColptrP1[i];
      while(ilow<ihigh){
        /* the following line is a sanity check...*/
        //if(j<AtARowidx[ilow] || j>AtARowidx[ihigh-1]) break; // not found

        ii=(ilow + ihigh) >> 1; //(ilow+ihigh)/2;
        //ii=ilow+((ihigh-ilow)>>1); /* ensures no index overflows */
        diff=j-AtARowidx[ii];
        if(diff<0)
          ihigh=ii;
        else if(diff>0)
          ilow=ii+1;
        else{ // found
          AtAVal[k]=AtAVal[ii];
          break;
        }
      }
    }
  }
}


/* CCS SpMV (sparse matrix-vector multiply) */

#if 0
#include <oski/oski.h>

/* compute y=A^T*x for a matrix A in CCS format and a vector x using OSKI */
void splm_calc_Atx_ccsA(struct splm_ccsm *A, double *x, double *y)
{
oski_matrix_t A_tunable;
oski_vecview_t x_view, y_view;
const double alpha=1.0, beta=0.0;

  oski_Init();
  A_tunable=oski_CreateMatCSC(A->colptr, A->rowidx, A->val, A->nr, A->nc, SHARE_INPUTMAT, 1, INDEX_ZERO_BASED);
  x_view=oski_CreateVecView(x, A->nr, STRIDE_UNIT);
  y_view=oski_CreateVecView(y, A->nc, STRIDE_UNIT);

  /*  y = a*A^T*x + b*y */
  oski_MatMult(A_tunable, OP_TRANS, alpha, x_view, beta, y_view);

  oski_DestroyMat(A_tunable);
  oski_DestroyVecView(x_view);
  oski_DestroyVecView(y_view);
  oski_Close();
}

#else
/* compute y=A^T*x for a matrix A in CCS format and a vector x */
void splm_calc_Atx_ccsA(struct splm_ccsm *A, double *const x, double *const y)
{
#ifndef USE_MKL_ATX
register int j, *ARowidx;
register double *AVal, sum0, sum1, sum2, sum3;
int i, nc, *AColptr, *AColptrP1, jlow;
const int blocksize=8, bpwr=3; /* 8=2^3 */
int n, blockn;

  nc=A->nc;
  AColptr=A->colptr; AColptrP1=AColptr+1; // for all i, AColptrP1[i]==AColptr[i+1] 

  /* compute y[i] as \sum_j A[j][i]*x[j] */
  for(i=nc; i-->0; ){ // A^T is nc x nr
    jlow=AColptr[i]; //jhigh=AColptrP1[i];
    AVal=A->val+jlow; ARowidx=A->rowidx+jlow;
    n=AColptrP1[i]-jlow;
#if 0
    /* straightforward implementation */
    for(j=n, sum0=0.0; j-->0;  )
      sum0+=AVal[j]*x[ARowidx[j]];
    y[i]=sum0;
#else
    /* use loop unrolling and blocking */
    blockn=(n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */
    /* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
    sum0=sum1=sum2=sum3=0.0;
    for(j=blockn-1; j>0; AVal+=blocksize, ARowidx+=blocksize, j-=blocksize){
      sum0+=AVal[0]*x[ARowidx[0]];
      sum1+=AVal[1]*x[ARowidx[1]];
      sum2+=AVal[2]*x[ARowidx[2]];
      sum3+=AVal[3]*x[ARowidx[3]];
      sum0+=AVal[4]*x[ARowidx[4]];
      sum1+=AVal[5]*x[ARowidx[5]];
      sum2+=AVal[6]*x[ARowidx[6]];
      sum3+=AVal[7]*x[ARowidx[7]];
    }
    /* see if there's any more work left */
    j=blockn;
    if(j<n){
      switch(n-j){
        case 7: sum0+=AVal[6]*x[ARowidx[6]];
        case 6: sum1+=AVal[5]*x[ARowidx[5]];
        case 5: sum2+=AVal[4]*x[ARowidx[4]];
        case 4: sum3+=AVal[3]*x[ARowidx[3]];
        case 3: sum0+=AVal[2]*x[ARowidx[2]];
        case 2: sum1+=AVal[1]*x[ARowidx[1]];
        case 1: sum2+=AVal[0]*x[ARowidx[0]];
      }
    }
    y[i]=(sum0+sum1)+(sum2+sum3);
#endif /* 0 */
  }
#else /* USE_MKL_ATX */
int nr, nc, nnz;
register int i;
double *AVal;
int *AColptr, *ARowidx;

const double alpha=1.0, beta=0.0;
extern void mkl_dcscmv(char *transa, int *m, int *k, double *alpha, char *matdescra, double *val, int *indx,
                       int *pntrb, int *pntre, double *x, double *beta, double *y);

  nr=A->nr; nc=A->nc; nnz=A->nnz;
  AColptr=A->colptr; ARowidx=A->rowidx; AVal=A->val;

  /* according to the documentation of MKL's mkl_dcscmv, it supports both
   * zero and one-based indexing. However specifying "GxxC" for the matdescra
   * argument did not work...
   */

  /* convert to one-based indexing */
  for(i=0; i<nnz; ++i) ARowidx[i]+=1;
  for(i=0; i<nc+1; ++i) AColptr[i]+=1;

  mkl_dcscmv("T", &nr, &nc, (double *)&alpha, "G", AVal, ARowidx, AColptr, AColptr+1, x, (double *)&beta, y);

  /* restore to zero-based */
  for(i=0; i<nnz; ++i) ARowidx[i]-=1;
  for(i=0; i<nc+1; ++i) AColptr[i]-=1;
#endif /* USE_MKL_ATX */
}
#endif

#ifdef USE_ATA_QUADS

/* NOTE: index compression functions below implement delta encoding, i.e.
 * compress the differences between successive indices rather than the
 * indices themselves
 */

#ifdef USE_VBYTE_COMPRESSION
/* 
 * byte-aligned index compression (vByte technique)
 *
 * Least significant bits come first, leftmost bit in each byte
 * signals if further blocks follow (1 if yes, 0 otherwise).
 * To improve compression, differences (i.e. deltas) of consecutive
 * integers are used. Two compression variants are supplied, one
 * assuming that the integer sequence to be compressed is 
 * non-decreasing (thus deltas are non-negative), and a more general
 * one that allows arbitrary integer sequences (i.e., deltas are of
 * arbitrary signs. The second variant needs to store an extra bit
 * for the sign in the second to leftmost bit of the first byte.  
 *
 * see Scholer et al. "Compression of Inverted Indexes for Fast Query Evaluation".
 * Proc. ACM SIGIR Conf. on Research and Development in Information Retrieval, 2002.
 * http://doi.acm.org/10.1145/564376.564416 
 * and
 * Williams and Zobel "Compressing Integers for Fast File Access".
 * Computer Journal, Vol. 42, Num. 3, 1999, pp. 193-201(9)
 */

/* non-decreasing sequences, i.e. uncompressed >= prev */
static int compressVB_ndecr(int uncompressed, int *prev, unsigned char *compressed)
{
register int delta, i=0;

  delta=uncompressed-*prev;
  if(delta<0){
    fprintf(stderr, "compressVB_ndecr(): cannot compress negative delta %d\n", delta);
    exit(1);
  }
  *prev=uncompressed;
  while(delta>=128){
    compressed[i]=(delta&127)|128; ++i;
    delta>>=7;
  }
  compressed[i]=delta; ++i;

  return i;
}

static int uncompressVB_ndecr(const unsigned char *compressed, int *prev, int *uncompressed)
{
register int p, shift, tmp, i=0;

  for(shift=0, p=*prev; 1; shift+=7){
    tmp=compressed[i]; ++i;
    p+=((tmp&127)<<shift);
    if(tmp<128) break;
  }
  *uncompressed=*prev=p;

  return i;
}

/* arbitrary relationship between uncompressed & prev */
static int compressVB(int uncompressed, int *prev, unsigned char *compressed)
{
register int delta, tmp, i=1;

  delta=uncompressed-*prev;
  *prev=uncompressed;
  if(delta>=0){
    tmp=(delta&63);
  }
  else{
    delta=-delta;
    tmp=(delta&63)|64; // set sign bit
  }
  *compressed=tmp;
  delta>>=6;
  if(delta>0){ // original delta > 64, continue to more bytes
    *compressed|=128; // mark first byte as non-final
   
    while(delta>=128){
      compressed[i]=(delta&127)|128; ++i;
      delta>>=7;
    }
    compressed[i]=delta; ++i;
  }

  return i;
}

static int uncompressVB(const unsigned char *compressed, int *prev, int *uncompressed)
{
register int shift, tmp, delta, i=1;
int sg;

  tmp=*compressed;
  delta=tmp&63;
  sg=tmp&64; // get sign bit
  for(shift=6; tmp>=128; shift+=7){
    tmp=compressed[i]; ++i;
    delta+=((tmp&127)<<shift);
  }
#if 1
  if(sg==0)
    *prev+=delta;
  else
    *prev-=delta;
#else
  sg=+1|-(sg>>6);
  *prev+=sg*delta;
#endif /* 0 */
  *uncompressed=*prev;

  return i;
}

/* special type of non-decreasing sequences: uncompressed == prev  || uncompressed == prev + 1 */
/* deltas are either 0 or 1, so each can be represented with one bit, i.e. 8 deltas per byte */
/* previous element passed in prev[0], length n in prev[1] */
static int compressVB_01(int uncompressed, int *prev, unsigned char *compressed)
{
register int delta, i;

  delta=uncompressed-*prev;
  *prev=uncompressed;
  i=(prev[1])&7; // n&7 == n%8
  ++prev[1];
  if(delta==1){
    *compressed|=1<<i;
  }
  else if(delta==0){
    *compressed&=~(1<<i);
  }else{
    fprintf(stderr, "invalid delta %d in compressVB_01()!\n", delta);
    exit(1);
  }
  if(i<7) return 0;

  return 1;
}

/* previous element passed in prev[0], length n in prev[1] */
inline static int uncompressVB_01(const unsigned char *compressed, int *prev, int *uncompressed)
{
register int i=0;

  i=(prev[1])&7; // n&7 == n%8
  ++prev[1];
  //*uncompressed=*prev=*prev + ((*compressed&(1<<i))>>i);
  *uncompressed=*prev=*prev + ((*compressed>>i)&1);
  if(i<7) return 0;

  return 1;
}
#endif /* USE_VBYTE_COMPRESSION */


#ifdef USE_S9_COMPRESSION

/* max unsigned int that can be represented with 28 bits */
#define __MAXUINT28   268435455 /* (2^28) - 1 */

/* Simple-9 integer compression
 *
 * Simple-9 packs several small integers into 28 bits of a 32 bit word,
 * leaving the remaining 4 bits as a selector. See 
 * p.93 of D. Salomon, G. Motta and D. Bryant: "Handbook of Data Compression"
 * or
 * V.N. Anh and A. Moffat: Inverted Index Compression Using Word-Aligned Binary Codes,
 * Information Retrieval(8):1, pp. 151-166, 2005
 * http://dx.doi.org/10.1023/B:INRT.0000048490.99518.5c
 *
 */

#define _MAXCOMP_S9NDEC 28 // max number of integers that can be compressed at once; equals the last element of codeno below
#define _UBUFSZ_S9NDEC (_MAXCOMP_S9NDEC+3) // 3: 1 for the unused empty slot + 2 for the start & end pointers
#define _CBUFSZ_S9NDEC (_UBUFSZ_S9NDEC-2)

#define _MAXCOMPS9  _MAXCOMP_S9NDEC //max(_MAXCOMP_S9NDEC, _MAXCOMP_S9ARB)
#define _UBUFSZS9   _UBUFSZ_S9NDEC //max(_UBUFSZ_S9NDEC, _UBUFSZ_S9ARB)

/* this is the original S-9 for compressing a nondecreasing, positive index sequence.
 * Since the number of indices that will be stored in 32 bits is
 * not known beforehand, a circular buffer is used to temporarily
 * store input indices for later examination. This circular buffer
 * is implemented in array uncompressed[], which for simplicity
 * stores the buffers' start and end pointers as its two last
 * elements.
 */
static int compressS9_ndecr(int idx, int *prev, int32_t uncompressed[_UBUFSZ_S9NDEC], int flush, uint32_t *compressed)
{
static int codelen[]={28, 14, 9, 7, 5, 4, 3, 2 , 1 };
static int codeno[] ={1 , 2 , 3, 4, 5, 7, 9, 14, 28};
register int j, sel;
int csz, cn, bl, delta;
uint32_t res, delta32, lim;
const int codenosz=sizeof(codeno)/sizeof(int);
/* circular buffer management */
int start, end, nidx;
register int cursor;

  delta=idx-*prev; *prev=idx;
  start=uncompressed[_UBUFSZ_S9NDEC-2];
  end=uncompressed[_UBUFSZ_S9NDEC-1];
  nidx=end-start; if(nidx<0) nidx+=_CBUFSZ_S9NDEC; // circular buffer size

  if(!flush){
    /* add to buffer */
    j=end++; if(end==_CBUFSZ_S9NDEC) end=0; // modulo increment
    if(end!=start){ // non full
      if((uint32_t)(delta) > (uint32_t)(__MAXUINT28)){ // delta<0 || delta>__MAXUINT28
        if(delta>0)
          fprintf(stderr, "compressS9_ndecr(): delta %d cannot be compressed with 28 bits; use vByte compression\n", delta);
        else
          fprintf(stderr, "compressS9_ndecr(): cannot compress negative delta %d\n", delta);
        exit(1);
      }
      uncompressed[j]=delta;
      uncompressed[_UBUFSZ_S9NDEC-1]=end; // update end
      ++nidx;
    }
    else{
      fprintf(stderr, "compressS9_ndecr(): Full circular buffer!\n");
      exit(1);
    }

    /* the maximum number of integers that can be compressed in a single word is
     * given by the last element of codeno[] (i.e., _MAXCOMP_S9NDEC). A premature request
     * to compress less than _MAXCOMP_S9NDEC elements might result in a compression worse than
     * that possible with all _MAXCOMP_S9NDEC elements. Thus, unless asked to "flush" the
     * uncompressed array, this function accumulates integers in it until they
     * become at least _MAXCOMP_S9NDEC, in which case it proceeds to compressing them
     */
    if(nidx<codeno[codenosz-1]) return 0;
  }

  csz=0;
  while(nidx>0){
    for(sel=codenosz; sel-->0; ){
      cn=codeno[sel];
      if(cn>nidx) continue;

      res=0;
      bl=codelen[sel];
      lim=1<<bl;
  
      /* iterate over elements in circular buffer */
      for(j=cn, cursor=start; j>0; --j){
#if 0
        if(cursor==end){ // empty buffer
          /* should not happen */
          fprintf(stderr, "compressS9_ndecr(): internal error, empty circular buffer!\n");
          exit(1);
        }
#endif
        delta32=uncompressed[cursor];
        ++cursor; if(cursor==_CBUFSZ_S9NDEC) cursor=0; // modulo increment

        if(delta32<lim)
          res=(res<<bl) | delta32;
        else break;
      }
      if(j==0){
        //printf("%d ints in 32 bits, code %d\n", cn, sel);
        res|=sel<<28;
        compressed[csz++]=res;
        nidx-=cn;
        /* remove cn elements from circular buffer
         * cursor==(start+cn)%_CBUFSZ_S9NDEC
         */
        uncompressed[_UBUFSZ_S9NDEC-2]=start=cursor; // update start
        break;
      }
    }
  }

  return csz;
}

/* uncompress a 32 bit int encoding several indices according to the S-9 scheme.
 * Produces at most _MAXCOMP_S9NDEC elements
 */
static int uncompressS9_ndecr(const uint32_t compressed, int *prev, int uncompressed[_MAXCOMP_S9NDEC])
{
register int i, tmp;
int j;
uint32_t sel;
register uint32_t low;

  sel=compressed>>28;
  low=compressed&0xfffffff;
  switch(sel){
    case 8:
      //28 codes of 1 bit
      /*
      for(i=j=0; i<28; ++i)
        uncompressed[j++]=1 & (low >> (27-i));
      */
      uncompressed[0 ]=1 & (low >> 27);
      uncompressed[1 ]=1 & (low >> 26);
      uncompressed[2 ]=1 & (low >> 25);
      uncompressed[3 ]=1 & (low >> 24);
      uncompressed[4 ]=1 & (low >> 23);
      uncompressed[5 ]=1 & (low >> 22);
      uncompressed[6 ]=1 & (low >> 21);
      uncompressed[7 ]=1 & (low >> 20);
      uncompressed[8 ]=1 & (low >> 19);
      uncompressed[9 ]=1 & (low >> 18);
      uncompressed[10]=1 & (low >> 17);
      uncompressed[11]=1 & (low >> 16);
      uncompressed[12]=1 & (low >> 15);
      uncompressed[13]=1 & (low >> 14);
      uncompressed[14]=1 & (low >> 13);
      uncompressed[15]=1 & (low >> 12);
      uncompressed[16]=1 & (low >> 11);
      uncompressed[17]=1 & (low >> 10);
      uncompressed[18]=1 & (low >> 9);
      uncompressed[19]=1 & (low >> 8);
      uncompressed[20]=1 & (low >> 7);
      uncompressed[21]=1 & (low >> 6);
      uncompressed[22]=1 & (low >> 5);
      uncompressed[23]=1 & (low >> 4);
      uncompressed[24]=1 & (low >> 3);
      uncompressed[25]=1 & (low >> 2);
      uncompressed[26]=1 & (low >> 1);
      uncompressed[27]=1 &  low;
      j=28;
      break;
    case 7:
      //14 codes of 2 bits
      uncompressed[0 ]=((low >> 26) & 0x3);
      uncompressed[1 ]=((low >> 24) & 0x3);
      uncompressed[2 ]=((low >> 22) & 0x3);
      uncompressed[3 ]=((low >> 20) & 0x3);
      uncompressed[4 ]=((low >> 18) & 0x3);
      uncompressed[5 ]=((low >> 16) & 0x3);
      uncompressed[6 ]=((low >> 14) & 0x3);
      uncompressed[7 ]=((low >> 12) & 0x3);
      uncompressed[8 ]=((low >> 10) & 0x3);
      uncompressed[9 ]=((low >> 8)  & 0x3);
      uncompressed[10]=((low >> 6)  & 0x3);
      uncompressed[11]=((low >> 4)  & 0x3);
      uncompressed[12]=((low >> 2)  & 0x3);
      uncompressed[13]=(low         & 0x3);
      j=14;
      break;
    case 6:
      //9 codes of 3 bits
      uncompressed[0]=((low >> 24) & 0x7);
      uncompressed[1]=((low >> 21) & 0x7);
      uncompressed[2]=((low >> 18) & 0x7);
      uncompressed[3]=((low >> 15) & 0x7);
      uncompressed[4]=((low >> 12) & 0x7);
      uncompressed[5]=((low >> 9)  & 0x7);
      uncompressed[6]=((low >> 6)  & 0x7);
      uncompressed[7]=((low >> 3)  & 0x7);
      uncompressed[8]=(low         & 0x7);
      j=9;
      break;
    case 5:
      //7 codes of 4 bits
      uncompressed[0]=((low >> 24) & 0xf);
      uncompressed[1]=((low >> 20) & 0xf);
      uncompressed[2]=((low >> 16) & 0xf);
      uncompressed[3]=((low >> 12) & 0xf);
      uncompressed[4]=((low >> 8)  & 0xf);
      uncompressed[5]=((low >> 4)  & 0xf);
      uncompressed[6]=(low         & 0xf);
      j=7;
      break;
    case 4:
      //5 codes of 5 bits
      uncompressed[0]=((low >> 20) & 0x1f);
      uncompressed[1]=((low >> 15) & 0x1f);
      uncompressed[2]=((low >> 10) & 0x1f);
      uncompressed[3]=((low >> 5)  & 0x1f);
      uncompressed[4]=(low         & 0x1f);
      j=5;
      break;
    case 3:
      //4 codes of 7 bits
      uncompressed[0]=((low >> 21) & 0x7f);
      uncompressed[1]=((low >> 14) & 0x7f);
      uncompressed[2]=((low >> 7)  & 0x7f);
      uncompressed[3]=(low         & 0x7f);
      j=4;
      break;
    case 2:
      //3 codes of 9 bits
      uncompressed[0]=((low >> 18) & 0x1ff);
      uncompressed[1]=((low >> 9)  & 0x1ff);
      uncompressed[2]=(low         & 0x1ff);
      j=3;
      break;
    case 1:
      //2 codes of 14 bits
      uncompressed[0]=((low >> 14) & 0x3fff);
      uncompressed[1]=(low         & 0x3fff);
      j=2;
      break;
    case 0:
      //1 code of 28 bits
      uncompressed[0]=low;
      j=1;
      break;
    default: /* should not reach this point */
      fprintf(stderr, "uncompressS9_ndecr(): internal error (unexpected selector %d)\n", sel);
      exit(1);
  }

  tmp=*prev;
  for(i=0; i<j; ++i){
    tmp+=uncompressed[i];
    uncompressed[i]=tmp;
  }
  *prev=tmp;

  return j;
}

#define _MAXCOMP_S9ARB 14 // max number of integers that can be compressed at once; equals the last element of codeno below
#define _UBUFSZ_S9ARB (_MAXCOMP_S9ARB+3) // 3: 1 for the unused empty slot + 2 for the start & end pointers
#define _CBUFSZ_S9ARB (_UBUFSZ_S9ARB-2)

/* variation of S-9 that supports arbitrary index sequences.
 * It is based on the idea of adding a sufficiently large number to
 * each index so that negative indices are mapped to positive values.
 * Since the number of indices that will be stored in 32 bits is
 * not known beforehand, a circular buffer is used to temporarily
 * store input indices for later examination. This circular buffer
 * is implemented in array uncompressed[], which for simplicity
 * stores the buffers' start and end pointers as its two last
 * elements.
 */
static int compressS9(int idx, int *prev, int32_t uncompressed[_UBUFSZ_S9ARB], int flush, uint32_t *compressed)
{
static int codelen[] ={28, 14  , 9  , 7 , 5 , 4, 3, 2 };
static int codeno[]  ={1 , 2   , 3  , 4 , 5 , 7, 9, 14};
static int codeoffs[]={0 , 8191, 255, 63, 15, 7, 3, 1 }; /* 2^{codelen-1}-1, excluding codelen=28 */
register int j, sel;
int csz, cn, bl, delta;
uint32_t res, tmpu, lim;
int32_t  delta32;
const int codenosz=sizeof(codeno)/sizeof(int);
/* circular buffer management */
int start, end, nidx;
register int cursor;

  delta=idx-*prev; *prev=idx;
  start=uncompressed[_UBUFSZ_S9ARB-2];
  end=uncompressed[_UBUFSZ_S9ARB-1];
  nidx=end-start; if(nidx<0) nidx+=_CBUFSZ_S9ARB; // circular buffer size

  if(!flush){
    /* add to buffer */
    j=end++; if(end==_CBUFSZ_S9ARB) end=0; // modulo increment
    if(end!=start){ // non full
      if((uint32_t)(delta+__MAXUINT28) > (uint32_t)(__MAXUINT28<<1)){ // delta < -__MAXUINT28 || delta > __MAXUINT28
        fprintf(stderr, "compressS9(): delta %d cannot be compressed with 28 bits; use vByte compression\n", delta);
        exit(1);
      }
      uncompressed[j]=delta;
      uncompressed[_UBUFSZ_S9ARB-1]=end; // update end
      ++nidx;
    }
    else{
      fprintf(stderr, "compressS9(): Full circular buffer!\n");
      exit(1);
    }

    /* the maximum number of integers that can be compressed in a single word is
     * given by the last element of codeno[] (i.e., _MAXCOMP_S9ARB). A premature request
     * to compress less than _MAXCOMP_S9ARB elements might result in a compression worse than
     * that possible with all _MAXCOMP_S9ARB elements. Thus, unless asked to "flush" the
     * uncompressed array, this function accumulates integers in it until they
     * become at least _MAXCOMP_S9ARB, in which case it proceeds to compressing them
     */
    if(nidx<codeno[codenosz-1]) return 0;
  }

  csz=0;
  while(nidx>0){
    for(sel=codenosz; sel-->0; ){
      cn=codeno[sel];
      if(cn>nidx) continue;

      res=0;
      bl=codelen[sel];
      lim=1<<bl;
  
      /* iterate over elements in circular buffer */
      for(j=cn, cursor=start; j>0; --j){
#if 0
        if(cursor==end){ // empty buffer
          /* should not happen */
          fprintf(stderr, "compressS9(): internal error, empty circular buffer!\n");
          exit(1);
        }
#endif
        delta32=uncompressed[cursor];
        ++cursor; if(cursor==_CBUFSZ_S9ARB) cursor=0; // modulo increment

        if(sel>0){
          tmpu=delta32+codeoffs[sel];
          if(tmpu>=lim) break; // (int)tmpu<0 || tmpu>=lim
        }
        else{
          if(delta32>=0)
            tmpu=delta32;
          else{
            tmpu=-delta32;
            sel=8;
          }
        }

        res=(res<<bl) | tmpu;
      }
      if(j==0){
        //printf("%d ints in 32 bits, code %d\n", cn, sel);
        res|=sel<<28;
        compressed[csz++]=res;
        nidx-=cn;
        /* remove cn elements from circular buffer
         * cursor==(start+cn)%_CBUFSZ_S9ARB
         */
        uncompressed[_UBUFSZ_S9ARB-2]=start=cursor; // update start
        break;
      }
    }
  }

  return csz;
}

/* uncompress a 32 bit int encoding several indices according to a modified S-9 scheme.
 * Produces at most _MAXCOMP_S9ARB elements
 */
static int uncompressS9(const uint32_t compressed, int *prev, int uncompressed[_MAXCOMP_S9ARB])
{
register int i, tmp;
int j;
uint32_t sel;
register uint32_t low;

  sel=compressed>>28;
  low=compressed&0xfffffff;
  switch(sel){
    case 8:
      //1 code of 28 bits
      *uncompressed=-((int)low);
      j=1;
      break;
    case 7:
      //14 codes of 2 bits
      uncompressed[0 ]=((low >> 26) & 0x3) - 1;
      uncompressed[1 ]=((low >> 24) & 0x3) - 1;
      uncompressed[2 ]=((low >> 22) & 0x3) - 1;
      uncompressed[3 ]=((low >> 20) & 0x3) - 1;
      uncompressed[4 ]=((low >> 18) & 0x3) - 1;
      uncompressed[5 ]=((low >> 16) & 0x3) - 1;
      uncompressed[6 ]=((low >> 14) & 0x3) - 1;
      uncompressed[7 ]=((low >> 12) & 0x3) - 1;
      uncompressed[8 ]=((low >> 10) & 0x3) - 1;
      uncompressed[9 ]=((low >> 8)  & 0x3) - 1;
      uncompressed[10]=((low >> 6)  & 0x3) - 1;
      uncompressed[11]=((low >> 4)  & 0x3) - 1;
      uncompressed[12]=((low >> 2)  & 0x3) - 1;
      uncompressed[13]=(low         & 0x3) - 1;
      j=14;
      break;
    case 6:
      //9 codes of 3 bits
      uncompressed[0]=((low >> 24) & 0x7) - 3;
      uncompressed[1]=((low >> 21) & 0x7) - 3;
      uncompressed[2]=((low >> 18) & 0x7) - 3;
      uncompressed[3]=((low >> 15) & 0x7) - 3;
      uncompressed[4]=((low >> 12) & 0x7) - 3;
      uncompressed[5]=((low >> 9)  & 0x7) - 3;
      uncompressed[6]=((low >> 6)  & 0x7) - 3;
      uncompressed[7]=((low >> 3)  & 0x7) - 3;
      uncompressed[8]=(low         & 0x7) - 3;
      j=9;
      break;
    case 5:
      //7 codes of 4 bits
      uncompressed[0]=((low >> 24) & 0xf) - 7;
      uncompressed[1]=((low >> 20) & 0xf) - 7;
      uncompressed[2]=((low >> 16) & 0xf) - 7;
      uncompressed[3]=((low >> 12) & 0xf) - 7;
      uncompressed[4]=((low >> 8)  & 0xf) - 7;
      uncompressed[5]=((low >> 4)  & 0xf) - 7;
      uncompressed[6]=(low         & 0xf) - 7;
      j=7;
      break;
    case 4:
      //5 codes of 5 bits
      uncompressed[0]=((low >> 20) & 0x1f) - 15;
      uncompressed[1]=((low >> 15) & 0x1f) - 15;
      uncompressed[2]=((low >> 10) & 0x1f) - 15;
      uncompressed[3]=((low >> 5)  & 0x1f) - 15;
      uncompressed[4]=(low         & 0x1f) - 15;
      j=5;
      break;
    case 3:
      //4 codes of 7 bits
      uncompressed[0]=((low >> 21) & 0x7f) - 63;
      uncompressed[1]=((low >> 14) & 0x7f) - 63;
      uncompressed[2]=((low >> 7)  & 0x7f) - 63;
      uncompressed[3]=(low         & 0x7f) - 63;
      j=4;
      break;
    case 2:
      //3 codes of 9 bits
      uncompressed[0]=((low >> 18) & 0x1ff) - 255;
      uncompressed[1]=((low >> 9)  & 0x1ff) - 255;
      uncompressed[2]=(low         & 0x1ff) - 255;
      j=3;
      break;
    case 1:
      //2 codes of 14 bits
      uncompressed[0]=((low >> 14) & 0x3fff) - 8191;
      uncompressed[1]=(low         & 0x3fff) - 8191;
      j=2;
      break;
    case 0:
      //1 code of 28 bits
      *uncompressed=low;
      j=1;
      break;
    default: /* should not reach this point */
      fprintf(stderr, "uncompressS9(): internal error (unexpected selector %d)\n", sel);
      exit(1);
  }

  tmp=*prev;
  for(i=0; i<j; ++i){
    tmp+=uncompressed[i];
    uncompressed[i]=tmp;
  }
  *prev=tmp;

  return j;
}
#endif /* USE_S9_COMPRESSION */


/* quads for CCS matrix */

/* instead of dynamically determining the indices of A's values needed to compute A^T*A in CCS
 * for a matrix A in CCS format (as splm_calc_AtA_ccsA() does), find and store groups of contiguous
 * nonzero indices for later use in multiplication. This strategy is faster at the expense of
 * using more memory. If a contiguous group of N index pairs starting at I & J in A and contributing
 * element at index K of AtA is found, then it is stored as a quadruple (N, I, J, K)
 *
 * It is assumed that the structure of AtA has been determined with a previous call to
 * splm_setup_AtA_ccsA(). Note also that most of the code is based on splm_calc_AtA_ccsA().
 *
 * Only the quadruples concerning the lower triangular & diagonal part of the product
 * are found. To save memory, the sequences formed by the N, I, J, K are independently
 * compressed with the vByte or S9 technique.
 */
int splm_get_AtA_ccsA_quads(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job, void *quads[4], int quadsz[4])
{
register int i, j, ii, jj, diff;
register int *AColptr, *AColptrP1, *ARowidx;
int *AtAColptr, *AtAColptrP1, *AtARowidx;
int nc, idx, k;
int nquads=0, curlen=0;
#ifdef USE_VBYTE_COMPRESSION
unsigned char *lenptr, *iiptr, *jjptr, *kptr; // len, ii0, jj0, k
int (*const compress_k)(int, int *, unsigned char *)=(job==1)? compressVB_01 : compressVB_ndecr;
//const int eltsz=sizeof(unsigned char);
int lastk[2]={0, 0};
#else
int32_t lenbuf[_UBUFSZS9]={0}, iibuf[_UBUFSZS9]={0}, jjbuf[_UBUFSZS9]={0}, kbuf[_UBUFSZS9]={0};
uint32_t *lenptr, *iiptr, *jjptr, *kptr; // len, ii0, jj0, k
//const int eltsz=sizeof(uint32_t);
int lastk=0;
#endif /* USE_VBYTE_COMPRESSION */
int lensz, iisz, jjsz, ksz; // allocated sizes for the above
int lenN=0, iiN=0, jjN=0, kN=0; // actual lengths for the above
int previi, prevjj, lastlen=0, lastii=0, lastjj=0; // note: previi, prevjj have a different role than lastii, lastjj!
int totlen=0;

  if(AtA->nnz<=0){
    fprintf(stderr, "product matrix must have been preallocated in splm_get_AtA_ccsA_quads()!\n");
    exit(1);
  }
  AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1] 
  AtARowidx=AtA->rowidx;

  nc=A->nc;
  AColptr=A->colptr; AColptrP1=AColptr+1; // for all i, AColptrP1[i]==AColptr[i+1] 
  ARowidx=A->rowidx;

  lensz=iisz=jjsz=ksz=AtA->nnz;
  /* in the worst case, a single compress invocation will need at most:
   *   a) sizeof(int) bytes for a compressed integer when using vByte
   *   b) _MAXCOMP 32 bit words for compressing _MAXCOMP large integers
   *      into separate words when using S-9
   * Thus, the extra sizes in the allocations below simplify checks by ensuring
   * that when xxN<xxsz, there exists enough space for calling the compression
   * function once more
   */
#ifdef USE_VBYTE_COMPRESSION
  lenptr=(unsigned char *)emalloc(lensz*sizeof(unsigned char)+sizeof(int));
  iiptr =(unsigned char *)emalloc(iisz*sizeof(unsigned char)+sizeof(int));
  jjptr =(unsigned char *)emalloc(jjsz*sizeof(unsigned char)+sizeof(int));
  kptr  =(unsigned char *)emalloc(ksz*sizeof(unsigned char)+sizeof(int));
#else
  lenptr=(uint32_t *)emalloc((lensz+_MAXCOMPS9)*sizeof(uint32_t));
  iiptr =(uint32_t *)emalloc((iisz+_MAXCOMPS9)*sizeof(uint32_t));
  jjptr =(uint32_t *)emalloc((jjsz+_MAXCOMPS9)*sizeof(uint32_t));
  kptr  =(uint32_t *)emalloc((ksz+_MAXCOMPS9)*sizeof(uint32_t));
#endif /* USE_VBYTE_COMPRESSION */

  /* find the quadruples corresponding to the lower triangle of AtA,
     which is computed (column by column) AtA[i][j] as \sum_k A[k][i]*A[k][j]
     NOTE: order of loops must remain as is, use of compressVB()/compressVB_ndecr() below
     depends upon it!
   */
  for(j=0; j<nc; ++j){

    if(job==1){
      /* assume that splm_setup_AtA_ccsA() has also been called with job==1,
       * so that the first element on column j is the one on the diagonal
       */
      //assert(AtARowidx[AtAColptr[j]]==j);
      k=AtAColptr[j];
    }
    else{
      /* ignore i<j, i.e. compute elements below or on the diagonal only */
      k=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
      k+=AtAColptr[j];
    }
    for(  ; k<AtAColptrP1[j]; ++k){ // loop over nonzero [A^T*A]_ij, i=0...nc-1 in the lower triangular part (i.e., i>=j)
      i=AtARowidx[k];
      previi=-2; prevjj=-2;

      for(ii=AColptr[i], jj=AColptr[j]; ii<AColptrP1[i] && jj<AColptrP1[j];  ){ // ii & jj are modified within the loop
        diff=ARowidx[ii]-ARowidx[jj];
        if(diff==0){
          /* the possibility of an inner product summing up to zero is ignored */
          if(ii==previi+1 && jj==prevjj+1){
            ++curlen;
          }else{ /* beginning of a new quadruple */
            /* check if more memory is needed and double current sizes if necessary */
            if(unlikely(lenN>=lensz)){
              lensz<<=1; //lensz*=2;
#ifdef USE_VBYTE_COMPRESSION
              lenptr=(unsigned char *)realloc(lenptr, lensz*sizeof(unsigned char)+sizeof(int));
#else
              lenptr=(uint32_t *)realloc(lenptr, (lensz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
              if(lenptr==NULL){
                fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_ccsA_quads() [len]\n", lensz);
                exit(1);
              }
            }
            if(unlikely(iiN>=iisz)){
              iisz<<=1; //iisz*=2;
#ifdef USE_VBYTE_COMPRESSION
              iiptr=(unsigned char *)realloc(iiptr, iisz*sizeof(unsigned char)+sizeof(int));
#else
              iiptr=(uint32_t *)realloc(iiptr, (iisz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
              if(iiptr==NULL){
                fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_ccsA_quads() [i]\n", iisz);
                exit(1);
              }
            }
            if(unlikely(jjN>=jjsz)){
              jjsz<<=1; //jjsz*=2;
#ifdef USE_VBYTE_COMPRESSION
              jjptr=(unsigned char *)realloc(jjptr, jjsz*sizeof(unsigned char)+sizeof(int));
#else
              jjptr=(uint32_t *)realloc(jjptr, (jjsz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
              if(jjptr==NULL){
                fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_ccsA_quads() [j]\n", jjsz);
                exit(1);
              }
            }
            if(unlikely(kN>=ksz)){
              ksz<<=1; //ksz*=2;
#ifdef USE_VBYTE_COMPRESSION
              kptr=(unsigned char *)realloc(kptr, ksz*sizeof(unsigned char)+sizeof(int));
#else
              kptr=(uint32_t *)realloc(kptr, (ksz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
              if(kptr==NULL){
                fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_ccsA_quads() [k]\n", ksz);
                exit(1);
              }
            }

            /* new quad of yet unknown length, starting at A's ii, jj and contributing to AtA's k */
            if(curlen>0){
#ifdef USE_VBYTE_COMPRESSION
              lenN+=compressVB(curlen, &lastlen, lenptr+lenN); // this stores the length of the previous quad!
#else
              lenN+=compressS9(curlen, &lastlen, lenbuf, 0, lenptr+lenN); // this stores the length of the previous quad!
#endif
              totlen+=curlen;
            }
            curlen=1;
#ifdef USE_VBYTE_COMPRESSION
            iiN+=compressVB(ii, &lastii, iiptr+iiN);
            jjN+=compressVB(jj, &lastjj, jjptr+jjN);
            kN +=(*compress_k)(k, lastk, kptr+kN);
#else
            iiN+=compressS9(ii, &lastii, iibuf, 0, iiptr+iiN);
            jjN+=compressS9(jj, &lastjj, jjbuf, 0, jjptr+jjN);
            kN +=compressS9_ndecr(k, &lastk, kbuf, 0, kptr+kN);
#endif
            ++nquads;
          }
          previi=ii; prevjj=jj;
          ++ii; ++jj;
          continue;
        } else if(diff<0){
          idx=binsearch_geq(ARowidx+ii, AColptrP1[i]-ii, ARowidx[jj]);
          if(idx<0) break;
          ii+=idx;
          continue;
        }
        else{ // diff>0
          idx=binsearch_geq(ARowidx+jj, AColptrP1[j]-jj, ARowidx[ii]);
          if(idx<0) break;
          jj+=idx;
          continue;
        }
      }
    }
  }
  totlen+=curlen;
#ifdef USE_VBYTE_COMPRESSION
  lenN+=compressVB(curlen, &lastlen, lenptr+lenN); // store the length of the last quad
  if(job==1 && nquads%8) ++kN; // count partially filled last byte when using compressVB_01()

  lenptr=(unsigned char *)realloc(lenptr, lenN*sizeof(unsigned char));
  iiptr=(unsigned char *)realloc(iiptr, iiN*sizeof(unsigned char));
  jjptr=(unsigned char *)realloc(jjptr, jjN*sizeof(unsigned char));
  kptr=(unsigned char *)realloc(kptr, kN*sizeof(unsigned char));
#else
  lenN+=compressS9(curlen, &lastlen, lenbuf, 0, lenptr+lenN); // store the length of the last quad

  /* flush circular buffers */
  lenN+=compressS9(0, &lastlen, lenbuf, 1, lenptr+lenN);
  iiN +=compressS9(0, &lastii, iibuf, 1, iiptr+iiN);
  jjN +=compressS9(0, &lastjj, jjbuf, 1, jjptr+jjN);
  kN +=compressS9_ndecr(0, &lastk, kbuf, 1, kptr+kN);

  lenptr=(uint32_t *)realloc(lenptr, lenN*sizeof(uint32_t));
  iiptr=(uint32_t *)realloc(iiptr, iiN*sizeof(uint32_t));
  jjptr=(uint32_t *)realloc(jjptr, jjN*sizeof(uint32_t));
  kptr=(uint32_t *)realloc(kptr, kN*sizeof(uint32_t));
#endif
  if(!lenptr || !iiptr || !jjptr || !kptr){
    fprintf(stderr, "one of the final memory reallocation requests failed in splm_get_AtA_ccsA_quads()\n");
    exit(1);
  }

#if 0
  fprintf(stderr, "splm_get_AtA_ccsA_quads(): found %d quads, avg. length %.2lf, %s compressed size %.2lf Mb, uncompressed %.2lf\n",
          nquads, totlen/(double)nquads, 
#ifdef USE_VBYTE_COMPRESSION
          "vByte",
#else
          "S9",
#endif
          (lenN+iiN+jjN+kN)*eltsz/1048576.0, nquads*sizeof(int)*4/1048576.0);
#endif

  quads[0]=(void *)lenptr; quads[1]=(void *)iiptr;
  quads[2]=(void *)jjptr; quads[3]=(void *)kptr;

  quadsz[0]=lenN; quadsz[1]=iiN;
  quadsz[2]=jjN; quadsz[3]=kN;

  return nquads;
}

void splm_calc_AtA_ccsA_quads(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job, void *quads[4], int nquads)
{
register int i, j;
int j0, ii0, jj0, k;
register double sum0, sum1, sum2, sum3;
double *AVal, *AtAVal;
register const double *pAii, *pAjj;
#ifdef USE_VBYTE_COMPRESSION
const unsigned char *lenptr, *iiptr, *jjptr, *kptr; // len, ii0, jj0, k
int (*const uncompress_k)(const unsigned char *, int *, int *)=(job==1)? uncompressVB_01 : uncompressVB_ndecr;
int lastk[2]={0, 0};
#else
const uint32_t *lenptr, *iiptr, *jjptr, *kptr; // len, ii0, jj0, k
int lenbuf[_MAXCOMPS9], iibuf[_MAXCOMPS9], jjbuf[_MAXCOMPS9], kbuf[_MAXCOMPS9];
int lenbufsz=0, iibufsz=0, jjbufsz=0, kbufsz=0;
int ilen=_MAXCOMPS9, iii=_MAXCOMPS9, ijj=_MAXCOMPS9, ik=_MAXCOMPS9;
int lastk=0;
#endif /* USE_VBYTE_COMPRESSION */
int lenN=0, iiN=0, jjN=0, kN=0;
int lastlen=0, lastii=0, lastjj=0;

  AVal=A->val;
  AtAVal=AtA->val;
#ifdef USE_VBYTE_COMPRESSION
  lenptr=(unsigned char *)quads[0]; iiptr=(unsigned char *)quads[1];
  jjptr=(unsigned char *)quads[2]; kptr=(unsigned char *)quads[3];
#else
  lenptr=(uint32_t *)quads[0]; iiptr=(uint32_t *)quads[1];
  jjptr=(uint32_t *)quads[2]; kptr=(uint32_t *)quads[3];
#endif

  /* clear AtA values */
  memset(AtAVal, 0, AtA->nnz*sizeof(double));

  for(i=nquads; i-->0;  ){
#ifdef USE_VBYTE_COMPRESSION
    lenN+=uncompressVB(lenptr+lenN, &lastlen, &j0);
    iiN +=uncompressVB(iiptr+iiN, &lastii, &ii0);
    jjN+=uncompressVB(jjptr+jjN, &lastjj, &jj0);
    kN +=(*uncompress_k)(kptr+kN, lastk, &k);
#else
    if(ilen>=lenbufsz){
      ilen=0; lenbufsz=uncompressS9(lenptr[lenN], &lastlen, lenbuf); ++lenN;
    }
    j0=lenbuf[ilen++];
    if(iii>=iibufsz){
      iii=0; iibufsz=uncompressS9(iiptr[iiN], &lastii, iibuf); ++iiN;
    }
    ii0=iibuf[iii++];
    if(ijj>=jjbufsz){
      ijj=0; jjbufsz=uncompressS9(jjptr[jjN], &lastjj, jjbuf); ++jjN;
    }
    jj0=jjbuf[ijj++];
    if(ik>=kbufsz){
      ik=0; kbufsz=uncompressS9_ndecr(kptr[kN], &lastk, kbuf); ++kN;
    }
    k=kbuf[ik++];
#endif /* USE_VBYTE_COMPRESSION */
    //printf("%d  %d %d  %d\n", j0, ii0, jj0, k);

    /* due to the small lengths of sequences to be multiplied,
     * loop unrolling and blocking didn't make any noticeable 
     * difference here. Instead, a scheme favoring sequences
     * up to eight elements long is used.
     */
    pAjj=AVal+jj0; pAii=AVal+ii0;
    sum0=sum1=sum2=sum3=0.0;
#if 0
    /* straightforward implementation */
    for(j=j0; j-->0; ++pAii, ++pAjj)
      sum0+=(*pAii)*(*pAjj);
    AtAVal[k]+=sum0;
#endif
    switch(j0){
      /* note the fallthrough and the lack of "+=" for cases 5-8 */
      case 8: sum0 =pAjj[7]*pAii[7];
      case 7: sum1 =pAjj[6]*pAii[6];
      case 6: sum2 =pAjj[5]*pAii[5];
      case 5: sum3 =pAjj[4]*pAii[4];
      case 4: sum0+=pAjj[3]*pAii[3];
      case 3: sum1+=pAjj[2]*pAii[2];
      case 2: sum2+=pAjj[1]*pAii[1];
      case 1: sum3+=pAjj[0]*pAii[0];
              AtAVal[k]+=(sum0+sum1)+(sum2+sum3);
              continue;

      default:{
        const double *lastp;

        /* depth 2 unrolling */
        if(j0&1){ // odd
          --j0;
          sum0=(*pAii)*(*pAjj);
          ++pAii; ++pAjj;
        }

        for(lastp=pAii+j0; pAii!=lastp; pAii+=2, pAjj+=2){
          sum0+=pAii[0]*pAjj[0];
          sum1+=pAii[1]*pAjj[1];
        }
        AtAVal[k]+=sum0+sum1;
      }
    }
  }

  if(job==1) return; /* no upper part */

{
int *AtAColptr, *AtAColptrP1, *AtARowidx;
register int ii, jj, diff;
int nc, k;

  nc=A->nc;
  AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1] 
  AtARowidx=AtA->rowidx;

  /* copy the elements below the diagonal to those above it */
  for(j=0; j<nc; ++j){
    int ilow, ihigh;

    /* find nonzero [A^T*A]_ij, i=0...nc-1 above the diagonal (i.e., i<j) */
    /* note that element _jj on the diagonal is guaranteed to be nonzero */
    jj=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
    jj+=AtAColptr[j];
    for(k=AtAColptr[j]; k<jj; ++k){
      i=AtARowidx[k];

      /* for all elements (i, j) copy the corresponding elements (j, i) to them */

      /* a binary search over all rows of column i is performed next, some speed
       * might be gained by searching below the diagonal only by setting ilow to
       *
       * ilow=binsearch_geq(AtARowidx+AtAColptr[i], AtAColptrP1[i]-AtAColptr[i], i) + AtAColptr[i];
       */
      ilow=AtAColptr[i];
      ihigh=AtAColptrP1[i];
      while(ilow<ihigh){
        /* the following line is a sanity check...*/
        //if(j<AtARowidx[ilow] || j>AtARowidx[ihigh-1]) break; // not found

        ii=(ilow + ihigh) >> 1; //(ilow+ihigh)/2;
        //ii=ilow+((ihigh-ilow)>>1); /* ensures no index overflows */
        diff=j-AtARowidx[ii];
        if(diff<0)
          ihigh=ii;
        else if(diff>0)
          ilow=ii+1;
        else{ // found
          AtAVal[k]=AtAVal[ii];
          break;
        }
      }
    }
  }
}
}

/* quads for CRS matrix */

/* instead of dynamically determining the indices of A's values needed to compute A^T*A in CCS
 * for a matrix A in CRS format (as splm_calc_AtAx_crsA() does), find and store groups of contiguous
 * nonzero indices for later use in multiplication. This strategy is faster at the expense of
 * using more memory. If a group of N index pairs I, J ... I, J+N in A that contribute to
 * elements at indices K ... K+N of AtA is found, then it is stored as a quadruple (N, I, J, K).
 * Note that the quadruples determined here are different from those in splm_get_AtA_ccsA_quads()
 * and compress less due to their increased irregularity
 *
 * It is assumed that the structure of AtA has been determined with a previous call to
 * splm_setup_AtA_crsA(). Note also that most of the code is based on splm_calc_AtAx_crsA().
 *
 * The first 'zerocols' columns of A are asumed to be equal to zero and
 * are therefore ignored
 *
 * Only the quadruples concerning the lower triangular & diagonal part of the product
 * are found. To save memory, the sequences formed by the N, I, J, K are independently
 * compressed with the vByte or S9 technique.
 */

int splm_get_AtA_crsA_quads(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, void *quads[4], int quadsz[4])
{
int nr, idx=0; /* -Wall */
register int i, j, k, r1, r2;
register int low, high;
int *ARowptr, *ARowptrP1, *AColidx, *AtAColptr, *AtAColptrP1, *AtARowidx, *AtARowidxP1;
int nquads=0, curlen=0;
#ifdef USE_VBYTE_COMPRESSION
unsigned char *lenptr, *r1ptr, *r2ptr, *idxptr;
//const int eltsz=sizeof(unsigned char);
#else
int32_t lenbuf[_UBUFSZS9]={0}, r1buf[_UBUFSZS9]={0}, r2buf[_UBUFSZS9]={0}, idxbuf[_UBUFSZS9]={0};
uint32_t *lenptr, *r1ptr, *r2ptr, *idxptr;
//const int eltsz=sizeof(uint32_t);
#endif /* USE_VBYTE_COMPRESSION */
int lensz, r1sz, r2sz, idxsz; // allocated sizes for the above
int lenN=0, r1N=0, r2N=0, idxN=0; // actual lengths for the above
int previdx, prevr2, lastlen=0, lastr1=0, lastr2=0, lastidx=0; // note: previdx, prevr2 have a different role than lastidx, lastr2!
int totlen=0;

  nr=A->nr;
  ARowptr=A->rowptr; ARowptrP1=A->rowptr+1; // for all i, ARowptrP1[i]==ARowptr[i+1]
  AColidx=A->colidx;

  AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1]
  AtARowidx=AtA->rowidx; AtARowidxP1=AtARowidx+1;

  lensz=r1sz=r2sz=idxsz=AtA->nnz;
  /* in the worst case, a single compress invocation will need at most:
   *   a) sizeof(int) bytes for a compressed integer when using vByte
   *   b) _MAXCOMP 32 bit words for compressing _MAXCOMP large integers
   *      into separate words when using S-9
   * Thus, the extra sizes in the allocations below simplify checks by ensuring
   * that when xxN<xxsz, there exists enough space for calling the compression
   * function once more
   */
#ifdef USE_VBYTE_COMPRESSION
  lenptr=(unsigned char *)emalloc(lensz*sizeof(unsigned char)+sizeof(int));
  r1ptr =(unsigned char *)emalloc(r1sz*sizeof(unsigned char)+sizeof(int));
  r2ptr =(unsigned char *)emalloc(r2sz*sizeof(unsigned char)+sizeof(int));
  idxptr=(unsigned char *)emalloc(idxsz*sizeof(unsigned char)+sizeof(int));
#else
  lenptr=(uint32_t *)emalloc((lensz+_MAXCOMPS9)*sizeof(uint32_t));
  r1ptr =(uint32_t *)emalloc((r1sz+_MAXCOMPS9)*sizeof(uint32_t));
  r2ptr =(uint32_t *)emalloc((r2sz+_MAXCOMPS9)*sizeof(uint32_t));
  idxptr=(uint32_t *)emalloc((idxsz+_MAXCOMPS9)*sizeof(uint32_t));
#endif /* USE_VBYTE_COMPRESSION */

  /* find the quadruples corresponding to the lower triangle of AtA,
     which is computed as A^T*A_ij = \sum_k A^T_ik * A_kj = \sum_k A_ki * A_kj.
     NOTE: order of loops must remain as is, use of compress below
     depends upon it!
   */
  for(k=nr; k-->0;  ){
    for(r1=binsearch_geq(AColidx+ARowptr[k], ARowptrP1[k]-ARowptr[k], zerocols) + ARowptr[k]; // ignore elements whose colidx < zerocols
                          r1<ARowptrP1[k]; ++r1){
      j=AColidx[r1]-zerocols;
      previdx=-2; prevr2=-2;
      /* j-th column of AtA += k-th row of A * AVal[r1] */
      /* ARowptr[k] instead of r1 below computes full product */
      for(r2=binsearch_geq(AColidx+r1, ARowptrP1[k]-r1, zerocols) + r1, // ignore elements whose colidx < zerocols
                          idx=-1; r2<ARowptrP1[k]; ++r2){
        i=AColidx[r2]-zerocols;
        /* find idx (in column j) s.t. AtARowidx[idx]==i */

        /* if we have previously found an element in a row of this column, check
         * to see if the element in the sought row is next to it. If not, search
         * for it using binary search: binsearch(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], i)
         */
        if(idx<=0 || AtARowidxP1[idx]!=i){
          low=(idx>=0)? idx : AtAColptr[j]; high=AtAColptrP1[j];

          while(low<high){
            idx=(low + high) >> 1; //(low+high)/2;
            //idx=low+((high-low)>>1); /* ensures no index overflows */
            if(i>AtARowidx[idx])
              low=idx+1;
            else
              high=idx;
          }
          idx=high;

#if 0
          if(high==AtAColptrP1[j] || AtARowidx[high]!=i){
            /* should not reach this point... */
            fprintf(stderr, "splm_get_AtA_crsA_quads() internal error: no element in column %d of A^t*A "
                            "with row index equal to %d found!\n\tmalformed Jacobian?\n", i+zerocols, j+zerocols);
            exit(1);
          }
#endif
        }
        else /* idx>0 && AtARowidxP1[idx]==i */
          ++idx;

        if(idx==previdx+1 && r2==prevr2+1){
          ++curlen;
        }else{ /* beginning of a new quadruple */
          /* check if more memory is needed and double current sizes if necessary */
          if(unlikely(lenN>=lensz)){
            lensz<<=1; //lensz*=2;
#ifdef USE_VBYTE_COMPRESSION
            lenptr=(unsigned char *)realloc(lenptr, lensz*sizeof(unsigned char)+sizeof(int));
#else
            lenptr=(uint32_t *)realloc(lenptr, (lensz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
            if(lenptr==NULL){
              fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_crsA_quads() [len]\n", lensz);
              exit(1);
            }
          }
          if(unlikely(r1N>=r1sz)){
            r1sz<<=1; //r1sz*=2;
#ifdef USE_VBYTE_COMPRESSION
            r1ptr=(unsigned char *)realloc(r1ptr, r1sz*sizeof(unsigned char)+sizeof(int));
#else
            r1ptr=(uint32_t *)realloc(r1ptr, (r1sz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
            if(r1ptr==NULL){
              fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_crsA_quads() [r1]\n", r1sz);
              exit(1);
            }
          }
          if(unlikely(r2N>=r2sz)){
            r2sz<<=1; //r2sz*=2;
#ifdef USE_VBYTE_COMPRESSION
            r2ptr=(unsigned char *)realloc(r2ptr, r2sz*sizeof(unsigned char)+sizeof(int));
#else
            r2ptr=(uint32_t *)realloc(r2ptr, (r2sz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
            if(r2ptr==NULL){
              fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_crsA_quads() [r2]\n", r2sz);
              exit(1);
            }
          }
          if(unlikely(idxN>=idxsz)){
            idxsz<<=1; //idxsz*=2;
#ifdef USE_VBYTE_COMPRESSION
            idxptr=(unsigned char *)realloc(idxptr, idxsz*sizeof(unsigned char)+sizeof(int));
#else
            idxptr=(uint32_t *)realloc(idxptr, (idxsz+_MAXCOMPS9)*sizeof(uint32_t));
#endif
            if(idxptr==NULL){
              fprintf(stderr, "memory reallocation request for %d bytes failed in splm_get_AtA_crsA_quads() [idx]\n", idxsz);
              exit(1);
            }
          }

          /* new quad of yet unknown length */
          if(curlen>0){
#ifdef USE_VBYTE_COMPRESSION
            lenN+=compressVB(curlen, &lastlen, lenptr+lenN); // this stores the length of the previous quad!
#else
            lenN+=compressS9(curlen, &lastlen, lenbuf, 0, lenptr+lenN); // this stores the length of the previous quad!
#endif
            totlen+=curlen;
          }
          curlen=1;
#ifdef USE_VBYTE_COMPRESSION
          r1N +=compressVB(r1, &lastr1, r1ptr+r1N);
          r2N +=compressVB(r2, &lastr2, r2ptr+r2N);
          idxN+=compressVB(idx, &lastidx, idxptr+idxN);
#else
          r1N +=compressS9(r1, &lastr1, r1buf, 0, r1ptr+r1N);
          r2N +=compressS9(r2, &lastr2, r2buf, 0, r2ptr+r2N);
          idxN+=compressS9(idx, &lastidx, idxbuf, 0, idxptr+idxN);
#endif
          ++nquads;
        }
        previdx=idx; prevr2=r2;
      }
    }
  }
  totlen+=curlen;
#ifdef USE_VBYTE_COMPRESSION
  lenN+=compressVB(curlen, &lastlen, lenptr+lenN); // store the length of the last quad

  lenptr=(unsigned char *)realloc(lenptr, lenN*sizeof(unsigned char));
  r1ptr=(unsigned char *)realloc(r1ptr, r1N*sizeof(unsigned char));
  r2ptr=(unsigned char *)realloc(r2ptr, r2N*sizeof(unsigned char));
  idxptr=(unsigned char *)realloc(idxptr, idxN*sizeof(unsigned char));
#else
  lenN+=compressS9(curlen, &lastlen, lenbuf, 0, lenptr+lenN); // store the length of the last quad

  /* flush circular buffers */
  lenN+=compressS9(0, &lastlen, lenbuf, 1, lenptr+lenN);
  r1N +=compressS9(0, &lastr1, r1buf, 1, r1ptr+r1N);
  r2N +=compressS9(0, &lastr2, r2buf, 1, r2ptr+r2N);
  idxN+=compressS9(0, &lastidx, idxbuf, 1, idxptr+idxN);

  lenptr=(uint32_t *)realloc(lenptr, lenN*sizeof(uint32_t));
  r1ptr=(uint32_t *)realloc(r1ptr, r1N*sizeof(uint32_t));
  r2ptr=(uint32_t *)realloc(r2ptr, r2N*sizeof(uint32_t));
  idxptr=(uint32_t *)realloc(idxptr, idxN*sizeof(uint32_t));
#endif /* USE_VBYTE_COMPRESSION */
  if(!lenptr || !r1ptr || !r2ptr || !idxptr){
    fprintf(stderr, "one of the final memory reallocation requests failed in splm_get_AtA_crsA_quads()\n");
    exit(1);
  }

#if 0
  fprintf(stderr, "splm_get_AtA_crsA_quads(): found %d quads, avg. length %.2lf, %s compressed size %.2lf Mb, uncompressed %.2lf\n",
          nquads, totlen/(double)nquads,
#ifdef USE_VBYTE_COMPRESSION
          "vByte",
#else
          "S9",
#endif
          (lenN+r1N+r2N+idxN)*eltsz/1048576.0, nquads*sizeof(int)*4/1048576.0);
#endif

  quads[0]=(void *)lenptr; quads[1]=(void *)r1ptr;
  quads[2]=(void *)r2ptr; quads[3]=(void *)idxptr;

  quadsz[0]=lenN; quadsz[1]=r1N;
  quadsz[2]=r2N; quadsz[3]=idxN;

  return nquads;
}


void splm_calc_AtA_crsA_quads(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, void *quads[4], int nquads)
{
register int i, j;
int j0, r1, r2, idx, nc;
double *AVal, *AtAVal;
#ifdef USE_VBYTE_COMPRESSION
const unsigned char *lenptr, *r1ptr, *r2ptr, *idxptr; // len, r1, r2, idx
#else
const uint32_t *lenptr, *r1ptr, *r2ptr, *idxptr; // len, r1, r2, idx
int lenbuf[_MAXCOMPS9], r1buf[_MAXCOMPS9], r2buf[_MAXCOMPS9], idxbuf[_MAXCOMPS9];
int lenbufsz=0, r1bufsz=0, r2bufsz=0, idxbufsz=0;
int ilen=_MAXCOMPS9, ir1=_MAXCOMPS9, ir2=_MAXCOMPS9, iidx=_MAXCOMPS9;
#endif /* USE_VBYTE_COMPRESSION */
int lenN=0, r1N=0, r2N=0, idxN=0;
int lastlen=0, lastr1=0, lastr2=0, lastidx=0;

  AVal=A->val;
  AtAVal=AtA->val;
#ifdef USE_VBYTE_COMPRESSION
  lenptr=(unsigned char *)quads[0]; r1ptr=(unsigned char *)quads[1];
  r2ptr=(unsigned char *)quads[2]; idxptr=(unsigned char *)quads[3];
#else
  lenptr=(uint32_t *)quads[0]; r1ptr=(uint32_t *)quads[1];
  r2ptr=(uint32_t *)quads[2]; idxptr=(uint32_t *)quads[3];
#endif /* USE_VBYTE_COMPRESSION */

  /* clear AtA values */
  memset(AtAVal, 0, AtA->nnz*sizeof(double));

  for(i=nquads; i-->0;  ){
    register const double *pA;
    register double alpha, *pAtA;

#ifdef USE_VBYTE_COMPRESSION
    lenN+=uncompressVB(lenptr+lenN, &lastlen, &j0);
    r1N +=uncompressVB(r1ptr+r1N, &lastr1, &r1);
    r2N +=uncompressVB(r2ptr+r2N, &lastr2, &r2);
    idxN+=uncompressVB(idxptr+idxN, &lastidx, &idx);
#else
    if(ilen>=lenbufsz){
      ilen=0; lenbufsz=uncompressS9(lenptr[lenN], &lastlen, lenbuf); ++lenN;
    }
    j0=lenbuf[ilen++];
    if(ir1>=r1bufsz){
      ir1=0; r1bufsz=uncompressS9(r1ptr[r1N], &lastr1, r1buf); ++r1N;
    }
    r1=r1buf[ir1++];
    if(ir2>=r2bufsz){
      ir2=0; r2bufsz=uncompressS9(r2ptr[r2N], &lastr2, r2buf); ++r2N;
    }
    r2=r2buf[ir2++];
    if(iidx>=idxbufsz){
      iidx=0; idxbufsz=uncompressS9(idxptr[idxN], &lastidx, idxbuf); ++idxN;
    }
    idx=idxbuf[iidx++];
#endif /* USE_VBYTE_COMPRESSION */

    //printf("%d  %d %d  %d\n", j0, r1, r2, idx);

    /* due to the small lengths of sequences to be multiplied,
     * loop unrolling and blocking didn't make any noticeable 
     * difference here. Instead, a scheme favoring sequences
     * up to eight elements long is used.
     */
    alpha=AVal[r1];
    pAtA=AtAVal+idx;
    pA=AVal+r2;
    switch(j0){
      /* note the fallthrough */
      case 8: pAtA[7]+=alpha*pA[7];
      case 7: pAtA[6]+=alpha*pA[6];
      case 6: pAtA[5]+=alpha*pA[5];
      case 5: pAtA[4]+=alpha*pA[4];
      case 4: pAtA[3]+=alpha*pA[3];
      case 3: pAtA[2]+=alpha*pA[2];
      case 2: pAtA[1]+=alpha*pA[1];
      case 1: pAtA[0]+=alpha*pA[0];
              continue;

      default:
        /* straightforward implementation */
        for(j=j0; j-->0; ++pAtA, ++pA)
          *pAtA+=alpha*(*pA);
    }
  }

  if(job==1) return; /* no upper part */

  /* see if strictly upper part is needed as well */
{
    int ilow, ihigh;
    register int ii, jj, k, diff;
    int *AtAColptr, *AtAColptrP1, *AtARowidx;

    AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1]
    AtARowidx=AtA->rowidx;
    nc=A->nc-zerocols;

    /* copy the elements below the diagonal to those above it */
    for(j=0; j<nc; ++j){
      /* find nonzero [A^T*A]_ij, i=0...nc-1 above the diagonal (i.e., i<j) */
      /* note that element _jj on the diagonal is guaranteed to be nonzero */
      jj=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
      jj+=AtAColptr[j];
      for(k=AtAColptr[j]; k<jj; ++k){
        i=AtARowidx[k];

        /* for all elements (i, j) copy the corresponding elements (j, i) to them */
        /* a binary search over all rows of column i is performed next, some speed
         * might be gained by searching below the diagonal only by setting ilow to
         *
         * ilow=binsearch_geq(AtARowidx+AtAColptr[i], AtAColptrP1[i]-AtAColptr[i], i) + AtAColptr[i];
         */

        ilow=AtAColptr[i];
        ihigh=AtAColptrP1[i];
        while(ilow<ihigh){
          /* the following line is a sanity check...*/
          //if(j<AtARowidx[ilow] || j>AtARowidx[ihigh-1]) break; // not found

          ii=(ilow + ihigh) >> 1; //(ilow+ihigh)/2;
          //ii=ilow+((ihigh-ilow)>>1); /* ensures no index overflows */
          diff=j-AtARowidx[ii];
          if(diff<0)
            ihigh=ii;
          else if(diff>0)
            ilow=ii+1;
          else{ // found
            AtAVal[k]=AtAVal[ii];
            break;
          }
        }
      }
    }
  }
}

/* compute y=A^T*x for a matrix A in CRS format and a vector x
 * y is computed as y_j += \sum_k A_kj * e_k
 */
void splm_calc_Atx_crsA(struct splm_crsm *A, int zerocols, double *x, double *y)
{
register int k, r1;
int j, nr, nc, *ARowptr, *ARowptrP1, *AColidx;
double *AVal;
register double xk;

  nr=A->nr; nc=A->nc-zerocols;
  ARowptr=A->rowptr; ARowptrP1=A->rowptr+1; // for all i, ARowptrP1[i]==ARowptr[i+1]
  AColidx=A->colidx;
  AVal=A->val;

  /* clear y values */
  memset(y, 0, nc*sizeof(double));

  for(k=nr; k-->0;  ){
    xk=x[k];
    for(r1=binsearch_geq(AColidx+ARowptr[k], ARowptrP1[k]-ARowptr[k], zerocols) + ARowptr[k]; // ignore elements whose colidx < zerocols
              r1<ARowptrP1[k]; ++r1){
      j=AColidx[r1]-zerocols;
      y[j]+=AVal[r1]*xk;
    }
  }
}
#endif /* USE_ATA_QUADS */

/************************** At*A and A^t*x for a CRS matrix A **************************/


/* compute the CCS structure of A^T*A for a matrix A in CRS format. if AtA->nnz>0, assumes that AtA
 * is large enough to hold the product. Otherwise, it computes its number of nonzeros and allocates
 * it accordingly.
 *
 * If job is 1, computes just the structure of the lower triangular & diagonal part of the product
 * If job is > 1, computes the structure of the full matrix product (i.e., all L, D & U parts)
 */
void splm_setup_AtA_crsA(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job)
{
register int i, j, k, jj;
register int *AtRowptr, *AColidx;
int *ARowptr, *AtColidx, *ARowptrP1;
int nr, nc, nnz;
int *colcounts; // counters for the number of nonzeros in each column
struct splm_ccsm ccsA;

  nr=A->nr; nc=A->nc;
  nnz=A->nnz;
  ARowptr=A->rowptr; AColidx=A->colidx;
  ARowptrP1=ARowptr+1; // for all i, ARowptrP1[i]==ARowptr[i+1]

  /* compute the nonzero structure of A^T in AtRowptr & AtColidx */

  colcounts=(int *)calloc(nc, sizeof(int)); /* init to zero */
  AtRowptr=(int *)malloc((nc+1+nnz)*sizeof(int));
  if(!colcounts || !AtRowptr){
    fprintf(stderr, "memory allocation request failed in splm_setup_AtA_crsA() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
    exit(1);
  }
  AtColidx=AtRowptr+nc+1;
  
  /* 1st pass: count #nonzeros in each column */
  for(j=ARowptr[nr]; j-->0;  )
    ++(colcounts[AColidx[j]]);

  /* 2nd pass: setup AtRowptr, AtColidx */
  for(j=k=0; j<nc; ++j){
    AtRowptr[j]=k;
    k+=colcounts[j];
    colcounts[j]=0;
  }
  AtRowptr[nc]=nnz;

  /* colcounts[j] will count the #nonzeros in col. j seen before the current row; note that it is cleared above */
  for(i=0; i<nr; ++i){
    for(j=ARowptr[i]; j<ARowptrP1[i]; ++j){
      jj=AColidx[j];
      k=AtRowptr[jj];
      k+=colcounts[jj]++;
      AtColidx[k]=i;
    }
  }

  free(colcounts); colcounts=NULL;
  ARowptr=AColidx=NULL; 

  /* AtRowptr & AtColidx now contain the nonzero structure of A^T */

  /* drop first "zerocols" rows from At (which correspond to A's first "zerocols" columns) */
  if(zerocols>0){
    int nrnnz;

    /* number of nonzeros in the first "zerocols" rows */
    nrnnz=AtRowptr[zerocols];

    /* adjust the contents of the row pointers */
    for(i=zerocols; i<nc+1; ++i)
      AtRowptr[i]-=nrnnz;

    nnz-=nrnnz;
    nc-=zerocols;

    /* advance pointers */
    AtRowptr+=zerocols;
    AtColidx+=nrnnz;
  }

  /* the CRS structure for A^T in AtRowptr, AtColidx is equivalent to
   * the CCS structure for A, hence splm_setup_AtA_ccsA() can be used
   */
  ccsA.nr=nr; ccsA.nc=nc; 
  ccsA.nnz=nnz;
  ccsA.rowidx=AtColidx;
  ccsA.colptr=AtRowptr;
  ccsA.val=NULL;

  splm_setup_AtA_ccsA(&ccsA, AtA, job);

  /* undo the action of dropping "zerocols" rows (which contained "nrnnz" nonzero elements) */
  if(zerocols>0){
    /* decrement AtRowptr, dont bother to adjust AtColidx and its contents */
    AtRowptr-=zerocols;
  }

  free(AtRowptr);
}


/* compute A^T*A in CCS for a matrix A in CRS format. It is assumed that the structure
 * of AtA has been determined with a previous call to splm_setup_AtA_crsA().
 *
 * Also, compute y=A^T*x
 *
 * The first 'zerocols' columns of A are asumed to be equal to zero and
 * are therefore ignored
 *
 * If job is 1, computes just the lower triangular & diagonal part of the product
 * If job is > 1, computes the full matrix product (i.e., all L, D & U parts)
 *
 * Algorithm:
 * A^T*A_ij = \sum_k A^T_ik * A_kj = \sum_k A_ki * A_kj.
 * Thus, the product can be implemented using an outer loop
 * for k that adds A_ki*A_kj to each element ij of the result.
 * If either A_ki or A_kj are zero, nothing is added to AtA_ij.
 * Note that with this scheme, the accesses to A and AtA are
 * always along rows, thus less likely to induce cache misses
 *
 * y=A^T*x is computed as y_j += \sum_k A_kj * e_k
 *
 */
void splm_calc_AtAx_crsA(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, double *const x, double *const y)
{
int nr, nc, idx=0; /* -Wall */
register int i, j, k, r1, r2;
register int low, high;
double *AVal, *AtAVal;
register double alpha;
int *ARowptr, *ARowptrP1, *AColidx, *AtAColptr, *AtAColptrP1, *AtARowidx, *AtARowidxP1;

  nr=A->nr; nc=A->nc-zerocols;
  ARowptr=A->rowptr; ARowptrP1=A->rowptr+1; // for all i, ARowptrP1[i]==ARowptr[i+1]
  AColidx=A->colidx;
  AVal=A->val;

  AtAColptr=AtA->colptr; AtAColptrP1=AtAColptr+1; // for all i, AtAColptrP1[i]==AtAColptr[i+1]
  AtARowidx=AtA->rowidx; AtARowidxP1=AtARowidx+1;
  AtAVal=AtA->val;

  /* clear AtA values */
  memset(AtAVal, 0, AtA->nnz*sizeof(double));

  /* clear y values */
  memset(y, 0, nc*sizeof(double));

  /* compute product's lower part */
  for(k=nr; k-->0;  ){
    for(r1=binsearch_geq(AColidx+ARowptr[k], ARowptrP1[k]-ARowptr[k], zerocols) + ARowptr[k]; // ignore elements whose colidx < zerocols
                          r1<ARowptrP1[k]; ++r1){
      j=AColidx[r1]-zerocols;
      alpha=AVal[r1];
      /* j-th column of AtA += k-th row of A * alpha */
      /* ARowptr[k] instead of r1 below computes full product */
      for(r2=binsearch_geq(AColidx+r1, ARowptrP1[k]-r1, zerocols) + r1, // ignore elements whose colidx < zerocols
                          idx=-1; r2<ARowptrP1[k]; ++r2){
        i=AColidx[r2]-zerocols;
        /* find idx (in column j) s.t. AtARowidx[idx]==i */

        /* if we have previously found an element in a row of this column, check
         * to see if the element in the sought row is next to it. If not, search
         * for it using binary search: binsearch(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], i)
         */
        if(idx<=0 || AtARowidxP1[idx]!=i){
          low=(idx>=0)? idx : AtAColptr[j]; high=AtAColptrP1[j];

          while(low<high){
            idx=(low + high) >> 1; //(low+high)/2;
            //idx=low+((high-low)>>1); /* ensures no index overflows */
            if(i>AtARowidx[idx])
              low=idx+1;
            else
              high=idx;
          }
          idx=high;

#if 0
          if(high==AtAColptrP1[j] || AtARowidx[high]!=i){
            /* should not reach this point... */
            fprintf(stderr, "splm_calc_AtAx_crsA() internal error: no element in column %d of A^t*A "
                            "with row index equal to %d found!\n\tmalformed Jacobian?\n", i+zerocols, j+zerocols);
            exit(1);
          }
#endif
        }
        else /* idx>0 && AtARowidxP1[idx]==i */
          ++idx;

        AtAVal[idx]+=AVal[r2]*alpha;
      }
      y[j]+=alpha*x[k];
    }
  }

  /* see if strictly upper part is needed as well */
  if(job>1){
    int ilow, ihigh;
    register int ii, jj, diff;

    /* copy the elements below the diagonal to those above it */
    for(j=0; j<nc; ++j){
      /* find nonzero [A^T*A]_ij, i=0...nc-1 above the diagonal (i.e., i<j) */
      /* note that element _jj on the diagonal is guaranteed to be nonzero */
      jj=binsearch_geq(AtARowidx+AtAColptr[j], AtAColptrP1[j]-AtAColptr[j], j);
      jj+=AtAColptr[j];
      for(k=AtAColptr[j]; k<jj; ++k){
        i=AtARowidx[k];

        /* for all elements (i, j) copy the corresponding elements (j, i) to them */
        /* a binary search over all rows of column i is performed next, some speed
         * might be gained by searching below the diagonal only by setting ilow to
         *
         * ilow=binsearch_geq(AtARowidx+AtAColptr[i], AtAColptrP1[i]-AtAColptr[i], i) + AtAColptr[i];
         */

        ilow=AtAColptr[i];
        ihigh=AtAColptrP1[i];
        while(ilow<ihigh){
          /* the following line is a sanity check...*/
          //if(j<AtARowidx[ilow] || j>AtARowidx[ihigh-1]) break; // not found

          ii=(ilow + ihigh) >> 1; //(ilow+ihigh)/2;
          //ii=ilow+((ihigh-ilow)>>1); /* ensures no index overflows */
          diff=j-AtARowidx[ii];
          if(diff<0)
            ihigh=ii;
          else if(diff>0)
            ilow=ii+1;
          else{ // found
            AtAVal[k]=AtAVal[ii];
            break;
          }
        }
      }
    }
  }
}

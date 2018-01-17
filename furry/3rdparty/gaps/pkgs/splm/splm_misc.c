/////////////////////////////////////////////////////////////////////////////////
////// 
//////  Miscelaneous functions for sparse Levenberg - Marquardt minimization  
//////  Copyright (C) 2008-2011  Manolis Lourakis (lourakis at ics.forth.gr)
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
#include <float.h>
#include <time.h>

#include "splm.h"
#include "splm_priv.h"

#define FABS(x)           (((x)>=0)? (x) : -(x))

/* auxiliary memory allocation routine with error checking */
void *splm_emalloc(char *file, int line, size_t sz)
{
void *ptr;

  ptr=(void *)malloc(sz);
  if(ptr==NULL){
    fprintf(stderr, "memory allocation request for %u bytes failed in file %s, line %d, exiting", (unsigned int) sz, file, line);
    exit(1);
  }

  return ptr;
}

/* auxiliary memory re-allocation routine with error checking */
void *splm_realloc(char *file, int line, void *oldptr, size_t sz)
{
void *ptr;

  ptr=(void *)realloc(oldptr, sz);
  if(ptr==NULL){
    fprintf(stderr, "memory re-allocation request for %u bytes failed in file %s, line %d, exiting\n", (unsigned int) sz, file, line);
    exit(1);
  }

  return ptr;
}

#define _HORZ_OFFSET 5 // horizontal offset from left edge of paper in PostScript points (~ 1 mm); >=1
#define _VERT_OFFSET 5 // vertical offset from lower edge of paper in PostScript points (~ 2 mm)

#define _A4_WIDTH    595 // A4 paper width in PostScript points (210 mm, 1 pt = 1/72 inch = 25.4/72 mm)
#define _A4_HEIGHT   842 // A4 paper height in PostScript points (297 mm)

/* 1.1x1.1 looks better than 1.0x1.0 */
#define _DOT_WIDTH   1.0
#define _DOT_HEIGHT  1.0


/* saves to a file the nonzero pattern of a sparse CRS/CCS matrix in
 * EPS (Encapsulated PostScript) format.
 *
 * regions specifies a set of rectangles (using ul_i ul_j, br_i, br_j) coordinates
 * that are to be drawn with the RGB triplets in colors. NOTE: it is the
 * responsibility of the caller to ensure that the specified regions 
 * cover the whole matrix!
 * Setting to NULL draws everything in black.
 */
static void splm_crcsm2eps(int isCRS, void *sm, int (*regions)[4], double (*colors)[3], int nregions, char *fname)
{
struct splm_crsm *crs=NULL;
struct splm_ccsm *ccs=NULL;
register int i, j, k;
int nrows, ncols, brect=1;
int dispw, disph; // size of displayed array in pts
FILE *fp;
time_t tim;
double scx, scy, scl;

  if(isCRS){
    crs=(struct splm_crsm *)sm;
    nrows=crs->nr;
    ncols=crs->nc;
  }
  else{
    ccs=(struct splm_ccsm *)sm;
    nrows=ccs->nr;
    ncols=ccs->nc;
  }

  if((fp=fopen(fname, "w"))==NULL){
    fprintf(stderr, "splm_crcsm2eps(): failed to open file %s for writing\n", fname);
    exit(1);
  }

  scx=((double)(_A4_WIDTH-2*_HORZ_OFFSET))/(double)(ncols);
  scy=((double)(_A4_HEIGHT-2*_VERT_OFFSET))/(double)(nrows);
  scl=(scx<=scy)? scx : scy; // minimum
  dispw=(int)(ncols*scl);
  disph=(int)(nrows*scl);

  /* print EPS preamble */
  time(&tim);
  fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n%%%%Title: %s\n%%%%Creator: splm ver. %s\n", fname, SPLM_VERSION);
  fprintf(fp, "%%%%Nonzero pattern for a %d x %d %s matrix\n", nrows, ncols, isCRS? "CRS" : "CCS");
  fprintf(fp, "%%%%CreationDate: %s", ctime(&tim));
  if(dispw>=disph) // "landscape"
    fprintf(fp, "%%%%BoundingBox: 0 %d %d %d\n", _A4_HEIGHT-2*_VERT_OFFSET-disph, _A4_WIDTH, _A4_HEIGHT);
  else // "portrait"
    fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", dispw+2*_HORZ_OFFSET, _A4_HEIGHT);
  fprintf(fp, "%%%%Magnification: 1.0\n%%%%Page: 1 1\n%%%%EndComments\n\n");
  fprintf(fp, "/origstate save def\ngsave\n0 setgray\n");
  if(brect) fprintf(fp, "%d %d %d %d rectstroke\n", _HORZ_OFFSET-1, _A4_HEIGHT-_VERT_OFFSET-disph-1, dispw+1, disph+1);

  /* move the coordinate system to the upper left corner with axes pointing as shown below:
   * +------> y
   * |
   * |
   * |
   * v x
   */
  fprintf(fp, "%d %d translate\n-90 rotate\n", _HORZ_OFFSET, _A4_HEIGHT-_VERT_OFFSET);
  /* scale the coordinate system so that the matrix's pattern fits in the narrowest page dimension */
  fprintf(fp, "%.4lf %.4lf scale\n", scl, scl);
  //if(brect) fprintf(fp, "%d %d %d %d rectstroke\n", _HORZ_OFFSET, _VERT_OFFSET, nrows, ncols);

  /* define dot dimensions */
  fprintf(fp, "/w %g def\n/h %g def\n", _DOT_WIDTH, _DOT_HEIGHT);
  /* define shorthand for rectfill */
  fprintf(fp, "/R { w h rectfill } bind def\n\n"); // x y assumed already in PS stack

  if(nregions<=0) regions=NULL;

  if(isCRS){
    if(!regions)
      for(i=0; i<nrows; ++i)
        for(j=crs->rowptr[i]; j<crs->rowptr[i+1]; ++j)
          fprintf(fp, "%d %d R\n", i, crs->colidx[j]);
    else
      for(k=0; k<nregions; ++k){
        fprintf(fp, "%g %g %g setrgbcolor\n", colors[k][0], colors[k][1], colors[k][2]);
        for(i=regions[k][0]; i<regions[k][2]; ++i)
          for(j=crs->rowptr[i]; j<crs->rowptr[i+1]; ++j)
            if(crs->colidx[j]>=regions[k][1] && crs->colidx[j]<regions[k][3])
              fprintf(fp, "%d %d R\n", i, crs->colidx[j]);
      }
  }
  else{
    if(!regions)
      for(j=0; j<ncols; ++j)
        for(i=ccs->colptr[j]; i<ccs->colptr[j+1]; ++i)
          fprintf(fp, "%d %d R\n", ccs->rowidx[i], j);
    else
      for(k=0; k<nregions; ++k){
        fprintf(fp, "%g %g %g setrgbcolor\n", colors[k][0], colors[k][1], colors[k][2]);
        for(j=regions[k][1]; j<regions[k][3]; ++j)
          for(i=ccs->colptr[j]; i<ccs->colptr[j+1]; ++i)
            if(ccs->rowidx[i]>=regions[k][0] && ccs->rowidx[i]<regions[k][2])
              fprintf(fp, "%d %d R\n", ccs->rowidx[i], j);
      }
  }

  fprintf(fp, "grestore\norigstate restore\n\n%%%%Trailer");
  fclose(fp);
}


/* interfaces to splm_crcsm2eps() */

void splm_crsm2eps(struct splm_crsm *sm, int (*regions)[4], double (*colors)[3], int nregions, char *fname)
{
  splm_crcsm2eps(1, (void *)sm, regions, colors, nregions, fname);
}

void splm_ccsm2eps(struct splm_ccsm *sm, int (*regions)[4], double (*colors)[3], int nregions, char *fname)
{
  splm_crcsm2eps(0, (void *)sm, regions, colors, nregions, fname);
}


/* 
 * Check the CCS Jacobian of a n-valued nonlinear function in m variables
 * evaluated at a point p, for consistency with the function itself.
 *
 * Based on fortran77 subroutine CHKDER by
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
 * Argonne National Laboratory. MINPACK project. March 1980.
 *
 *
 * func points to a function from R^m --> R^n: Given a p in R^m it yields hx in R^n
 * fjac is the Jacobian of func at p, whose correctness is to be tested.
 *     Note that row i of fjac corresponds to the gradient of
 *     the i-th component of func, evaluated at p.
 * p is an input array of length m containing the point of evaluation.
 * m is the number of variables
 * n is the number of functions
 * adata points to possible additional data and is passed uninterpreted
 *     to func.
 * err is an array of length n. On output, err contains measures
 *     of correctness of the respective gradients. if there is
 *     no severe loss of significance, then if err[i] is 1.0 the
 *     i-th gradient is correct, while if err[i] is 0.0 the i-th
 *     gradient is incorrect. For values of err between 0.0 and 1.0,
 *     the categorization is less certain. In general, a value of
 *     err[i] greater than 0.5 indicates that the i-th gradient is
 *     probably correct, while a value of err[i] less than 0.5
 *     indicates that the i-th gradient is probably incorrect.
 *
 *
 * The function does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the components
 * of p should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 */

static void sparselm_chkjac_core(
    void (*func)(double *p, double *hx, int m, int n, void *adata),
    struct splm_ccsm *fjac,
    double *p, int m, int n, int jnnz, void *adata, double *err)
{
double factor=100.0;
double one=1.0;
double zero=0.0;
double *fvec, *pp, *fvecp, *buf;

register int i, j;
double eps, epsf, temp, epsmch;
double epslog;
int fvec_sz=n, pp_sz=m, fvecp_sz=n;

  epsmch=DBL_EPSILON;
  eps=(double)sqrt(epsmch);

  buf=(double *)malloc((fvec_sz + pp_sz + fvecp_sz)*sizeof(double));
  if(!buf){
    fprintf(stderr, "sparselm_chkjac_core(): memory allocation request failed\n");
    exit(1);
  }
  fvec=buf;
  pp=fvec+fvec_sz;
  fvecp=pp+pp_sz;

  /* compute fvec=func(p) */
  (*func)(p, fvec, m, n, adata);

  /* compute pp */
  for(j=0; j<m; ++j){
    temp=eps*FABS(p[j]);
    if(temp==zero) temp=eps;
    pp[j]=p[j]+temp;
  }

  /* compute fvecp=func(pp) */
  (*func)(pp, fvecp, m, n, adata);

  epsf=factor*epsmch;
  epslog=(double)log10(eps);

  for(i=0; i<n; ++i)
    err[i]=zero;

  for(j=0; j<m; ++j){
    temp=FABS(p[j]);
    if(temp==zero) temp=one;

    for(i=fjac->colptr[j]; i<fjac->colptr[j+1]; ++i) // nonzero elements in column j
      err[fjac->rowidx[i]]+=temp*fjac->val[i];
  }

  for(i=0; i<n; ++i){
    temp=one;
    if(fvec[i]!=zero && fvecp[i]!=zero && FABS(fvecp[i]-fvec[i])>=epsf*FABS(fvec[i]))
        temp=eps*FABS((fvecp[i]-fvec[i])/eps - err[i])/(FABS(fvec[i])+FABS(fvecp[i]));
    err[i]=one;
    if(temp>epsmch && temp<eps)
        err[i]=((double)log10(temp) - epslog)/epslog;
    if(temp>=eps) err[i]=zero;
  }

  free(buf);

  return;
}

/* error checking for a CCS Jacobian */
void sparselm_chkjacccs(
    void (*func)(double *p, double *hx, int m, int n, void *adata),
    void (*jacf)(double *p, struct splm_ccsm *jac, int m, int n, void *adata),
    double *p, int m, int n, int jnnz, void *adata, double *err)
{
struct splm_ccsm fjac;

  splm_ccsm_alloc(&fjac, n, m, jnnz);

  /* compute the Jacobian at p */
  (*jacf)(p, &fjac, m, n, adata);

  sparselm_chkjac_core(func, &fjac, p, m, n, jnnz, adata, err);

  splm_ccsm_free(&fjac);
}

/* as above for a CRS Jacobian */
void sparselm_chkjaccrs(
    void (*func)(double *p, double *hx, int m, int n, void *adata),
    void (*jacf)(double *p, struct splm_crsm *jac, int m, int n, void *adata),
    double *p, int m, int n, int jnnz, void *adata, double *err)
{
struct splm_crsm crsjac;
struct splm_ccsm ccsjac;

  splm_crsm_alloc(&crsjac, n, m, jnnz);
  splm_ccsm_alloc(&ccsjac, n, m, jnnz);

  /* compute the Jacobian at p */
  (*jacf)(p, &crsjac, m, n, adata);

  splm_crsm2ccsm(&crsjac, &ccsjac);
  sparselm_chkjac_core(func, &ccsjac, p, m, n, jnnz, adata, err);

  splm_crsm_free(&crsjac);
  splm_ccsm_free(&ccsjac);
}

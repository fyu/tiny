/////////////////////////////////////////////////////////////////////////////////
//// 
////  Sparse Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2004-2010  Manolis Lourakis (lourakis at ics.forth.gr)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "splm.h"
#include "splm_priv.h"
#include "compiler.h"

#define SPLM_EPSILON       1E-12
#define SPLM_EPSILON_SQ    ( (SPLM_EPSILON)*(SPLM_EPSILON) )

#define SPLM_ONE_THIRD     0.3333333334 /* 1.0/3.0 */

#define FABS(x)           (((x)>=0)? (x) : -(x))

#define HAVE_CSPARSE      // Funk added this

/* attempt to guess the Jacobian's non-zero pattern 
 * The idea is to add a small value to each parameter in turn
 * and identify the observations that are influenced.
 *
 * This function should be used with caution as it cannot guarantee
 * that the true non-zero pattern will be found. Furthermore, it can
 * give rise to domain errors.
 *
 * Returns the number of nonzero elements found
 */
static int splm_jacpatguess(void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                             double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata,
                             double *hx, double delta
                            )
{
register int i, j, k;
register double d, tmp;
double *hxx;
int *colptr, *rowidx;

  colptr=jac->colptr;
  rowidx=jac->rowidx;

  (*func)(p, hx, nvars, nobs, adata); // hx=f(p)
  hxx=(double *)emalloc(nobs*sizeof(double));
  
  for(j=k=0; j<nvars; ++j){
    colptr[j]=k;

    /* determine d=max(SPLM_DELTA_SCALE*|p[j]|, delta), see HZ */
    d=SPLM_DELTA_SCALE*p[j]; // force evaluation
    d=FABS(d);
    if(d<delta) d=delta;

    tmp=p[j]; p[j]+=d;
    (*func)(p, hxx, nvars, nobs, adata); // hxx=f(p+d)
    p[j]=tmp; /* restore */

    for(i=0; i<nobs; ++i){
      tmp=FABS(hxx[i]-hx[i]);
      if(tmp>0.0){ // tmp>1E-07*d
        /* element (i, j) of jac != 0 */
        //printf("observation %d depends on variable %d [%.10lf]\n", i, j, tmp);
        if(k>=jac->nnz){ // more memory needed, double current size
          splm_ccsm_realloc_novalues(jac, nobs, nvars, jac->nnz<<1); // 2*jac->nnz
          rowidx=jac->rowidx; // re-init
        }
        rowidx[k++]=i;
      }
    }
  }
  colptr[nvars]=k;
  splm_ccsm_realloc_novalues(jac, nobs, nvars, k); // adjust to actual size...

  free(hxx);
  return k;
}

/* ptrs to functions computing a sparse CRS or CCS Jacobian */
struct splm_fjacfuncs{
  void (*ccsfjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata);
  void (*crsfjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata);
};

/* contains information necessary for computing a finite difference approximation to a Jacobian,
 * e.g. function to differentiate, problem dimensions and pointers to working memory buffers
 */
struct fdj_data_ {
  double delta;     /* parameter increment for computing Jacobian */
  double *hx;       /* func evaluated at p, i.e. hx=func(p), nobsx1 */
  double *hxx;      /* work array for evaluating func(p+delta), nobsx1 */
  void (*func)(double *p, double *hx, int nvars, int nobs, void *adata); /* function to differentiate */
  int forw;         /* 1: forward differences, 0: central differences */
  int nfeval;       /* number of func evaluations for computing Jacobian */
  void *adata;
};

/* finite difference approximation to the Jacobian of func
 * using either forward or central differences.
 * The structure of the Jacobian is assumed known.
 * Uses the strategy described in Nocedal-Wright, ch. 7, pp. 169.
 */
static void splm_fdif_jac(
    double *p,              /* I: current parameter estimate, nvarsx1 */
    struct splm_ccsm *jac,  /* O: CCS array storing approximate Jacobian, nobsxnvars */
    int nvars,
    int nobs,
    void *dat)
{
register int i, j, jj, k;
register double d;
int ii, m, *jcol, *varlist, *coldone, forw;
int *vidxs, *ridxs;
double *tmpd;
struct fdj_data_ *fdjd=(struct fdj_data_ *)dat;
void (*func)(double *p, double *hx, int m, int n, void *adata);
double *hx, *hxx, delta;

  /* retrieve problem-specific information passed in *dat */
  func=fdjd->func;
  hx=fdjd->hx;
  hxx=fdjd->hxx;
  forw=fdjd->forw;
  delta=fdjd->delta;

  //(*func)(p, hx, nvars, nobs, fdjd->adata); ++(fdjd->nfeval);

  jcol=(int *)emalloc(nobs*sizeof(int)); /* keeps track of measurements influenced by the set of variables currently in "varlist" below */
  for(i=0; i<nobs; ++i) jcol[i]=-1;

  k=splm_ccsm_col_maxnelms(jac);
  vidxs=(int *)emalloc(2*k*sizeof(int));
  ridxs=vidxs+k;

  varlist=(int *)emalloc(nvars*sizeof(int)); /* stores indices of J's columns which are computed with the same "func" call */
  coldone=(int *)emalloc(nvars*sizeof(int)); /* keeps track of J's columns which have been already computed */
  memset(coldone, 0, nvars*sizeof(int)); /* initialize to zero */

  tmpd=(double *)emalloc(nvars*sizeof(double));

  for(j=0; j<nvars; ++j){
    double scl;

    if(coldone[j]) continue; /* column j already computed */

    //for(i=0; i<nobs; ++i) jcol[i]=-1;
    k=splm_ccsm_col_elmidxs(jac, j, vidxs, ridxs);
    for(i=0; i<k; ++i) jcol[ridxs[i]]=j;
    varlist[0]=j; m=1; coldone[j]=1;

    for(jj=j+1; jj<nvars; ++jj){
      if(coldone[jj]) continue; /* column jj already computed */

      k=splm_ccsm_col_elmidxs(jac, jj, vidxs, ridxs);
      for(i=0; i<k; ++i)
        if(jcol[ridxs[i]]!=-1) goto nextjj;

      if(k==0) { coldone[jj]=1; continue; } /* all zeros column, ignore */

      /* column jj does not clash with previously considered ones, mark it */
      for(i=0; i<k; ++i) jcol[ridxs[i]]=jj;
      varlist[m++]=jj; coldone[jj]=1;

nextjj:
      continue;
    }

//printf("Jacobian for colunn %d computed with %d other columns\n", j, m-1);

    for(k=0; k<m; ++k){
      /* determine d=max(SPLM_DELTA_SCALE*|p[varlist[k]]|, delta), see HZ */
      d=SPLM_DELTA_SCALE*p[varlist[k]]; // force evaluation
      d=FABS(d);
      if(d<delta) d=delta;

      tmpd[varlist[k]]=d;
      p[varlist[k]]+=d;
    }

    (*func)(p, hxx, nvars, nobs, fdjd->adata); ++(fdjd->nfeval); // hxx=f(p+d)

    if(forw){
      for(k=0; k<m; ++k)
        p[varlist[k]]-=tmpd[varlist[k]]; /* restore */

      scl=1.0;
    }
    else{ // central
      for(k=0; k<m; ++k)
        p[varlist[k]]-=2*tmpd[varlist[k]];

      (*func)(p, hx, nvars, nobs, fdjd->adata); ++(fdjd->nfeval); // hx=f(p-d)

      for(k=0; k<m; ++k)
        p[varlist[k]]+=tmpd[varlist[k]]; /* restore */

      scl=0.5; // 1./2.
    }

    for(k=0; k<m; ++k){
      d=tmpd[varlist[k]];
      d=scl/d; /* invert so that divisions can be carried out faster as multiplications */

      jj=splm_ccsm_col_elmidxs(jac, varlist[k], vidxs, ridxs);
      for(i=0; i<jj; ++i){
        ii=ridxs[i];
        jac->val[vidxs[i]]=(hxx[ii]-hx[ii])*d;
        jcol[ii]=-1; /* restore */
      }
    }
  }

  free(tmpd);
  free(coldone);
  free(varlist);
  free(vidxs);
  free(jcol);

//printf("Function evaluations for Jacobian computation: %d\n", fdjd->nfeval);
}


/* Compute e=x-y for two n-vectors x and y and return the squared L2 norm of e.
 * e can coincide with either x or y. 
 * Uses loop unrolling and blocking to reduce bookkeeping overhead & pipeline
 * stalls and increase instruction-level parallelism; see http://www.abarnett.demon.co.uk/tutorial.html
 */
static double nrmL2xmy(double *e, const double *x, const double *y, const int n)
{
const int blocksize=8, bpwr=3; /* 8=2^3 */
register int i;
int blockn;
register double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0;

  /* n may not be divisible by blocksize, 
   * go as near as we can first, then tidy up.
   */
  blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */

  /* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
  for(i=blockn; i>0; e+=blocksize, x+=blocksize, y+=blocksize, i-=blocksize){
    e[0]=x[0]-y[0]; sum0+=e[0]*e[0];
    e[1]=x[1]-y[1]; sum1+=e[1]*e[1];
    e[2]=x[2]-y[2]; sum2+=e[2]*e[2];
    e[3]=x[3]-y[3]; sum3+=e[3]*e[3];
    e[4]=x[4]-y[4]; sum0+=e[4]*e[4];
    e[5]=x[5]-y[5]; sum1+=e[5]*e[5];
    e[6]=x[6]-y[6]; sum2+=e[6]*e[6];
    e[7]=x[7]-y[7]; sum3+=e[7]*e[7];
  }

  /*
   * There may be some left to do.
   * This could be done as a simple for() loop, 
   * but a switch is faster (and more interesting) 
   */

  i=blockn;
  if(i<n){ 
  /* Jump into the case at the place that will allow
   * us to finish off the appropriate number of items. 
   */
    switch(n - i){
      case 7 : e[6]=x[6]-y[6]; sum0+=e[6]*e[6];
      case 6 : e[5]=x[5]-y[5]; sum1+=e[5]*e[5];
      case 5 : e[4]=x[4]-y[4]; sum2+=e[4]*e[4];
      case 4 : e[3]=x[3]-y[3]; sum3+=e[3]*e[3];
      case 3 : e[2]=x[2]-y[2]; sum0+=e[2]*e[2];
      case 2 : e[1]=x[1]-y[1]; sum1+=e[1]*e[1];
      case 1 : e[0]=x[0]-y[0]; sum2+=e[0]*e[0];
    }
  }

  return (sum0+sum1)+(sum2+sum3);
}

/* As above for a zero x, i.e. compute e=-y and return its squared L2 norm.
 * Note that arg x is unused below!
 */
static double nrmL2my(double *e, const double *x, const double *y, const int n)
{
const int blocksize=8, bpwr=3; /* 8=2^3 */
register int i;
int blockn;
register double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0;

  /* n may not be divisible by blocksize, 
   * go as near as we can first, then tidy up.
   */
  blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */

  /* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
  for(i=blockn; i>0; e+=blocksize, y+=blocksize, i-=blocksize){
    e[0]=-y[0]; sum0+=e[0]*e[0];
    e[1]=-y[1]; sum1+=e[1]*e[1];
    e[2]=-y[2]; sum2+=e[2]*e[2];
    e[3]=-y[3]; sum3+=e[3]*e[3];
    e[4]=-y[4]; sum0+=e[4]*e[4];
    e[5]=-y[5]; sum1+=e[5]*e[5];
    e[6]=-y[6]; sum2+=e[6]*e[6];
    e[7]=-y[7]; sum3+=e[7]*e[7];
  }

  /*
   * There may be some left to do.
   * This could be done as a simple for() loop, 
   * but a switch is faster (and more interesting) 
   */

  i=blockn;
  if(i<n){ 
  /* Jump into the case at the place that will allow
   * us to finish off the appropriate number of items. 
   */
    switch(n - i){
      case 7 : e[6]=-y[6]; sum0+=e[6]*e[6];
      case 6 : e[5]=-y[5]; sum1+=e[5]*e[5];
      case 5 : e[4]=-y[4]; sum2+=e[4]*e[4];
      case 4 : e[3]=-y[3]; sum3+=e[3]*e[3];
      case 3 : e[2]=-y[2]; sum0+=e[2]*e[2];
      case 2 : e[1]=-y[1]; sum1+=e[1]*e[1];
      case 1 : e[0]=-y[0]; sum2+=e[0]*e[0];
    }
  }

  return (sum0+sum1)+(sum2+sum3);
}

/* main sparse Levenberg-Marquardt function; accepts either a CRS or a CCS Jacobian */
/* 
 * This is the main sparse Levenberg-Marquardt function.
 * It seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function accepts either a) an analytic Jacobian in CRS or CCS format or
 * b) the nonzero pattern of the Jacobian in CRS or CCS and then approximates its
 * elements via finite differencing.
 *
 * Returns the number of iterations (>=0) if successful, SPLM_ERROR if failed
 */

static int splm_core(
    void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                                              /* functional relation describing measurements. Given a parameter vector p,
                                               * computes a prediction of the measurements \hat{x}. p is nvarsx1,
                                               * \hat{x} is nobsx1, maximum
                                               */
    struct splm_fjacfuncs *const jacfs,       /* contains pointer to appropriate function for evaluating the sparse
                                               * Jacobian dX/dp or initialize its structure with its nonzero pattern,
                                               * see jactyp below. The computed Jacobian is in one of CCS/CRS formats
                                               */
    int jactyp,                               /* if equal to SPLM_ANJAC, jacfs->crsfjac/jacfs->ccsfjac computes
                                               * analytically the sparse Jacobian dX/da in CRS/CCS format.
                                               *
                                               * if equal to SPLM_ZPJAC, jacfs->crsfjac/jacfs->ccsfjac only
                                               * initializes the Jacobian structure according to its nonzero pattern;
                                               * the Jacobian is then approximated with repetitive func calls and
                                               * finite differences.
                                               * This can be computationally inefficient and is thus NOT recommended.
                                               *
                                               * if equal to SPLM_NOJAC, splm_core() attempts to detect the non-zero
                                               * pattern of the Jacobian. The correctness of this procedure cannot be
                                               * guaranteed and this option should be used with caution.
                                               */
    double *const p,    /* I/O: initial parameter estimates. On output has the estimated solution. size nvars */
    double *const x,    /* measurement vector. size nobs. NULL implies a zero vector */
    const int nvars,    /* I: parameter vector dimension (i.e. #unknowns) [m] */
    const int nconvars, /* I: number of parameters (starting from the 1st) whose values should not be modified. >=0 */
    const int nobs,     /* I: measurement vector dimension [n] */
    int Jnnz,           /* I: number of nonzeros for the Jacobian J */
    int JtJnnz,         /* I: number of nonzeros for the product J^t*J, -1 if unknown */
    const int itmax,    /* I: maximum number of iterations */
    double opts[SPLM_OPTS_SZ],
	                      /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, delta, spsolver]. Respectively the scale factor
                         * for initial \mu, stopping thresholds for ||J^T e||_inf, ||dp||_2 and ||e||_2, the step used in difference
                         * approximation to the Jacobian and the sparse direct solver to employ. Set to NULL for defaults to be used.
                         * If \delta<0, the Jacobian is approximated with central differences which are more accurate (but more
                         * expensive to compute!) compared to the forward differences employed by default.
                         */
    double info[SPLM_INFO_SZ],
	                     /* O: information regarding the minimization. Set to NULL if don't care
                        * info[0]=||e||_2 at initial p.
                        * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                        * info[5]= # iterations,
                        * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                        *                                 2 - stopped by small dp
                        *                                 3 - stopped by itmax
                        *                                 4 - singular matrix. Restart from current p with increased mu 
                        *                                 5 - too many failed attempts to increase damping. Restart with increased mu
                        *                                 6 - stopped by small ||e||_2
                        *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. User error
                        * info[7]= # function evaluations
                        * info[8]= # Jacobian evaluations
			                  * info[9]= # linear systems solved, i.e. # attempts	for reducing error
                        */
    void *adata       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 
)
{
register int i, j, ii=0, jj; /* -Wall */
const int verbose=0;

/* variables related to the Jacobian; jac_ccs & jac_crs are *not* used simultaneously */
struct splm_ccsm jac_ccs; /* CCS structure for the Jacobian J */
struct splm_crsm jac_crs; /* CRS structure for the Jacobian J */
const int jacisccs=jacfs->ccsfjac || jactyp==SPLM_ZPJAC || jactyp==SPLM_NOJAC, /* is J in CCS? note that finite differentiation always uses CCS */
          jacisan=jactyp==SPLM_ANJAC; /* is J analytic? */
struct splm_ccsm jacTjac; /* CCS structure for the product J^T J */

/* The following are work arrays that are dynamically allocated by splm_core() */
double *e;    /* work array for storing the errors, size nobsx1 */
double *jacTe; /* work array for storing J^T e, size nvars */

double *dp;   /* work array for storing the parameter vector updates, size nvars */

/* The following variables are needed by the LM algorithm */
register int itno;  /* iteration counter */
int issolved;
/* temporary work arrays that are dynamically allocated */
double *hx,         /* \hat{x}_i, max. size m*n*mnp */
       **diagHess,  /* pointers to Hessian diagonal, size nvars */
       *pdp;        /* p + dp, size nvars */

double mu;                /* damping constant */
register double muincr,   /* damping increment */
                tmp;      /* mainly used in matrix & vector multiplications */
double p_eL2, jacTe_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
double p_L2, dp_L2=DBL_MAX, dF, dL;
double tau, eps1, eps2, eps2_sq, eps3, delta;
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0, nlss=0;
double (*const L2xmy)(double *e, const double *x, const double *y, const int n)=(x!=NULL)? nrmL2xmy : nrmL2my;
struct fdj_data_ fdj_data;
void *jac_adata;

/* direct sparse solver */
int (*solver)(struct splm_ccsm *A, double *B, void **state, int what, double *x)=NULL;
void *solvstate=NULL;
int spsolver, hess_howto=2;
//double timea, timeb;

#ifdef USE_ATA_QUADS
int nquads=0, quadsz[4];
void *quads[4]={NULL, NULL, NULL, NULL};
#endif

  /* error checking */
  if(nobs<nvars){
    fprintf(stderr, "splm_core(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
    return SPLM_ERROR;
  }

  if(nconvars<0){
    fprintf(stderr, "splm_core(): number of variables to be kept constant must be non negative ([%d])\n", nconvars);
    return SPLM_ERROR;
  }

  if(!jacfs->crsfjac && !jacfs->ccsfjac){
    if(jacisan){
      fprintf(stderr, "splm_core(): a function computing the Jacobian should be supplied!\n");
      return SPLM_ERROR;
    }
    else
      fprintf(stderr, "\n*** warning: requested attempt to discover the Jacobian's nonzero pattern might fail quietly in splm_core()!\n");
  }
  else if(Jnnz<=0){ /* a function computing either the Jacobian or its nonzero pattern has been supplied */
    fprintf(stderr, "splm_core(): invalid Jacobian size speficied, number of nonzero elements must be positive! [%d]\n", Jnnz);
    return SPLM_ERROR;
  }

  /* Initialization */
  mu=jacTe_inf=tmp=0.0; /* -Wall */
  if(opts){
    tau=FABS(opts[0]);
    eps1=FABS(opts[1]);
    eps2=FABS(opts[2]);
    eps2_sq=opts[2]*opts[2];
    eps3=FABS(opts[3]);
    delta=opts[4];
    spsolver=(int)opts[5];
  }
  else{ // use default values
    tau=SPLM_INIT_MU;
    eps1=SPLM_STOP_THRESH;
    eps2=SPLM_STOP_THRESH;
    eps2_sq=SPLM_STOP_THRESH*SPLM_STOP_THRESH;
    eps3=SPLM_STOP_THRESH;
    delta=SPLM_DIFF_DELTA;
    spsolver=SPLM_CSPARSE;
    // spsolver=SPLM_CHOLMOD;
  }

  /* Note that when the nonzero pattern of an approximate Jacobian is supplied in
   * CRS, the finite differences are computed with the aid of splm_fdif_jac() that
   * accepts a CCS Jacobian!
   */

  /* allocate work arrays */
  e=(double *)emalloc(nobs*sizeof(double));
  jacTe=(double *)emalloc(nvars*sizeof(double));
  dp=(double *)emalloc(nvars*sizeof(double));

  hx=(double *)emalloc(nobs*sizeof(double));
  diagHess=(double **)emalloc(nvars*sizeof(double *));
  pdp=(double *)emalloc(nvars*sizeof(double));

  if(jacisccs){
    if(jactyp!=SPLM_NOJAC)
      splm_ccsm_alloc(&jac_ccs, nobs, nvars, Jnnz);
    else{
      splm_ccsm_alloc_novalues(&jac_ccs, nobs, nvars, Jnnz>0? Jnnz : 128); // no better estimate of Jnnz yet...
      Jnnz=splm_jacpatguess(func, p, &jac_ccs, nvars, nobs, adata, hx, delta);
      splm_ccsm_alloc_values(&jac_ccs);
      JtJnnz=-1; // do not trust any estimates regarding the size of JtJ
    }

    jac_crs.rowptr=jac_crs.colidx=NULL;
    jac_crs.val=NULL;
    jac_crs.nr=jac_crs.nc=0;
  }
  else{
    splm_crsm_alloc(&jac_crs, nobs, nvars, Jnnz);

    jac_ccs.colptr=jac_ccs.rowidx=NULL;
    jac_ccs.val=NULL;
    jac_ccs.nr=jac_ccs.nc=0;
  }

  /* 
   * hess_howto 1: compute just the lower triangular & diagonal part of the J^T*J matrix product
   * hess_howto > 1: compute the full matrix product
   */
  /* choose sparse solver */
  switch(spsolver){
#ifdef HAVE_CHOLMOD
    case SPLM_CHOLMOD:
      solver=splm_Axb_CHOLMOD; hess_howto=1;
      break;
#endif

#ifdef HAVE_CSPARSE
    case SPLM_CSPARSE:
      solver=splm_Axb_CSPARSE; hess_howto=2;
      break;
#endif

#ifdef HAVE_LDL
    case SPLM_LDL:
      solver=splm_Axb_LDLP; hess_howto=2; // LDL actually uses the upper triangular part!
      break;
#endif

#ifdef HAVE_UMFPACK
    case SPLM_UMFPACK:
      solver=splm_Axb_UMFPACK; hess_howto=2;
      break;
#endif

#if 0
#ifdef HAVE_MA77
    case SPLM_MA77:
      solver=splm_Axb_MA77; hess_howto=1;
      break;
#endif
#endif

#ifdef HAVE_MA57
    case SPLM_MA57:
      solver=splm_Axb_MA57; hess_howto=1;
      break;
#endif

#ifdef HAVE_MA47
    case SPLM_MA47:
      solver=splm_Axb_MA47; hess_howto=1;
      break;
#endif

#ifdef HAVE_MA27
    case SPLM_MA27:
      solver=splm_Axb_MA27; hess_howto=1;
      break;
#endif

#ifdef HAVE_PARDISO
    case SPLM_PARDISO:
      solver=splm_Axb_PARDISO; hess_howto=1;
      break;
#endif

#ifdef HAVE_MKL
    case SPLM_DSS:
      solver=splm_Axb_DSS; hess_howto=1;
      break;
#endif

#ifdef HAVE_SuperLU
    case SPLM_SuperLU:
      solver=splm_Axb_SuperLU; hess_howto=2;
      break;
#endif

#ifdef HAVE_TAUCS
    case SPLM_TAUCS:
      solver=splm_Axb_TAUCS; hess_howto=1;
      break;
#endif

#ifdef HAVE_SPOOLES
    case SPLM_SPOOLES:
      solver=splm_Axb_SPOOLES; hess_howto=1;
      break;
#endif

#ifdef HAVE_MUMPS
    case SPLM_MUMPS:
      solver=splm_Axb_MUMPS; hess_howto=1;
      break;
#endif

#if 0
    case SPLM_sparseQR:
      solver=splm_Axb_sparseQR; hess_howto=2;
      break;
#endif

    default: {
      char *solvname[]={
        "**unused**",
        "CHOLMOD", "CSPARSE", "LDL", "UMFPACK",
        "MA77", "MA57", "MA47", "MA27",
        "PARDISO", "DSS", "SuperLU", "TAUCS", "SPOOLES", "MUMPS",
      };

      if(spsolver>0 && spsolver<sizeof(solvname)/sizeof(char *))
        fprintf(stderr, "\n*** unsupported sparse direct solver \"%s\" specified to splm_core(), exiting!\n\n", solvname[spsolver]);
      else
        fprintf(stderr, "\n*** unknown sparse direct solver \"%d\" specified to splm_core(), exiting!\n\n", spsolver);
      exit(1);
    }
  }

  /* hess_howto should not be modified beyond this point */

  if(JtJnnz>0){
    if(hess_howto==1) // do not allocate space for the strictly upper part
      splm_ccsm_alloc(&jacTjac, nvars-nconvars, nvars-nconvars, ((JtJnnz-(nvars-nconvars))>>1) + nvars-nconvars);
    else
      splm_ccsm_alloc(&jacTjac, nvars-nconvars, nvars-nconvars, JtJnnz);
  }
  else
    jacTjac.nnz=-1; // #nonzeros in J^t*J is unknown, request its computation

  if(!jacisan){
    fprintf(stderr, "\n*** warning: no analytic Jacobian function supplied to splm_core(), performance will suffer!\n\n");

    fdj_data.nfeval=0;
    /* init non zero structure in jac_ccs */
    if(jacfs->ccsfjac)
      (*(jacfs->ccsfjac))(p, &jac_ccs, nvars, nobs, adata);
    else if(jacfs->crsfjac){ /* user's nonzero pattern Jacobian function expects a CRS matrix */
      struct splm_crsm tmp;

      /* splm_crsm_alloc_novalues() allocates space just for the structure
       * of the Jacobian and not for the unused values
       */
      splm_crsm_alloc_novalues(&tmp, nobs, nvars, Jnnz);
      (*(jacfs->crsfjac))(p, &tmp, nvars, nobs, adata);
      splm_crsm2ccsm(&tmp, &jac_ccs);
      splm_crsm_free(&tmp);
      jacfs->crsfjac=NULL; // not really necessary...
    }
    else /* both jacfs->ccsfjac, jacfs->crsfjac are NULL, J pattern jas been detected above */
      fdj_data.nfeval=nvars+1; /* pattern discovery involves n+1 calls */

    fdj_data.func=func;
    if(delta>0.0)
      fdj_data.forw=1; /* use forward differencing */
    else{
      delta=-delta; /* make positive */
      fdj_data.forw=0; /* use central differencing */
    }
    fdj_data.delta=delta;
    fdj_data.hx=hx;
    fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
    fdj_data.adata=adata;
    jacfs->ccsfjac=splm_fdif_jac;
    jac_adata=(void *)&fdj_data;
  }
  else{
    fdj_data.func=NULL;
    fdj_data.hxx=fdj_data.hx=NULL;
    jac_adata=adata;
  }

  /* compute the error vectors e_ij in hx */
  (*func)(p, hx, nvars, nobs, adata); nfev=1;
  /* ### compute e=x - f(p) and its L2 norm */
  p_eL2=(*L2xmy)(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

  if(verbose) printf("sparseLM: initial error %g\n", p_eL2);
  init_p_eL2=p_eL2;
  if(!SPLM_FINITE(p_eL2)) stop=7;

  for(itno=stop=0; itno<itmax && !stop; ++itno){
    /* Note that p, e and ||e||_2 have been updated at the previous iteration */

    if(p_eL2<=eps3){ /* error is small */
      stop=6;
      break;
    }

    if(jacisccs){
      /* compute the Jacobian at p in CCS */
      (*(jacfs->ccsfjac))(p, &jac_ccs, nvars, nobs, jac_adata); ++njev;

      if(nconvars>0) ii=splm_ccsm_drop_cols(&jac_ccs, nconvars); // get rid of the constant parameters in the Jacobian
      /* compute J^T J and J^T e. */
      if(itno==0){
        splm_setup_AtA_ccsA(&jac_ccs, &jacTjac, hess_howto); // compute product's structure
#ifdef USE_ATA_QUADS
        nquads=splm_get_AtA_ccsA_quads(&jac_ccs, &jacTjac, hess_howto, quads, quadsz);
#endif
      }
//timea=splm_gettime();
#ifdef USE_ATA_QUADS
      splm_calc_AtA_ccsA_quads(&jac_ccs, &jacTjac, hess_howto, quads, nquads);
#else
      splm_calc_AtA_ccsA(&jac_ccs, &jacTjac, hess_howto);
#endif /* USE_ATA_QUADS */
      splm_calc_Atx_ccsA(&jac_ccs, e, jacTe+nconvars); // Note: all of e is needed, not just the e+nconvars part!
//timeb=splm_gettime();
//printf("%d: %.2f\n", itno, timeb-timea);
      if(nconvars>0) splm_ccsm_restore_cols(&jac_ccs, nconvars, ii); // undo drop

      if(verbose && itno==0){
        printf("Jacobian and approximate Hessian nonzeros: %d  %d\n", jac_ccs.colptr[nvars], jacTjac.colptr[nvars-nconvars]);
        printf("Densities: %.3g%%  %.3g%%\n", 100*((double)jac_ccs.nnz)/(jac_ccs.nr*jac_ccs.nc),
                                              100*((double)jacTjac.nnz)/(jacTjac.nr*jacTjac.nc));
      }
    }
    else{
      /* compute the Jacobian at p in CRS */
      (*(jacfs->crsfjac))(p, &jac_crs, nvars, nobs, jac_adata); ++njev;

      if(itno==0){
        splm_setup_AtA_crsA(&jac_crs, nconvars, &jacTjac, hess_howto); // compute product's structure
#ifdef USE_ATA_QUADS
        nquads=splm_get_AtA_crsA_quads(&jac_crs, nconvars, &jacTjac, hess_howto, quads, quadsz);
#endif
      }
      /* compute J^T J and J^T e. */
//timea=splm_gettime();
#ifdef USE_ATA_QUADS
      splm_calc_AtA_crsA_quads(&jac_crs, nconvars, &jacTjac, hess_howto, quads, nquads);
      splm_calc_Atx_crsA(&jac_crs, nconvars, e, jacTe+nconvars);
#else
      splm_calc_AtAx_crsA(&jac_crs, nconvars, &jacTjac, hess_howto, e, jacTe+nconvars);
#endif /* USE_ATA_QUADS */
//timeb=splm_gettime();
//printf("%d: %.2f\n", itno, timeb-timea);

      if(verbose && itno==0){
        printf("Jacobian and approximate Hessian nonzeros: %d  %d\n", jac_crs.rowptr[nobs], jacTjac.colptr[nvars-nconvars]);
        printf("Densities: %.3g%%  %.3g%%\n", 100*((double)jac_crs.nnz)/(jac_crs.nr*jac_crs.nc),
                                              100*((double)jacTjac.nnz)/(jacTjac.nr*jacTjac.nc));
      }
    }

#if 0
    {
    const int nreg=6;
    int regbnds[6][4];
    double cols[6][3]={
      {1, 0.17, 0},
      {0, 0.83, 0},
      {.08, .08, .92},
      {.08, .08, .92},
      //{1, .83, 0},
      {0.17, 0.83, .83},
      {0.17, 0.83, .83},
    };
    int nc, np, cnp, pnp, knp;

    nc=11-1; np=305; cnp=6; pnp=3; knp=5; // bt
    //nc=2; np=114; cnp=12; pnp=3; knp=0; // tensor

    regbnds[0][0]=0;    regbnds[0][1]=0;    regbnds[0][2]=nc*cnp;       regbnds[0][3]=nc*cnp;
    regbnds[1][0]=nc*cnp; regbnds[1][1]=nc*cnp; regbnds[1][2]=nc*cnp+np*pnp; regbnds[1][3]=nc*cnp+np*pnp;
    regbnds[2][0]=0;    regbnds[2][1]=nc*cnp; regbnds[2][2]=nc*cnp;       regbnds[2][3]=nc*cnp+np*pnp;
    regbnds[3][0]=nc*cnp; regbnds[3][1]=0;    regbnds[3][2]=nc*cnp+np*pnp; regbnds[3][3]=nc*cnp;

    regbnds[4][0]=0; regbnds[4][1]=nc*cnp+np*pnp;    regbnds[4][2]=nc*cnp+np*pnp+knp; regbnds[4][3]=nc*cnp+np*pnp+knp;
    regbnds[5][0]=nc*cnp+np*pnp; regbnds[5][1]=0;    regbnds[5][2]=nc*cnp+np*pnp+knp; regbnds[5][3]=nc*cnp+np*pnp;

    if(jacisccs)
      splm_ccsm2eps(&jac_ccs, NULL, NULL, 0, "/tmp/jac.eps");
    else
      splm_crsm2eps(&jac_crs, NULL, NULL, 0, "/tmp/jac.eps");

    splm_crcsm2eps(0, (void *)&jacTjac, regbnds, cols, nreg, "/tmp/jacTjac.eps");
    }
#endif

    /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=nconvars, p_L2=jacTe_inf=0.0; i<nvars; ++i){
      if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;
      p_L2+=p[i]*p[i];
    }

#if 0
if(!(itno%10)){
  printf("Current estimate: ");
  for(i=0; i<nvars; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g\n", jacTe_inf, p_eL2);
}
#endif

    /* check for convergence */
    if((jacTe_inf <= eps1)){
      dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* set up pointers to the diagonal of the Hessian & compute initial damping factor */
    if(itno==0){
      for(j=0; j<nconvars; ++j){ // just in case...
        diagHess[j]=NULL;
        jacTe[j]=0.0;
      }

      for(j=0; j<nvars-nconvars; ++j){ // CHECKME nconvars!
        ii=jacTjac.colptr[j];
        jj=jacTjac.colptr[j+1]-1;
        if(ii>jj){
          fprintf(stderr, "splm_core(): Hessian's column %d is empty, the corresponding Jacobian column is also empty?\n", j+nconvars);
          exit(1);
        }

        /* binary search for finding the element at row j */
        while(ii<=jj){
          if(j<jacTjac.rowidx[ii] || j>jacTjac.rowidx[jj]) break; // not found

          i=((unsigned int)ii + (unsigned int)jj) >> 1; //(ii+jj)/2;
          if(j<jacTjac.rowidx[i]){
            jj=i-1;
            continue;
          }
          else if(j>jacTjac.rowidx[i]){
            ii=i+1;
            continue;
          }
          else{
            diagHess[j+nconvars]=jacTjac.val+i; // add nconvars CHECKME
            goto nextj;
          }

          /* should never get here! */
          fprintf(stderr, "splm_core() internal error: could not find the location corresponding to var %d in the Hessian's diagonal\n", j+nconvars);
          exit(1);
        }
nextj:  continue;
      }

      for(i=nconvars, tmp=DBL_MIN; i<nvars; ++i)
        if(*(diagHess[i])>tmp) tmp=*(diagHess[i]); /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */
    /* NOTE: during the following loop, e might contain the error corresponding to *unaccepted* dampings! */
    muincr=mu;
    while(1){
      /* solve the linear system (J^T J + mu I) dp = J^t e to compute dp. */

      /* augment Hessian */
      for(i=nconvars; i<nvars; ++i)
        *(diagHess[i])+=muincr;

      /* symbolic processing for the first time only; assume identical nonzero pattern for subsequent invocations */
      issolved=(*solver)(&jacTjac, jacTe+nconvars, &solvstate, (nlss==0)? 1 : 2, dp+nconvars); ++nlss;

      if(issolved){
        /* parameter vector updates are now in dp */

        for(i=nconvars; i-->0;  ) dp[i]=0.0; /* no change for the first nconvars params */

        /* compute p's new estimate and ||dp||^2 */
        for(i=0, dp_L2=0.0; i<nvars; ++i){
          pdp[i]=p[i] + (tmp=dp[i]);
          dp_L2+=tmp*tmp;
        }

        if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(dp_L2>=(p_L2+eps2)/SPLM_EPSILON_SQ){ /* almost singular */
         stop=4;
         break;
       }

        (*func)(pdp, hx, nvars, nobs, adata); ++nfev; /* evaluate function at p + dp */
        /* ### compute ||e(pdp)||_2 */
        pdp_eL2=(*L2xmy)(e, x, hx, nobs); /* e=x-hx, pdp_eL2=||e|| */
        if(!SPLM_FINITE(pdp_eL2)){
          if(verbose) /* identify the offending prediction */
            for(i=0; i<nobs; ++i)
              if(!SPLM_FINITE(hx[i]))
                printf("sparseLM: component %d of the predicted measurement vector is invalid!\n", i);

          stop=7;
          break;
        }

        for(i=0, dL=0.0; i<nvars; ++i)
          dL+=dp[i]*(mu*dp[i]+jacTe[i]);

        dF=p_eL2-pdp_eL2;

        if(verbose>1)
          printf("\ndamping term %12g, gain ratio %10g, errors %10g %10g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2, pdp_eL2);

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(2.0*dF/dL-1.0);
          tmp=1.0-tmp*tmp*tmp;
          mu=mu*( (tmp>=SPLM_ONE_THIRD)? tmp : SPLM_ONE_THIRD );
          nu=2;

          for(i=0; i<nvars; ++i) /* update p's estimate */
            p[i]=pdp[i];

          p_eL2=pdp_eL2; /* update ||e||_2 */
          break;
        }
      } /* issolved */

      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

#if 0
      /* restore diagonal entries */
      for(i=nconvars; i<nvars; ++i)
        *(diagHess[i])-=mu;
#endif
      muincr=-mu;
      mu*=nu;
      muincr+=mu; // muincr:=new_mu - old_mu
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown) */
        /* fprintf(stderr, "Too many failed attempts to increase the damping factor in splm_core()! Singular Hessian matrix?\n"); */
        stop=5;
        break;
      }
      nu=nu2;

    } /* inner loop */
  }

  if(itno>=itmax) stop=3;

  /* restore diagonal entries */
  if(itno) /* itno==0 implies that the LM outter loop has not been executed at all */
    for(i=nconvars; i<nvars; ++i)
      *(diagHess[i])-=mu;

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=dp_L2;

    if(itno)
      for(i=nconvars, tmp=DBL_MIN; i<nvars; ++i)
        if(tmp<*(diagHess[i])) tmp=*(diagHess[i]);

    info[4]=mu/tmp;
    info[5]=itno;
    info[6]=stop;
    info[7]=nfev + (jacisan? 0 : fdj_data.nfeval);
    info[8]=njev;
    info[9]=nlss;
  }
                                                               
  if(verbose){
    double density;

    density=jacfs->ccsfjac? ((double)jac_ccs.nnz/(jac_ccs.nr*jac_ccs.nc)) : ((double)jac_crs.nnz/(jac_crs.nr*jac_crs.nc));
    fflush(stdout);
    fprintf(stdout, "sparseLM using %d parameters (%d constant) and %d observations\n", nvars, nconvars, nobs);
    fprintf(stdout, "%s %s Jacobian, density %.1f%%\n\n", jacisan? "analytic" : "approximate", jacfs->ccsfjac? "CCS" : "CRS", density*100.0);
    fprintf(stdout, "sparseLM returning after %g iter, reason %g, error %g [initial %g], %d/%d func/fjac evals, %d lin. systems\n", info[5], info[6], info[1], info[0], (int)info[7], (int)info[8], (int)info[9]);
    fflush(stdout);
  }

   /* free whatever was allocated */
  if(jacisccs) splm_ccsm_free(&jac_ccs);
  else splm_crsm_free(&jac_crs);
  if(jacTjac.nnz>0) splm_ccsm_free(&jacTjac);
  free(e); free(jacTe);  
  free(dp);               

  free(hx); free(diagHess); free(pdp);
  if(fdj_data.hxx){
    free(fdj_data.hxx);
  }

  /* free solver's state */
  if(itno) (*solver)(NULL, NULL, &solvstate, 3, NULL);

#ifdef USE_ATA_QUADS
  if(quads[0]){
    free(quads[0]); free(quads[1]);
    free(quads[2]); free(quads[3]);
  }
#endif

  return (stop!=4 && stop!=7)?  itno : SPLM_ERROR;
}

/* wrappers around splm_core() depending on the format of the Jacobian.
 * For explanations of their arguments, please consult the comments referring
 * to the homonymous arguments of splm_core() above
 */

int sparselm_dercrs(void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                    void (*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
                    double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
                    int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
struct splm_fjacfuncs fjf;

  fjf.crsfjac=fjac;
  fjf.ccsfjac=NULL;

  return splm_core(func, &fjf, SPLM_ANJAC, p, x, nvars, nconvars, nobs,
                        Jnnz, JtJnnz, itmax, opts, info, adata);
}

int sparselm_derccs(void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                    void (*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
                    double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
                    int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
struct splm_fjacfuncs fjf;

  fjf.ccsfjac=fjac;
  fjf.crsfjac=NULL;

  return splm_core(func, &fjf, SPLM_ANJAC, p, x, nvars, nconvars, nobs,
                        Jnnz, JtJnnz, itmax, opts, info, adata);
}

int sparselm_difcrs(void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                    void (*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
                    double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
                    int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
struct splm_fjacfuncs fjf;

  fjf.crsfjac=fjac;
  fjf.ccsfjac=NULL;

  return splm_core(func, &fjf, fjac? SPLM_ZPJAC : SPLM_NOJAC, p, x, nvars, nconvars, nobs,
                        Jnnz, JtJnnz, itmax, opts, info, adata);
}

int sparselm_difccs(void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
                    void (*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
                    double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
                    int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
struct splm_fjacfuncs fjf;

  fjf.ccsfjac=fjac;
  fjf.crsfjac=NULL;

  return splm_core(func, &fjf, fjac? SPLM_ZPJAC : SPLM_NOJAC, p, x, nvars, nconvars, nobs,
                        Jnnz, JtJnnz, itmax, opts, info, adata);
}


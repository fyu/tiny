/*
/////////////////////////////////////////////////////////////////////////////////////////////
//// 
////  Prototypes and definitions for the sparse Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2005-2011  Manolis Lourakis (lourakis at ics.forth.gr)
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
///////////////////////////////////////////////////////////////////////////////////////////////
*/

#ifndef _SPLM_H_
#define _SPLM_H_

#ifdef __cplusplus
extern "C" {
#endif

/* sparse direct solvers */
#define SPLM_CHOLMOD        1
#define SPLM_CSPARSE        2
#define SPLM_LDL            3
#define SPLM_UMFPACK        4
#define SPLM_MA77           5
#define SPLM_MA57           6
#define SPLM_MA47           7
#define SPLM_MA27           8
#define SPLM_PARDISO        9
#define SPLM_DSS            10
#define SPLM_SuperLU        11
#define SPLM_TAUCS          12
#define SPLM_SPOOLES        13
#define SPLM_MUMPS          14
/*#define SPLM_sparseQR     -1 */

#define SPLM_DIFF_DELTA     1E-06 /* finite differentiation minimum delta */
#define SPLM_DELTA_SCALE    1E-04 /* finite differentiation delta scale   */

#define SPLM_OPTS_SZ        6 /* max(5, 6) */
#define SPLM_INFO_SZ        10
#define SPLM_ERROR          -1
#define SPLM_INIT_MU        1E-03
#define SPLM_STOP_THRESH    1E-12
#define SPLM_VERSION        "1.3 (December 2011)"

#define SPLM_ANJAC          0 /* analytic Jacobian */
#define SPLM_ZPJAC          1 /* approximate Jacobian, only zero pattern supplied */ 
#define SPLM_NOJAC          2 /* approximate Jacobian, zero pattern to be guessed, use cautiously */ 


/* Sparse matrix representation using Compressed Row Storage (CRS) format.
 * See http://www.netlib.org/linalg/html_templates/node91.html
 */

struct splm_crsm{
    int nr, nc;   /* #rows, #cols for the sparse matrix */
    int nnz;      /* number of nonzero array elements */
    double *val;  /* storage for nonzero array elements. size: nnz */
    int *colidx;  /* column indexes of nonzero elements. size: nnz */
    int *rowptr;  /* locations in val that start a row. size: nr+1.
                   * By convention, rowptr[nr]=nnz
                   */
};

/* Sparse matrix representation using Compressed Column Storage (CCS) format.
 * See http://www.netlib.org/linalg/html_templates/node92.html
 */

struct splm_ccsm{
    int nr, nc;   /* #rows, #cols for the sparse matrix */
    int nnz;      /* number of nonzero array elements */
    double *val;  /* storage for nonzero array elements. size: nnz */
    int *rowidx;  /* row indexes of nonzero elements. size: nnz */
    int *colptr;  /* locations in val that start a column. size: nc+1.
                   * By convention, colptr[nc]=nnz
                   */
};


/* sparse LM for functions with CRS and CCS Jacobians */

extern int sparselm_dercrs(
       void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
       void (*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
       double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
       int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

extern int sparselm_derccs(
       void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
       void (*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
       double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
       int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

extern int sparselm_difcrs(
       void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
       void (*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
       double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
       int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

extern int sparselm_difccs(
       void (*func)(double *p, double *hx, int nvars, int nobs, void *adata),
       void (*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
       double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
       int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata);

/* error checking for CRS and CCS Jacobians */
extern void sparselm_chkjaccrs(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        void (*jacf)(double *p, struct splm_crsm *jac, int m, int n, void *adata),
        double *p, int m, int n, int jnnz, void *adata, double *err);

extern void sparselm_chkjacccs(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        void (*jacf)(double *p, struct splm_ccsm *jac, int m, int n, void *adata),
        double *p, int m, int n, int jnnz, void *adata, double *err);

/* CRS sparse matrices manipulation routines */
extern void splm_crsm_alloc(struct splm_crsm *sm, int nr, int nc, int nnz);
extern void splm_crsm_alloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz);
extern void splm_crsm_alloc_values(struct splm_crsm *sm);
extern void splm_crsm_realloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz);
extern void splm_crsm_free(struct splm_crsm *sm);
extern int splm_crsm_elmidx(struct splm_crsm *sm, int i, int j);
extern int splm_crsm_elmrow(struct splm_crsm *sm, int idx);
extern int splm_crsm_row_elmidxs(struct splm_crsm *sm, int i, int *vidxs, int *jidxs);
extern int splm_crsm_row_maxnelms(struct splm_crsm *sm);
extern int splm_crsm_col_elmidxs(struct splm_crsm *sm, int j, int *vidxs, int *iidxs);
extern void splm_crsm2ccsm(struct splm_crsm *crs, struct splm_ccsm *ccs);
extern void splm_crsm_row_sort(struct splm_crsm *sm);

/* CCS sparse matrices manipulation routines */
extern void splm_ccsm_alloc(struct splm_ccsm *sm, int nr, int nc, int nnz);
extern void splm_ccsm_alloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz);
extern void splm_ccsm_alloc_values(struct splm_ccsm *sm);
extern void splm_ccsm_realloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz);
extern void splm_ccsm_free(struct splm_ccsm *sm);
extern int splm_ccsm_elmidx(struct splm_ccsm *sm, int i, int j);
extern int splm_crsm_elmcol(struct splm_ccsm *sm, int idx);
extern int splm_ccsm_row_elmidxs(struct splm_ccsm *sm, int i, int *vidxs, int *jidxs);
extern int splm_ccsm_col_elmidxs(struct splm_ccsm *sm, int j, int *vidxs, int *iidxs);
extern int splm_ccsm_col_maxnelms(struct splm_ccsm *sm);
extern void splm_ccsm2crsm(struct splm_ccsm *ccs, struct splm_crsm *crs);
extern int splm_ccsm_drop_cols(struct splm_ccsm *A, int ncols);
extern void splm_ccsm_restore_cols(struct splm_ccsm *A, int ncols, int ncnnz);
extern void splm_ccsm_col_sort(struct splm_ccsm *sm);

extern double splm_gettime(void);


/* Sparse matrix representation using Sparse Triplet (ST) format.
 * Note that the matrix might have an allocated size (maxnnz) larger
 * than the number of elements it actually holds (nnz). 
 * Primarily intended to be used for setting up the structure
 * of a corresponding CCS matrix.
 * See http://people.sc.fsu.edu/~jburkardt/data/st/st.html
 */

struct splm_stm{
  int nr, nc;   /* #rows, #cols for the sparse matrix */
  int nnz;      /* number of nonzero array elements */
  int maxnnz;   /* maximum number of nonzero array elements the array can hold */
  int *rowidx;  /* row indexes of nonzero elements. size: nnz */
  int *colidx;  /* column indexes of nonzero elements. size: nnz */
  double *val;  /* nonzero elements, NULL when only structure is stored */
};

extern void splm_stm_alloc(struct splm_stm *sm, int nr, int nc, int maxnnz);
extern void splm_stm_allocval(struct splm_stm *sm, int nr, int nc, int maxnnz);
extern void splm_stm_free(struct splm_stm *sm);
extern int splm_stm_nonzero(struct splm_stm *sm, int i, int j);
extern int splm_stm_nonzeroval(struct splm_stm *sm, int i, int j, double val);
extern void splm_stm2ccsm(struct splm_stm *st, struct splm_ccsm *ccs);
extern void splm_tri2ccsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_ccsm *ccs);
extern void splm_stm2crsm(struct splm_stm *st, struct splm_crsm *crs);
extern void splm_tri2crsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_crsm *crs);


#ifdef __cplusplus
}
#endif

#endif /* _SPLM_H_ */

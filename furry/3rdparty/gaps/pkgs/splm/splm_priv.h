/*
/////////////////////////////////////////////////////////////////////////////////////////////
//// 
////  Private protos and defs for the sparse Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2005-2008  Manolis Lourakis (lourakis at ics.forth.gr)
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

#ifndef __SPLM_PRIV_H_
#define __SPLM_PRIV_H_

#ifdef __cplusplus
extern "C" {
#endif

#define emalloc(sz)         splm_emalloc(__FILE__, __LINE__, sz)
#define erealloc(p, sz)     splm_realloc(__FILE__, __LINE__, p, sz)
extern void *splm_emalloc(char *file, int line, size_t sz);
extern void *splm_realloc(char *file, int line, void *oldptr, size_t sz);


/* A^t*A & A^t*x computation */
extern void splm_setup_AtA_ccsA(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job);
extern void splm_calc_AtA_ccsA(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job);
extern void splm_calc_Atx_ccsA(struct splm_ccsm *A, double *const x, double *const y);
extern int splm_get_AtA_ccsA_quads(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job, void *quads[4], int quadsz[4]);
extern void splm_calc_AtA_ccsA_quads(struct splm_ccsm *A, struct splm_ccsm *AtA, const int job, void *quads[4], int nquads);

extern void splm_setup_AtA_crsA(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job);
extern void splm_calc_AtAx_crsA(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, double *const x, double *const y);
extern int splm_get_AtA_crsA_quads(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, void *quads[4], int quadsz[4]);
extern void splm_calc_AtA_crsA_quads(struct splm_crsm *A, int zerocols, struct splm_ccsm *AtA, const int job, void *quads[4], int nquads);
extern void splm_calc_Atx_crsA(struct splm_crsm *A, int zerocols, double *x, double *y);


/* LDL-based solution of linear systems */
extern int splm_Axb_LDLP(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* HSL MA47-based solution of linear systems */
extern int splm_Axb_MA47(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* HSL MA27-based solution of linear systems */
extern int splm_Axb_MA27(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* PARDISO-based solution of linear systems */
extern int splm_Axb_PARDISO(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* SuperLU-based solution of linear systems */
extern int splm_Axb_SuperLU(struct splm_ccsm *inA, double *inB, void **state, int what, double *x);

/* SparseQR-based solution of linear systems */
extern int splm_Axb_sparseQR(struct splm_ccsm *inA, double *inB, void **state, int what, double *x);

/* TAUCS-based solution of linear systems */
extern int splm_Axb_TAUCS(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* CSparse sparse Cholesky solver */
extern int splm_Axb_CSPARSE(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* HSL MA57-based solution of linear systems */
extern int splm_Axb_MA57(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* CHOLMOD-based solution of linear systems */
extern int splm_Axb_CHOLMOD(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* UMFPACK-based solution of linear systems */
extern int splm_Axb_UMFPACK(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* SPOOLES-based solution of linear systems */
extern int splm_Axb_SPOOLES(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* MUMPS-based solution of linear systems */
extern int splm_Axb_MUMPS(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* DSS-based solution of linear systems */
extern int splm_Axb_DSS(struct splm_ccsm *A, double *B, void **state, int what, double *x);

/* HSL MA77-based solution of linear systems */
extern int splm_Axb_MA77(struct splm_ccsm *A, double *B, void **state, int what, double *x);

#ifdef __cplusplus
}
#endif

#endif /* __SPLM_PRIV_H_ */

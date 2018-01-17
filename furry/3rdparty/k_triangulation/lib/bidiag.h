#ifndef	__BIDIAG_H__
#define	__BIDIAG_H__

bool	gsl_linalg_bidiag_decomp	(Matrix&	A,
					 Vector&	tau_U,
					 Vector&	tau_V);

bool	gsl_linalg_bidiag_unpack2 	(Matrix&	A, 
					 Vector&	tau_U, 
					 Vector&	tau_V,
					 Matrix&	V);

#endif	/* __BIFDIAG_H__ */

#ifndef	__HOUSEHOLDER_H__
#define	__HOUSEHOLDER_H__

double	gsl_linalg_householder_transform	(Vector&	v);
bool	gsl_linalg_householder_hm 		(double		tau,
						 const Vector&	v,
						 Matrix&	A);
bool	gsl_linalg_householder_mh 		(double		tau,
						 const Vector&	v,
						 Matrix&	A);
bool	gsl_linalg_householder_hv 		(double		tau,
						 const Vector&	v,
						 Vector&	w);
bool	gsl_linalg_householder_hm1 		(double		tau,
						 Matrix&	A);

#endif	/* __HOUSEHOLDER_H__ */

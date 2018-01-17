#ifndef	__SVDSTEP_H__
#define	__SVDSTEP_H__

void	chop_small_elements 	(const Vector&	d,
				 Vector&	f);
double	trailing_eigenvalue 	(const Vector&	d,
				 Vector&	f);
void	create_schur 		(double 	d0,
				 double 	f0,
				 double 	d1,
				 double 	*c,
				 double 	*s);
void	svd2 			(Vector&	d,
				 Vector&	f,
				 Matrix&	U,
				 Matrix&	V);
void	chase_out_intermediate_zero (Vector&	d,
				     Vector&	f,
				     Matrix&	U,
				     int	k0);
void	chase_out_trailing_zero	(Vector&	d,
				 Vector&	f,
				 Matrix&	V);
void	qrstep 			(Vector&	d,
				 Vector& 	f,
				 Matrix& 	U,
				 Matrix& 	V);
#endif	/* __SVDSTEP_H__ */

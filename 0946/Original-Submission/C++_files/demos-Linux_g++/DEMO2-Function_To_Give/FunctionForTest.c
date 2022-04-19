/*
	===========================================================================
	PURPOSE
	===========================================================================

	Here a Laplace Transform
	fzE
	and its Inverse function
	gzE
	are provided to test RELIADIFF function on them, by mean of the 
	TEST_FunctionToGive.c program.
	
	*******
	RELIADIFF routine computes an approximate value of the Inverse Laplace 
	Transform by means of its Laguerre polynomial expansion.
	The method which this software is based on uses the Automatic 
	Differentiation to compute the coefficients of the series expansion.
	The summation of successive terms is truncated when the truncation 
	error is less or equal than a tolerance, or when the provided maximum 
	number of terms is reached.
	The inverse function is computed as:

	f(x)=exp((sigma0-b)*x [ c0 L0(2bx)+c1 L1(2bx)+...+ cn Ln(2bx)]
  
	where
	    - SIGMA0 (>=0) is (an estimate of) the convergence
	      abscissa of the Laplace Transform. 
	    - b is a parameter. 
	    - Lk (0<=k<=n) is the k-th degree Laguerre Polynomial.  
   
	===========================================================================
	ARGUMENTS
	===========================================================================
	
	fzE requires in input:
	- z: 		(TADIFF) double precision: the evaluation point of the Transform
	
	gzE requires in input: 
	- t: 		double precision: the evaluation point of the Inverse function

	===========================================================================
	RETURN VALUE
	===========================================================================

	fzE provides in output:
	- (TADIFF) double precision: the evaluation of the Transform in z
	
	gzE requires in input: 
	- double precision: the evaluation of the Inverse function in t 

	===========================================================================
	AUTHORS
	===========================================================================
	
		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI
   	
	===========================================================================	

*/
#include "../utility/Util.h"

/*==================================================================*/
/*WRITE HERE THE LAPLACE TRASFORM FOR THE TEST*/
/*eventually using the GSL and the provided utility functions (see DEMOS_userguide.pdf)*/

/*Laplace Transform*/
/* F(z)=1/((z+a)^5)  a=3/5 */
T<double> fzE(T<double> z) {
   int n=5;      
   return 1./pow_int_T((z+3./5.),n);            
}

/*==================================================================*/
/*WRITE HERE THE RELATED INVERSE LAPLACE TRASFORM TO COMPARE WITH THE RELIADIFF RESULT*/
/*eventually using the GSL and the provided utility functions (see DEMOS_userguide.pdf)*/

/*Related Inverse function*/
/* f(t) = t^(n-1)/((n-1)!exp(a*t))    a=3/5 , n=5*/
double gzE(double t){
    double a=3./5.;
    return gsl_sf_pow_int(t,4)/(gsl_sf_fact(4)*exp(a*t));
}
/*==================================================================*/

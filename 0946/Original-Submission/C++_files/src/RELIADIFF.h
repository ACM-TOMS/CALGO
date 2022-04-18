/*
	==================================================================
	PURPOSE
	==================================================================

	Header file for RELIADIFF software.
   
	==================================================================
	AUTHORS
	==================================================================
	
		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI
   	
	==================================================================	

*/


#ifndef RELIADIFF_H
	#define RELIADIFF_H
	#define MaxLength 5000
	#include "fadbad/tadiff.h"	
	#include <math.h>
	#include <float.h>
	#include <string.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <iostream>
	
	using namespace std;
	using namespace fadbad;	
	
	/*
	==================================================================                                              
	PURPOSE
	==================================================================                                               

	Calculate the Phi function 
	
	Phi (z) = (2*b) /(1-z) * F( (2*b)/(1-z)  + sigma – b  )
	
	where F is a Laplace Transform

	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	z  :		(TADIFF) double precision
				it contains value at which the Phi function
				is required. 
					       
	b       	double precision 
				it contains a parameter of the Phi function
	
	sigma   	double precision
				it contains a parameter of the Phi function that should 
				be greater than the convergence abscissa of F function
	
	sigma0  	double precision    	
				it contains the convergence abscissa of F function
				IT SHOULD BE GREATER OR EQUAL THAN ZERO
	
	
	fz :        (TADIFF) double precision function pointer
				it contains the name of the Laplace Transform function F.

	szero:		integer
				it can contain a parameter to say if the F function has
				a singularity at zero.
				IT CAN BE 0 (there isn't) OR GREATER (there is).
	
	==================================================================
	RETURN VALUE
	==================================================================

	(TADIFF) DOUBLE PRECISION	numerical value of Phi function at the point z  
									 
	==================================================================
	*/
	T<double>  phi (T<double> , double ,double, double, T<double> (*)(T<double>), int);

	
	/* 
	===========================================================================
	PURPOSE
	===========================================================================
	
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
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	x  :	double precision
		it contains value at which the Inverse Laplace Transform 
		is required. Each component of x has to have a value  
		greater or  equal to zero.
					       
	fz :    (TADIFF) double precision  function pointer
		it contains the name of the Laplace Transform function.
		
	==================================================================
	INPUT/OUTPUT PARAMETERS
	==================================================================

	szero:	string
		on entry it can contain 
		-	a parameter, interpreted as integer, to say if the 
			transform has a singularity at zero.
			IT CAN BE 0 (there isn't) OR GREATER (there is).
		or
		-	the character "n", to mean that the user doesn't 
			want to provide it.
			It will be posed to the default value (szero=0).
		
		on exit it will contain the value actually used by the 
			software in computation, to interpret as integer.
	
	pcoeff:	string
		on entry it can contain 
		-	parameter, interpreted as integer, to print or not 
			all the used coefficients in a file at the end of 
			work. IT CAN BE 0 (don't print) OR GREATER (print).
		or
		-	the character "n", to mean that the user doesn't 
			want to provide it.
			It will be posed to the default value (pcoeff=0).
		
		on exit it will contain the value actually used by the 
			software in computation, to interpret as integer.
				
	sigma0: string
		on entry it can contain 
		-	the abscissa of convergence of f: it will be 
			interpreted as double precision.
		-	If it is less than zero, it will be posed to zero.
		or
		-	the character "n", to mean that the user doesn't 
			want to provide it.
			It will be posed to the default value (sigma0=1).
		
		on exit it will contain the value actually used by the 
			software in computation, to interpret as double.
				
	tol :   string
		on entry it can contain 
		-	the required accuracy on f, interpreted as double 
			precision.
			+	if it is less or equal to 0 it will be posed 
				to the default value (tol=10^(-3)).
			+	if it is greater than 1 it will be posed to 
				the default value (tol=10^(-3)).
			+	if it is less than the machine precision 
				(single precision), it will be posed to the 
				value of the machine precision.
		or
		-	the character "n", to mean that the user doesn't 
			want to provide it.
			It will be posed to the default value (tol=10^(-3)).
		
		on exit it will contain the value actually used by the 
			software in computation, to interpret as double.
				
	nmax:	string
		on entry it can contain 
		-	the required maximum number of Laguerre series terms, 
			interpreted as integer.
			+	if nmax < 8, it is posed at a default value 
				(nmax = 2000)
			+	if nmax > MaxLength=5000, it is posed at 
				nmax=MaxLength.
		or
		-	the character "n", to mean that the user doesn't 
			want to provide it.
			It will be posed to the default value (nmax = 2000).
		
		on exit it will contain the value actually used by the 
			software in computation, to interpret as integer.
				
	==================================================================
	OUTPUT PARAMETERS
	==================================================================
			  
	ILf:		double precision 
			It will contain the computed value of f at x.

	absesterr: 	double precision
			It will contain the estimate of the absolute error on f.
                             			
	relesterr:  	double precision 
			It will contain the estimate of the relative error on f.

	flag:       	integer  (DIAGNOSTIC) 
			It will contain an information on the obtained result 
			accuracy.
					   			   
	nopt:       	integer (DIAGNOSTIC) 
			It will contain the found optimal number of Laguerre 
			series terms.    

	ncalc		integer (DIAGNOSTIC) 
			It will contain the total number of terms calculated by 
			the software to find the optimal one.

	==================================================================
	RETURN VALUE
	==================================================================

	INTEGER: 	It will contain a diagnostic on the input data.
			if it is equal to 1, it means that the routine ended 
			without any output because the input x was less than 0.

	==================================================================
	*/
	int RELIADIFF( double, T<double> (*)(T<double>), char *, char *, char *, char *, char *, 
					double *, double *, double *, int *, int *, int *);
        
#endif


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

	(TADIFF) DOUBLE PRECISION	numerical value of Phi function the 
								point z  
									 
	==================================================================
	AUTHORS
	==================================================================

		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI

	==================================================================
		
*/

#include "RELIADIFF.h"

T<double> phi(T<double> z,double b, double sigma,double sigma0,T<double> (*fz)(T<double>), int szero){
	T<double> y1,x,y,pt;  
	double shift=-1.;
	x=(2*b)/(1-z);							
	y=x+sigma-b;	
	if(szero!=0)
		pt=y+sigma0+shift; /*shift in case F has a singularity at zero*/
	else
		pt=y+sigma0;
	y1=x*(fz(pt)); 
	return y1;
}




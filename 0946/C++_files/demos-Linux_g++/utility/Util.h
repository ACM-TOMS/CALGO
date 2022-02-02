/*
	==================================================================
	PURPOSE
	==================================================================

	Header file for auxiliraty routines, useful for the 
	building of test drivers and demos for RELIADIFF.   
	
	==================================================================
	AUTHORS
	==================================================================
	
		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI
   	
	==================================================================	

*/


#ifndef UTIL_H 
	#define UTIL_H 
	#include "src/RELIADIFF.h"
	#include <sys/time.h>
	#include <string.h>
	
	#include <gsl/gsl_errno.h>
	#include "gsl/gsl_math.h"
	#include "gsl/gsl_sf_pow_int.h"
	#include "gsl/gsl_sf_dawson.h"
	#include "gsl/gsl_sf_bessel.h"
	#include "gsl/gsl_sf_expint.h"
	#include "gsl/gsl_sf_gamma.h"
	#include "gsl/gsl_sf_erf.h"
	#include "gsl/gsl_integration.h"
	#include <gsl/gsl_roots.h>
	
	
	/*AUXILIARY ROUTINES*/

	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine computes the integer power of a T<double> value.
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	a  :	(TADIFF) double precision
			it contains the base
					       
	n :     integer
			it contains the exponent
				
	==================================================================
	RETURN VALUE
	==================================================================

	(TADIFF) double precision 	
	It will contain the n-power of a
	
	==================================================================
	*/
	T<double> pow_int_T( T<double>, int );

	
	
	
	
	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine computes the Laguerre Polynomial in x of degree K 
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	K  :	integer: degree of the Laguerre polynomial to calculate 
					       
	x :     double precision: point of evaluation
				
	==================================================================
	RETURN VALUE
	==================================================================

	double precision: Laguerre polynomial in x
	
	==================================================================
	*/
	double Laguerre(int , double );
	
	
	
	
	
	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine computes the Hermite Polynomial in x of degree K 
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	K  :	integer: degree of the Hermite polynomial to calculate 
					       
	x :     double precision: point of evaluation
				
	==================================================================
	RETURN VALUE
	==================================================================

	double precision: Hermite polynomial in x
	
	==================================================================
	*/
	double Hermite(int , double );

	
	
	
	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine prints in file the flag (diagnostics parameter) 
	interpretation for RELIADIFF routine
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	fp:		file handler: handler of the file where to write the 
			flag meaning in RELIADIFF
					      				
	==================================================================
	*/
	void print_flags_file(FILE *);
	
	
	
	
	
	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine prints in file the ncalc and nopt (diagnostics 
	parameters) interpretation for RELIADIFF routine
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	fp:		file handler: handler of the file where to write Nopt/Ncalc 
			meaning in RELIADIFF
			
	nmax:	integer: required maximum number of Taylor coefficients to 
			calculate in RELIADIFF
					      				
	==================================================================
	*/
	void print_N_file(FILE *, int );
	
	
	
	/*
	==================================================================
	PURPOSE
	==================================================================

	This routine prints in file a warning about the interpretation for 
	RELIADIFF output
	
	==================================================================
	INPUT PARAMETERS
	==================================================================
		
	fp:		file handler: handler of the file where to write Nopt/Ncalc 
			meaning in RELIADIFF
			
	==================================================================
	*/
	void print_warning(FILE *);
        
#endif

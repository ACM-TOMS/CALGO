/*
	==================================================================
	PURPOSE
	==================================================================

	In this file there are some auxiliraty routines, useful for the 
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


#include "Util.h"

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
T<double> pow_int_T(T<double> a, int n){
	int i=0;
	T<double> pot=a;
	if (n==0) return 1;
	for (i=2;i<=n;i++) pot=pot*a;
	return pot;
}



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
double Laguerre(int K, double x){
	int k=0;
	double temp1=0, temp2=0, LK=0;
	double L0=1;
	double L1=1-x;
	for(k=2;k<=K;k++){
		temp1=(-x+2*((double) k)-1)/(double)(k);
		temp2=(((double) k)-1)/(double)(k); 
		LK=temp1*L1-temp2*L0;
		L0=L1;
		L1=LK;
	}	
	return LK;
}


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
double Hermite(int K, double x){
	int k=0;
	double temp1=0, temp2=0, HK=0;
	double H0=1;
	double H1=2*x;
	for(k=2;k<=K;k++){
		temp1=2*x*H1;
		temp2=2*(k-1)*H0; 
		HK=temp1-temp2;
		H0=H1;
		H1=HK;
	}	
	return HK;
}


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
void print_flags_file(FILE *fp){

	fprintf(fp,"\n******************************************************************************************************\n");
	fprintf(fp,"\n                 	                       ERROR ESTIMATE DIAGNOSTIC\n");
	fprintf(fp,"\na. flag =0: corresponds to the case of output:\n");
	fprintf(fp,"\t- absesterr < tol\n");
	fprintf(fp,"\t- relesterr < tol\n");
	fprintf(fp,"\tBoth the absolute and the relative error estimates are smaller than tolerance,\n");
	fprintf(fp,"\tso the software fully satisfies the required accuracy.\n ");
	fprintf(fp,"\tThis means that the software can obtain more accurate result\n ");
	fprintf(fp,"\tif user requires a smaller tolerance.\n ");	
	fprintf(fp,"b. flag =1: corresponds to the case of output:\n");
	fprintf(fp,"\t- absesterr is the minimum obtained value, but greater than tolerance.\n"); 
	fprintf(fp,"\t- relesterr is the minimum obtained value, but greater than tolerance.\n");
	fprintf(fp,"\tThis means that within nmax terms of the series expansion, the algorithm cannot satisfy\n");
	fprintf(fp,"\tthe required accuracy. So, it provides the numerical result within the maximum\n");
	fprintf(fp,"\tattainable accuracy with no more than nmax terms, and nopt will be the  number\n");
	fprintf(fp,"\tof terms  at which such minimum is reached (eventually different from ncalc,\n");
	fprintf(fp,"\tthat is the total number of calculated coefficients).\n"); 
	fprintf(fp,"\tMoreover, this means that the series seems to converge too slowly or to diverge.\n");
	fprintf(fp,"\tSo, user can try to obtain a better result tuning some of the optional parameters:\n");
	fprintf(fp,"\tsinf, sigma0, nmax.\n"); 
	fprintf(fp,"\tUser is also invited  to verify if the Transform satisfies algorithm's requirements.\n");	
	fprintf(fp,"c. flag =2: corresponds to the case of output:\n");
	fprintf(fp,"\t- absesterr < tol.\n");
	fprintf(fp,"\t- relesterr is the minimum obtained value, but greater than tolerance.\n");
	fprintf(fp,"\tOnly the absolute error estimate is smaller than the user's required accuracy.\n");
	fprintf(fp,"\tThis means that within nmax terms of the series expansion, the algorithm cannot\n");
	fprintf(fp,"\tsatisfy the required accuracy. So, it provides the numerical result within the\n");
	fprintf(fp,"\tmaximum attainable accuracy with no more than nmax terms, and nopt will be the\n");
	fprintf(fp,"\tnumber of terms  at which such minimum is reached (eventually different from ncalc,\n");
	fprintf(fp,"\tthat is the total number of calculated coefficients).\n"); 
	fprintf(fp,"\tMoreover, this means that the inverse function f rapidly decreases towards zero.\n");
	fprintf(fp,"\tUser can try to obtain a better result tuning some of the optional parameters:\n");
	fprintf(fp,"\tsinf, sigma0, nmax.\n");
	fprintf(fp,"d. flag =3: corresponds to the case of output:\n");
	fprintf(fp,"\t- absesterr is the minimum obtained value, but greater than tolerance.\n");
	fprintf(fp,"\t- relesterr < tol.\n");
	fprintf(fp,"\tOnly the relative error estimate is smaller than the user's required accuracy.\n");
	fprintf(fp,"\tThis means that within nmax terms of the series expansion, the algorithm can\n");
	fprintf(fp,"\tsatisfy the required accuracy, but not for the relative error. This means also\n");
	fprintf(fp,"\tthat the inverse function f increases rapidly.\n");	 
	fprintf(fp,"\n*******************************************************************************************************\n");
	fprintf(fp,"\n\n\n");
	return;
}


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
void print_N_file(FILE *fp, int nmax){
	fprintf(fp,"\n************************************************************************************************************\n");
	fprintf(fp,"\n                 	                       nopt and ncalc\n");
	fprintf(fp,"\nncalc is the maximum number of terms of the Laguerre series expansion calculated by the algorithm.\n");
	fprintf(fp,"nopt is the number of terms of the Laguerre series expansion that gives the numerical result within\n");
	fprintf(fp,"the maximum attainable accuracy, less or equal to nmax (the required maximum number of terms).\n");
	fprintf(fp,"You can find one of three cases:\n");
	fprintf(fp,"a. nopt=ncalc<nmax\n");
	fprintf(fp,"\tThe computed  value of the inverse Laplace function agrees with the true one within log(tol)\n");
	fprintf(fp,"\tsignificant and decimal digits. This value is obtained calculating nopt terms of the Laguerre\n");
	fprintf(fp,"\tseries expansion. It should correspond to flag = 0 or flag = 3.\n");
	fprintf(fp,"b. nopt<ncalc=nmax\n");
	fprintf(fp,"\tWithin nmax terms of the Laguerre series expansion, the algorithm cannot satisfy the user's\n");
	fprintf(fp,"\trequired accuracy, and the series seems to diverge. The algorithm provides a numerical result\n");
	fprintf(fp,"\twithin  the maximum attainable accuracy, and nopt is the  number of terms at which such maximum\n");
	fprintf(fp,"\tis reached. It should correspond to flag = 1 or flag = 2.\n");
	fprintf(fp,"c. nopt=ncalc=nmax\n");
	fprintf(fp,"\tThis occurs if:\n");
	fprintf(fp,"\t\tWithin nmax terms of the Laguerre series expansion, the algorithm cannot satisfy the user's\n");
	fprintf(fp,"\t\trequired accuracy, but the series could converge, even if very slowly, or diverge: the\n");
	fprintf(fp,"\t\talgorithm provides numerical result with the maximum attainable accuracy reached within nmax\n");
	fprintf(fp,"\t\tterms of the Laguerre series expansion. If the series diverges, nopt accidentally corresponds to nmax.\n");
	fprintf(fp,"\t\tIt should correspond to flag = 1 or flag = 2.\n");
	fprintf(fp,"\tor if:\n");
	fprintf(fp,"\t\tThe algorithm can satisfy the user's required accuracy, within exactly nmax terms of the Laguerre\n");
	fprintf(fp,"\t\texpansion, so the series seems to converge, even if quite slowly: the algorithm provides numerical\n");
	fprintf(fp,"\t\tresult within the required accuracy, reached within nmax terms of the Laguerre series expansion.\n");
	fprintf(fp,"\t\tIt should be flag = 0 or flag = 3. \n");
	fprintf(fp,"\n**************************************************************************************************************\n");
	fprintf(fp,"\n\n\n");
	return;
}


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

void print_warning(FILE *fp){
	fprintf(fp,"\n************************************************************************************************************\n");
	fprintf(fp,"                 	                             WARNING\n");
	fprintf(fp,"-------------------Calculated errors are just errors estimate, they are not errors bounds-------------------\n");
	fprintf(fp,"************************************************************************************************************\n");
	return;
}
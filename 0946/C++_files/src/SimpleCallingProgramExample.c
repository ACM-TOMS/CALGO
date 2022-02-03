/*
	==================================================================
	RELIADIFF Calling Example
	==================================================================
 
	The program computes the numerical inversion of a Laplace transform.
	The Transform is

	F(z)=1/(z^4-a^4)  with a=3/5

	at point x=2.5
    
	================================================================== 
	Description of RELIADIFF
	==================================================================
 
	RELIADIFF computes an approximate value of the Inverse Laplace 
	Transform by means of its Laguerre polynomial expansion, using the 
	Automatic Differentiation to computes the Taylor coefficients ot 
	the series.

	The inverse function is computed as:
	f(x)=exp((sigma0-b)*x [ c0 L0(2bx)+c1 L1(2bx)+...+ cn Ln(2bx)]
	where
	- SIGMA0 (>=0) is (an estimate of) the convergence
	abscissa of the Laplace Transform. 
	- b is a parameter. 
	- Lk (0<=k<=n) is the k-th degree Laguerre Polynomial.   
   
	==================================================================
	Example Program Results.
	==================================================================
	
	Used Tolerance on accuracy: 1.000000e-03;
	Used abscissa of convergence on F: 1.000000e+00;
	Used maximum number of Laguerre series coefficients: 2000;
	Singularity in zero:  no;
	Used Taylor coefficients printed:  no;

				    TABLE
	----------------------------------------------------------------------------
	|     x    |     f_comp      | estabserr | estrelerr |  Nopt |  Ncal |FLAG|
	---------------------------------------------------------------------------
	| 2.50e+00 | 2.61987232e+00 | 3.210e-04 | 1.225e-04 |    12 |    12 |  0 |

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

T<double> fzTest(T<double> z) {
	return 1./(pow(z,4)-pow(3./5.,4));                    
}


int main(int argc,char **argv){
/*RELIADIFF INPUT*/
	double		x=2.5;			//Inverse Function evaluation point(s)
	T<double>	(*fz)(T<double>);   	//Function F Pointer
	
/*RELIADIFF OPTIONAL INPUT: we don't give any*/
	char        	szero[10]="n"; 		//optional: there's or not a singularity at zero for F
	char		pcoeff[10]="n";		//optional: if user doesn't say yes, the driver doesn't print coefficients
	char		sigma0[10]="n";    	//optional: abscissa of convergence of F 
	char 		tol[10]="n";		//optional: tolerance on accuracy
	char		nmax[10]="n"; 		//optional: maximum number of series coefficients
	
/*RELIADIFF OUTPUT*/	
	double		absesterr;	//absolute error estimate
	double		relesterr;	//relative error estimate	
	double 		ILf;		//Inverse Function f computed
	int		flag;		//diagnostics on the result
	int		nopt;		//diagnostics on the software work
	int		ncalc;		//diagnostics on the software work
	int		ierr;		//diagnostics on the software work	
	
/*AUXILIARY VARIABLES*/
	char	name[20]="";
	
	fz = fzTest; 
							
	/************************************** CALLING RELIADIFF *********************************************/
	ierr=RELIADIFF(x,fz,szero,pcoeff,sigma0,tol,nmax,&ILf,&absesterr,&relesterr,&flag,&nopt,&ncalc); 	
	/******************************************************************************************************/
			   
	switch(ierr){
	case 1:
		printf("\nx=%f - RELIADIFF stopped: x<0.0!\n It must be x>=0\n\n",x);
		break;
	default: 		
		printf("Used Tolerance on accuracy: %e;\n",atof(tol));
		printf("Used abscissa of convergence on F: %e;\n",atof(sigma0));
		printf("Used maximum number of Laguerre series coefficients: %d;\n",atoi(nmax));
		printf("Singularity in zero: ");
		if(atoi(szero)) printf(" yes;\n");
		else printf(" no;\n");
		printf("Used Taylor coefficients printed: ");
		if(atoi(pcoeff)){
			sprintf(name,"coeff_x%.1f.txt",x);
			printf(" yes, they are in file %s;\n",name);
		}
		else printf(" no;\n");
		printf("\n                             TABLE\n");
		printf("----------------------------------------------------------------------------\n");
		printf("|     x    |     f_comp     | estabserr | estrelerr |  Nopt |  Ncal |FLAG|\n");
		printf("----------------------------------------------------------------------------\n");			
		printf("| %4.2e | %14.8e | %9.3e | %9.3e | %5d | %5d | %2d |\n", x, ILf, absesterr, relesterr, nopt, ncalc, flag);
		break;
	}/*endswitch*/
	printf("\n************************************************************************************************************\n");
	printf("                 	                             WARNING\n");
	printf("-------------------Calculated errors are just errors estimate, they are not errors bounds-------------------\n");
	printf("************************************************************************************************************\n");
	printf("\n\nProgram termineted. Please press a key to exit.");
	getchar();
	return 0;
}
/*END OF MAIN*/


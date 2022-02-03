/*
	==================================================================
	PURPOSE
	==================================================================

	This routine computes an approximate value of the Inverse Laplace 
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
	PROCEDURE NEEDED
	==================================================================                   
				 
	T<double> phi(T<double> z,double b, double sigma, double sigma0,
				T<double> (*fz)(T<double>), int szero)

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

int RELIADIFF( double x, T<double> (*fz)(T<double>),  //input
	       char *szeroc, char *pcoeffc,char *sigma0c, char *tolc, char *nmaxc,   //input-output
	       double *ILf, double *absesterr, double *relesterr, int *flag, int *nopt, int *ncalc) //output

{
    /* ...Local Variables specifications... */
	int 	k,n,nstart=4,i,pcoeff=0,szero=0,nmax=2000;
	double  b,expsx,esterr,esterr_new,sigma;
	double  fx,Ntol;
	double  esp,L1,L2,LK,y;
	double  temp1,temp2;
	double	sigma0=1.;
	double	defTol=1.0e-03;
	double	tol=defTol;
	char	defau[5]="n";
	double	shift=-1;
	
	/*... Output File Pointer and Name .....................*/
	FILE	*fp;
	char	name[20]="";
	
	T<double> y1,z=0; /*for the evaluation of the phi function*/
        
    /*********************************************************************************\
	 INPUT DATA CHECK 
    \********************************************************************************/
      
	if( x<0.0) return 1; 
	
	if (strcmp(tolc, defau)!=0){
		tol=atof(tolc);	
		if (tol>1) tol=defTol;
		else if (tol<FLT_EPSILON) tol=FLT_EPSILON;
	}
	if (strcmp(nmaxc, defau)!=0){
		nmax=atoi(nmaxc);
		if(nmax<nstart) nmax=2000;
		else if(nmax>MaxLength) nmax=MaxLength;
	}
	if (strcmp(pcoeffc, defau)!=0){
		pcoeff=atoi(pcoeffc);	
	}
	if (strcmp(sigma0c, defau)!=0) {
		sigma0=atof(sigma0c);
		if(sigma0<0.0) sigma0=0.0;
	}	
	if (strcmp(szeroc, defau)!=0){
		szero=atoi(szeroc);
		if(szero) sigma0=sigma0-shift;	/*shift in case F has a singularity at zero*/
	}
	
	/****************************************************************************\
	 SPECIFICATION OF THE PARAMETERS EXPSX, SIGMA, B AND SCALING OF THE TOLERANCE TOL. 
	 TOL IS SCALED AS FOLLOWS:  Ntol=tol/exp(sigma*x)/exp(sigma0*x)	  
	\**************************************************************************/
				
	if(x>0.0){		
		sigma=2.5/x;  		/*  optimal value of sigma*x = 2.5 */
		b=2.5*sigma;    	/*  optimal value of b=2.5*sigma */		
		expsx=exp(sigma*x);		
		Ntol=(tol/expsx)/exp(sigma0*x);
	}
	else{	        
		/* if x=0 sigma =4, b=1, ntol=tol, expsx=1 */		
		sigma=4.;												
		b=1;
		expsx=1;			
		Ntol=tol;			
	}	
	y=2*b*x;
	esp=exp((sigma-b)*x);	
	n=nstart;
	*nopt=n;	
	z[1]=1;
	/* PHI FUNCTION CALL FOR THE COMPUTATION OF COEFFICIENTS CK */
	y1=phi(z,b,sigma,sigma0,fz,szero);	
	do{	
		/* COMPUTATION OF COEFFICIENTS CK */		
		y1.eval(n);										
		esterr_new=(fabs(y1[n-2])+fabs(y1[n-1])+fabs(y1[n]))/3;
		
		if (n==nstart){ esterr=esterr_new; }
		else if(esterr_new<esterr){
			esterr=esterr_new;      
			*nopt=n; 
		}
		n+=2;		 		
	}while((esterr>Ntol)&&(n<=nmax));     
	
	*ncalc=n-2;
	
	/***************************************************************************************
		COMPUTATION OF  THE LAGUERRE SERIES EXPANSION APPROXIMATING THE INVERSE FUNCTION    
		[ c0*L0(2bx)+c1*L1(2bx)+...+ cn*Lnopt(2bx)]
	*************************************************************************************/
	
	if(pcoeff){
		sprintf(name,"coeff_x%.1f.txt",x);
		if(!(fp=fopen(name,"w"))){ printf("RELIADIFF, x=%.1f. File %s to print coefficients not open!",x,name); exit(1); }	
		fprintf(fp,"x=%e\n",x);
		for(i=0;i<=*nopt;i++)
			fprintf(fp,"coeff[%d]\t%e\n",i,y1[i]);
		fclose(fp);
	}
	
	L1=1;
	L2=1-y;	
	fx=y1[0]+y1[1]*L2;
	
	for(k=2;k<=(*nopt);k++){   
		temp1=(-y+2*((double) k)-1)/(double)(k);
		temp2=(((double) k)-1)/(double)(k); 
		LK=temp1*L2-temp2*L1;
		fx=fx+y1[k]*LK;
		L1=L2;
		L2=LK;
	}	
	/* ILf IS THE APPROXIMATED VALUE OF THE INVERSE LAPLACE FUNCTION*/							
	*ILf=fx*esp;								
	*ILf=*ILf*exp(sigma0*x);
	
	if(szero!=0)
		/*damping in case F has a singularity at zero*/
		*ILf=*ILf*exp(shift*x);       
		

	/************************************************************************************
		ABSOLUTE AND RELATIVE ERRORS ESTIMATE
	*************************************************************************************/
	/************************************************************************************
	 * 			           WARNING
	--------Calculated errors are just errors estimate, they are not errors bounds-------	
	
	*************************************************************************************/
		
	/* ABSOLUTE ERROR ESTIMATE ON ILF*/
	*absesterr=(esterr*expsx)*exp(sigma0*x); 
	
	/* RELATIVE ERROR ESTIMATE ON ILF*/
	if (fx<=FLT_EPSILON) *relesterr=(esterr)/fabs((fx+1)*exp(-b*x));
	else *relesterr=(esterr)/fabs(fx*exp(-b*x));			
	

	/*IF relesterr<tol  AND  absesterr<tol THEN IFLAG = 0*/
	*flag=0;						 
	if(*relesterr>tol){
		/*IF relesterr>tol absesterr>tol, n>=nmax THEN IFLAG= 1*/
		if(*absesterr>tol) *flag=1; 
		/*IF relesterr>tol AND  absesterr<tol THEN IFLAG = 2*/
		else *flag=2;				 			
	}
	else{
		/*IF relesterr<tol AND absesterr>tol THEN IFLAG = 3 */
		if(*absesterr>tol) *flag=3; 
	}
	
	
	/***************************************************************************************
		USED VALUES FOR OPTIONAL PARAMETERS
	*************************************************************************************/
	sprintf(sigma0c,"%f",sigma0);
	sprintf(pcoeffc,"%d",pcoeff);
	sprintf(szeroc,"%d",szero);
	sprintf(nmaxc,"%d",nmax);
	sprintf(tolc,"%f",tol);	
	
	return 0;
}
/*END OF RELIADIFF()*/


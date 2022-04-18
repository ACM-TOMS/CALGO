/*
	===========================================================================
	PURPOSE
	===========================================================================

	This program is a driver to test the behavior of RELIADIFF on the
	provided Laplace Transform
	fzE
	with provided Inverse function
	gzE
	The function are defined in the file FunctionForTest.c
	
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
	User can give:

	NO ARGUMENTS :   
	    The driver will use all default values:
	    tol= depends on the function        required accuracy
	    sigma0= 1				convergence abscissa for F
	    nmax= 2000				maximum number of series coefficients
	    pcoeff= 0				don't print coefficients
	    szero=0				the Transform has not a singularity at zero
            x= 1, 5, 10, 15       		evaluation points

	1 ARGUMENT: tol (that could be a string "n" to give it the default value)
	    The driver will use other default values:
	    sigma0= 1							
	    nmax= 2000							
	    pcoeff= 0		
	    szero=0	
            x= 1, 5, 10, 15       				     				

	2 ARGUMENTS: tol,sigma0 (each one could be a string "n" to give it the
		default value)
	    The driver will use other default values:				 				
	    nmax= 2000							
	    pcoeff= 0		
	    szero=0	
            x= 1, 5, 10, 15       			

	3 ARGUMENTS: tol,sigma0,nmax (each one could be a string "n" to give it the 
		default value)
	    The driver will use other default values:				 					    
	    pcoeff= 0		
	    szero=0	
            x= 1, 5, 10, 15    	
				
	4 ARGUMENTS: tol,sigma0,nmax,pcoeff (each one could be a string "n" to give 
			      it the default value)
	    The driver will use other default values:				 					    	
	    szero=0	
            x= 1, 5, 10, 15       	

	5 ARGUMENTS: tol,sigma0,nmax,pcoeff,szero (each one could be a string "n" to give 
		it the default value)
	    The driver will use other default values:				 					    
            x= 1, 5, 10, 15       	
   
	>8 ARGUMENTS: tol,sigma0,nmax,pcoeff,szero (each one could be a string "n" to give 
		it the default value)
	    and:
	    6-TH ARGUMENT: range	that should be 0 or greater
	      IF range>0 the inverse functions are computed on
		      a set of equispaced points belonging to
		      the interval [a,b] with step size "step".
		  In this case:
		  7-TH ARGUMENT: a,
		  8-TH ARGUMENT: b,
		  9-TH ARGUMENT: step

	      IF range=0 the inverse are computed on a given set
		      of values.
		  In this case
		  7-TH ARGUMENT: dim = number of evaluation
					points
		  8-TH TO (7+dim)-TH ARGUMENT: the evaluation 
					points. 	

	===========================================================================
	RETURN VALUE
	===========================================================================

	INTEGER: 	It will contain a diagnostic on work of the driver.
			If it is 0, everything works fine.

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
#include "./FunctionForTest.c"


int main(int argc,char **argv){

/*RELIADIFF INPUT*/
	double		*x=NULL;		//Inverse Function evaluation point(s)
	T<double>	(*fz)(T<double>);   	//Function F Pointer
		
/*RELIADIFF OPTIONAL INPUT*/
	char        szeroc[10]="n"; 		//optional: there's or not a singularity at zero for F
	char		pcoeff[10]="n";		//optional: if user doesn't say yes, the driver doesn't print coefficients
	char		sigma0c[10]="n";   	//optional: abscissa of convergence of F  	
	char 		tolc[10]="n";		//optional: tolerance on accuracy
	char		nmaxc[10]="n"; 		//optional: maximum number of series coefficients
	
/*RELIADIFF OUTPUT*/
	double		*absesterr=NULL;	//absolute error estimate
	double		*relesterr=NULL;	//relative error estimate	
	double 		ILf=0;			//Inverse Function f computed
	int		*flag=NULL;		//diagnostics on the result
	int		nopt;			//diagnostics on the software work
	int		ncalc;			//diagnostics on the software work
	int		ierr=0;			//diagnostics on the software work
	
/*AUXILIARY VARIABLES*/	
	int 		szero=0;	
	int		dim=0;
	double		sigma0;
	int		i; 
	char		defau[5]="n";
	int         	range=0;
	double      	a=0,b=0,step=0,dimd=0.;
	int 		nmax=0;
	double		g;		
/*.... true absolute error ............................ */
	double		*abserr=NULL;				
/*... true relative error...............................*/
	double		*pabserr=NULL;		
/*... Function Pointers	................................*/
	double		(*gz)(double);		
/*... Output File Pointer and Name .....................*/
	FILE		*fp;
	char		name[20]="";
/*... Timing variables .................................*/
	struct 		timeval startsingle, endsingle;
	double		elapsedsingle, elapsedtotal,maxsingle,maxsingleold,xmax=-1;
	double		tstart,tend;

	printf("\n++++++++++++++++++++++++++++++\n");
	printf("  EXAMPLE PROGRAM STARTED!      \n");
	printf("++++++++++++++++++++++++++++++\n\n\n");

	/* .............. input reading ..........................*/

	if(argc>1){
		strcpy(tolc,argv[1]);			/*first argument: tol*/
		if(argc>2){
			strcpy(sigma0c,argv[2]);	/*second argoments: sigma0*/
			if(argc>3){
				strcpy(nmaxc,argv[3]); /*third argument: nmax*/
				if(argc>4){
					strcpy(pcoeff,argv[4]); /*fourth argument: pcoeff*/
					if(argc>5){
						strcpy(szeroc,argv[5]); /*fifth argument: szero*/
						if(strcmp(szeroc, defau)!=0) szero=atoi(szeroc);
						if(argc>6){
							range=atoi(argv[6]); /*sixth argument: range*/
							if(range){				/*if range=1 then sixth, seventh and eighth arguments 
														give interval and step for evaluation points*/
								if(argc>=10){
									a=(double)atof(argv[7]);
									b=(double)atof(argv[8]);
									step=(double)atof(argv[9]);
									dimd=(b-a)/step;
									dim=(int)(dimd+1);
									x=(double*)calloc(dim,sizeof(double));
									i=0;
									while(i<dim){
										x[i]=a;
										i++;
										a=a+step;
									}
								}/*endif*/
								else{
									printf("\nToo few arguments!! Please repeat the command.\n");
									printf("Note: 6th argument > 0 ----> #arguments >= 9 (>9 ignored)\n");
									exit(1);
								}/*endelse*/

							}/*endif*/
							else{/*if range=0 then seventh argument is the number of element of x*/
								if((argc>7)&&(strcmp(argv[7], defau)!=0)){
									dim=atoi(argv[7]);	/*seventh argument is the number of element of x*/						
									if(argc>=dim+8){
										x=(double*)calloc(dim,sizeof(double));
										for(i=0;i<dim;i++){
											x[i]=atof(argv[8+i]);		/*... and next dim arguments are 
																			the evaluation points*/
										}/*endfor*/
									}/*endif*/
									else{
										printf("\nToo few arguments!! Please repeat the command.\n");
										printf("Note: 6th argument = 0 ----> #arguments >= %d (>%d ignored)\n",dim+7, dim+7);
										exit(1);
									}/*endelse*/
								}/*endif*/
								else{
									printf("\nToo few arguments!! Please repeat the command.\n");
									printf("Note: 6th argument = 0 ----> #arguments > 7\n");
									exit(1);
								}/*endelse*/
							}/*endelse*/
						}/*endif 6*/
					}/*endif 5*/
				}/*endif 4*/
			}/*endif 3*/
		}/*endif 2*/
	} /*endif 1*/
	
	printf("\n--------------------------------------------------------------------");
	printf("\nuser required accuracy:\ttol=%s",tolc);
	if(strcmp(tolc, defau)==0) printf("\n----RELIADIFF will use tol default value");
	printf("\nuser required maximum number of coefficients:\tnmax=%s",nmaxc);
	if(strcmp(nmaxc, defau)==0) printf("\n----RELIADIFF will use nmax default value");
	else {
	      nmax=atoi(nmaxc);
	      if(nmax>MaxLength) {
			printf("\nmaximum number of coefficients allowed = TADIFF MaxLength: %d", MaxLength);
			printf("\n----RELIADIFF will pose nmax = TADIFF MaxLength");
	      }
	      else if(nmax<8) printf("\n----RELIADIFF will use nmax default value");
	}
	printf("\nuser required convergence abscissa:\tsigma0=%s",sigma0c);
	if(strcmp(sigma0c, defau)==0) printf("\n----RELIADIFF will use sigma0 default value\n");
	else sigma0=atof(sigma0c);
	if(x==NULL) {
		dim=4;
		x=(double*)calloc(dim,sizeof(double));
		x[0]=1;
		x[1]=5;
		x[2]=10;
		x[3]=15;
	}/*endif*/
	printf("\nnumber of evaluation points=\t%d\n",dim);
	printf("--------------------------------------------------------------------\n");
/* .............. memory space allocation of local variables..........................*/
	flag=(int*)calloc(dim,sizeof(int));
	relesterr=(double*)calloc(dim,sizeof(double));
	abserr=(double*)calloc(dim,sizeof(double));
	absesterr=(double*)calloc(dim,sizeof(double));
	pabserr=(double*)calloc(dim,sizeof(double));
			
/* ......................... TEST STARTS..........................................................*/
	fz = fzE; gz = gzE; 
	
	if (szero){
		printf("\nThe Transform has a singularity at zero.\n");	
	}
	
	/*... Managing output file ...............................................................*/
	sprintf(name,"output_demo2.txt");
	printf("\nTable File: %s\n",name);
	if(!(fp=fopen(name,"w"))){ printf("File %s not open!",name); exit(1); }

	fprintf(fp,"\n");
	fprintf(fp,"*******************************************************************************************************\n");
	fprintf(fp,"*                                       Test Function                                                 *\n");
	fprintf(fp,"*******************************************************************************************************\n");
	fprintf(fp,"\n");
	
	fprintf(fp,"user required Tolerance on accuracy: %s;\n",tolc);
	if (strcmp(tolc, defau)==0) fprintf(fp,"---- RELIADIFF will give to the tolerance a default value.\n");
	fprintf(fp,"user given abscissa of convergence on F: %s;\n",sigma0c);
	if (strcmp(sigma0c, defau)==0) fprintf(fp,"---- RELIADIFF will give to the convergence abscissa a default value.\n");
	fprintf(fp,"user required maximum number of Laguerre series coefficients: %s;\n",nmaxc);
	if (strcmp(nmaxc, defau)==0) fprintf(fp,"---- RELIADIFF will give to maximum number of coefficients a default value.\n");

	if (szero){
		fprintf(fp,"\n*******************************************************************************************************");
		fprintf(fp,"\nRequired TRASLATION ON z in Laplace Trasform fx (szero=%d): z-1\n",szero);
		fprintf(fp,"TRASLATION ON sigma0 : sigma0+1\n");
		fprintf(fp,"DAMPING ON THE COMPUTED INVERSE: f=f*exp(-1*x)\n");
		fprintf(fp,"*******************************************************************************************************\n");
		if (strcmp(sigma0c, defau)!=0) fprintf(fp,"If sigma0<0 RELIADIFF will pose sigma0=0: here sigma0=%f\n",sigma0+1); 
	}
	else if (strcmp(sigma0c, defau)!=0) fprintf(fp,"If sigma0<0 RELIADIFF will pose sigma0=0: here sigma0=%f\n",sigma0); 
	
	fprintf(fp,"\n---------------------------------------------------------------\n");
			
	/*... timing variables initializing for the next function .....................................................*/
	elapsedsingle=0;
	elapsedtotal=0;
	maxsingle=0;
	maxsingleold=0;
	xmax=-1;
			
	for(i=0; i<dim; ++i){					
		gettimeofday(&startsingle, NULL);		
		/************************************** CALLING RELIADIFF ***************************************************/
		ierr=RELIADIFF(x[i],fz,szeroc,pcoeff,sigma0c,tolc,nmaxc,&ILf,&absesterr[i],&relesterr[i],&flag[i],&nopt,&ncalc); 	
		/************************************************************************************************************/
		gettimeofday(&endsingle, NULL);
					
		/*... Timing ..................................................*/
		tstart=startsingle.tv_sec+(startsingle.tv_usec/1000000.0);
		tend=endsingle.tv_sec+(endsingle.tv_usec/1000000.0);
		elapsedsingle=tend-tstart;
		maxsingleold=maxsingle;
		maxsingle=max(elapsedsingle, maxsingle);
		if (maxsingle!=maxsingleold) xmax=x[i];
		elapsedtotal=elapsedtotal+elapsedsingle;
		
		switch(ierr){
		case 1:
			printf("\nx=%f - RELIADIFF stopped: x<0.0!\n It must be x>=0\n\n\n",x[i]);
			break;
		default: 
			if (i==0){
				fprintf(fp,"\n\n\n\nUsed Tolerance on accuracy: %s;\n",tolc);
				fprintf(fp,"Used abscissa of convergence on F: %s;\n",sigma0c);
				fprintf(fp,"Used maximum number of Laguerre series coefficients: %s;\n\n",nmaxc);
				fprintf(fp,"\n                                               TABLE\n");
				fprintf(fp,"--------------------------------------------------------------------------------------------------------------\n");
				fprintf(fp,"|     x    |      f_comp      | trueabserr | estabserr | truerelerr | estrelerr |   Nopt   |   Ncal   | FLAG |\n");
				fprintf(fp,"--------------------------------------------------------------------------------------------------------------\n");
			}
			/******************************************************************************************
			*                                  COMPUTED ERROR ESTIMATES                               *
			*******************************************************************************************
			REMARK:
			AT EACH STEP OF THE DO-LOOP RELIADIFF RETURNS:
			  1.absesterr    ABSOLUTE ERROR ESTIMATE
			  2.relesterr    RELATIVE ERROR ESTIMATE
			  3.flag         ERROR DIAGNOSTIC
			  4.ILf          INVERSE LAPLACE FUNCTION
			  5.nopt         OPTIMAL VALUE OF n
			  5.ncalc        REACHED VALUE OF n								
			*******************************************************************************************\
			See Doc/RELIADIFF_userguide.pdf, section 8 to know about the meaning of flag and nopt/ncalc
			\*******************************************************************************************/
			/*... g is the true INVERSE LAPLACE FUNCTION ...............................................................*/
			g=gz(x[i]);
			/*... Managing true errors ...............................................................*/
			if(g==0){ 
				abserr[i]=fabs(ILf)+DBL_EPSILON; 
				pabserr[i]=abserr[i];
			}
			else{
				/* abserr   absolute error       	---> trueabserr	*/
				abserr[i]=fabs(ILf-g)+DBL_EPSILON; 																
				/* pabserr	relative error       	---> truerelerr	*/
				pabserr[i]=fabs(ILf-g)/fabs(g)+DBL_EPSILON; 
				/* absesterr                     	-------> estabserr	*/
				/* relesterr                     	-------> estrelerr	*/
			}					 
			/*... Managing output file ...............................................................*/
			fprintf(fp,"| %4.2e | %14.8e | %10.3e | %9.3e | %10.3e | %9.3e | %8d | %8d | %4d |\n", x[i], ILf, abserr[i], absesterr[i], pabserr[i], relesterr[i], nopt, ncalc, flag[i]);
						
			break;
		}/*endswitch*/
		/* UPDATE OF i (and x)*/
	}/*endfor on x*/

	printf("\n*****************************************************************************\n");
	fprintf(fp,"\n\n\nelapsed time (all values) in sec = %e;",elapsedtotal); 
	fprintf(fp,"\nmax elapsed time (x=%4.2e) in sec = %e;\n\n\n",xmax,maxsingle);
	print_warning(fp);
	print_flags_file(fp);
	print_N_file(fp,nmax);
	fclose(fp);
						   
	free(flag);
	free(relesterr);
	free(abserr);
	free(absesterr);
	free(pabserr);

	printf("\n\n\n++++++++++++++++++++++++++++++\n");
	printf("EXAMPLE PROGRAM TERMINED!\n");
	printf("++++++++++++++++++++++++++++++\n");
			
	/**************************************************************************/		
	/*                         END OF TEST                                    */
	/**************************************************************************/
	return 0;
}/*endmain*/
/*END OF MAIN*/


/*
	===========================================================================
	PURPOSE
	===========================================================================

	This program is a driver to test the behavior of RELIADIFF on the
	provided database of Laplace Transform/Inverse functions.
	
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
	    ntf= 1				number of functions to test RELIADIFF with
	    nmax= 2000				maximum number of series coefficients
	    pcoeff= 0				don't print coefficients
        x= 1, 5, 10, 15       			evaluation points

	1 ARGUMENT: tol (that could be a string "n" to give it the default
			 value)
	    The driver will use other default values:
	    ntf= 1								
	    nmax= 2000							
	    pcoeff= 0							
            x= 1, 5, 10, 15       				   				

	2 ARGUMENTS: tol,ntf (each one could be a string "n" to give it the 
			      default value, ntf could be "a" to mean all the 
			      database functions)
	    The driver will use other default values:				 				
	    nmax= 2000							
	    pcoeff= 0							
            x= 1, 5, 10, 15       				       			

	3 ARGUMENTS: tol,ntf,nmax (each one could be a string "n" to give it the 
			      default value, ntf could be "a" to mean all the 
			      database functions)
	    The driver will use other default values:				 					    
	    pcoeff= 0							
            x= 1, 5, 10, 15       				      	
				
	4 ARGUMENTS: tol,ntf,nmax,pcoeff (each one could be a string "n" to give 
			      it the default value, ntf could be "a" to mean all 
			      the database functions)
	    The driver will use other default values:				 					    
            x= 1, 5, 10, 15       	
   
	>7 ARGUMENTS: tol,ntf,nmax,pcoeff (each one could be a string "n" to give 
			      it the default value, ntf could be "a" to mean all 
			      the database functions)
	    and:
	    5-TH ARGUMENT: range	that should be 0 or greater
	      IF range>0 the inverse functions are computed on
		      a set of equispaced points belonging to
		      the interval [a,b] with step size "step".
		  In this case:
		  6-TH ARGUMENT: a,
		  7-TH ARGUMENT: b,
		  8-TH ARGUMENT: step

	      IF range=0 the inverse are computed on a given set
		      of values.
		  In this case
		  6-TH ARGUMENT: dim = number of evaluation
					points
		  7-TH TO (6+dim)-TH ARGUMENT: the evaluation 
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

#include "./database/dbL.h"
#include <sys/stat.h>

int main(int argc,char **argv){
/*... Local Variables specification................................. */
/*RELIADIFF INPUT*/
	double		*x=NULL;		//Inverse Function evaluation point(s)
	T<double>	(*fz)(T<double>);   	//Function F Pointer
	
	
/*RELIADIFF OPTIONAL INPUT*/
	char        	szeroc[10]="n"; 	//optional: there's or not a singularity at zero for F
	char		pcoeffc[10]="n";		//optional: if user doesn't say yes, the driver doesn't print coefficients
	char		sigma0c[10]="n";    //optional: abscissa of convergence of F  	
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
/*....ntf is the number of the functions to test together. By default it is equal to 1......*/
	int		dim=0;
	int 		pcoeff=0;
	int		tf,i,j,ntf=1,ntfmax=89; 
	char		defau[5]="n",all[5]="a";
	int         	range=0;
	double      	a1=3./5., b1=5./7., c1=-9./7.;
	double      	a=0,b=0,step=0,dimd=0.;
	double		sigma0=0.0;
	double 		tol=0;		
	int 		szero=0;
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
	char		path[100]="";
	char		cmd[25]="";
/*... Timing variables .................................*/
	struct 		timeval startsingle, endsingle;
	double		elapsedsingle, elapsedtotal,maxsingle,maxsingleold,xmax=-1;
	double		tstart,tend;

	
/************************************************************************
THE TEST-DRIVER MAIN PROGRAM IS DESIGNED IN ORDER TO MANAGE THE TEST 
OF THE RELATIVE FUNCTION ON THE DATABASE FUNCTIONS
************************************************************************/
	
	printf("\n++++++++++++++++++++++++++++++\n");
	printf("   DRIVER PROGRAM STARTED!      \n");
	printf("++++++++++++++++++++++++++++++\n\n\n");

	if(argc>1){
			strcpy(tolc,argv[1]);			/*first argument: tol*/
		if(argc>2){
			if (strcmp(argv[2], defau)!=0){
				if (strcmp(argv[2], all)!=0)
					ntf=atoi(argv[2]);	/*second argument: ntf*/
				else ntf=ntfmax;
			}
			if(argc>3){
				strcpy(nmaxc,argv[3]); 	/*third argument: nmax*/
				if(argc>4){
					strcpy(pcoeffc,argv[4]);	/*fourth argument: pcoeff*/
					if(argc>5){
						range=atoi(argv[5]); 	/*fifth argument: range*/
						if(range){	      	/*if range=1 then fifth, sixth, sevent arguments 
									  give interval and step for evaluation points*/
							if(argc>=9){
								a=(double)atof(argv[6]);
								b=(double)atof(argv[7]);
								step=(double)atof(argv[8]);
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
								printf("Note: 5th argument > 0 ----> #arguments >= 8 (>8 ignored)\n");
								exit(1);
							}/*endelse*/
						}/*endif*/
						else{			/*if range=0 then sixth argument is the number of element of x*/
							if((argc>6)&&(strcmp(argv[6], defau)!=0)){
								dim=atoi(argv[6]);	/*sixth argument is the number of element of x*/		
								if(argc>=dim+7){
									x=(double*)calloc(dim,sizeof(double));
									for(i=0;i<dim;i++){
										x[i]=atof(argv[7+i]);	/*... and next dim arguments are 
													  the evaluation points*/
									}/*endfor*/
								}/*endif*/
								else{
									printf("\nToo few arguments!! Please repeat the command.\n");
									printf("Note: 5th argument = 0 ----> #arguments >= %d (>%d ignored)\n",dim+6, dim+6);
									exit(1);
								}/*endelse*/
							}/*endif*/
							else{
								printf("\nToo few arguments!! Please repeat the command.\n");
								printf("Note: 5th argument = 0 ----> #arguments > 6\n");
								exit(1);
							}/*endelse*/
						}/*endelse*/
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
	if(strcmp(pcoeffc, defau)!=0) {
	    printf("\nRELIADIFF will print the used Taylot coefficients");
	    pcoeff=atoi(pcoeffc);
	}
	printf("\nnumber of database functions to test:\tntf=%d",ntf);
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
	
/*..............TESTS START: THE OUTPUT IS PUT ON FILES NUMBERED AS THE FUNCTION NUMBER........ */
	for(j=0;j<ntf;j++){
		if(ntf<ntfmax){ printf("\nInsert the next Function to test (number):\n"); scanf("%d",&tf); }
		else tf=j+1;
		if ((tf<1)||(tf>ntfmax)){ printf("Number of function not valid! It must be in [1,%d]\n",ntfmax); j=j-1; }
		else{			
/* ......................... TEST STARTS..........................................................*/
			szero=0;
			printf("\n************************* Function number: %d ******************************\n",tf);

			switch(tf){
			case 1:			
				fz=fz1; gz=gz1; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-03; break;
			case 2:
				fz=fz2; gz=gz2; sigma0=0; tol= 1.0e-01; break; 
			case 3:
				fz=fz3; gz=gz3; sigma0=a1; tol= 1.0e-01; break; 
			case 4:
				fz=fz4; gz=gz4; sigma0=a1; tol= 1.0e-03; break; 
			case 5:
				fz=fz5; gz=gz5; sigma0=0; tol= 1.0e-03; break;
			case 6:
				fz = fz6; gz=gz6; sigma0=-a1; tol= 1.0e-03; break;  
			case 7:
				fz=fz7; gz=gz7; sigma0=0; tol= 1.0e-04; break;
			case 8:
				fz=fz8; gz=gz8; sigma0=0; tol= 1.0e-02; break; 
			case 9:
				fz = fz9; gz = gz9; sigma0=0; tol= 1.0e-05; break;
			case 10:
				fz = fz10; gz = gz10; sigma0=0; tol= 1.0e-04; break;
			case 11:
				fz = fz11; gz = gz11; sigma0=0; tol= 1.0e-05; break;
			case 12:
				fz = fz12; gz = gz12; sigma0=0.5; tol= 1.0e-06; break;
			case 13:
				fz = fz13; gz = gz13; sigma0=a1; tol= 1.0e-05; break;
			case 14:
				fz=fz14; gz=gz14; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-06; break;
			case 15:
				fz=fz15; gz=gz15; sigma0=GSL_MAX_DBL(-a1,0); tol= 1.0e-05; break;
			case 16:
				fz = fz16; gz = gz16; sigma0=0; tol= 1.0e-06; break;
			case 17:
				fz = fz17; gz = gz17; sigma0=0; tol= 1.0e-06; break;
			case 18:
				fz = fz18; gz = gz18; sigma0=a1; tol= 1.0e-06; break;
			case 19:
				fz=fz19; gz=gz19; sigma0=GSL_MAX_DBL(-a1,GSL_MAX_DBL(-b1,-c1)); tol= 1.0e-06; break;
			case 20:
				fz = fz20; gz = gz20; sigma0=0; tol= 1.0e-06; break;
			case 21:
				fz = fz22; gz = gz22; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-04; break;
			case 22:
				fz = fz22; gz = gz22; sigma0=0; tol= 1.0e-06; break;
			case 23:
				fz = fz23; gz = gz23; sigma0=0; tol= 1.0e-06; break;
			case 24:
				fz = fz24; gz = gz24; sigma0=0; tol= 1.0e-04; break;
			case 25:
				fz = fz25; gz = gz25; sigma0=0; tol= 1.0e-06; break;
			case 26:
				fz = fz26; gz = gz26; sigma0=a1; tol= 1.0e-06; break;
			case 27:
				fz=fz27; gz=gz27; sigma0=0; tol= 1.0e-06; break;
			case 28:
				fz=fz28; gz=gz28; sigma0=0; tol= 1.0e-06; break;
			case 29:
				fz=fz29; gz=gz29; sigma0=-a1; tol= 1.0e-03; break;
			case 30:
				fz=fz30; gz=gz30; sigma0=-a1; tol= 1.0e-04; break;
			case 31:
				fz=fz31; gz=gz31; sigma0=GSL_MAX_DBL(gsl_sf_pow_int(GSL_MAX_DBL(-b1,0),2),gsl_sf_pow_int(a1,2)); szero=1; tol= 1.0e-01; break; 
			case 32:
				fz=fz32; gz=gz32; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-04; break;
			case 33:
				fz=fz33; gz=gz33; sigma0=0; szero=1;  tol= 1.0e-06; break; 
			case 34:
				fz=fz34; gz=gz34; sigma0=0; szero=1; tol= 1.0e-01; break;  
			case 35:
				fz=fz35; gz=gz35; sigma0=GSL_MAX_DBL(gsl_sf_pow_int(a1,2),gsl_sf_pow_int(GSL_MAX_DBL(-b1,0),2)); szero=1; tol= 1.0e-01; break;
			case 36:
				fz=fz36; gz=gz36; sigma0=gsl_sf_pow_int(a1,2); szero=1; tol= 1.0e-01; break;  
			case 37:
				fz=fz37; gz=gz37; sigma0=GSL_MAX_DBL(-a1,-b1); szero=1; tol= 1.0e-02; break;
			case 38:
				fz=fz38; gz=gz38; sigma0= 0; tol= 1.0e-02; break; 
			case 39:
				fz = fz39; gz = gz39; sigma0=0; tol= 1.0e-06; break; 
			case 40:
				fz = fz40; gz = gz40; sigma0=-0.5; tol= 1.0e-04; break;
			case 41:
				fz = fz41; gz = gz41; sigma0=1; tol= 1.0e-06; break; 
			case 42:
				fz = fz42; gz = gz42; sigma0=2; tol= 1.0e-06; break;
			case 43:
				fz = fz43; gz = gz43; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-04; break;   
			case 44:
				fz = fz44; gz = gz44; sigma0=0; tol= 1.0e-05; break;
			case 45:
				fz = fz45; gz = gz45; sigma0=0; tol= 1.0e-04; break;
			case 46:
				fz = fz46; gz = gz46; sigma0=a1; tol= 1.0e-06; break;			
			case 47:
				fz=fz47; gz=gz47; sigma0=0; tol= 1.0e-05; break;
			case 48:
				fz=fz48; gz=gz48; sigma0=0; tol= 1.0e-05; break;
			case 49:
				fz=fz49; gz=gz49; sigma0=gsl_sf_pow_int(GSL_MAX_DBL(-a1,0),2); szero=1; tol= 1.0e-01; break;
			case 50:
				fz=fz50; gz=gz50; sigma0=GSL_MAX_DBL(-2*a1,0); tol= 1.0e-04; break;
			case 51:
				fz=fz51; gz=gz51; sigma0=0; szero=1; tol= 1.0e-06; break; 
			case 52:
				fz=fz52; gz=gz52; sigma0=0; szero=1; tol= 1.0e-01; break; 
			case 53:
				fz=fz53; gz=gz53; sigma0=-a1; szero=1; tol= 1.0e-02; break;
			case 54:
				fz=fz54; gz=gz54; sigma0=0; szero=1; tol= 1.0e-01; break; 
			case 55:
				fz=fz55; gz=gz55; sigma0=-a1; szero=1; tol= 1.0e-01; break; 
			case 56:
				fz=fz56; gz=gz56; sigma0=GSL_MAX_DBL(-a1,-b1); szero=1; tol= 1.0e-01; break; 
			case 57:
				fz=fz57; gz=gz57; sigma0=gsl_sf_pow_int(GSL_MAX_DBL(-a1,0),2); szero=1; tol= 1.0e-01; break; 
			case 58:
				fz=fz58; gz=gz58; sigma0=0; tol= 1.0e-01; break; 
			case 59:
				fz=fz59; gz=gz59; sigma0=-1; tol= 1.0e-01; break; 
			case 60:
				fz=fz60; gz=gz60; sigma0=0; tol= 1.0e-01; break; 
			case 61:
				fz=fz61; gz=gz61; sigma0=0; tol= 1.0e-02; break; 
			case 62:
				fz=fz62; gz=gz62; sigma0=0; tol= 1.0e-02; break; 
			case 63:
				fz = fz63; gz = gz63; sigma0=GSL_MAX_DBL(-a1,-b1); tol= 1.0e-04; break;
			case 64:
				fz = fz64; gz = gz64; sigma0=0; tol= 1.0e-04; break;
			case 65:
				fz = fz65; gz = gz65; sigma0=GSL_MAX_DBL(a1,0); tol= 1.0e-06; break;	
			case 66:
				fz=fz66; gz=gz66; sigma0=GSL_MAX_DBL(0,-2*a1); tol= 1.0e-05; break;	
			case 67:
				fz=fz67; gz=gz67; sigma0=0; tol= 1.0e-05; break; 
			case 68:
				fz=fz68; gz=gz68; sigma0=0; tol= 1.0e-04; break;
			case 69:
				fz=fz69; gz=gz69; sigma0=0; szero=1; tol= 1.0e-01; break; 
			case 70:			
				fz=fz70; gz=gz70; sigma0=0; szero=1; tol= 1.0e-03; break;
			case 71:
				fz=fz71; gz=gz71; sigma0= 0; tol= 1.0e-06; break; 
			case 72:
				fz=fz72; gz=gz72; sigma0=0; tol= 1.0e-06; break;
			case 73:
				fz=fz73; gz=gz73; sigma0=gsl_sf_pow_int(GSL_MAX_DBL(-a1,0),2); szero=1; tol= 1.0e-04; break;
			case 74:
				fz=fz74; gz=gz74; sigma0=0; szero=1; tol= 1.0e-04; break;
			case 75:
				fz=fz75; gz=gz75; sigma0=0; tol= 1.0e-05; break; 
			case 76:
				fz=fz76; gz=gz76; sigma0=gsl_sf_pow_int(GSL_MAX_DBL(-a1,0),2); tol= 1.0e-05; break; 
			case 77:
				fz=fz77; gz=gz77; sigma0=0; tol= 1.0e-05; break; 
			case 78:
				fz=fz78; gz=gz78; sigma0=0; szero=1; tol= 1.0e-04; break;
			case 79:
				fz=fz79; gz=gz79; sigma0=gsl_sf_pow_int(GSL_MAX_DBL(-a1,0),2); szero=1; tol= 1.0e-03; break;
			case 80:
				fz=fz80; gz=gz80; sigma0=0; tol= 1.0e-01; break; 
			case 81:
				fz=fz81; gz=gz81; sigma0=0; tol= 1.0e-06; break; 
			case 82:
				fz=fz82; gz=gz82; sigma0=0; szero=1; tol= 1.0e-01; break; 
			case 83:
				fz=fz83; gz=gz83; sigma0=0; tol= 1.0e-04; break;
			case 84:
				fz=fz84; gz=gz84; sigma0=0; szero=1; tol= 1.0e-01; break;
			case 85:
				fz=fz85; gz=gz85; sigma0=0; tol= 1.0e-04; break; 
			case 86:
				fz=fz86; gz=gz86; sigma0=0; tol= 1.0e-03; break; 				
			case 87:	
				fz=fz87; gz=gz87; sigma0=0; tol= 1.0e-06; break;
			case 88:	
				fz=fz88; gz=gz88; sigma0=0; tol= 1.0e-06; break; 
			case 89:
				fz=fz89; gz=gz89; sigma0=1; tol= 1.0e-05; break; 				
			default:
				printf("\n\n***********\nFunction number must be in [1,%d].\nWe test on a default function: 01\n***********\n\n",ntfmax);
				fz = fz1; gz = gz1; sigma0=0;  break;
			} 
			if (szero){
				printf("\nThe Transform has a singularity at zero.\n");	
			}
			/*... Managing output file ...............................................................*/

			sprintf(path,"fz%02dfiles",tf);
			if(mkdir(path,S_IRWXU|S_IRWXG|S_IRWXO) == -1){ 
				printf("\nFunction n. %d. Directory for output files not created!\n",tf);
				return 1;
			}
			sprintf(name,"%s/fz%02dtable.txt",path,tf);
			printf("\nTable File: %s\n",name);
			if(!(fp=fopen(name,"w"))){ printf("File %s not open!",name); exit(1); }	
			
			
			fprintf(fp,"\n");
			fprintf(fp,"*******************************************************************************************************\n");
			fprintf(fp,"*                                     Test Function n.%02d                                             *\n", tf);
			fprintf(fp,"*******************************************************************************************************\n");
			fprintf(fp,"\n");
			fprintf(fp,"user required Tolerance on accuracy: %s;\n",tolc);
			if(strcmp(tolc, defau)==0) fprintf(fp,"---- using database default tolerance for this function: %e;\n",tol);
			fprintf(fp,"default abscissa of convergence on this transform: %f;\n",sigma0);
			fprintf(fp,"user required maximum number of Laguerre series coefficients: %s;\n",nmaxc);
			if (strcmp(nmaxc, defau)==0) fprintf(fp,"---- RELIADIFF will give to maximum number of coefficients a default value.\n");

			sprintf(sigma0c,"%f",sigma0);
			if(strcmp(tolc, defau)==0) sprintf(tolc,"%f",tol);

			if (szero){
				fprintf(fp,"\n*******************************************************************************************************");
				fprintf(fp,"\nRequired TRASLATION ON z in Laplace Trasform fx (szero=%d): z-1\n",szero);
				fprintf(fp,"TRASLATION ON sigma0 : sigma0+1=%f\n", sigma0+1);
				fprintf(fp,"DAMPING ON THE COMPUTED INVERSE: f=f*exp(-1*x)\n");
				fprintf(fp,"*******************************************************************************************************\n");
				fprintf(fp,"If sigma0<0 RELIADIFF will pose sigma0=0: here sigma0=%f\n",sigma0+1); 
			}
			else fprintf(fp,"If sigma0<0 RELIADIFF will pose sigma0=0: here sigma0=%f\n",sigma0); 
			sprintf(szeroc,"%d",szero);
			
			
			fprintf(fp,"\n---------------------------------------------------------------\n");
			
			/*... timing variables initializing for the next function .....................................................*/
			elapsedsingle=0;
			elapsedtotal=0;
			maxsingle=0;
			maxsingleold=0;
			xmax=-1;
				
			/*... Evaluating The Inverse Function on all the points in x[] .................................................*/
			for(i=0; i<dim; ++i){					
				gettimeofday(&startsingle, NULL);		
				/************************************** CALLING RELIADIFF ***************************************************/
				ierr=RELIADIFF(x[i],fz,szeroc,pcoeffc,sigma0c,tolc,nmaxc,&ILf,&absesterr[i],&relesterr[i],&flag[i],&nopt,&ncalc); 	
				/************************************************************************************************************/
				gettimeofday(&endsingle, NULL);
				
				/*... Timing ..................................................*/
				tstart=startsingle.tv_sec+(startsingle.tv_usec/1000000.0);
				tend=endsingle.tv_sec+(endsingle.tv_usec/1000000.0);
				elapsedsingle=tend-tstart;
				maxsingleold=maxsingle;
				maxsingle=GSL_MAX_DBL(elapsedsingle, maxsingle);
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
						fprintf(fp,"\n                                                          TABLE\n");
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
					See Doc/RELIADIFF-userguide.pdf, section 8 to know about the meaning of flag and nopt/ncalc
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
				}
				/* UPDATE OF i (and x)*/
			}/*endfor on x*/
			if(pcoeff){
				sprintf(cmd,"mv coeff_x* %s",path);
				system(cmd);
			}
			printf("\n*****************************************************************************\n");
			fprintf(fp,"\n\n\nelapsed time (all values) in sec = %e;",elapsedtotal); 
			fprintf(fp,"\nmax elapsed time (x=%4.2e) in sec = %e;\n\n\n",xmax,maxsingle);
			print_warning(fp);
			print_flags_file(fp);
			print_N_file(fp,nmax);
			fclose(fp);                    
		}/*endelse (check on tf)*/               
	} /*endfor on ntf*/

	free(flag);
	free(relesterr);
	free(abserr);
	free(absesterr);
	free(pabserr);

    printf("\n\n\n++++++++++++++++++++++++++++++\n");
	printf("DRIVER PROGRAM TERMINED!\n");
	printf("++++++++++++++++++++++++++++++\n");
	/*                         END OF TESTS                                  */
	/**************************************************************************/
	
	return 0;
}
/*END OF MAIN*/

   
   
 





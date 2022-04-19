(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`Texts`*)


(* ::Text:: *)
(*MathTIDES. Version 1.3.0.*)
(*This file is part of TIDES.*)
(*  *)


(* ::Text:: *)
(*Alberto Abad, Roberto Barrio, Fernando Blesa, Marcos Rodriguez*)
(*Grupo de Mec\[AAcute]nica Espacial.  IUMA.*)
(*University of Zaragoza*)
(*50009 Zaragoza. Spain.*)
(**)
(* http://gme.unizar.es/software/tides*)
(* Contact: <tides@unizar.es>*)
(**)


(* ::Title:: *)
(*Texts*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection:: *)
(*Comienzo*)


BeginPackage["MathTIDES`Texts`"]


(* ::Subsection:: *)
(*S\[IAcute]mbolos*)


{
mincgen,
minhgen,
minfgen,
headstdDP,
headstdMP
}


(* ::Subsection:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`Texts`*"]
Clear @@ Names["MathTIDES`Texts`*"]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*C*)


mincgen =
"
/****************************************************************************
    
    minc_tides: kernel of the C Minimal version of TIDES
    This file is part of TIDES.

 
    Contributors:
 
    A. Abad, R. Barrio, F. Blesa, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

*****************************************************************************/

#include \"minc_tides.h\" 


double    *v, *p, **XVAR, **XVAR2, TOL_ABS, TOL_REL, tzero, deltat; 
int       VARS, PARS, MAX_ORDER, tflag;
double    fac1=0.95e0,fac2=10.e0,fac3=0.8e0 ; 
double    rmaxstep=1.e2,rminstep=1.e-2; 
int       nitermax=5, nordinc=5, minord=6, maxord=26; 
int       dense_output = -1, defect_error_control = 0;
int       accepted_steps = 0, rejected_steps = 0;
int       ipos = 1;
FILE      *fd;



/************************************************************************/ 

void minc_tides(double *var, int nvar, double *par, int npar,  double tini, double tend, double dt,
                double tolrel, double tolabs) 
{ 
	double    tol, tolo;
	volatile  double  temp, extra;
	double    t0,t1,t2, step, nstep;
	int       i, ORDER; 

	
	VARS = nvar;
	PARS = npar; 
	MAX_ORDER = maxord;
	TOL_REL = tolrel;
	TOL_ABS = tolabs;
	tzero = tini;
	deltat = dt;
	v = var; 
	p = par;
	
	declare_matrix_coefs_mc();

	t0 = tini; 
	t1= 0.e0;
	t2 = 0.e0; 

	extra = 0.e0; 
	
	
	if(t0 < tend) {
		tflag = 0;
		if(deltat < 0) deltat = -deltat;
	} else {
		tflag = 1; 
		if(deltat > 0) deltat = -deltat;
	}

	if(dense_output) {
		fprintf(fd, \"%25.15le \", tini ); 
		for(i = 0; i < VARS; i++) fprintf(fd, \"%25.15le \" , v[i]); 
		fprintf(fd, \"\\n\"); 	
	}

	while (((t0 < tend) && (tflag == 0)) ||
		   ((t0 > tend) && (tflag == 1))  ) { 
		tolerances_mc(&tol, &tolo, &ORDER);
		mincseries(t0, v, p, XVAR, ORDER, MAX_ORDER); 
		step = steps_mc(ORDER, tol);
		if(defect_error_control) steps_DEC_mc(t0, tol, ORDER, &step);
		temp = t0; 
		nstep = step + extra; 
		t0 = temp + nstep; 
		extra = (temp-t0) +nstep; 
		if(((t0 > tend) && (tflag == 0)) ||
		   ((t0 < tend) && (tflag == 1)))  nstep = (tend-temp); 
		accepted_steps +=1; 
		if(dense_output) dense_output_mc(temp, nstep, ORDER);
		horner_mc(v, nstep, ORDER) ; 
	} 
    ipos = 1;

} 

/************************************************************************/ 
void dense_output_mc(double t0, double step, int ORDER)
{
	int i;
	double tend, ti, tit, vh[VARS]; 
	tend = t0 + step;
	ti = tzero + ipos*deltat;
	tit = ti - t0;
	while(((tit <= step) && (tflag == 0)) ||
		   ((tit >= step) && (tflag == 1))  ) {
		horner_mc(vh, tit, ORDER); 
		fprintf(fd, \"%25.15le \", ti ); 
		for(i = 0; i < VARS; i++) fprintf(fd, \"%25.15le \" , vh[i]); 
		fprintf(fd, \"\\n\");
		ipos++;
		ti = tzero + ipos*deltat;
		tit = ti -t0;
	}
}

/************************************************************************/ 

void horner_mc(double *v, double t, int ORDER) 
{ 
	int i,j; 
	double temp; 
	for(j = 1; j <= VARS; j++) { 
		temp = XVAR[ORDER][j] * t; 
		for(i = ORDER-1; i >= 1; i--) { 
			temp = t * (temp + XVAR[i][j]);
		} 
		v[j-1] = temp + XVAR[0][j]; 
	} 
} 

void hornerd_mc(double *v, double t, int ORDER) 
{ 
	int i,j; 
	double temp; 
	for(j = 1; j <= VARS; j++) { 
		temp = ORDER * XVAR[ORDER][j] * t; 
		for(i = ORDER-1; i > 1; i--)	
			temp = t * (temp + i * XVAR[i][j]); 
		v[j-1]= temp + XVAR[1][j]; 
	} 
} 
/************************************************************************/ 
double norm_inf_vec_mv()
{
	int i; 
	double norm;
	norm = fabs(v[0]);
	for(i = 1; i < VARS; i++) 
		norm = max_d(norm, fabs(v[i]));
	return norm;
}

double norm_inf_mat_mc(int ord)
{
	int i; 
	double norm;
	norm = 0.e0;
	for(i = 1; i <= VARS; i++) {
		norm = max_d(norm, fabs(XVAR[ord][i]));
	}
	return norm;
}
/************************************************************************/ 

void tolerances_mc(double *tol, double *tolo, int *ORDER)
{
	static double ynb =0.e0;
	double yna, miny;
	yna = norm_inf_vec_mv();
	*tol = TOL_ABS + max_d(yna,ynb) * TOL_REL;
	miny = min_d(yna,ynb);
	if(miny > 0.) 
		*tolo = min_d(TOL_ABS/miny, TOL_REL); 
	else
		*tolo = min_d(TOL_ABS, TOL_REL);	
	*ORDER = min_i(MAX_ORDER, floor(-log(*tolo)/2)+nordinc); 
	*ORDER = max_i(minord, *ORDER); 
	ynb = yna;
}

double steps_mc(int ORDER, double tol)
{
	static double stepant =0.e0;
	double ynu, ynp, tu, tp, step;
	int ord, orda, ordp;
	double  dord, dorda, dordp, rstep;
	
	ord = ORDER+1;
	do {
		ord--;
		ynu = norm_inf_mat_mc(ord);
	} while (ynu == 0.e0 && ord > 0) ;
	
	if(ord !=0 ) {
		orda = ord-1;
		ordp = ord+1;
		dord = 1.e0/ord;
		dorda = 1.e0/orda;
		dordp = 1.e0/ordp;
		ynp = norm_inf_mat_mc(orda);		
		if(ynp == 0.e0 ) 
			step = pow(tol, dordp) * pow(1.e0/ynu, dord); 
		else {	
			tp = pow(tol, dord) * pow(1.e0/ynp, dorda);
			tu = pow(tol, dordp) * pow(1.e0/ynu, dord); 
			step = min_d(tp,tu); 
		} 
		if(stepant != 0.e0) {
			rstep = step/stepant;
			if( rstep > rmaxstep )
				step = rmaxstep * stepant;
			else if( rstep <rminstep ) 
				step = rmaxstep * stepant;
	    }
		step *= fac1;
	} else {
		printf(\"*********Error*********\");
		abort();
	}
	if(tflag ==1) step = -step;
	return step;	
}

void steps_DEC_mc(double t0, double tol, int ORDER, double *step)
{
	int iter, i;
	double norma, vh[VARS], vdh[VARS], t; 
	iter = 1; 
	norma = 0.e150; 
	t = *step;
	while((norma> fac2*tol) && (iter < nitermax)) { 
		horner_mc(vh, t, ORDER);
		hornerd_mc(vdh, t, ORDER);
		mincseries(t0+t, vh, p, XVAR2, 1, 1); 
		norma = fabs(XVAR2[1][1]-vdh[0]); 
		for(i=2; i <=VARS; i++) 
			norma +=max_d(fabs(XVAR2[1][i]-vdh[i-1]), norma); 
		if(iter >1) {
			t = fac3 *t; 
			rejected_steps +=1 ;
		}
		iter++; 
	} 
	*step = t;
}
/************************************************************************/ 
void	declare_matrix_coefs_mc()
{
	int i,j;
	XVAR  = (double **) calloc(MAX_ORDER+1,sizeof(double *));
	XVAR2 = (double **) calloc(2,sizeof(double *));
	
	for( i=0; i<=MAX_ORDER; i++ )
		XVAR[i] = (double *) calloc(VARS+1, sizeof(double));
	for( i=0; i< 2; i++ )
		XVAR2[i] = (double *) calloc(VARS+1, sizeof(double));
	
	for( i=0; i<=MAX_ORDER; i++ )
		for( j=0; j<=VARS; j++ )
			XVAR[i][j] = 0; 
	for( i=0; i< 2; i++ )
		for( j=0; j<=VARS; j++ )
			XVAR2[i][j] = 0; 
}

/************************************************************************/ 

double mul_mc(double* u, double* v, int k) 

{ 
	int j; 
	double w = 0.e0; 
	for(j = 0; j <= k; j++) 
		w += (u[j] * v[k-j]); 
	return w; 
} 

double div_mc(double* u, double* v, double* w, int k) 
{ 

	int j; 
	double ww; 
	if(v[0] == 0.e0) { 
		printf(\"*** Function div_mc found division by zero ***\");
		exit(EXIT_FAILURE); 
	} 
	ww = u[k]; 
	for(j = 1; j <= k; j++) ww -= (v[j] *w[k-j]);
	ww /= v[0]; 
	return ww; 
} 

double inv_mc(double p, double* u, double* w, int k) 
{ 
	
	int j; 
	double ww = 0.e0; 
	if(u[0] == 0.e0) { 
		printf(\"*** Function inv_mc found division by zero ***\");
		exit(EXIT_FAILURE); 
	} 
	if(k == 0) 
		ww = 1.e0; 
	else  
		for(j = 0; j < k; j++) ww -= (u[k-j] *w[j]);
	ww /= u[0]; 
	return ww*p; 
} 

double exp_mc(double* u, double* v, int k) 

{ 
	int j; 
	double w; 
	if(k == 0) 
		w = exp(u[0]); 
	else { 
		w = k *v[0]*u[k]; 
		for(j = 1; j < k; j++) 
			w += ( (k-j) * v[j] * u[k-j] );
		w /= k; 
	} 
	return w; 
} 

double pow_mc_c(double* u, double c, double* w, int k) 
{ 
	int j; 
	double ww; 
	if(k == 0) { 
		if(u[0] == 0.e0) { 
			printf(\"*** Function pow_mc _c found division by zero ***\");
			exit(EXIT_FAILURE); 
		} 
		ww = pow(u[0], c);  
	} else { 
		if(u[0] != 0.e0) {
			ww = c * k * w[0] * u[k]; 
			for(j = 1; j < k; j++) 
				ww += ( (c * (k - j) -j ) * w[j] * u[k-j]); 
			ww /= (k * u[0]); 
		} else ww = 0.e0;   
	} 
	return ww; 
} 

double log_mc(double* u, double* w, int k) 
{ 
	int j; 
	double ww; 
	if (k == 0) { 
		if(u[0] == 0.e0) { 
			printf(\"*** Function log_mc found division by zero ***\");
			exit(EXIT_FAILURE); 
		} 
		ww = log(u[0]); 
	} else { 
		ww=k*u[k]; 
		for(j = 1; j < k; j++) 
			ww -= ((k - j) * u[j] * w[k-j]); 
		ww /= (k * u[0]); 
	} 
	return ww; 
} 

double sin_mc(double* u, double* v, int k) 
{
	int j; 
	double ww; 
	if (k == 0) { 
		ww = sin(u[0]); 
	} else { 
		ww = u[1]*v[k-1];  
		for(j = 2; j <= k; j++) 
			ww += ( j * u[j] * v[k-j] ); 
		ww /= k; 
	} 
	return ww; 
} 

double cos_mc(double* u, double* v, int k) 

{ 

	int j; 
	double ww; 
	if(k == 0) { 
		ww = cos(u[0]); 
	} else { 
		ww = -u[1]*v[k-1]; 
		for(j = 2; j <= k; j++) 
			ww -= ( j * u [j] * v [k-j] ); 
		ww /=   k;   
	}     
	return ww;     
} 


    
";


(* ::Subsection::Closed:: *)
(*H*)


minhgen = 
"
/****************************************************************************
    
    minc_tides: kernel of the C Minimal version of TIDES
    This file is part of TIDES.

 
    Contributors:
 
    A. Abad, R. Barrio, F. Blesa, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

*****************************************************************************/

#ifndef minc_tides_HeadFile 
#define minc_tides_HeadFile 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

double mul_mc(double* u,   double* v,   int k);     
double div_mc(double* u,   double* v,   double* w,int k);   
double inv_mc(double p, double* u, double* w, int k); 
double exp_mc(double* u,   double* v,   int k);     
double pow_mc_c(double* u,   double e,   double* w,   int k);   
double log_mc(double* u,   double* w,   int k);     
double sin_mc(double* u,   double* v,   int k);   
double cos_mc(double* u,   double* v,   int k);

void	horner_mc( double *v, double t, int ORDER) ;
void	hornerd_mc(double *v, double t, int ORDER) ; 
void	dense_output_mc(double t0, double step, int ORDER);

double	norm_inf_vec_mv();
double	norm_inf_mat_mc(int ord);
void	declare_matrix_coefs_mc();
void	tolerances_mc(double *tol, double *tolo, int *ORDER);
double	steps_mc(int ORDER, double tol);
void	steps_DEC_mc(double t0, double tol, int ORDER, double *step);


void	mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO); 
void	minc_tides(double *var, int nvar, double *par, int npar,  double tini, double tend, double dt,
                double tol_rel, double tol_abs); 


static inline int min_i(int a, int b) { 
	return a < b ? a: b; 
} 

static  inline int max_i(int a, int b) { 
	return a > b ? a: b; 
} 

static  inline double min_d(double a, double b) { 
	return a < b ? a: b; 
} 

static  inline double max_d(double a, double b) { 
	return a > b ? a: b; 
} 


#endif 



";


(* ::Subsection::Closed:: *)
(*Fortran*)


minfgen = 
"
C****************************************************************************
C       
C     minf_tides: kernel of the Fortran Minimal version of TIDES
C     This file is part of TIDES.
C        
C   
C     Contributors:            
C   
C     A. Abad, R. Barrio, F. Blesa, M. Rodriguez
C     Grupo de Mecanica Espacial
C     University of Zaragoza
C     SPAIN
C       
C     http://gme.unizar.es/software/tides
C     Contact: <tides@unizar.es>
C        
C*****************************************************************************

      BLOCKDATA CONSTMETHOD
      REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
      INTEGER nitermax,nordinc,minord,maxord
      INTEGER accepted_steps, rejected_steps
      INTEGER IPOS
      LOGICAL dense_output, defect_error_control
      COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
      COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
      COMMON /OPT/ dense_output, defect_error_control
      COMMON /ARS/ accepted_steps, rejected_steps
      COMMON /DOPOS/IPOS
	  DATA fac1,fac2,fac3/0.95d0,10.d0,0.8d0/
	  DATA rminstep,rmaxstep/1.0d2,1.0d-2/
	  DATA nitermax,nordinc/5,5/
	  DATA minord,maxord/6,26/
      DATA dense_output, defect_error_control/.TRUE.,.FALSE./
      DATA accepted_steps, rejected_steps/0,0/
      DATA IPOS/1/
      END


      SUBROUTINE minf_tides(v,numvar,p,numpar,tini,tend,dt,
     &   tolrel,tolabs)

        IMPLICIT NONE
        LOGICAL dense_output, defect_error_control
        CHARACTER fname*20
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        INTEGER FL,IPOS
        INTEGER numvar,numpar
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /OPT/ dense_output, defect_error_control
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /DOPOS/IPOS
        COMMON /FILE/FL
        REAL*8 v(numvar),p(numpar),XVAR(0:maxord,0:numvar)
        REAL*8 t0,tini,tend,dt,tolrel,tolabs
        REAL*8 tol,tolo, step, temp, nstep, extra
        INTEGER ORDER
        INTEGER tflag

        TZERO = tini
        t0 = tini
        DELTAT = dt
        TOL_REL = tolrel
        TOL_ABS = tolabs
        NVAR = numvar
        NPAR = numpar
        extra = 0.d0

        IF (t0 .LT. tend) THEN
            tflag = 0
            IF (DELTAT .LT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        ELSE
            tflag = 1
            IF (DELTAT .GT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        END IF

        IF(dense_output) THEN
       		WRITE(FL,'(90E25.16)') TZERO, v
        END IF
           
        DO WHILE(((t0 .LT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .GT. tend).AND.(tflag.EQ.1)))
            CALL tolerances_mf(v, tol,tolo, ORDER)
            CALL minfseries(t0,v,NVAR,p,NPAR,XVAR,ORDER,maxord)
            CALL steps_mf(tol, XVAR, ORDER, step, tflag)

            IF (defect_error_control) THEN
                CALL steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)
            END IF
 
            temp = t0
            nstep = step + extra
            t0 = temp + nstep
            extra = (temp-t0)+nstep
            IF(((t0 .GT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .LT. tend).AND.(tflag.EQ.1))) THEN
                nstep  = (tend-temp)
            END IF
            accepted_steps =  accepted_steps + 1
            IF(dense_output) THEN
                CALL dense_output_mf(temp,nstep,XVAR,ORDER,tflag)
            END IF
            CALL horner_mf(v,XVAR,ORDER,nstep)
                
        END DO
           
        IPOS = 1

        RETURN
      END SUBROUTINE
           
C--------------------------------------------------------------
C--------------------------------------------------------------
C     tolerances_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE tolerances_mf(v, tol, tolo, ORDER)       
        IMPLICIT NONE
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        INTEGER ORDER
        REAL*8 yna,ynb,tol,tolo,v(NVAR),miny
           DATA ynb /0.0d0/
           SAVE ynb
          
        CALL norm_inf_vec_mf(v, yna)
        tol = TOL_ABS + MAX(yna,ynb)*TOL_REL
        miny = MIN(yna,ynb)
        IF(miny .gt. 0.d0) THEN
             tolo = MIN(TOL_ABS/miny, TOL_REL)
        ELSE
             tolo = MIN(TOL_ABS, TOL_REL)
        END IF


        ORDER = MIN(maxord, int(-log(tolo)/2)+nordinc)
        ORDER = MAX(minord,ORDER)
        
        ynb = yna
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_mf(tol, XVAR, ORDER, step, tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 XVAR(0:maxord,0:NVAR)
        REAL*8 tol, ynu,ynp,dord,dorda,dordp,sp,su
        INTEGER ORDER, ord, orda, ordp
        REAL*8 step, stepant, rstep
        DATA stepant /0.0d0/
        SAVE stepant
        INTEGER tflag
          
        ord = ORDER+1
        ynu = 0.d0
        DO WHILE((ynu .EQ. 0.d0) .AND. (ord .GT. 0))
           ord = ord-1
           CALL norm_inf_mat_mf(XVAR,ord, ynu)
        END DO
         
        IF(ord .NE. 0)  THEN
             orda = ORDER-1
             ordp = ORDER+1
             dord  = 1.d0/DFLOAT(ord)
             dorda = 1.d0/DFLOAT(orda)
             dordp = 1.d0/DFLOAT(ordp)
             CALL norm_inf_mat_mf(XVAR,orda, ynp)
             IF(ynp .eq. 0.d0) THEN
               step = (tol**dordp) * ((1.d00/ynu)**dord)
             ELSE
              sp = (tol**dord) * ((1.d00/ynp)**dorda)
              su = (tol**dordp) * ((1.d00/ynu)**dord)
              step = MIN(sp,su)
            END IF 
            IF(stepant. NE. 0.0d0) THEN
              rstep = step/stepant
              IF(rstep .GT. rmaxstep) THEN
                 step = rmaxstep*stepant
              ELSE IF(rstep .LT. rminstep) THEN
                 step = rminstep*stepant
              END IF
            END IF
            step = fac1*step
        ELSE
            WRITE(*,*) '*********Error*********'
            STOP
        END IF
        IF (tflag .EQ. 1) THEN
            step = -step
        END IF
        RETURN
      END SUBROUTINE

 

C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_DEC_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)      
        IMPLICIT NONE
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 vh(NVAR),vdh(NVAR), NORM
        REAL*8 t0,step,tolo, p(NPAR)
        REAL*8 XVAR(0:maxord,0:NVAR), XVAR2(0:maxord,0:NVAR)
        INTEGER ITER, ORDER,I
c
        ITER = 1
        NORM = 1.d99
           
        DO WHILE ((NORM .GT. fac2*tolo) .AND. (ITER .LT. nitermax))
             
            CALL horner_mf(vh,XVAR,ORDER,step)
            CALL hornerd_mf(vdh,XVAR,ORDER,step)
            CALL minfseries(t0+step,vh,NVAR,p,NPAR,XVAR2,1,1)
            NORM = ABS(XVAR2(1,1)-vdh(1))
            DO I=2, NVAR
                    NORM = NORM + ABS(XVAR2(1,i)-vdh(i))
            END DO
            IF(ITER .GT.1) THEN
                rejected_steps = rejected_steps+1
                step = fac3*step
            END IF
            ITER = ITER+1
               
        END DO
c
        RETURN
      END SUBROUTINE

         
C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_vec_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE norm_inf_vec_mf(v, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        REAL*8 norm
        REAL*8 v(NVAR)
        INTEGER I
        norm = 0.d0
        DO I =1,NVAR
             norm =MAX(ABS(v(I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_mat_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE norm_inf_mat_mf(XVAR,ord, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 norm
        REAL*8 XVAR(0:maxord,0:NVAR)
        INTEGER I,ord
        norm = 0.e0
        DO I =1,NVAR
             norm =MAX(abs(XVAR(ord,I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     dense_output_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE dense_output_mf(t0,step,XVAR,ORDER,tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER FL,IPOS
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /FILE/ FL
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /DOPOS/IPOS
        REAL*8 v(NVAR), XVAR(0:maxord,0:NVAR)
        REAL*8 t0,ti,tit,tend,step
        INTEGER  ORDER
        INTEGER tflag

        tend = t0 + step
        ti = TZERO + IPOS*DELTAT
        tit= ti    - t0
        DO WHILE(((tit .LE. step).AND.(tflag.EQ.0)) .OR.
     &           ((tit .GE. step).AND.(tflag.EQ.1)))
            CALL horner_mf(v,XVAR,ORDER,tit)
            WRITE(FL,'(90E25.16)') ti, v
            IPOS =IPOS + 1
            ti = TZERO + IPOS*DELTAT
            tit= ti    - t0
        END DO
        RETURN
      END SUBROUTINE
C--------------------------------------------------------------
C--------------------------------------------------------------
C     horner_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE horner_mf(v,XVAR,ORDER,t)       
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = XVAR(ORDER,nd)*t
          DO i = ORDER-1, 1, -1
            au = (au+XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(0,nd)
        END DO
        RETURN
      END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     hornerd_mf 
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE hornerd_mf(v,XVAR,ORDER,t)     
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = ORDER*XVAR(ORDER,nd)*t
          DO i = ORDER-1, 2, -1
            au = (au+i*XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(1,nd)
        END DO
        RETURN
       END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     mul_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION mul_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 mul_mf
        au = XX(0,n1)*XX(i,n2)
        DO m = 1, i
          au = au + XX(m,n1)*XX(i-m,n2)
        END DO
        mul_mf = au
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     div_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION div_mf(n1,n2,n3,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,n3,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 div_mf
        IF(XX(0,n2) .EQ. 0.d0)THEN
          WRITE(*,*) ' Function div_mf found division by zero'
          STOP
        ELSE
          au = XX(i,n1)
          DO m = 1, i
            au = au - XX(m,n2)*XX(i-m,n3)
          END DO
          div_mf = au/XX(0,n2)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     inv_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION inv_mf(p,n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au,p
        REAL*8 inv_mf
        inv_mf = 0.d0
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .EQ. 0.d0) THEN
            WRITE(*,*) ' Function inv_mf found division by zero'
            STOP
          ELSE
            inv_mf = p/XX(0,n1)
          END IF
        ELSE
          au = 0.d0
          DO m = 0, i -1
            au = au - (XX(i-m,n1)*XX(m,n2))
          END DO
          inv_mf = au*p/XX(0,n1)
        END IF
        RETURN
      END


C--------------------------------------------------------------
C--------------------------------------------------------------
C     pow_mf_c
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION pow_mf_c(n1,ex,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 ex,au
        REAL*8 pow_mf_c
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .EQ. 0.d0) THEN
            WRITE(*,*) ' Function pow_mf_c found division by zero'
            STOP
          ELSE
            pow_mf_c = XX(0,n1)**ex
          END IF
        ELSE
          IF(XX(0,n1) .NE. 0.d0) THEN
            au = ex*i*XX(0,n2)*XX(i,n1)
            DO m = 1, i - 1
              au = au+(ex*(i-m)-m)*XX(m,n2)*XX(i-m,n1)
            END DO
            pow_mf_c = au/(i*XX(0,n1))
          ELSE
            pow_mf_c = 0.D0
          END IF
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     exp_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION exp_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 exp_mf
        IF(i.EQ.0)THEN
          exp_mf = EXP(XX(0,n1))
        ELSE
          au = i*XX(0,n2)*XX(i,n1)
          DO m = 1, i - 1
            au = au+(i-m)*XX(m,n2)*XX(i-m,n1)
          END DO
          exp_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     log_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION log_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 log_mf
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .LE. 0.d0)THEN
            WRITE(*,*) 'Function log_mf found log of a negative value'
            STOP
          ELSE
            log_mf = LOG(XX(0,n1))
          END IF
        ELSE
          au = i*XX(i,n1)
          DO m = 1, i - 1
            au = au-(i-m)*XX(m,n1)*XX(i-m,n2)
          END DO
          log_mf = au/i/XX(0,n1)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     sin_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION sin_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 sin_mf
        IF(i.EQ.0)THEN
          sin_mf = SIN(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          sin_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     cos_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION cos_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 cos_mf
        IF(i.EQ.0)THEN
          cos_mf = COS(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          cos_mf = -au/i
        END IF
        RETURN
      END




";


(* ::Subsection::Closed:: *)
(*Header Standard DP*)


headstdDP = 
"
#ifndef Header_DP_TIDES_h
#define Header_DP_TIDES_h

#define real_Double

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*********************************************************************************************/

#define set_iterations() \\
set_iteration_parameters( \\
NUM_DERIVATIVES,VARIABLES,PARAMETERS,FUNCTIONS, \\
LINKS, PARTIALS_VARS, ORDER);\\
set_iteration_lists(POS__PARTIALS, POS_FUNCTIONS, \\
POS_ACCUM , POS_COEFS , POS_PREVI , POS_PREIV,\\
POS_ACCUM_S, POS_COEFS_S, POS_PREVI_S, POS_PREIV_S);


typedef long (*position_derivative)(char *der);

long position_variable(int v, position_derivative posder, char* der);
long position_function(int f, position_derivative posder, char* der);

int  is_variable(int num);

void set_iteration_parameters(long nvd, int v, int p, int f, int l, int prt, int ord);
void set_max_order(int ord);
void set_iteration_lists(int *prt, int *flst, 
	long *pra, long *prvi, long *priv, long *prc,
	long *prsa, long *prsvi, long *prsiv, long * prcs);

/*********************************************************************************************/


extern	int		MAX_ORDER;
extern	long	NDER;
extern	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
extern	int		*PARTIAL_LIST, *FUNCTION_LIST;
extern	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
extern	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;



void	varDB_init(double var[NVARS+1][NDER][MAX_ORDER+1], 
		double v[], double t);

void	parDB_init(double par[NPARS][NDER][MAX_ORDER+1], double p[]);

void	linkDB_init(double lk[NLINKS][NDER][MAX_ORDER+1]);

void	derDB_init(double var[NVARS+1][NDER][MAX_ORDER+1],
		double par[NPARS][NDER][MAX_ORDER+1], double v[]);

void	write_solution_DB(double cvf[][MAX_ORDER+1],
		double var[NVARS+1][NDER][MAX_ORDER+1], 
		double link[NLINKS][NDER][MAX_ORDER+1]);

void	double_htilde(double h[NDER][MAX_ORDER+1], long j, long v, long i, 
		double *ht, int ORDER_INDEX);

void	double_var_t(double f[NDER][MAX_ORDER+1], 
		double u[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_var_t_c(char* cs, double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_var_t_cc(double c, double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_add_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_add_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_add_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);


void	double_mul_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_mul_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_mul_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t_vc(double u[NDER][MAX_ORDER+1], double c, 
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t_cv(double c, double u[NDER][MAX_ORDER+1],
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_inv_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	double_exp_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	double_pow_t_c(double u[NDER][MAX_ORDER+1], char* cs, 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_pow_t_cc(double u[NDER][MAX_ORDER+1], double c, 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_sct_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i);
void	double_sct_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	double_sin_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_cos_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    double_sin_cos_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	double_sinh_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_cosh_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    double_sinh_cosh_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);

void	double_fgt_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], long i);
void	double_fgt_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	double_log_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	double_asin_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_acos_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_atan_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_asinh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_acosh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_atanh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);

/*********************************************************************************************/



void	double_init(double *rop); 
void	double_set_d(double *rop, double op); 
void	double_set_str(double *rop, char *op); 
void	double_set(double *rop, double op); 
double	double_get_d(double op); 
void	double_set_prec(int dig); 
void	double_clear (double op);
void	double_add(double *rop, double op1, double op2); 
void	double_sub(double *rop, double op1, double op2); 
void	double_mul(double *rop, double op1, double op2); 
void	double_div(double *rop, double op1, double op2); 
void	double_pow(double *rop, double op1, double op2); 
void	double_abs(double *rop, double op);  
void	double_add_i (double *rop, double op1, long   op2); 
void	double_sub_i (double *rop, double op1, long   op2); 
void	double_i_sub (double *rop, long   op1, double op2);
void	double_mul_i (double *rop, double op1, long   op2); 
void	double_div_i (double *rop, double op1, long  op2); 
void 	double_i_div (double *rop, long   op1, double op2); 
void 	double_pow_i (double *rop, double op1, long   op2);
void	double_i_pow (double *rop, unsigned long  op1, double op2);
int		double_greater(double op1, double op2); 
int		double_greaterequal(double op1, double op2);
int		double_less(double op1, double op2); 
int		double_lessequal(double op1, double op2);
int		double_equal(double op1, double op2); 
void	double_log(double *rop, double op); 
void	double_log10(double *rop, double op); 
void	double_exp(double *rop, double op); 
void	double_exp2(double *rop, double op); 
void	double_exp10(double *rop, double op); 
void	double_cos(double *rop, double op); 
void	double_sin(double *rop, double op);
void	double_sin_cos(double *rsin, double *rcos, double op); 
void	double_tan(double *rop, double op);
void	double_sec(double *rop, double op);
void	double_csc(double *rop, double op);
void	double_cot(double *rop, double op);
void	double_acos(double *rop, double op);
void	double_asin(double *rop, double op);
void	double_atan(double *rop, double op);
void	double_atan2(double *rop, double op1, double op2); 
void	double_cosh(double *rop, double op);
void	double_sinh(double *rop, double op);
void	double_tanh(double *rop, double op);
void	double_sech(double *rop, double op);
void	double_csch(double *rop, double op);
void	double_coth(double *rop, double op);
void	double_acosh(double *rop, double op);
void	double_asinh(double *rop, double op);
void	double_atanh(double *rop, double op);

typedef  double*	Array1DB;
typedef  double**	Array2DB;

void	Array1DB_init(Array1DB *vec, long dim); 
void	Array2DB_init(Array2DB *vec, long rows, long columns); 
void	Array3DB_init(Array2DB *vec, long dim, long rows, long columns); 
void	Array1DB_set(Array1DB rop, Array1DB op, long dim);
void	Array2DB_set(Array2DB rop, Array2DB op, long rows, long columns);
void	Array2DB_column_set(Array2DB rop, Array1DB op, long c, long dim);
void	Array2DB_row_set(Array2DB rop, Array1DB op, long r, long dim);

void	double_write (char *c, double op);

/*********************************************************************************************/
typedef double		realNUM;
typedef double* 	realVEC;
typedef double**	realMAT;


#define  realMAT_init		Array2DB_init
#define  realVEC_init		Array1DB_init
#define  variables_init		varDB_init
#define  parameters_init	parDB_init
#define  links_init			linkDB_init
#define  derivatives_init	derDB_init
#define  write_solution		write_solution_DB
#define  variables_free		varDB_free
#define  parameters_free	parDB_free
#define  links_free			linkDB_free
#define  set_precision_digits	double_set_prec
#define  var_t				double_var_t
#define  var_t_c			double_var_t_c
#define  var_t_cc		    double_var_t_cc
#define  add_t				double_add_t
#define  add_t_c			double_add_t_c
#define  add_t_cc			double_add_t_cc
#define  sub_t				double_sub_t
#define  sub_t_c			double_sub_t_c
#define  sub_t_cc			double_sub_t_cc
#define  mul_t				double_mul_t
#define  mul_t_c			double_mul_t_c
#define  mul_t_cc			double_mul_t_cc
#define  divide_t			double_div_t
#define  divide_t_cv		double_div_t_cv
#define  divide_t_vc		double_div_t_vc
#define  inv_t				double_inv_t
#define  exp_t				double_exp_t
#define  pow_t_c			double_pow_t_c
#define  pow_t_cc			double_pow_t_cc
#define  sin_t				double_sin_t
#define  sincos_t			double_sin_cos_t
#define  sincosh_t			double_sinh_cosh_t
#define  cos_t				double_cos_t
#define  sinh_t				double_sinh_t
#define  cosh_t				double_cosh_t
#define  asin_t				double_asin_t
#define  acos_t				double_acos_t
#define  atan_t				double_atan_t
#define  asinh_t			double_asinh_t
#define  acosh_t			double_acosh_t
#define  atanh_t			double_atanh_t
#define  log_t				double_log_t
#define tides_write			double_write

/*********************************************************************************************/



typedef long (*LinkedFunction)(double t, double v[], 
							   double p[], int orden, double cvfd[][orden+1]);

typedef void (*EvalFGFunction)(double t, double *x, double *p, double *f, double *grad);

void use_default_step_estimator ();
void set_info_taylor();
void unset_info_taylor();
void str_info_taylor();
void add_info_step(realNUM tstep);


int  taylor_order(double eps);
void norm_inf(double *rop, int n, int k, double coef[][k+1]);
void compute_step(double *rop, double tol, int n, int ord, double coef[][ord+1]);
void compute_step0(double *rop, double tol, int n, int ord, double coef[][ord+1]);
void compute_step1(double *rop, double tol, int n, int ord, double coef[][ord+1]);
void compute_step2(double *rop, double tol, int n, int ord, double coef[][ord+1]);

void compute_tol (double *tol, double tolrel, double tolabs, int n, int ord, double coef[][ord+1]);
void taylor_horner(int n, int ord, double coef[][ord+1], double t, double x[]);
void taylor_horner_der(int n, int ord, double coef[][ord+1], double t, double x[]);
void write_taylor_solution(int n, int k, int j, double tini, 
		double x[], double y[], double **mat, FILE* fileout);

void dp_tides(LinkedFunction fcn, 
		int nvar, int npar, int nfun, 
		double x[], double p[],
		double lt[], int ntes, 	
		double tolrel, double tolabs,   
		double **mat, FILE* fileout);

void dp_tides_point(LinkedFunction fcn, 
		int nvar, int npar, int nfun, 
		double x[], double p[],
		double t0, double tf, double dt, 	
		double tolrel, double tolabs,   
		double **mat, FILE* fileout);

void dp_tides_poc_end(LinkedFunction fcn, 
		int nvar, int npar, 
		double *x, double *p,
		double tini, double tend, 	
		double tolrel, double tolabs,   
		double* derini, double *derfin, double *partials) ;

void dp_tides_end(LinkedFunction fcn, 
		int nvar, int npar, 
		double *x, double *p,
		double tini, double tend, 	
		double tolrel, double tolabs,   
		double *partials) ;

void dp_tides_find_zeros(LinkedFunction fcn, 
		int nvar, int npar,double *x, double *p,
		double tini, double tend, double tol, 
		int *numevents, double **events,   
		double** mat, FILE* fileout) ;

void dp_tides_find_extrema(LinkedFunction fcn, 
		int nvar, int npar,double *x, double *p,
		double tini, double tend, double tol, 
		int *numevents, double **events,   
		double** mat, FILE* fileout) ;

void dp_tides_find_minimum(LinkedFunction fcn, 
		int nvar, int npar,double *x, double *p,
		double tini, double tend, double tol, 
		int *numevents, double **events,   
		double** mat, FILE* fileout) ;

void dp_tides_find_maximum(LinkedFunction fcn, 
		int nvar, int npar,double *x, double *p,
		double tini, double tend, double tol, 
		int *numevents, double **events,   
		double** mat, FILE* fileout);

void dp_tides_events(LinkedFunction fcn, 
		int nvar, int npar,double *x, double *p,
		double tini, double tend, double tol, 
		int *numevents, double **events,  int evcase, 
		double** mat, FILE* fileout) ;

int dp_tides_find_period(LinkedFunction fcn, 
		int nvar, int npar, double *x, double *p, double tini, double tend, 	
		double tol, double *period, double *distance);

void dp_distance(int nvar, double *xf, double *xi, double *distance);

#endif
";


(* ::Subsection::Closed:: *)
(*Header Standard MP*)


headstdMP = 
"
#ifndef Header_MP_TIDES_h
#define Header_MP_TIDES_h

#define real_MP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include \"mpfr.h\"


/*********************************************************************************************/

#define set_iterations() \\
set_iteration_parameters( \\
NUM_DERIVATIVES,VARIABLES,PARAMETERS,FUNCTIONS, \\
LINKS, PARTIALS_VARS, ORDER);\\
set_iteration_lists(POS__PARTIALS, POS_FUNCTIONS, \\
POS_ACCUM , POS_COEFS , POS_PREVI , POS_PREIV,\\
POS_ACCUM_S, POS_COEFS_S, POS_PREVI_S, POS_PREIV_S);


typedef long (*position_derivative)(char *der);

long position_variable(int v, position_derivative posder, char* der);
long position_function(int f, position_derivative posder, char* der);

int  is_variable(int num);

void set_iteration_parameters(long nvd, int v, int p, int f, int l, int prt, int ord);

void set_max_order(int ord);

void set_iteration_lists(int *prt, int *flst, 
		long *pra, long *prvi, long *priv, long *prc,  
		long *prsa, long *prsvi, long *prsiv, long * prcs);

/*********************************************************************************************/


extern	int		MAX_ORDER;
extern	long	NDER;
extern	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
extern	int		*PARTIAL_LIST, *FUNCTION_LIST;
extern	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
extern	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;


		
void	varMP_init(mpfr_t var[NVARS+1][NDER][MAX_ORDER+1], 
		mpfr_t v[], mpfr_t t);

void	parMP_init(mpfr_t par[NPARS][NDER][MAX_ORDER+1], mpfr_t p[]);

void	linkMP_init(mpfr_t lk[NLINKS][NDER][MAX_ORDER+1]);

void	derMP_init(mpfr_t var[NVARS+1][NDER][MAX_ORDER+1],
		mpfr_t par[NPARS][NDER][MAX_ORDER+1], mpfr_t v[]);

void 	clear (mpfr_t var[NVARS+1][NDER][MAX_ORDER+1], 
		mpfr_t par[NPARS][NDER][MAX_ORDER+1],
		mpfr_t link[NLINKS][NDER][MAX_ORDER+1]);

void	write_solution_MP(mpfr_t cvf[][MAX_ORDER+1],
		mpfr_t var[NVARS+1][NDER][MAX_ORDER+1], 
		mpfr_t link[NLINKS][NDER][MAX_ORDER+1]);

void	mpfrts_htilde(mpfr_t h[NDER][MAX_ORDER+1], long j, long v, long i, 
		mpfr_t *ht, int ORDER_INDEX);

void	mpfrts_var_t(mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t u[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_var_t_c(char* cs, mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_var_t_cc(mpfr_t c, mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	mpfrts_add_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_sub_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_add_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_sub_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_add_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1], 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_sub_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1], 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);


void	mpfrts_mul_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_mul_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_mul_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1], 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void    mpfrts_div_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
					 mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_div_t_vc(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t c, 
						mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_div_t_cv(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1],
						mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	mpfrts_inv_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	mpfrts_exp_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	mpfrts_pow_t_c(mpfr_t u[NDER][MAX_ORDER+1], char* cs, 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_pow_t_cc(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t c, 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	mpfrts_sct_0(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], long i);
void	mpfrts_sct_i(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	mpfrts_sin_t(mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_cos_t(mpfr_t c[NDER][MAX_ORDER+1], mpfr_t s[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    mpfrts_sin_cos_t (mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	mpfrts_sinh_t(mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_cosh_t(mpfr_t c[NDER][MAX_ORDER+1], mpfr_t s[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    mpfrts_sinh_cosh_t (mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);


void	mpfrts_fgt_0(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], long i);
void	mpfrts_fgt_i(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	mpfrts_log_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	mpfrts_asin_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_acos_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_atan_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_asinh_t(mpfr_t f[NDER][MAX_ORDER+1],mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_acosh_t(mpfr_t f[NDER][MAX_ORDER+1],mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_atanh_t(mpfr_t f[NDER][MAX_ORDER+1],mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX);


/*********************************************************************************************/

#define GMP_RND GMP_RNDN

int		binary_precision(int prec) ;

void	mpfrts_init (mpfr_t *rop); /* INITIALIZE ROP */
void	mpfrts_set_i(mpfr_t  *rop, long op); /* SET INTEGER NUMBER */
void	mpfrts_set_d (mpfr_t *rop, double op); /* SET DOUBLE NUMBER */
void	mpfrts_set_str (mpfr_t *rop, char *op); /* SET STRING AS NUMBER */
void	mpfrts_set (mpfr_t *rop, mpfr_t op); /* SET REALNUM AS NUMBER */
double	mpfrts_get_d (mpfr_t op); /* TRANSFORM TO DOUBLE */
long    mpfrts_get_i(mpfr_t op);

int		mpfrts_get_prec (); 
void 	mpfrts_set_prec (int dig);

void	mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 + OP2 */
void	mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 - OP2 */
void 	mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2); 
void 	mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 x OP2 */
void	mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 / OP2 */
void 	mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 ^ OP2 */
void	mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2);
void 	mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2);
void	mpfrts_abs(mpfr_t  *rop, mpfr_t op); /* ROP = |OP| */
void	mpfrts_neg(mpfr_t  *rop, mpfr_t op);

int	mpfrts_greater(mpfr_t op1, mpfr_t op2); /* 1 IF OP1 > OP2; 0 OTHERWHISE */
int	mpfrts_greaterequal(mpfr_t op1, mpfr_t op2);
int	mpfrts_less(mpfr_t op1, mpfr_t op2); 
int	mpfrts_lessequal(mpfr_t op1, mpfr_t op2);
int	mpfrts_equal(mpfr_t op1, mpfr_t op2); /* 1 IF OP1 == OP2; 0 OTHERWHISE */

void	mpfrts_sqrt(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log(mpfr_t  *rop, mpfr_t op); /* ROP = LOG (OP) */
void	mpfrts_log2(mpfr_t  *rop, mpfr_t op); /* ROP = LOG_2 (OP) */
void	mpfrts_log10(mpfr_t  *rop, mpfr_t op); /* ROP = LOG_1(OP) */
void	mpfrts_exp(mpfr_t  *rop, mpfr_t op); /* ROP = e ^ OP */
void	mpfrts_exp2(mpfr_t  *rop, mpfr_t op); /* ROP = 2 ^ OP */
void	mpfrts_exp10(mpfr_t  *rop, mpfr_t op); /* ROP = 10 ^ OP */
void	mpfrts_cos(mpfr_t  *rop, mpfr_t op); /* ROP = COS (OP) */
void	mpfrts_sin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op); /* RSIN = SIN (OP); RCOS = COS (OP) */
void	mpfrts_tan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sec(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csc(mpfr_t  *rop, mpfr_t op);
void	mpfrts_cot(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acos(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = ATAN2 (OP1, OP2), THE SAME AS IN DOUBLE */
void	mpfrts_cosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_tanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sech(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csch(mpfr_t  *rop, mpfr_t op);
void	mpfrts_coth(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atanh(mpfr_t  *rop, mpfr_t op);


void mpfrts_write_var(mpfr_t op);
void mpfrts_write (char *c, mpfr_t op);
void mpfrts_fread (FILE *file, mpfr_t rop);
void mpfrts_fwrite (FILE *file, mpfr_t op, int prec); 

typedef  mpfr_t*	Array1MP;
typedef  mpfr_t**	Array2MP;

void	Array1MP_init(Array1MP *vec, long dim);
void	Array2MP_init(Array2MP *vec, long rows, long columns);
void	Array1MP_set(Array1MP rop, Array1MP op, long dim);
void	Array2MP_set(Array2MP rop, Array2MP op, long rows, long columns);
void	Array2MP_column_set(Array2MP rop, Array1MP op, long c, long dim);
void	Array2MP_row_set(Array2MP rop, Array1MP op, long r, long dim);

/*********************************************************************************************/
typedef mpfr_t 		realNUM;
typedef mpfr_t* 	realVEC;
typedef mpfr_t**	realMAT;

#define  realMAT_init		Array2MP_init
#define	 realVEC_init		Array1MP_init
#define  variables_init		varMP_init
#define  parameters_init	parMP_init
#define  links_init			linkMP_init
#define  derivatives_init	derMP_init
#define  variables_free		varMP_free
#define  parameters_free	parMP_free
#define  links_free			linkMP_free
#define  write_solution		write_solution_MP
#define  set_precision_digits	mpfrts_set_prec
#define  var_t				mpfrts_var_t
#define  var_t_c			mpfrts_var_t_c
#define  var_t_cc			mpfrts_var_t_cc
#define  add_t				mpfrts_add_t
#define  add_t_c			mpfrts_add_t_c
#define  add_t_cc			mpfrts_add_t_cc
#define  sub_t				mpfrts_sub_t
#define  sub_t_c			mpfrts_sub_t_c
#define  sub_t_cc			mpfrts_sub_t_cc
#define  mul_t				mpfrts_mul_t
#define  mul_t_c			mpfrts_mul_t_c
#define  mul_t_cc			mpfrts_mul_t_cc
#define  divide_t			mpfrts_div_t
#define  divide_t_cv		mpfrts_div_t_cv
#define  divide_t_vc		mpfrts_div_t_vc
#define  inv_t				mpfrts_inv_t
#define  exp_t				mpfrts_exp_t
#define  pow_t_c			mpfrts_pow_t_c
#define  pow_t_cc			mpfrts_pow_t_cc
#define  sincos_t			mpfrts_sin_cos_t
#define  sincosh_t			mpfrts_sinh_cosh_t
#define  sin_t				mpfrts_sin_t
#define  cos_t				mpfrts_cos_t
#define  sinh_t				mpfrts_sinh_t
#define  cosh_t				mpfrts_cosh_t
#define  asin_t				mpfrts_asin_t
#define  acos_t				mpfrts_acos_t
#define  atan_t				mpfrts_atan_t
#define  asinh_t			mpfrts_asinh_t
#define  acosh_t			mpfrts_acosh_t
#define  atanh_t			mpfrts_atanh_t
#define  log_t				mpfrts_log_t


/*********************************************************************************************/

typedef long (*LinkedFunction)(mpfr_t t, mpfr_t v[], 
				mpfr_t p[], int orden, mpfr_t cvfd[][orden+1]);

typedef void (*EvalFGFunction)(mpfr_t t, mpfr_t *x, mpfr_t *p, mpfr_t *f, mpfr_t *grad);

int getOrder ();
int getNsteps ();

void mpfrts_use_default_step_estimator ();
void mpfrts_set_info_taylor();
void mpfrts_unset_info_taylor();
void mpfrts_str_info_taylor();
void mpfrts_add_info_step(mpfr_t tstep);

int mpfrts_taylor_order(mpfr_t eps);
void mpfrts_norm_inf(mpfr_t *rop, int n, int k, mpfr_t coef[][k+1]);
void mpfrts_compute_step(mpfr_t *rop, mpfr_t tol, int n, int ord, 
		mpfr_t coef[][ord+1]);
void mpfrts_compute_step0(mpfr_t *rop, mpfr_t tol, int n, int ord, 
		mpfr_t coef[][ord+1]);
void mpfrts_compute_step1(mpfr_t *rop, mpfr_t tol, int n, int ord, 
		mpfr_t coef[][ord+1]);
void mpfrts_compute_step2(mpfr_t *rop, mpfr_t tol, int n, int ord,
		mpfr_t coef[][ord+1]);
void mpfrts_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, 
		int n, int ord, mpfr_t coef[][ord+1]);
void mpfrts_taylor_horner(int n, int ord, mpfr_t coef[][ord+1],mpfr_t t, mpfr_t x[]);
void mpfrts_taylor_horner_der(int n, int ord, mpfr_t coef[][ord+1],mpfr_t t, mpfr_t x[]);
void mpfrts_write_taylor_solution( int n, int k, int j, mpfr_t tini, 
		mpfr_t x[], mpfr_t y[], mpfr_t** mat, FILE* fileout);
int mpfrts_valid_step (LinkedFunction fcn, mpfr_t *step, mpfr_t tip, 
		mpfr_t eps, int nvar, int ncol, int order,
		mpfr_t cvfd[][order+1], mpfr_t p[]);

void mp_tides(LinkedFunction fcn, 
		int nvar, int npar, int nfun, 
		mpfr_t x[], mpfr_t p[],
		mpfr_t lt[], int ntes, 	
		mpfr_t tolrel, mpfr_t tolabs, 
		mpfr_t** mat, FILE* fileout);

void mp_tides_point(LinkedFunction fcn, 
		int nvar, int npar, int nfun, 
		mpfr_t x[], mpfr_t p[],
		mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
		mpfr_t tolrel, mpfr_t tolabs,   
		mpfr_t** mat, FILE* fileout);

void mp_tides_poc_end(LinkedFunction fcn, 
		int nvar, int npar,  
		mpfr_t x[], mpfr_t p[],
		mpfr_t tini, mpfr_t tend, 	
		mpfr_t tolrel, mpfr_t tolabs, 
		mpfr_t* derini, mpfr_t *derfin, mpfr_t *partials) ;

void mp_tides_end(LinkedFunction fcn, 
		int nvar, int npar,  
		mpfr_t x[], mpfr_t p[],
		mpfr_t tini, mpfr_t tend, 	
		mpfr_t tolrel, mpfr_t tolabs, 
		mpfr_t *partials) ;

void mp_tides_find_zeros(LinkedFunction fcn, 
		int nvar, int npar,mpfr_t *x, mpfr_t *p,
		mpfr_t tini, mpfr_t tend, mpfr_t tol, 
		int *numevents, mpfr_t **events,   
		mpfr_t** mat, FILE* fileout) ;

void mp_tides_find_extrema(LinkedFunction fcn, 
		int nvar, int npar,mpfr_t *x, mpfr_t *p,
		mpfr_t tini, mpfr_t tend, mpfr_t tol, 
		int *numevents, mpfr_t **events,   
		mpfr_t** mat, FILE* fileout);

void mp_tides_find_minimum(LinkedFunction fcn, 
		int nvar, int npar,mpfr_t *x, mpfr_t *p,
		mpfr_t tini, mpfr_t tend, mpfr_t tol, 
		int *numevents, mpfr_t **events,   
		mpfr_t** mat, FILE* fileout) ;

void mp_tides_find_maximum(LinkedFunction fcn, 
		int nvar, int npar,mpfr_t *x, mpfr_t *p,
		mpfr_t tini, mpfr_t tend, mpfr_t tol, 
		int *numevents, mpfr_t **events,   
		mpfr_t** mat, FILE* fileout);

void mp_tides_events(LinkedFunction fcn, 
		int nvar, int npar,mpfr_t *x, mpfr_t *p,
		mpfr_t tini, mpfr_t tend, mpfr_t tol, 
		int *numevents, mpfr_t **events,  int evcase, 
		mpfr_t** mat, FILE* fileout) ;

int mp_tides_find_period(LinkedFunction fcn, 
		int nvar, int npar, mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, 	
		mpfr_t tol, mpfr_t *period, mpfr_t *distance);

void mp_distance(int nvar, mpfr_t *xf, mpfr_t *xi, mpfr_t *distance);

#endif
";


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`Texts`"]

EndPackage[]

Null

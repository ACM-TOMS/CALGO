/****************************************************************************
  This file is part of TIDES.
 
 Contributors:
 
 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 *****************************************************************************/

#include "minc_tides.h" 


double    *v, *p, **XVAR, **XVAR2, TOL_ABS, TOL_REL, tzero, deltat; 
int       VARS, PARS, MAX_ORDER, tflag;
double    fac1=0.95e0,fac2=10.e0,fac3=0.8e0 ; 
double    rmaxstep=1.e2,rminstep=1.e-2; 
int       nitermax=5, nordinc=5, minord=6, maxord=26; 
int       dense_output = -1, defect_error_control = 0;
int       accepted_steps = 0, rejected_steps = 0;
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
		fprintf(fd, "%25.15le ", tini ); 
		for(i = 0; i < VARS; i++) fprintf(fd, "%25.15le " , v[i]); 
		fprintf(fd, "\n"); 	
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


} 

/************************************************************************/ 
void dense_output_mc(double t0, double step, int ORDER)
{
	static int ipos = 1;
	int i;
	double tend, ti, tit, vh[VARS]; 
	tend = t0 + step;
	ti = tzero + ipos*deltat;
	tit = ti - t0;
	while(((tit <= step) && (tflag == 0)) ||
		   ((tit >= step) && (tflag == 1))  ) {
		horner_mc(vh, tit, ORDER); 
		fprintf(fd, "%25.15le ", ti ); 
		for(i = 0; i < VARS; i++) fprintf(fd, "%25.15le " , vh[i]); 
		fprintf(fd, "\n");
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
		printf("*********Error*********");
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
		printf("*** Function div_mc found division by zero ***");
		exit(EXIT_FAILURE); 
	} 
	ww = u[k]; 
	for(j = 1; j <= k; j++) ww -= (v[j] *w[k-j]);
	ww /= v[0]; 
	return ww; 
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
			printf("*** Function pow_mc _c found division by zero ***");
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
			printf("*** Function log_mc found division by zero ***");
			exit(EXIT_FAILURE); 
		} 
		ww = log(u[0]); 
	} else { 
		ww=k*u[k]; 
		for(j = 1; j < k; j++) 
			ww += ((k - j) * u[j] * w[k-j]); 
		ww /= (k * u[0]); 
	} 
	return ww; 
} 

double sin_mc(double* u, double* v, double* w, int k) 
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

double cos_mc(double* u, double* v, double* w, int k) 

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


    

/****************************************************************************
 libTIDES. Version 1.3.0.
 This file is part of TIDES.
 
 Contributors:
 
 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 *****************************************************************************/

#include "taylorODE.h"

int		_info_steps_taylor_ = 0;
double	etapa_minima;
double	etapa_maxima;
double	etapa_total;
int		num_etapas = 0;
int		order_series = 0;
int 	clearPartials = 1;
int		order_estimator = 0;

static double hant = 0.;

double 	fac1 = 0.95;
double	fac2 = 10.;
double	fac3 = 0.8;
double	rmaxstep = 1.e2;
double	rminstep = 1.e-2;
double	nitermax = 5;
int		nordinc = 4;
int		minord  = 6;
int		defect_error_control = 0;



void use_default_step_estimator () {
	defect_error_control = 1;
}
	
void set_info_taylor() {  
	_info_steps_taylor_ = 1;

	etapa_minima = 1.e30;
	etapa_maxima = 0.;
	etapa_total = 0.;
	num_etapas = 0; 
	order_series = 0;
}

void unset_info_taylor() {  
	_info_steps_taylor_ = 0;
	etapa_minima = 1.e30;
	etapa_maxima = 0.;
	etapa_total = 0.;
	num_etapas = 0; 
	order_series = 0;
}

void str_info_taylor() {
	double tot, num, average;

	tot = etapa_total;
	num = num_etapas;
	average = tot / num;

	if(_info_steps_taylor_ == 1) {
		printf("============================================================\n");
		printf ("Number of integration steps in Taylor method: %d\n", num_etapas);
		tides_write ("Minimum   integration step:", etapa_minima);
		tides_write ("Maximum   integration step:", etapa_maxima);
		tides_write ("Averaged  integration step:",average);
		printf("Order of the Taylor series  : %d\n",  order_series);
		printf("============================================================\n\n");
	}
}


void add_info_step(double tstep) {
	if (tstep < etapa_minima) 
		etapa_minima = tstep;
	if (tstep > etapa_maxima)
		etapa_maxima = tstep;
	num_etapas++;
	etapa_total = etapa_total + tstep;
}


int taylor_order(double tol) {
	double rop;
	int		nrop;
	rop = nordinc - log(tol)/2.0;
	nrop = (int) ceil (rop);
	if(nrop < minord) nrop = minord;
	return nrop;
}

void norm_inf(double *rop, int n, int k, double coef[][k+1]) {
	int i;
	double ck[n], ninf;

	ninf = -1.;
	for (i=0; i<n; i++) 
		ck[i] = fabs(coef[i][k]); 
	
	for (i=0; i<n; i++) 
		if (ck[i] > ninf)
			ninf = ck[i];
	tides_set (rop, ninf);
	
}

void compute_step(double *rop, realNUM tol, int n, int ord, double coef[][ord+1]) {
	double tult, tpen, nor, aux;
	int i;
	norm_inf (&nor, n, ord, coef);

	if (nor == 0.) tult = 0.;
	else tult = pow(tol, 1./(ord+1)) * pow(nor, -1./ord);
	
	norm_inf (&nor, n, ord-1, coef);
	if (nor == 0.) tpen = 0.;
	else tpen = pow(tol, 1./ord) * pow(nor, -1./(ord-1));

	if ((tult == 0.) && (tpen == 0.)) {
		i = ord -1;
		*rop = 0.;

		while ((i > 0) && (nor == 0.)) {
			norm_inf (&nor, n, i-1, coef);
			tides_set_d (&aux, 1./i);
			tides_pow (&aux, tol, aux);
			tides_set_d (rop, -1./(i-1));
			tides_pow (rop, tol, *rop);
			i--;
		}
	} else if (tult == 0.) *rop = tpen;
	else if (tpen == 0.) *rop = tult;
	else if (tult < tpen) *rop = tult;
	else *rop = tpen;

	*rop = (*rop) * fac1;
	if (hant != 0.) {
		if ((*rop)/hant > rmaxstep) *rop = hant * rmaxstep;
		else if ((*rop)/hant < rminstep) *rop = hant * rminstep;	
	}
	hant = *rop;
	
	if(_info_steps_taylor_ == 1) add_info_step(*rop);
}


void compute_tol (double *tol, double tolrel, double tolabs, int n, int ord, double coef[][ord+1]) {
	double norm_y0, norm_y1, max_norm;

	norm_inf(&norm_y0, n, 0, coef);
	norm_inf(&norm_y1, n, 1, coef);

	if (norm_y0 > norm_y1) max_norm = norm_y0;
	else max_norm = norm_y1;

	*tol = tolrel * max_norm;

	if (*tol > tolabs) *tol = tolabs;

}



void taylor_horner (int n, int ord, double coef[][ord+1], double t, double x[]) {
	int i,j;
	if(n > 0 ) {
		double total;
		for(i = 0; i < n; i++ ){
			total = 0.;
			for (j=ord; j>=0; j--) 
				total = total*t + coef[i][j];
			tides_set (&x[i], total);
		}
	}
}

void taylor_horner_der (int n, int ord, double coef[][ord+1], double t, double x[]) 
{
	int i,j;
	double total;
	for(i = 0; i < n; i++ ){
		total = 0.;
		for (j=ord; j>0; j--) 
			total = total*t + j*coef[i][j];
		x[i] = total;
	}
}


void write_taylor_solution( int n, int j, double tini, 
	double x[], realNUM** mat, FILE* fileout) {

	int i;
	char formato[] = "%.16le ";

	if(fileout != NULL) {
		fprintf (fileout, formato, tini);
		for(i = 0; i < n; i++) {
			if(x[i] >= 0) fprintf(fileout, " ");
			fprintf (fileout, formato, x[i]);
		}
		fprintf (fileout, "\n");
	}
	if(mat != NULL) {
		mat[j][0] = tini;
		for(i = 0; i < n; i++)
			mat[j][i+1] = x[i];
	}
}




int valid_step (LinkedFunction fcn, double *step, double tip, 
		double tol, int nvar, int ncol, int order, double cvfd[][order+1],
		 double p[]) {
	int i, j;
	int accepted = 1;
	double b[nvar][order], y[nvar], yp[nvar], t, cn[ncol][order+1];
	double dif[nvar], nor, aux, up;


	nor = -1.;
	up  = fac2*tol; 

	for (i=0; i<nvar; i++) for (j=0; j<=order-1; j++)
		b[i][j] = cvfd[i][j+1] * (j+1);

	t = tip + *step;
	taylor_horner (nvar, order, cvfd, *step, y);
	taylor_horner (nvar, order-1, b, *step, yp);
	fcn (t, y, p, 1, cn);

	for (i=0; i<nvar; i++) {
		tides_sub (&dif[i], yp[i], cn[i][1]);
		aux = fabs(dif[i]);
		if (tides_greater (aux, nor)) tides_set (&nor, aux);
	}
	accepted = 0;
	for (i=0; i<nitermax; i++) 
		if (nor > up) {
			*step = *step * fac3;
		} else {
			accepted = 1;
			return accepted;
		}
	return accepted;
}


void dp_tides_point(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	double x[], double p[],
	double t0, double tf, double dt, 	
	double tolrel, double tolabs,   
	double** mat, FILE* fileout) {

	int ntot = (int) floor ((tf - t0)/dt) + 1;
	int i;
	double lt[ntot];
	
	for (i=0; i<ntot; i++) tides_init (&lt[i]);

		
	for(i = 0; i < ntot; i++) 
		lt[i] = t0 + i * dt;

	dp_tides(fcn, nvar, npar, nfun, x, p, lt, ntot, tolrel, tolabs, 
			mat, fileout);
}


void dp_tides(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	double x[], double p[],
	double lt[], int ntes, 	
	double tolrel, double tolabs,   
	double** mat, FILE* fileout) {

	int order, i, ilt, signo;
	long ncol;

	double tini, tfin, tip, tipant, tstep, ststep, dlt, eps;
	double utipant, utip, utend, absdlt;

	double kahan_c = 0., kahan_t, kahan_y;


	if(_info_steps_taylor_ == 1) set_info_taylor();
	else unset_info_taylor();
	
	if(ntes < 0) {
		printf("Error: ntes debe ser mayor o igual que 2"); 
	} else {
		tini = lt[0]; tip = tini;
		tfin = lt[ntes-1];

		order = taylor_order (tolabs); 
		order_series = order;

		if (lt[1] > lt[0]) signo = 1;
		else signo = -1;
		ncol = fcn (0., NULL, NULL, -1, NULL);
		double cvfd[ncol][order+1];
		double y[ncol], yn[ncol];

		fcn (tip, x, p, order, cvfd);

		for (i=0; i<nvar; i++) {
			yn[i] = x[i];
		}
		
		taylor_horner (ncol, order, cvfd, 0., y);

		write_taylor_solution(ncol,  0, tini, y, mat, fileout);
		

		compute_tol (&eps, tolrel, tolabs, ncol, order, cvfd);
		compute_step(&tstep, eps, nvar, order, cvfd);
		ststep = tstep * signo;
		if (defect_error_control)
			valid_step (fcn, &ststep, tip, eps, nvar, ncol, order, cvfd, p);

		taylor_horner(ncol, order, cvfd, ststep, y);

		for (i=0; i<nvar; i++) {
			x[i] = y[i];
		}
		tipant = tip;
		tip = tip + ststep;
		ilt = 1;
		utipant = tipant * signo;
		utend = tfin * signo;
		while (utipant < utend) {
			if(ilt < ntes) {
				dlt = lt[ilt] - tipant;
				absdlt = fabs(dlt);
				while (ilt<ntes && (absdlt < tstep)) {
					taylor_horner(ncol,order,cvfd,dlt,yn);
					write_taylor_solution(ncol,ilt,lt[ilt],yn, mat, fileout);
					ilt++;
					dlt = lt[ilt] - tipant;
					absdlt = fabs(dlt);
				}
			} 
			utip = tip * signo;
			if (utip < utend) {
				fcn (tip, y, p, order, cvfd);


				compute_tol (&eps, tolrel, tolabs, ncol, 
					order, cvfd);
				compute_step(&tstep, eps, nvar, order, cvfd);
				
				ststep = tstep * signo;
				
				if (defect_error_control)
					valid_step (fcn, &ststep, tip, eps, nvar, ncol, order, cvfd, p);
				taylor_horner(ncol, order, cvfd, ststep, y);


				for (i=0; i<nvar; i++) 
					x[i] = y[i];
			}
			tipant = tip;
			utipant = tipant * signo;

			kahan_y = ststep - kahan_c;
			kahan_t = tip + kahan_y;
			kahan_c = (kahan_t - tip) - kahan_y;
			tip = kahan_t;
		}
		if(_info_steps_taylor_==1)  str_info_taylor(); 
	clearPartials = 1;
	for (i=0; i<nvar; i++) x[i] = yn[i];

	}

}

void dp_tides_poaux(LinkedFunction fcn, 
					int nvar, int npar, 
					double *x, double *p,
					double tini, double tend, 	
					double tolrel, double tolabs,   
					double* derini, double *derfin, double *partials) 
{
	
	int order, i, first = 1;
	long ncol;
	
	double tip, tipant, tstep, dlt, eps;	
	double kahan_c = 0., kahan_t, kahan_y;
	
	order = taylor_order (tolabs); 
	
	ncol = fcn (0., NULL, NULL, -1, NULL);
	double cvfd[ncol][order+1];
	double y[ncol];
	
	tip = tini;
	for(i=0; i<nvar; i++) y[i] = x[i];
	
	
	while (tip < tend) {
		
		fcn (tip, y, p, order, cvfd);
		if(first) {
			for(i=0; i<nvar; i++) derini[i] = cvfd[i][1];
			first = 0;
		}
		compute_tol (&eps, tolrel, tolabs, ncol, order, cvfd);
		compute_step(&tstep, eps, nvar, order, cvfd);
		if (defect_error_control)
			valid_step (fcn, &tstep, tip, eps, nvar, ncol, order, cvfd, p);
		taylor_horner(ncol, order, cvfd, tstep, y);		
		
		dlt = tend - tip;
		if (dlt < tstep) {
			taylor_horner(ncol,order,cvfd,dlt,y);
			taylor_horner_der(nvar,order,cvfd,dlt,derfin);
			for(i=0; i<nvar; i++) x[i] = y[i];
			for(i=nvar; i<ncol; i++) partials[i-nvar] = y[i];
			clearPartials = 1;
			return; 
		}
		
		
		tipant = tip;
		
		kahan_y = tstep - kahan_c;
		kahan_t = tip + kahan_y;
		kahan_c = (kahan_t - tip) - kahan_y;
		tip = kahan_t;
	}
	
}


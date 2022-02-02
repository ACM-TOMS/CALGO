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

#include "mpfrNUM.h"
#include "mpfrtaylorODE.h"


extern int	_info_steps_taylor_;
mpfr_t	etapa_minima;
mpfr_t	etapa_maxima;
mpfr_t	etapa_total;
extern int	num_etapas;
extern int	order_series;
extern int 	clearPartials;

extern int		order_estimator;
extern double 	fac1;
extern double	fac2;
extern double	fac3;
extern double	rmaxstep;
extern double	rminstep;
extern double	nitermax;
extern int		nordinc;
extern int		minord;
extern int		defect_error_control;

static double hant = 0.;


	
void mpfrts_set_info_taylor() {  
	_info_steps_taylor_ = 1;
	mpfrts_init (&etapa_minima);
	mpfrts_init (&etapa_maxima);
	mpfrts_init (&etapa_total);

	mpfrts_set_str (&etapa_minima, "1.e30");
	mpfrts_set_str (&etapa_maxima, "0");
	mpfrts_set_str (&etapa_total, "0.");
	num_etapas = 0; 
	order_series = 0;
}

void mpfrts_unset_info_taylor() { 
	_info_steps_taylor_ = 0;
	mpfrts_init (&etapa_minima);
	mpfrts_init (&etapa_maxima);
	mpfrts_init (&etapa_total);
	mpfrts_set_str (&etapa_minima, "1.e30");
	mpfrts_set_str (&etapa_maxima, "0");
	mpfrts_set_str (&etapa_total, "0.");
	num_etapas = 0; 
	order_series = 0;
}

void mpfrts_str_info_taylor() {
	double aux;
	mpfr_t tot, num, average;

	mpfrts_init (&tot); mpfrts_init (&num); mpfrts_init (&average);
	mpfrts_set (&tot, etapa_total);
	aux = num_etapas;
	mpfrts_set_d (&num, aux);
	mpfrts_div (&average, tot, num);

	if(_info_steps_taylor_ == 1) {
		printf("============================================================\n");
		printf ("Number of integration steps in Taylor method: %d\n", num_etapas);
		mpfrts_write ("Minimum   integration step:", etapa_minima);
		mpfrts_write ("Maximum   integration step:", etapa_maxima);
		mpfrts_write ("Averaged  integration step:",average);
		printf("Order of the Taylor series  : %d\n",  order_series);
		printf("============================================================\n\n");
	}
}


void mpfrts_add_info_step(mpfr_t tstep) {
	if (mpfrts_less (tstep, etapa_minima)) 
		mpfrts_set (&etapa_minima, tstep);
	if (mpfrts_greater (tstep, etapa_maxima)) 
		mpfrts_set (&etapa_maxima, tstep);
	num_etapas++;
	mpfrts_add (&etapa_total, etapa_total, tstep);
}


int mpfrts_taylor_order(mpfr_t eps) {
	mpfr_t rop, coc;
	double aux, nrop;

	mpfrts_init (&rop);
	mpfrts_log (&rop, eps);
	mpfrts_init (&coc);
	mpfrts_set_str(&coc, "2.0");
	mpfrts_div (&rop, rop, coc);
	mpfrts_i_sub (&rop, nordinc, rop);
	aux = mpfrts_get_d (rop);
	mpfr_clear (rop);
	mpfr_free_cache ();
	nrop = (int) ceil (aux);
	if(nrop < minord) nrop = minord;
	return nrop;
}


void mpfrts_norm_inf(mpfr_t *rop, int n, int k, mpfr_t coef[][k+1]) {
	int i;
	mpfr_t ck[n], ninf;

	mpfrts_init (&ninf); mpfrts_set_str (&ninf, "-1.e10000");
	for (i=0; i<n; i++) {
		mpfrts_init (&ck[i]);
		mpfrts_abs (&ck[i], coef[i][k]); 
	}
	for (i=0; i<n; i++) 
		if (mpfrts_greater (ck[i], ninf)) 
			mpfrts_set (&ninf, ck[i]);
	mpfrts_set (rop, ninf);
	for (i=0; i<n; i++)
		mpfr_clear (ck[i]);
	mpfr_clear (ninf);
	mpfr_free_cache ();
	
}

void mpfrts_compute_step(mpfr_t *rop, mpfr_t tol, int n, int ord,
		mpfr_t coef[][ord+1]) {
	mpfr_t tult, tpen, nor, zero, aux;
	int i;
	mpfrts_init (&zero); mpfrts_init (&aux);
	mpfrts_init (&tult); mpfrts_init (&tpen);
	mpfrts_init (&nor);
	mpfrts_norm_inf (&nor, n, ord, coef);
	if(mpfrts_equal(nor, zero)) mpfrts_set_d (&tult, 0.);
	else {
		mpfrts_set_d (&aux, 1./(ord+1));
		mpfrts_pow (&aux, tol, aux);
		mpfrts_set_d (&tult, -1./ord);
		mpfrts_pow (&tult, nor, tult);
		mpfrts_mul (&tult, aux, tult);	
	}


	mpfrts_norm_inf (&nor, n, ord-1, coef);
	if(mpfrts_equal(nor, zero)) mpfrts_set_d (&tpen, 0.);
	else {
		mpfrts_set_d (&aux, 1./(ord));
		mpfrts_pow (&aux, tol, aux);
		mpfrts_set_d (&tpen, -1./(ord-1));
		mpfrts_pow (&tpen, nor, tpen);
		mpfrts_mul (&tpen, aux, tpen);	
	} 

	if (mpfrts_equal (tult, zero) && (mpfrts_equal (tpen, zero))) {
		i = ord -1;
		mpfrts_set_d (rop, 0);

		while ((i > 0) && mpfrts_equal (nor, zero)) {
			mpfrts_norm_inf (&nor, n, i-1, coef);
			mpfrts_set_d (&aux, 1./i);
			mpfrts_pow (&aux, tol, aux);
			mpfrts_set_d (rop, -1./(i-1));
			mpfrts_pow (rop, tol, *rop);
			i--;
		}
	} else if (mpfrts_equal (tult, zero)) mpfrts_set (rop, tpen);
	else if (mpfrts_equal (tpen, zero)) mpfrts_set (rop, tult);
	else {
		if (mpfrts_less (tult, tpen)) mpfrts_set (rop, tult);
		else mpfrts_set (rop, tpen);
	}
	mpfr_t mp_fac; mpfrts_init (&mp_fac); mpfrts_set_d (&mp_fac, fac1);
	mpfrts_mul (rop, *rop, mp_fac);

	mpfr_t mp_hant; mpfrts_init (&mp_hant); mpfrts_set_d (&mp_hant, hant);	
	mpfr_t mp_rmax; mpfrts_init (&mp_rmax); 
		mpfrts_set_d (&mp_rmax, rmaxstep);	
	mpfr_t mp_rmin; mpfrts_init (&mp_rmin); 
		mpfrts_set_d (&mp_rmin, rminstep);	

	if (hant != 0.) {
		mpfrts_div (&aux, *rop, mp_hant);
		if (mpfrts_greater(aux, mp_rmax))
			mpfrts_mul (rop, mp_hant, mp_rmax);
		if (mpfrts_less (aux, mp_rmin))
			mpfrts_mul (rop, mp_hant, mp_rmin);

	}

	hant = mpfrts_get_d (*rop);

	if(_info_steps_taylor_ == 1) mpfrts_add_info_step(*rop);
	mpfr_clear (tult);
	mpfr_clear (tpen);
	mpfr_clear (nor);
	mpfr_clear (zero);
	mpfr_clear (aux);
	mpfr_clear (mp_hant);
	mpfr_clear (mp_rmax);
	mpfr_clear (mp_rmin);
	mpfr_clear (mp_fac);
	mpfr_free_cache ();	
}



void mpfrts_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int n, int ord, mpfr_t coef[][ord+1]) {
	mpfr_t norm_y0, norm_y1, max_norm;

	mpfrts_init (&norm_y0); mpfrts_init (&norm_y1); mpfrts_init (&max_norm);


	mpfrts_norm_inf(&norm_y0, n, 0, coef);
	mpfrts_norm_inf(&norm_y1, n, 1, coef);

	if (mpfrts_greater (norm_y0, norm_y1)) 
		mpfrts_set (&max_norm, norm_y0);
	else 
		mpfrts_set (&max_norm, norm_y1);


	mpfrts_mul (tol, tolrel, max_norm);
	if (mpfrts_greater (*tol, tolabs)) mpfrts_set (tol, tolabs);
	
	mpfr_clear (norm_y0);
	mpfr_clear (norm_y1);
	mpfr_clear (max_norm);

	mpfr_free_cache ();

}


void mpfrts_taylor_horner (int n, int ord, mpfr_t coef[][ord+1], 
		mpfr_t t, mpfr_t x[]) {
	int i, j;
	if(n > 0 ) {
		mpfr_t total, aux;
		mpfrts_init (&total); mpfrts_init (&aux);
		for(i = 0; i < n; i++ ){
			mpfrts_set_d (&total,0.);
			for (j=ord; j>=0; j--) {
				mpfrts_mul (&aux, t, total);
				mpfrts_add (&total, aux, coef[i][j]);
			}
			mpfrts_set (&x[i], total);
		}
		mpfr_clear (total);
		mpfr_clear (aux);
		mpfr_free_cache ();
	}
}
void mpfrts_taylor_horner_der(int n, int ord, mpfr_t coef[][ord+1], 
						   mpfr_t t, mpfr_t x[]) {
	int i, j;
	if(n > 0 ) {
		mpfr_t total, aux1, aux2;
		mpfrts_init (&total); 
		mpfrts_init (&aux1);
		mpfrts_init (&aux2);
		for(i = 0; i < n; i++ ){
			mpfrts_set_d (&total,0.);
			for (j=ord; j>0; j--) {
				mpfrts_mul (&aux1, t, total);
				mpfrts_mul_i (&aux2, coef[i][j], j);
				mpfrts_add (&total, aux1, aux2);
			}
			mpfrts_set (&x[i], total);
		}
		mpfr_clear (total);
		mpfr_clear (aux1);
		mpfr_clear (aux2);
		mpfr_free_cache ();
	}
}


void mpfrts_write_taylor_solution (int n, int j, mpfr_t tini, 
	mpfr_t x[], mpfr_t** mat, FILE* fileout) {

	int i, precision;
	mpfr_t zero;
	mpfrts_init (&zero);
	precision = mpfrts_get_prec ();
	if (fileout != NULL) {
		mpfr_out_str (fileout, 10, precision, tini, GMP_RND);
		for (i=0; i<n; i++){ 
			fprintf (fileout, "   ");
			if (mpfrts_greaterequal (x[i], zero))
				fprintf (fileout, " ");
			mpfr_out_str (fileout, 10, precision, x[i], GMP_RND);
		}
		fprintf (fileout, "\n");
	}
	if (mat != NULL) {
		mpfrts_set (&mat[j][0], tini);

		for (i=0; i<n; i++)
			mpfrts_set (&mat[j][i+1], x[i]);
	}
	mpfr_clear (zero);
	mpfr_free_cache ();	
}




int mpfrts_valid_step (LinkedFunction fcn, mpfr_t *step, mpfr_t tip, 
		mpfr_t tol, int nvar, int ncol, int order, 
		mpfr_t cvfd[][order], mpfr_t p[]) {
	int i, j, k;
	int accepted = 1;
	mpfr_t b[nvar][order], y[nvar], yp[nvar], t, cn[ncol][order+1],
		dif[nvar], nor, aux, up;

	mpfrts_init (&t); mpfrts_init (&nor); mpfrts_init (&aux);
	mpfrts_init (&up); 

	mpfrts_set_str (&nor, "-1e10000");
	mpfrts_set_d (&up, fac2); mpfrts_mul (&up, up, tol);
	for (i=0; i<nvar; i++) {
		mpfrts_init (&y[i]);
		mpfrts_init (&dif[i]);
		mpfrts_init (&yp[i]);
		for (j=0; j<order; j++)
			mpfrts_init (&b[i][j]);
	}
	for (i=0; i<ncol; i++) for (j=0; j<=order; j++)
		mpfrts_init (&cn[i][j]);
	for (i=0; i<nvar; i++) for (j=0; j<order; j++) mpfrts_init (&b[i][j]);
	for (i=0; i<nvar; i++) mpfrts_init (&y[i]);
	for (i=0; i<nvar; i++) mpfrts_init (&dif[i]);
	for (i=0; i<nvar; i++) mpfrts_init (&yp[i]);
	for (i=0; i<ncol; i++) for(j=0; j<=order; j++) mpfrts_init(&cn[i][j]);

	for (i=0; i<nvar; i++) for (j=0; j<=order-1; j++)
		mpfrts_mul_i (&b[i][j], cvfd[i][j+1], j+1);

	mpfrts_add (&t, tip, *step);
	mpfrts_taylor_horner (nvar, order, cvfd, *step, y);
	mpfrts_taylor_horner (nvar, order-1, b, *step, yp);
	fcn (t, y, p, 1, cn);
	for (i=0; i<nvar; i++) {
		mpfrts_sub (&dif[i], yp[i], cn[i][1]);
		mpfrts_abs (&aux, dif[i]);
		if (mpfrts_greater (aux, nor)) mpfrts_set (&nor, aux);
	}
	accepted = 0;
	mpfr_t fac; mpfrts_init (&fac); mpfrts_set_d (&fac, fac3);
	for (k=0; k<nitermax; k++)
		if (mpfrts_greater (nor, up)){
			mpfrts_mul (step, *step, fac);
		}
		else {
			accepted = 1;
			mpfr_clear (t);
			mpfr_clear (nor);
			mpfr_clear (aux);
			mpfr_clear (up);
			mpfr_clear (fac);
			for (i=0; i<nvar; i++) {
				for (j=0; j<order; j++) 
					mpfr_clear (b[i][j]);
				mpfr_clear (y[i]);
				mpfr_clear (yp[i]);
				mpfr_clear (dif[i]);
			}
			for (i=0; i<ncol; i++) {
				for (j=0; j<=order; j++)
					mpfr_clear (cn[i][j]);
			}

			mpfr_free_cache ();
			return accepted;
		}
	
	mpfr_clear (t);
	mpfr_clear (nor);
	mpfr_clear (aux);
	mpfr_clear (up);
	mpfr_clear (fac);
	for (i=0; i<nvar; i++) {
		for (j=0; j<order; j++) 
			mpfr_clear (b[i][j]);
		mpfr_clear (y[i]);
		mpfr_clear (yp[i]);
		mpfr_clear (dif[i]);
	}
	for (i=0; i<ncol; i++) {
		for (j=0; j<=order; j++)
			mpfr_clear (cn[i][j]);
	}

	mpfr_free_cache ();	
	return accepted;
}


void mp_tides_point(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	mpfr_t x[], mpfr_t p[],
	mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
	mpfr_t tolrel, mpfr_t tolabs,   
	mpfr_t** mat, FILE* fileout) {

	mpfr_t aux; mpfrts_init (&aux); mpfrts_sub (&aux, tf, t0);
		mpfrts_div (&aux, aux, dt);
	mpfr_floor (aux, aux);
	int ntot = mpfr_get_si (aux, GMP_RNDN)+1;
	int i;
	mpfr_t lt[ntot];
	for (i=0; i<ntot; i++) mpfrts_init (&lt[i]);

	for(i = 0; i < ntot; i++) {
		//lt[i] = t0 + i * inctime;
		mpfrts_mul_i (&lt[i], dt, i);
		mpfrts_add (&lt[i], t0, lt[i]);
	}
	mp_tides(fcn,nvar,npar,nfun,x,p,lt,ntot,tolrel, tolabs,
			mat,fileout);
}


void mp_tides(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	mpfr_t x[], mpfr_t p[],
	mpfr_t lt[], int ntes, 	
	mpfr_t tolrel, mpfr_t tolabs,   
	mpfr_t** mat, FILE* fileout) {

	int order, i,  j, ilt, signo;
	long ncol;

	mpfr_t tini, tfin, tip, tipant, tstep, ststep, dlt, zero, eps;
	mpfr_t utipant, utip, utend, absdlt;

	mpfrts_init (&tini); mpfrts_init (&tfin); mpfrts_init (&tip);
	mpfrts_init (&tipant); mpfrts_init (&tstep); mpfrts_init (&ststep);
	mpfrts_init (&dlt); mpfrts_init (&zero);
	mpfrts_init (&utipant); mpfrts_init (&utend); mpfrts_init (&utip);
	mpfrts_init (&absdlt); mpfrts_init (&eps);

	mpfr_t kahan_c, kahan_y, kahan_t;
		mpfrts_init (&kahan_c); mpfrts_init (&kahan_y);
		mpfrts_init (&kahan_t);
		mpfrts_set_str (&kahan_c, "0.");


	if(_info_steps_taylor_ == 1) mpfrts_set_info_taylor();
	else mpfrts_unset_info_taylor();
		
	if(ntes < 2) {
		printf("Error: ntes debe ser mayor o igual que 2"); 
	} else {
		mpfrts_set (&tini, lt[0]); mpfrts_set (&tip, tini);
		mpfrts_set (&tfin, lt[ntes-1]);
		order = mpfrts_taylor_order (tolabs);
		order_series = order;
		
		if (mpfrts_greater(lt[1], lt[0])) signo = 1;
		else signo = -1;
		ncol = fcn (zero, NULL, NULL, -1, NULL);
		mpfr_t cvfd[ncol][order+1];
		mpfr_t y[ncol], yn[ncol];
		
		for (i=0; i<ncol; i++) {
			mpfrts_init (&y[i]);
			mpfrts_init (&yn[i]);
			for (j=0; j<=order; j++) mpfrts_init (&cvfd[i][j]);
		}

		
		fcn (tip, x, p, order, cvfd);

	
		for (i=0; i<nvar; i++) {
			mpfrts_set (&yn[i], x[i]);
		}
		mpfrts_taylor_horner (ncol, order, cvfd, zero, y);
		mpfrts_write_taylor_solution(ncol, 0, tini, y, mat, fileout);

		mpfrts_compute_tol (&eps, tolrel, tolabs, ncol, order, cvfd);
		mpfrts_compute_step(&tstep, eps, nvar, order, cvfd);
		mpfrts_mul_i (&ststep, tstep, signo);			

		if (defect_error_control)
			mpfrts_valid_step (fcn, &ststep, tip, eps, nvar, ncol, order, cvfd, p);
		mpfrts_taylor_horner(ncol, order, cvfd, ststep, y);

		mpfrts_set (&tipant, tip);
		mpfrts_add (&tip, tip, ststep);
		ilt = 1;
		mpfrts_mul_i (&utipant, tipant, signo);
		mpfrts_mul_i (&utend, tfin, signo);

		while (mpfrts_less (utipant, utend)) {
			if(ilt < ntes) {
				mpfrts_sub (&dlt, lt[ilt], tipant);
				mpfrts_abs (&absdlt, dlt);
				while (ilt<ntes && mpfrts_less(absdlt, tstep)) {
					mpfrts_taylor_horner(ncol,order,cvfd,dlt,yn);
					mpfrts_write_taylor_solution(ncol,ilt,lt[ilt],yn, mat, fileout);
					ilt++;
					if (ilt > ntes-1) {
						break;
					}
					mpfrts_sub (&dlt, lt[ilt], tipant);
					mpfrts_abs (&absdlt, dlt);
					mpfr_free_cache ();				}
			}
			mpfrts_mul_i (&utip, tip, signo);
			if (mpfrts_less (utip, utend)) {
				fcn (tip, y, p, order, cvfd);
		
				mpfrts_compute_tol (&eps, tolrel, tolabs, 
					ncol, order, cvfd);
				mpfrts_compute_step(&tstep, eps, nvar, order, cvfd);
				mpfrts_mul_i (&ststep, tstep, signo);
					if (defect_error_control)
						mpfrts_valid_step (fcn, &ststep, tip, eps, nvar, ncol, order, cvfd, p);
				mpfrts_taylor_horner(ncol, order, cvfd, ststep, y);

				for (i=0; i<nvar; i++) {	
					mpfrts_set (&x[i], y[i]);						}
			}

			mpfrts_set (&tipant, tip);
			mpfrts_mul_i (&utipant, tipant, signo);

			mpfrts_sub (&kahan_y, ststep, kahan_c);
			mpfrts_add (&kahan_t, tip, kahan_y);
			mpfrts_sub (&kahan_c, kahan_t,tip);
				mpfrts_sub (&kahan_c, kahan_c, kahan_y);
			mpfrts_set (&tip, kahan_t);



			mpfr_free_cache ();	
		}
		if(_info_steps_taylor_==1)  mpfrts_str_info_taylor();
		for (i=0; i<nvar; i++) {	
			mpfrts_set (&x[i], yn[i]);						}

	clearPartials=1;
	}
}

void mp_tides_poaux(LinkedFunction fcn, 
					int nvar, int npar, 
					mpfr_t x[], mpfr_t p[],
					mpfr_t tini, mpfr_t tend,  	
					mpfr_t tolrel, mpfr_t tolabs,   
					mpfr_t *derini, mpfr_t *derfin, mpfr_t *partials)
{
	
	int order, i,j, first = 1;
	long ncol;
	
	mpfr_t tip, tipant, tstep, dlt, eps,zero;
	
	mpfrts_init (&tip); mpfrts_init (&tipant);
	mpfrts_init (&tstep); mpfrts_init (&dlt);
	mpfrts_init (&eps); 
	mpfrts_init (&zero); mpfrts_set_str (&zero, "0.");
	
	mpfr_t kahan_c, kahan_y, kahan_t;
	mpfrts_init (&kahan_c); mpfrts_init (&kahan_y);
	mpfrts_init (&kahan_t);
	mpfrts_set_str (&kahan_c, "0.");
	
	order = mpfrts_taylor_order (tolabs);
	
	ncol = fcn (zero, NULL, NULL, -1, NULL);
	mpfr_t cvfd[ncol][order+1];
	mpfr_t y[ncol];
	
	for (i=0; i<ncol; i++) {
		mpfrts_init (&y[i]);
		for (j=0; j<=order; j++) mpfrts_init (&cvfd[i][j]);
	}
	
	mpfrts_set (&tip, tini);
	for (i=0; i<nvar; i++) mpfrts_set (&y[i], x[i]);
	
	
	while (mpfrts_less (tip, tend) ) {
		
		fcn (tip, y, p, order, cvfd);
		if(first) {
			for(i=0; i<nvar; i++) mpfrts_set (&derini[i], cvfd[i][1]);
			first = 0;
		}
		mpfrts_compute_tol (&eps, tolrel, tolabs, ncol, order, cvfd);
		mpfrts_compute_step(&tstep, eps, nvar, order, cvfd);
		if (defect_error_control)
			mpfrts_valid_step (fcn, &tstep, tip, eps, nvar, ncol, order, cvfd, p);
		mpfrts_taylor_horner(ncol, order, cvfd, tstep, y);		
		
		mpfrts_sub (&dlt, tend, tip);
		if (mpfrts_less (dlt, tstep)) {
			mpfrts_taylor_horner(ncol,order,cvfd,dlt,y);
			mpfrts_taylor_horner_der(nvar,order,cvfd,dlt,derfin);
			for(i=0; i<nvar; i++) mpfrts_set (&x[i], y[i]);  
			for(i=nvar; i<ncol; i++) mpfrts_set (&partials[i-nvar], y[i]); 
			clearPartials=1;
			return; 
		}
		
		mpfrts_set (&tipant, tip);
		
		mpfrts_sub (&kahan_y, tstep, kahan_c);
		mpfrts_add (&kahan_t, tip, kahan_y);
		mpfrts_sub (&kahan_c, kahan_t,tip);
		mpfrts_sub (&kahan_c, kahan_c, kahan_y);
		mpfrts_set (&tip, kahan_t);
	}
}


int getNsteps () {return num_etapas;}
int getOrder () {return order_series;}

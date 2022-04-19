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

#include "commonITER.h"
#include "mpfrITER.h"


void varMP_init(mpfr_t var[NVARS+1][NDER][MAX_ORDER+1],  mpfr_t v[], mpfr_t t)
{
	int i, j, k;
	for (i=0; i<=NVARS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfrts_init (&var[i][j][k]);
	mpfrts_set (&var[0][0][0], t);
	mpfrts_set_str (&var[0][0][1], "1.");
	for(i=1; i<=NVARS; i++) 
		mpfrts_set (&var[i][0][0], v[i-1]);
}
void parMP_init(mpfr_t par[NPARS][NDER][MAX_ORDER+1], mpfr_t p[])
{
	int i, j, k;
	for (i=0; i<NPARS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfrts_init (&par[i][j][k]);
	for(i=0; i<NPARS; i++) mpfrts_set (&par[i][0][0], p[i]);
	for(i=0; i <NPARTIALS; i++)
		if(!is_variable(PARTIAL_LIST[i]))
			mpfrts_set_d (&par[PARTIAL_LIST[i]-NVARS-1][i+1][0], 1);	
}
void linkMP_init(mpfr_t lk[NLINKS][NDER][MAX_ORDER+1])
{
	int i, j, k;
	for (i=0; i<NLINKS; i++) 
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfrts_init (&lk[i][j][k]);
}


void derMP_init(mpfr_t var[NVARS+1][NDER][MAX_ORDER+1],
		mpfr_t par[NPARS][NDER][MAX_ORDER+1], mpfr_t v[])
{
	int i;
	extern int clearPartials;
	if (clearPartials)
	for(i=0; i <NPARTIALS; i++) {
		if(is_variable(PARTIAL_LIST[i])) 
			mpfrts_set_d (&var[PARTIAL_LIST[i]][i+1][0], 1.);
	}
	else {
		int j, nvf;
		nvf = NVARS+NFUNS;
		for( i = 0; i < NVARS; i++) 
			for(j = 1; j < NDER; j++) 
				mpfrts_set (&var[i+1][j][0], v[i+(j*nvf)]);
	}



	clearPartials = 0;
}


void clear (mpfr_t var[NVARS+1][NDER][MAX_ORDER+1], 
		mpfr_t par[NPARS][NDER][MAX_ORDER+1],
		mpfr_t link[NLINKS][NDER][MAX_ORDER+1]) {

	int i, j, k;
	for (i=0; i<NVARS+1; i++)
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfr_clear (var[i][j][k]);
	for (i=0; i<NPARS; i++)
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfr_clear (par[i][j][k]);
	for (i=0; i<NLINKS; i++)
		for (j=0; j<NDER; j++)
			for (k=0; k<=MAX_ORDER; k++)
				mpfr_clear (link[i][j][k]);

}


void	write_solution_MP(mpfr_t cvf[][MAX_ORDER+1],
		mpfr_t var[NVARS+1][NDER][MAX_ORDER+1],
		mpfr_t link[NLINKS][NDER][MAX_ORDER+1])
{
	int i,j,k, nvf;
	nvf = NVARS+NFUNS;
	mpfr_t vfun;
	mpfrts_init(&vfun);
	for( i = 0; i < NVARS; i++) 
		for(j = 0; j < NDER; j++) 
			for(k = 0; k <= MAX_ORDER; k++)
				mpfrts_set (&cvf[i+(j*nvf)][k], var[i+1][j][k]);
	for( i = 0; i < NFUNS; i++) 
		for(j = 0; j < NDER; j++) 
			for(k = 0; k <= MAX_ORDER; k++) {
				if (FUNCTION_LIST[i] > 0) mpfrts_set(&vfun, link[FUNCTION_LIST[i]][j][k]);
				else mpfrts_set(&vfun, var[-FUNCTION_LIST[i]][j][k]);
				mpfrts_set (&cvf[NVARS+i+(j*nvf)][k], vfun);
			}
	mpfr_clear (vfun);
}

/************************************************************************/

void mpfrts_htilde(mpfr_t h[NDER][MAX_ORDER+1], long j, long v, long i, mpfr_t *ht, int ORDER_INDEX)
{
	mpfr_t cero;
	mpfrts_init (&cero);
	if(ORDER_INDEX >= 0) {
		if( j == ORDER_INDEX && i==v) mpfrts_set(ht, cero);
		else mpfrts_set(ht, h[v][j]);
	}
	mpfr_clear (cero);
	mpfr_free_cache ();
}

/************************************************************************/

void	mpfrts_var_t(mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	mpfr_t val; mpfrts_init (&val);
	int i;
	if(ORDER_INDEX > 0) {
		for(i = 0; i < NDER; i++) {
			mpfrts_div_i(&val, f[i][ORDER_INDEX-1], ORDER_INDEX);
			mpfrts_set(&w[i][ORDER_INDEX], val);
		}
	}
	mpfr_clear (val);
	mpfr_free_cache ();

}

void	mpfrts_var_t_c(char* cs, mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	mpfr_t c;
	if(ORDER_INDEX == 1) {
		mpfrts_init(&c);
		mpfrts_set_str(&c, cs);
		mpfrts_set(&w[0][1], c);
		mpfr_clear (c);
	}
}
void	mpfrts_var_t_cc(mpfr_t c, mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		mpfrts_set(&w[0][1], c);
	}
}


void	mpfrts_add_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1],
		 mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			mpfrts_add(&w[i][ORDER_INDEX], u[i][ORDER_INDEX], v[i][ORDER_INDEX]);
	}
}

void	mpfrts_sub_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1],
		 mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			mpfrts_sub(&w[i][ORDER_INDEX], u[i][ORDER_INDEX], v[i][ORDER_INDEX]);
	}
}

void	mpfrts_add_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			mpfrts_add(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			mpfrts_set(&w[0][ORDER_INDEX], u[0][ORDER_INDEX]);
		for(i = 1; i < NDER; i++) 
			mpfrts_set(&w[i][ORDER_INDEX], u[i][ORDER_INDEX]);
	}
	mpfr_clear (c);
	mpfr_free_cache ();
}
void	mpfrts_add_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1], 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			mpfrts_add(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			mpfrts_set(&w[0][ORDER_INDEX], u[0][ORDER_INDEX]);
		for(i = 1; i < NDER; i++) 
			mpfrts_set(&w[i][ORDER_INDEX], u[i][ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void	mpfrts_sub_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			mpfrts_sub(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			mpfrts_mul_i(&w[0][ORDER_INDEX], u[0][ORDER_INDEX],-1);
		for(i = 1; i < NDER; i++) 
			mpfrts_mul_i(&w[i][ORDER_INDEX], u[i][ORDER_INDEX],-1);
	}
	mpfr_clear (c);
	mpfr_free_cache ();
}
void	mpfrts_sub_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1], 
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			mpfrts_sub(&w[0][ORDER_INDEX], c, u[0][ORDER_INDEX]);
		else
			mpfrts_mul_i(&w[0][ORDER_INDEX], u[0][ORDER_INDEX],-1);
		for(i = 1; i < NDER; i++) 
			mpfrts_mul_i(&w[i][ORDER_INDEX], u[i][ORDER_INDEX],-1);
	}
	mpfr_free_cache ();
}


void	mpfrts_mul_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1],
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, vi, j;
	mpfr_t sumint, sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);

		for(i = 0; i < NDER; i++) {
			mpfrts_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {

					mpfrts_mul(&partial, u[PREV_VI[vi]][ORDER_INDEX-j], v[PREV_IV[vi]][j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_set(&w[i][ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_free_cache ();

}


void	mpfrts_mul_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1],
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		for(i = 0; i < NDER; i++) 
			mpfrts_mul(&w[i][ORDER_INDEX], c, u[i][ORDER_INDEX]);
	}
	mpfr_clear (c);
	mpfr_free_cache ();

}
void	mpfrts_mul_t_cc(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1],
					   mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			mpfrts_mul(&w[i][ORDER_INDEX], c, u[i][ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void mpfrts_div_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, ht;
	if(g[0][0] == 0 ) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_str (&g[0][0], "1");
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);
		mpfrts_init(&ht);
		for(i = 0; i < NDER; i++) {
			mpfrts_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(h, ORDER_INDEX-j, PREV_VI[vi], i, &ht, ORDER_INDEX);
					mpfrts_mul(&partial, ht, g[PREV_IV[vi]][j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f[i][ORDER_INDEX], sumext);
			mpfrts_div(&sumext, sumext, g[0][0]);
			mpfrts_set(&h[i][ORDER_INDEX], sumext);
		}
		mpfr_clear (sumint);
		mpfr_clear (sumext);
		mpfr_clear (partial);
		mpfr_clear (zero);
		mpfr_clear (ht);
	}
	
}

void	mpfrts_div_t_vc(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t c, 
						mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i;
	mpfr_t ci;
	mpfrts_init(&ci);
	mpfrts_i_div(&ci, 1, c);
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < NDER; i++) 
			mpfrts_mul(&w[i][ORDER_INDEX], ci, u[i][ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void	mpfrts_div_t_cv(mpfr_t c, mpfr_t u[NDER][MAX_ORDER+1],
						mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, f, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&u[0][0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&f);
		mpfrts_init(&wt);
		for(i = 0; i < NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) mpfrts_set_d(&f, 1.);
			else mpfrts_set(&f, zero);
			mpfrts_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(w, ORDER_INDEX-j, PREV_VI[vi], i, &wt, ORDER_INDEX);
					mpfrts_div(&wt, wt, c);
					mpfrts_mul(&partial, wt, u[PREV_IV[vi]][j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f, sumext);
			mpfrts_div(&sumext, sumext, u[0][0]);
			mpfrts_mul(&sumext, sumext, c);
			mpfrts_set(&w[i][ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_clear (f);
	mpfr_clear (wt);
	mpfr_free_cache ();
	
}


void mpfrts_inv_t(mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, f, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&u[0][0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&f);
		mpfrts_init(&wt);
		for(i = 0; i < NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) mpfrts_set_d(&f, 1.);
			else mpfrts_set(&f, zero);
			mpfrts_set(&sumext, zero);
			for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(w, ORDER_INDEX-j, PREV_VI[vi], i, &wt, ORDER_INDEX);
					mpfrts_mul(&partial, wt, u[PREV_IV[vi]][j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f, sumext);
			mpfrts_div(&sumext, sumext, u[0][0]);
			mpfrts_set(&w[i][ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_clear (f);
	mpfr_clear (wt);
	mpfr_free_cache ();

}

void mpfrts_exp_t(mpfr_t u[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint,sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_exp(&sumext, u[0][0]);
					mpfrts_set(&w[0][0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						mpfrts_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j < ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[PREV_VI[vi]][j], u[PREV_IV[vi]][ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_mul_i(&sumint, sumint, ORDER_INDEX-j);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
	}
}


void mpfrts_pow_t_c(mpfr_t u[NDER][MAX_ORDER+1], char* cs, 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, partialb, zero, wt, c;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&u[0][0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&partialb);
		mpfrts_init(&wt);
		mpfrts_init(&c);
		mpfrts_set_str(&c, cs);
		mpfrts_set (&sumext, zero);
	
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_pow(&sumext, u[0][0], c);
					mpfrts_set(&w[0][0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						mpfrts_mul(&partial, partial, c);
						mpfrts_htilde(w, 0, PREVSTAR_IV[vi], i, &wt, 0);
						mpfrts_mul(&partialb, wt, u[PREVSTAR_VI[vi]][0]);
						mpfrts_sub(&partial, partial, partialb);
						mpfrts_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_div(&sumint, sumint, u[0][0]);
					mpfrts_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						mpfrts_htilde(w, j, PREV_VI[vi], i, &wt, ORDER_INDEX);
						mpfrts_mul(&partial, wt, u[PREV_IV[vi]][ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_add_i(&partial, c, 1);
					mpfrts_mul_i(&partial, partial, j);
					mpfrts_mul_i(&partialb, c, ORDER_INDEX);
					mpfrts_sub(&partial, partialb, partial);
					mpfrts_mul(&sumint, sumint, partial);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div(&sumext, sumext, u[0][0]);
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
		
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (partialb);
	mpfr_clear (zero);
	mpfr_clear (c);
	mpfr_clear (wt);
	mpfr_free_cache ();

	
	
}

void mpfrts_pow_t_cc(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t c, 
					mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, partialb, zero, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&u[0][0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&partialb);
		mpfrts_init(&wt);
		mpfrts_set (&sumext, zero);
		
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_pow(&sumext, u[0][0], c);
					mpfrts_set(&w[0][0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[PREVSTAR_VI[vi]][0], u[PREVSTAR_IV[vi]][0]);
						mpfrts_mul(&partial, partial, c);
						mpfrts_htilde(w, 0, PREVSTAR_IV[vi], i, &wt, 0);
						mpfrts_mul(&partialb, wt, u[PREVSTAR_VI[vi]][0]);
						mpfrts_sub(&partial, partial, partialb);
						mpfrts_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_div(&sumint, sumint, u[0][0]);
					mpfrts_set(&w[i][0], sumint);
				}
			}
		} else {		
			for(i = 0; i < NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
						mpfrts_htilde(w, j, PREV_VI[vi], i, &wt, ORDER_INDEX);
						mpfrts_mul(&partial, wt, u[PREV_IV[vi]][ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_add_i(&partial, c, 1);
					mpfrts_mul_i(&partial, partial, j);
					mpfrts_mul_i(&partialb, c, ORDER_INDEX);
					mpfrts_sub(&partial, partialb, partial);
					mpfrts_mul(&sumint, sumint, partial);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div(&sumext, sumext, u[0][0]);
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i][ORDER_INDEX], sumext);
			}
			
		}
		
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (partialb);
	mpfr_clear (zero);
	mpfr_clear (wt);
	mpfr_free_cache ();
}

/************************************************************************/

void mpfrts_sct_0(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
		mpfr_t h[NDER][MAX_ORDER+1], long i)
{
	long vi;
	mpfr_t sumint, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_set(&sumint, zero);					
	for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
		mpfrts_mul(&partial, g[PREVSTAR_VI[vi]][0], f[PREVSTAR_IV[vi]][0]);
		mpfrts_mul_i(&partial, partial, PREVSTAR_COEF[vi]);		
		mpfrts_add(&sumint, sumint, partial);
	}
	mpfrts_set(&h[i][0], sumint);
	mpfr_clear (sumint);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_free_cache ();
}

void mpfrts_sct_i(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX)
{
	long vi,j;
	mpfr_t sumint, sumext, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_init(&sumext);
	mpfrts_set(&sumext, zero);
	for(j = 1; j <= ORDER_INDEX; j++) {
		mpfrts_set(&sumint, zero);
		for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
			mpfrts_mul(&partial, g[PREV_VI[vi]][ORDER_INDEX-j], f[PREV_IV[vi]][j]);
			mpfrts_mul_i(&partial, partial, PREV_COEF[vi]);		
			mpfrts_add(&sumint, sumint, partial);
		}
		mpfrts_mul_i(&sumint, sumint, j);		
		mpfrts_add(&sumext, sumext, sumint);
	}
	mpfrts_div_i(&sumext, sumext, ORDER_INDEX);		
	mpfrts_set(&h[i][ORDER_INDEX], sumext);
	mpfr_clear (partial);
	mpfr_clear (sumint);
	mpfr_clear (zero);
	mpfr_clear (sumext);
}

void    mpfrts_sin_cos_t (mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_sin (&s[0][0], f[0][0]);
					mpfrts_cos (&c[0][0], f[0][0]);
				} else {
					mpfrts_sct_0 (f,c,s,i);
					mpfrts_sct_0 (f,s,c,i);
					mpfrts_mul_i (&c[i][0], c[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(f,s,c,i,ORDER_INDEX);
				mpfrts_sct_i(f,c,s,i,ORDER_INDEX);
				mpfrts_mul_i (&c[i][ORDER_INDEX],
					 c[i][ORDER_INDEX], -1);
			}
		}
	}
}


void mpfrts_sin_t(mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_sin(&w[0][0], s[0][0]);
				} else {
					
					mpfrts_sct_0(s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(s,c,w,i,ORDER_INDEX);
			}
		}
	}
}

void mpfrts_cos_t(mpfr_t c[NDER][MAX_ORDER+1], mpfr_t s[NDER][MAX_ORDER+1],
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_cos(&w[0][0], c[0][0]);
				} else {
					mpfrts_sct_0(c,s,w,i);
					mpfrts_mul_i(&w[i][0], w[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(c,s,w,i,ORDER_INDEX);
				mpfrts_mul_i(&w[i][ORDER_INDEX], w[i][ORDER_INDEX], -1);
			}
		}
	}
}


void    mpfrts_sinh_cosh_t (mpfr_t f[NDER][MAX_ORDER+1], 
		mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_sinh (&s[0][0], f[0][0]);
					mpfrts_cosh (&c[0][0], f[0][0]);
				} else {
					mpfrts_sct_0 (f,c,s,i);
					mpfrts_sct_0 (f,s,c,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(f,s,c,i,ORDER_INDEX);
				mpfrts_sct_i(f,c,s,i,ORDER_INDEX);
			}
		}
	}
}



void mpfrts_sinh_t(mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_sinh(&w[0][0], s[0][0]);
				} else {
					mpfrts_sct_0(s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(s,c,w,i,ORDER_INDEX);
			}
		}
	}
}
void mpfrts_cosh_t(mpfr_t c[NDER][MAX_ORDER+1], mpfr_t s[NDER][MAX_ORDER+1], 
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_cosh(&w[0][0], c[0][0]);
				} else {
					mpfrts_sct_0(c,s,w,i);
					mpfrts_mul_i(&w[i][0], w[i][0], -1);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_sct_i(c,s,w,i,ORDER_INDEX);
			}
		}
	}
}

/***************************************************************/

void mpfrts_fgt_0(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], long i)
{
	long vi;
	mpfr_t sumint, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_set(&sumint, zero);					
	if(mpfrts_equal (g[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&g[0][0], 1.);
	}
	
	for( vi = PREVSTAR_ACCUM[i]; vi < PREVSTAR_ACCUM[i+1]; vi++) {
		if(PREVSTAR_VI[vi]>0) {
			mpfrts_mul(&partial, h[PREVSTAR_IV[vi]][0], g[PREVSTAR_VI[vi]][0]);
			mpfrts_mul_i(&partial, partial, PREVSTAR_COEF[vi]);
			mpfrts_add(&sumint, sumint, partial);
		}
	}
	mpfrts_sub(&sumint, f[i][0], sumint);
	mpfrts_div(&sumint, sumint, g[0][0]);
	mpfrts_set(&h[i][0], sumint);
	mpfr_clear (zero);
	mpfr_clear (partial);
	mpfr_clear (sumint);
}
void mpfrts_fgt_i(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX)
{
	long vi,j;
	mpfr_t sumint, sumext, partial, zero, ht;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_init(&sumext);
	mpfrts_init(&ht);
	mpfrts_set(&sumext, zero);
	if(mpfrts_equal (g[0][0], zero)) {
		printf("**********  Divide by cero  ***********\n");
		printf("**********   bad   result   ***********\n");
		mpfrts_set_d (&g[0][0], 1.);
	}
	for(j = 0; j < ORDER_INDEX; j++) {
		mpfrts_set(&sumint, zero);
		for( vi = PREV_ACCUM[i]; vi < PREV_ACCUM[i+1]; vi++) {
			mpfrts_htilde(h, ORDER_INDEX-j, PREV_IV[vi], i, &ht, ORDER_INDEX);
			mpfrts_mul(&partial, ht, g[PREV_VI[vi]][j]);
			mpfrts_mul_i(&partial, partial, PREV_COEF[vi]);
			mpfrts_add(&sumint, sumint, partial);
		}
		mpfrts_mul_i(&sumint, sumint, ORDER_INDEX-j);
		mpfrts_add(&sumext, sumext, sumint);
	}
	mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
	mpfrts_sub(&sumext, f[i][ORDER_INDEX], sumext);
	mpfrts_div(&sumext, sumext, g[0][0]);
	mpfrts_set(&h[i][ORDER_INDEX], sumext);
	mpfr_clear (zero);
	mpfr_clear (partial);
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (ht);
}


void	mpfrts_log_t(mpfr_t u[NDER][MAX_ORDER+1],
		mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_log(&w[0][0], u[0][0]);
				} else {
					mpfrts_fgt_0(u,u,w,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(u,u,w,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_asin_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_asin(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_acos_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_acos(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_atan_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_atan(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_asinh_t(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_asinh(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_acosh_t(mpfr_t f[NDER][MAX_ORDER+1],mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_acosh(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_atanh_t(mpfr_t f[NDER][MAX_ORDER+1],mpfr_t g[NDER][MAX_ORDER+1], 
		mpfr_t h[NDER][MAX_ORDER+1], int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < NDER; i++) {
				if(i == 0) {
					mpfrts_atanh(&h[0][0], f[0][0]);
				} else {
					mpfrts_fgt_0(f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < NDER; i++) {
				mpfrts_fgt_i(f,g,h,i,ORDER_INDEX);
			}
		}
	}
}


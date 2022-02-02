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

#ifndef _mpfrITER_H
#define _mpfrITER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpfr.h"
#include "mpfrNUM.h"


extern	int		MAX_ORDER;
extern	long	NDER;
extern	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
extern	int		*PARTIAL_LIST, *FUNCTION_LIST;
extern	char 	**STR_DER;
extern	long	*POS_DER;
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


#endif


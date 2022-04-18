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


#ifndef Header_MP_TIDES_h
#define Header_MP_TIDES_h

#define real_MP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpfr.h"


/*********************************************************************************************/

#define set_iterations() \
set_iteration_parameters( \
NUM_DERIVATIVES,VARIABLES,PARAMETERS,FUNCTIONS, \
LINKS, PARTIALS_VARS, ORDER);set_iteration_lists(POS__PARTIALS, POS_FUNCTIONS, \
POS_ACCUM , POS_COEFS , POS_PREVI , POS_PREIV,\
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



void	varMP_init(mpfr_t var[NVARS+1][NDER][MAX_ORDER+1],mpfr_t v[], mpfr_t t);
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
void	mpfrts_add_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_sub_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_add_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_sub_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_mul_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t v[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_mul_t_c(char* cs, mpfr_t u[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_inv_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1],
			int ORDER_INDEX);
void	mpfrts_exp_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1],
			int ORDER_INDEX);
void	mpfrts_pow_t_c(mpfr_t u[NDER][MAX_ORDER+1], char* cs, 
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
			mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1],int ORDER_INDEX);
void	mpfrts_sinh_t(mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	mpfrts_cosh_t(mpfr_t c[NDER][MAX_ORDER+1], mpfr_t s[NDER][MAX_ORDER+1], 
			mpfr_t w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    mpfrts_sinh_cosh_t (mpfr_t f[NDER][MAX_ORDER+1], 
			mpfr_t s[NDER][MAX_ORDER+1], mpfr_t c[NDER][MAX_ORDER+1],int ORDER_INDEX);
void	mpfrts_fgt_0(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
			mpfr_t h[NDER][MAX_ORDER+1], long i);
void	mpfrts_fgt_i(mpfr_t f[NDER][MAX_ORDER+1], mpfr_t g[NDER][MAX_ORDER+1],
			mpfr_t h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	mpfrts_log_t(mpfr_t u[NDER][MAX_ORDER+1], mpfr_t w[NDER][MAX_ORDER+1],int ORDER_INDEX);
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

void	mpfrts_init (mpfr_t *rop); 
void	mpfrts_set_d (mpfr_t *rop, double op); 
void	mpfrts_set_str (mpfr_t *rop, char *op); 
void	mpfrts_set (mpfr_t *rop, mpfr_t op); 
double	mpfrts_get_d (mpfr_t op); 
int		mpfrts_get_prec (); 
void 	mpfrts_set_prec (int dig);
void	mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2); 
void 	mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2);
void 	mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2);
void	mpfrts_abs(mpfr_t  *rop, mpfr_t op); 
int		mpfrts_greater(mpfr_t op1, mpfr_t op2); 
int		mpfrts_greaterequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_less(mpfr_t op1, mpfr_t op2); 
int		mpfrts_lessequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_equal(mpfr_t op1, mpfr_t op2); 
void	mpfrts_log(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_cos(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_sin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op); 
void	mpfrts_tan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sec(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csc(mpfr_t  *rop, mpfr_t op);
void	mpfrts_cot(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acos(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_cosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_tanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sech(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csch(mpfr_t  *rop, mpfr_t op);
void	mpfrts_coth(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_write (char *c, mpfr_t op);
void	mpfrts_fread (FILE *file, mpfr_t rop);
void	mpfrts_fwrite (FILE *file, mpfr_t op, int prec); 

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
#define  add_t				mpfrts_add_t
#define  add_t_c			mpfrts_add_t_c
#define  sub_t				mpfrts_sub_t
#define  sub_t_c			mpfrts_sub_t_c
#define  mul_t				mpfrts_mul_t
#define  mul_t_c			mpfrts_mul_t_c
#define  divide_t			mpfrts_div_t
#define  inv_t				mpfrts_inv_t
#define  exp_t				mpfrts_exp_t
#define  pow_t_c			mpfrts_pow_t_c
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

int getOrder ();
int getNsteps ();

#endif


/*********************************************************************************************/

#ifndef mpfr_lorenz_tides_h
#define mpfr_lorenz_tides_h

long  mpfr_lorenz_columns();
long  mpfr_lorenz(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1]);
long  mpfr_lorenz_pos_der(char *der);
long  mpfr_lorenz_variable_column(int v, char *der);
long  mpfr_lorenz_function_column(int v, char *der);

#endif



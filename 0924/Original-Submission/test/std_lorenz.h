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


#ifndef Header_DP_TIDES_h
#define Header_DP_TIDES_h

#define real_Double

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*********************************************************************************************/

#define set_iterations() \
set_iteration_parameters( \
NUM_DERIVATIVES,VARIABLES,PARAMETERS,FUNCTIONS, \
LINKS, PARTIALS_VARS, ORDER);\
set_iteration_lists(POS__PARTIALS, POS_FUNCTIONS, \
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



void	varDB_init(double var[NVARS+1][NDER][MAX_ORDER+1],double v[], double t);

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

void	double_add_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_add_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);


void	double_mul_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_mul_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
	double h[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_inv_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
	int ORDER_INDEX);
void	double_exp_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
	int ORDER_INDEX);
void	double_pow_t_c(double u[NDER][MAX_ORDER+1], char* cs, 
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
	double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1],int ORDER_INDEX);
void	double_sinh_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_cosh_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
	double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    double_sinh_cosh_t (double f[NDER][MAX_ORDER+1], 
	double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1],int ORDER_INDEX);

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
#define  add_t				double_add_t
#define  add_t_c			double_add_t_c
#define  sub_t				double_sub_t
#define  sub_t_c			double_sub_t_c
#define  mul_t				double_mul_t
#define  mul_t_c			double_mul_t_c
#define  divide_t			double_div_t
#define  inv_t				double_inv_t
#define  exp_t				double_exp_t
#define  pow_t_c			double_pow_t_c
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


#endif


/*********************************************************************************************/

#ifndef std_lorenz_tides_h
#define std_lorenz_tides_h

long  std_lorenz_columns();
long  std_lorenz(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1]);
long  std_lorenz_pos_der(char *der);
long  std_lorenz_variable_column(int v, char *der);
long  std_lorenz_function_column(int v, char *der);

#endif



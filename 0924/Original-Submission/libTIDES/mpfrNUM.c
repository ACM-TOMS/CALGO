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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpfrNUM.h"


int PRECISION = 16;  
int BINARY_PRECISION = 53;  

int binary_precision(int prec) {
	return ceil(prec * 3.3219);
}

void mpfrts_init(mpfr_t *rop) {
	BINARY_PRECISION = binary_precision(PRECISION);
	mpfr_init2(*rop, BINARY_PRECISION);
	mpfrts_set_i(rop, 0);
}
void mpfrts_set_i(mpfr_t  *rop, long op) {
	mpfr_set_si(*rop, op, GMP_RND);
}
void mpfrts_set_d(mpfr_t  *rop, double op) {
	mpfr_set_d(*rop, op, GMP_RND);
}
void mpfrts_set_str(mpfr_t  *rop, char *op) {
	mpfr_set_str(*rop, op, 10, GMP_RND);
}
void mpfrts_set(mpfr_t  *rop, mpfr_t op) {
	mpfr_set(*rop, op, GMP_RND);
}
double mpfrts_get_d(mpfr_t op) {
	return mpfr_get_d(op, GMP_RND);
}
long    mpfrts_get_i(mpfr_t op) {
	return mpfr_get_si(op, GMP_RND);
}


int mpfrts_get_prec () {
	return PRECISION;
}

void mpfrts_set_prec (int dig) {
	PRECISION = dig;
}


void mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_add(*rop, op1, op2, GMP_RND);
}
void mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2) {
	mpfr_add_si (*rop, op1, op2, GMP_RND);
}
void mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_sub(*rop, op1, op2, GMP_RND);
}
void mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2) {
	mpfr_sub_si (*rop, op1, op2, GMP_RND);
}
void mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2) {
	mpfr_si_sub (*rop, op1, op2, GMP_RND);
}
void mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_mul(*rop, op1, op2, GMP_RND);
}
void mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2) {
	mpfr_mul_si (*rop, op1, op2, GMP_RND);
}
void mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_div(*rop, op1, op2, GMP_RND);
}
void mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2) {
	mpfr_div_si (*rop, op1, op2, GMP_RND);
}
void mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2) {
	mpfr_si_div (*rop, op1, op2, GMP_RND);
}
void mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_pow(*rop, op1, op2, GMP_RND);
}
void mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2) {
	mpfr_pow_si (*rop, op1, op2, GMP_RND);
}
void mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2) {
	mpfr_ui_pow (*rop, op1, op2, GMP_RND);
}
void mpfrts_abs (mpfr_t  *rop, mpfr_t op) {
	mpfr_abs(*rop, op, GMP_RND);
}
void mpfrts_neg(mpfr_t  *rop, mpfr_t op){
	mpfr_neg(*rop, op, GMP_RND);
}

int	mpfrts_greater(mpfr_t op1, mpfr_t op2) {
	return mpfr_greater_p(op1, op2);
}
int	mpfrts_greaterequal(mpfr_t op1, mpfr_t op2) {
	return mpfr_greaterequal_p(op1, op2);
}
int	mpfrts_less(mpfr_t op1, mpfr_t op2) {
	return mpfr_less_p(op1, op2);
} 
int	mpfrts_lessequal(mpfr_t op1, mpfr_t op2) {
	return mpfr_lessequal_p(op1, op2);
}
int	mpfrts_equal(mpfr_t op1, mpfr_t op2) {
	return mpfr_equal_p(op1, op2);
}


void	mpfrts_sqrt(mpfr_t  *rop, mpfr_t op){
	mpfr_sqrt(*rop, op, GMP_RND);
}
void mpfrts_log(mpfr_t  *rop, mpfr_t op) {
	mpfr_log(*rop, op, GMP_RND);
}
void mpfrts_log2(mpfr_t  *rop, mpfr_t op) {
	mpfr_log2(*rop, op, GMP_RND);
}
void mpfrts_log10(mpfr_t  *rop, mpfr_t op) {
	mpfr_log10(*rop, op, GMP_RND);
}
void mpfrts_exp(mpfr_t  *rop, mpfr_t op) {
	mpfr_exp(*rop, op, GMP_RND);
}
void mpfrts_exp2(mpfr_t  *rop, mpfr_t op) {
	mpfr_exp2(*rop, op, GMP_RND);	
}
void mpfrts_exp10(mpfr_t  *rop, mpfr_t op) {
	mpfr_exp10(*rop, op, GMP_RND);
}
void mpfrts_cos(mpfr_t  *rop, mpfr_t op) {
	mpfr_cos(*rop, op, GMP_RND);
}
void mpfrts_sin(mpfr_t  *rop, mpfr_t op) {
	mpfr_sin(*rop, op, GMP_RND);
}
void mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op) {
	mpfr_sin_cos(*rsin, *rcos, op, GMP_RND);
}
void mpfrts_tan(mpfr_t  *rop, mpfr_t op) {
	mpfr_tan(*rop, op, GMP_RND);
}
void mpfrts_sec(mpfr_t  *rop, mpfr_t op) {
	mpfr_sec(*rop, op, GMP_RND);
}
void mpfrts_csc(mpfr_t  *rop, mpfr_t op) {
	mpfr_csc(*rop, op, GMP_RND);
}
void mpfrts_cot(mpfr_t  *rop, mpfr_t op) {
	mpfr_cot(*rop, op, GMP_RND);
}
void mpfrts_acos(mpfr_t  *rop, mpfr_t op) {
	mpfr_acos(*rop, op, GMP_RND);
}
void mpfrts_asin(mpfr_t  *rop, mpfr_t op) {
	mpfr_asin(*rop, op, GMP_RND);
}
void mpfrts_atan(mpfr_t  *rop, mpfr_t op) {
	mpfr_atan(*rop, op, GMP_RND);
}
void mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2) {
	mpfr_atan2(*rop, op1, op2, GMP_RND);
}
void mpfrts_cosh(mpfr_t  *rop, mpfr_t op) {
	mpfr_cosh(*rop, op, GMP_RND);
}
void mpfrts_sinh(mpfr_t  *rop, mpfr_t op) {
	mpfr_sinh(*rop, op, GMP_RND);
}
void mpfrts_tanh(mpfr_t  *rop, mpfr_t op) {
	mpfr_tanh(*rop, op, GMP_RND);
}
void mpfrts_sech(mpfr_t  *rop, mpfr_t op) {
	mpfr_sech(*rop, op, GMP_RND);
}
void mpfrts_csch(mpfr_t  *rop, mpfr_t op) {
	mpfr_csch(*rop, op, GMP_RND);
}
void mpfrts_coth(mpfr_t  *rop, mpfr_t op) {
	mpfr_coth(*rop, op, GMP_RND);
}
void mpfrts_acosh(mpfr_t  *rop, mpfr_t op) {
	mpfr_acosh(*rop, op, GMP_RND);
}
void mpfrts_asinh(mpfr_t  *rop, mpfr_t op) {
	mpfr_asinh(*rop, op, GMP_RND);
}
void mpfrts_atanh(mpfr_t  *rop, mpfr_t op) {
	mpfr_atanh(*rop, op, GMP_RND);
}



void mpfrts_write_var(mpfr_t op) {
	mpfr_out_str(stdout, 10, PRECISION, op, GMP_RND);
}
void mpfrts_write(char *c, mpfr_t op) {
	printf("%s", c);
	printf(" = ");
	mpfr_out_str(stdout, 10, PRECISION, op, GMP_RND);
	putchar('\n');
}
void mpfrts_fread (FILE *file, mpfr_t rop) {
	mpfr_inp_str (rop, file, 10, GMP_RND);
}
void mpfrts_fwrite (FILE *file, mpfr_t op, int prec) {
	mpfr_out_str (file, 10, prec, op, GMP_RND);
	fprintf (file, "  ");
}


void	Array1MP_init(Array1MP *vec, long dim)
{
	long i;
	*vec = (Array1MP)calloc(dim, sizeof(mpfr_t));
	for(i=0; i <dim; i++) mpfrts_init((*vec)+i);;
}

void	Array2MP_init(Array2MP *vec, long rows, long columns)
{
	long i,j;
	(*vec) = (Array2MP)calloc(rows,sizeof(Array1MP));
	for( i=0; i<rows; i++ )
		(*vec)[i] = (Array1MP) calloc( columns, sizeof(mpfr_t));
	for( i=0; i<rows; i++ )
		for( j=0; j<columns; j++ )
			mpfrts_init((*vec)[i]+j);
}


void	Array1MP_set(Array1MP rop, Array1MP op, long dim)
{
	long i;
	for( i=0; i<dim; i++ ) mpfrts_set(&rop[i], op[i]);
}
void	Array2MP_set(Array2MP rop, Array2MP op, long rows, long columns)
{
	long i,j;
	for( i=0; i<rows; i++ )
		for( j=0; j<columns; j++ )
				mpfrts_set(&rop[i][j], op[i][j]);
}


void	Array2MP_column_set(Array2MP rop, Array1MP op, long c, long rows)
{
	long i; 
	for(i=0; i<rows; i++)
		mpfrts_set(&rop[i][c], op[i]);
}
void	Array2MP_row_set(Array2MP rop, Array1MP op, long r, long columns)
{
	long i; 
	for(i=0; i<columns; i++)
		mpfrts_set(&rop[r][i], op[i]);
}


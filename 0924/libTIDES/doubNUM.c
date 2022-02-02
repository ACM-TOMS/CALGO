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

#include "doubNUM.h"




void double_init(double *rop) {
	*rop = 0.;
}
void double_set_d(double *rop, double op) {
	*rop = op;
}
void double_set_str(double *rop, char *op) {
	*rop = strtod(op, NULL);
}
void double_set(double *rop, double op) {
	*rop = op;
}
double double_get_d(double op) {
	return op;
}
void double_get_prec (int dig) {}
void double_set_prec (int dig) {}
void double_clear (double op) {}



void double_add(double *rop, double	op1, double op2) {
	*rop = op1 + op2;
}
void double_sub(double *rop, double op1, double op2) {
	*rop = op1 - op2;
}
void double_mul(double *rop, double op1, double op2) {
	*rop = op1 * op2;
}
void double_div(double *rop, double op1, double op2) {
	*rop =  op1 / op2 ;
}
void double_pow(double *rop, double op1, double op2) {
	*rop = pow(op1, op2);
}
void double_abs(double *rop, double op) {
	*rop = fabs(op);
}

void double_add_i (double *rop, double op1, long  op2) {
	*rop = op1 + (double) op2;
}
void double_sub_i (double *rop, double op1, long  op2) {
	*rop = op1 - (double) op2;
}
void double_i_sub (double *rop, long  op1, double op2) {
	*rop = (double) op1 - op2;
}
void double_mul_i (double *rop, double op1, long  op2) {
	*rop = op1 * (double) op2;
}
void double_div_i (double *rop, double op1, long  op2) {
	*rop = op1 / (double) op2;
}
void double_i_div (double *rop, long  op1, double op2) {
	*rop = (double) op1 / op2;
}
void double_pow_i (double *rop, double op1, long  op2) {
	*rop = pow (op1, (double) op2);
}
void double_i_pow (double *rop, unsigned long  op1, double op2) {
	*rop = pow ((double) op1, op2);
}

int double_greater(double op1, double op2) {
	return op1 > op2;
}
int double_greaterequal(double op1, double op2) {
	return op1 >= op2;
}
int double_less(double op1, double op2) {
	return op1 < op2;
}
int double_lessequal(double op1, double op2) {
	return op1 <= op2;
}
int double_equal(double op1, double op2) {
	return op1 == op2;
}

void double_log(double *rop, double op) {
	*rop = log(op);
}
void double_log10(double *rop, double op) {
	*rop = log10(op);
}
void double_exp(double *rop, double op) {
	*rop = exp(op);
}
void double_exp2(double *rop, double op) {
	double_pow(rop, 2., op);
}
void double_exp10(double *rop, double op) {

	double_pow(rop, 10., op);
}
void double_cos(double *rop, double op) {
	*rop = cos(op);
}
void double_sin(double *rop, double op) {
	*rop = sin(op);
}
void double_sin_cos(double *rsin, double *rcos, double op) {
	*rsin = sin (op); *rcos = cos (op);
}
void double_tan(double *rop, double op) {
	*rop = tan(op);
}
void double_sec(double *rop, double op) {
	double_cos(rop, op) ;
	*rop = 1./(*rop);
}
void double_csc(double *rop, double op) {
	double_sin(rop, op);
	*rop = 1./(*rop);
}
void double_cot(double *rop, double op) {
	double_tan(rop, op);
	*rop = 1./(*rop);

}
void double_acos(double *rop, double op) {
	*rop = acos(op);
}
void double_asin(double *rop, double op) {
	*rop = asin(op);
}
void double_atan(double *rop, double op) {
	*rop = atan(op);
}
void double_atan2(double *rop, double op1, double op2) {
	*rop = atan2(op1, op2);
}
void double_cosh(double *rop, double op) {
	*rop = cosh(op);
}
void double_sinh(double *rop, double op) {
	*rop = sinh(op);
}
void double_tanh(double *rop, double op) {
	*rop = tanh(op);
}
void double_sech(double *rop, double op) {
	double_cosh(rop, op);
	*rop = 1./(*rop);

}
void double_csch(double *rop, double op) {
	double_sinh(rop, op);
	*rop = 1./(*rop);
}
void double_coth(double *rop, double op) {
	double_coth(rop,op);
	*rop = 1./(*rop);

}
void double_acosh(double *rop, double op) {
	*rop = acosh(op);
}
void double_asinh(double *rop, double op) {
	*rop = asinh(op);
}
void double_atanh(double *rop, double op) {
	*rop = atanh(op);
}

void	Array1DB_init(Array1DB *vec, long dim)
{
	long i;
	*vec = (Array1DB)calloc(dim, sizeof(double));
	for(i=0; i <dim; i++) double_init((*vec)+i);;
}

void	Array2DB_init(Array2DB *vec, long rows, long columns)
{
	long i,j;
	(*vec) = (Array2DB) calloc(rows,sizeof(Array1DB));
	for( i=0; i<rows; i++ )
		(*vec)[i] = (Array1DB) calloc( columns, sizeof(double));
	for( i=0; i<rows; i++ )
		for( j=0; j<columns; j++ )
			double_init((*vec)[i]+j);
}

void	Array3DB_init(Array2DB *rop, long dim, long rows, long columns)
{
	long i; 
	for( i=0; i<dim; i++ ) Array2DB_init(rop+i, rows, columns);
}


void	Array1DB_set(Array1DB rop, Array1DB op, long dim)
{
	long i;
	for( i=0; i<dim; i++ ) double_set(&rop[i], op[i]);
}
void	Array2DB_set(Array2DB rop, Array2DB op, long rows, long columns)
{
	long i,j;
	for( i=0; i<rows; i++ )
		for( j=0; j<columns; j++ )
				double_set(&rop[i][j], op[i][j]);
}

void	Array2DB_column_set(Array2DB rop, Array1DB op, long c, long rows)
{
	long i; 
	for(i=0; i<rows; i++)
		double_set(&rop[i][c], op[i]);
}
void	Array2DB_row_set(Array2DB rop, Array1DB op, long r, long columns)
{
	long i; 
	for(i=0; i<columns; i++)
		double_set(&rop[r][i], op[i]);
}


void double_write (char *c, double op) {
	printf ("%s", c); printf (" = %.16le\n", op);
}

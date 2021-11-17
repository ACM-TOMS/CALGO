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

#include "mpfr_kepler.h"


long  mpfr_kepler(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])
{
	static int   VARIABLES        = 4;
	static int   PARAMETERS       = 1;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 10;
	static int   PARTIALS_VARS    = 0;
	static long  NUM_DERIVATIVES  = 1;
	static long  NUM_COLUMNS      = 4;

	static int   POS__PARTIALS[1] = {0};
	static int   POS_FUNCTIONS[1] = {0};

	static long  POS_ACCUM[2] = {0,1};
	static long  POS_COEFS[1] = {1};
	static long  POS_PREVI[1] = {0};
	static long  POS_PREIV[1] = {0};

	static long  POS_ACCUM_S[2] = {0,1};
	static long  POS_COEFS_S[1] = {1};
	static long  POS_PREVI_S[1] = {0};
	static long  POS_PREIV_S[1] = {0};


	if(ORDER < 0) return NUM_COLUMNS;

	static int  NOT_INITIALIZED = 1;
	if(NOT_INITIALIZED)
	{
		set_iterations();
		NOT_INITIALIZED = 0; 
	}
	set_max_order(ORDER);

	realNUM var[VARIABLES+1][NUM_DERIVATIVES][ORDER+1];
	realNUM par[PARAMETERS][NUM_DERIVATIVES][ORDER+1];
	realNUM link[LINKS][NUM_DERIVATIVES][ORDER+1];
	variables_init(var,v,t);
	parameters_init(par,p);
	links_init(link);
	derivatives_init(var,par,v);

	int i;
	for(i=0;  i<=ORDER; i++) {
		var_t(var[3],var[1], i);
		var_t(var[4],var[2], i);
		var_t(link[8],var[3], i);
		var_t(link[9],var[4], i);
		mul_t_c("-1.",var[1],link[0],i);
		mul_t_c("-1.",var[2],link[1],i);
		mul_t(var[1],var[1],link[2],i);
		mul_t(var[2],var[2],link[3],i);
		add_t(link[2],link[3],link[4],i);
		mul_t(link[0],par[0],link[5],i);
		mul_t(link[1],par[0],link[6],i);
		pow_t_c(link[4],"-1.5",link[7],i);
		mul_t(link[5],link[7],link[8],i);
		mul_t(link[6],link[7],link[9],i);
	}

	write_solution(cvfd,var,link);
	clear(var,par,link);

	return NUM_COLUMNS;
}

long  mpfr_kepler_columns()
{
	 return 5;
}

long  mpfr_kepler_pos_der(char *der)
{
	static char* STR_DER[1] = {"0000"};
	long i;
	for(i=0; i < 1; i++)
		if(strcmp(der,STR_DER[i]) == 0) return i;
	return -1;
}

long  mpfr_kepler_variable_column(int v, char *der)
{
	 return position_variable(v, mpfr_kepler_pos_der, der);
}

long  mpfr_kepler_function_column(int f, char *der)
{
	 return position_function(f, mpfr_kepler_pos_der, der);
}



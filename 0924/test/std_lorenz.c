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

#include "std_lorenz.h"


long  std_lorenz(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])
{
	static int   VARIABLES        = 3;
	static int   PARAMETERS       = 3;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 11;
	static int   PARTIALS_VARS    = 0;
	static long  NUM_DERIVATIVES  = 1;
	static long  NUM_COLUMNS      = 3;

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
		var_t(link[6],var[1], i);
		var_t(link[10],var[2], i);
		var_t(link[9],var[3], i);
		sub_t(var[2],var[1],link[0],i);
		mul_t_c("-1.",par[2],link[1],i);
		mul_t_c("-1.",var[1],link[2],i);
		mul_t(par[1],var[1],link[3],i);
		mul_t(var[1],var[2],link[4],i);
		sub_t(link[3],var[2],link[5],i);
		mul_t(link[0],par[0],link[6],i);
		mul_t(link[1],var[3],link[7],i);
		mul_t(link[2],var[3],link[8],i);
		add_t(link[4],link[7],link[9],i);
		add_t(link[5],link[8],link[10],i);
	}

	write_solution(cvfd,var,link);

	return NUM_COLUMNS;
}

long  std_lorenz_columns()
{
	 return 4;
}

long  std_lorenz_pos_der(char *der)
{
	static char* STR_DER[1] = {"000"};
	long i;
	for(i=0; i < 1; i++)
		if(strcmp(der,STR_DER[i]) == 0) return i;
	return -1;
}

long  std_lorenz_variable_column(int v, char *der)
{
	 return position_variable(v, std_lorenz_pos_der, der);
}

long  std_lorenz_function_column(int f, char *der)
{
	 return position_function(f, std_lorenz_pos_der, der);
}



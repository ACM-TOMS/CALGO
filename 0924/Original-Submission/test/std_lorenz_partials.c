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

#include "std_lorenz_partials.h"


long  std_lorenz_partials(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])
{
	static int   VARIABLES        = 3;
	static int   PARAMETERS       = 3;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 11;
	static int   PARTIALS_VARS    = 2;
	static long  NUM_DERIVATIVES  = 6;
	static long  NUM_COLUMNS      = 18;

	static int   POS__PARTIALS[2] = {4,6};
	static int   POS_FUNCTIONS[1] = {0};

	static long  POS_ACCUM[7] = {0,1,3,5,8,12,15};
	static long  POS_COEFS[15] = {1,1,1,1,1,1,2,1,1,1,1,1,1,2,1};
	static long  POS_PREVI[15] = {0,0,1,0,2,0,1,3,0,2,1,4,0,2,5};
	static long  POS_PREIV[15] = {0,1,0,2,0,3,1,0,4,1,2,0,5,2,0};

	static long  POS_ACCUM_S[7] = {0,1,2,3,5,7,9};
	static long  POS_COEFS_S[9] = {1,1,1,1,1,1,1,1,1};
	static long  POS_PREVI_S[9] = {0,0,0,0,1,0,2,0,2};
	static long  POS_PREIV_S[9] = {0,1,2,3,1,4,1,5,2};


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

long  std_lorenz_partials_columns()
{
	 return 19;
}

long  std_lorenz_partials_pos_der(char *der)
{
	static char* STR_DER[6] = {"00","10","01","20","11","02"};
	long i;
	for(i=0; i < 6; i++)
		if(strcmp(der,STR_DER[i]) == 0) return i;
	return -1;
}

long  std_lorenz_partials_variable_column(int v, char *der)
{
	 return position_variable(v, std_lorenz_partials_pos_der, der);
}

long  std_lorenz_partials_function_column(int f, char *der)
{
	 return position_function(f, std_lorenz_partials_pos_der, der);
}



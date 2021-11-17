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

#ifndef _commonITER_H
#define _commonITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

int	 is_variable(int num);


void set_iteration_parameters(
		long nvd, int v, int p, int f, int l, int prt, int ord);

void set_max_order(int ord);

void set_iteration_lists( 
		int *prt, int *flst,  
		long *pra, long *prvi, long *priv, long *prc,  
		long *prsa, long *prsvi, long *prsiv, long * prcs);

#endif


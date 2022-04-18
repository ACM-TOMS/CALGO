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


int		MAX_ORDER;
long	NDER;
int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
int		*PARTIAL_LIST, *FUNCTION_LIST;
long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;



int is_variable(int num)
{
	return (num <= NVARS);
}

long position_variable(int v, position_derivative posder, char* der)
{
	long nvf, pd;
	nvf = NVARS+NFUNS;
	pd  = posder(der);
	if(pd <= -1 || v >= NVARS || v < 0 ) return -1;
	return (v+(pd*nvf)+1);
}
long position_function(int f, position_derivative posder, char* der)
{
	long nvf, pd;
	nvf = NVARS+NFUNS;
	pd  = posder(der);
	if(pd <= -1 || f >= NFUNS || f < 0 ) return -1;
	return (NVARS+f+(pd*nvf)+1);
}

/*************************************************************/
void set_iteration_parameters(
		long nvd, int v, int p, int f, int l, int prt, int ord)
{
	NDER			= nvd;
	NVARS			= v;
	NPARS			= p;
	NFUNS			= f;
	NLINKS			= l;
	NPARTIALS		= prt;
	MAX_ORDER		= ord;
}

void set_iteration_lists( int *prt, int *flst, 
		long *pra, long *prc, long *prvi, long *priv,   
		long *prsa,  long *prcs, long *prsvi, long *prsiv)
{
	PARTIAL_LIST	= prt;
	FUNCTION_LIST	= flst;
	PREV_ACCUM		= pra;
	PREV_COEF		= prc;
	PREV_VI			= prvi;
	PREV_IV			= priv;
	PREVSTAR_ACCUM	= prsa;
	PREVSTAR_COEF	= prcs;
	PREVSTAR_VI		= prsvi;
	PREVSTAR_IV		= prsiv;
}


void set_max_order (int ord) {
	MAX_ORDER = ord;
}
/************************************************************************/

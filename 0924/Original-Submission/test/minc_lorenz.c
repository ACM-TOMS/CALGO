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


#include "minc_tides.h"

void    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)
{
	int VAR,PAR,TT,i,j, inext;
	VAR = 3;
	PAR = 3;
	TT = 13;
	double XX[TT+1][MO+1];

	double pr[4];
	for(i=0; i<PAR; i++) pr[i] = p[i];
	pr[3] = -0.1e1 * pr[2];

	for(j=0; j<=TT; j++)
		for(i=0; i<=ORDER; i++)
			XX[j][i] = 0.e0;
	XX[0][0] = t;
	XX[0][1] = 1.e0;
	for(i=1;i<=VAR;i++) {
		XX[i][0] = v[i-1];
	}

	for(i=0;i<ORDER;i++) {
		XX[4][i] = XX[2][i]-XX[1][i];
		XX[5][i] = -0.1e1*XX[1][i];
		XX[6][i] = pr[1]*XX[1][i];
		XX[7][i] = mul_mc(XX[1],XX[2],i);
		XX[8][i] = XX[6][i]-XX[2][i];
		XX[9][i] = pr[0]*XX[4][i];
		XX[10][i] = pr[3]*XX[3][i];
		XX[11][i] = mul_mc(XX[5],XX[3],i);
		XX[12][i] = XX[7][i]+XX[10][i];
		XX[13][i] = XX[8][i]+XX[11][i];
		inext = i + 1;
		XX[1][inext] =  XX[9][i]/inext;
		XX[2][inext] =  XX[13][i]/inext;
		XX[3][inext] =  XX[12][i]/inext;
	}

	for(j=0; j<=VAR; j++)
		for(i=0; i<=ORDER; i++)
			XVAR[i][j] = XX[j][i];
}


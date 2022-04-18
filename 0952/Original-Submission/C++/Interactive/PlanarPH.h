#pragma once
#include <complex>
#include <math.h>
#include "PointPH.h"
using namespace std;

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

struct PHquintic
{
	complex<double> p[6] ;		//Bezier control points of PH quintic
	complex<double> w[3] ;		//Bernstein coefficients of w(t) polynomial
	double sigma[5] ;			//parametric speed Bernstein coefficients
	double s[6] ;				//arc length Bernstein coefficients
};

class PlanarPH
{
protected:
	void Spline(BOOL closed);
	void Hermite(complex<double> p0, complex<double> p1, complex<double> p4, complex<double> p5);

public:
	PHquintic spline[Max];
	complex<double> q[Max+1];
	int num, iter;

	PlanarPH();
	PlanarPH(const CPointPH m_pt[], int index);
	PlanarPH(const CPointPH m_pt[]);
	virtual ~PlanarPH();
	friend complex<double> beval(int n, const complex<double> b[], double t);
	friend double beval(int n, const double b[], double t);
	friend void tridiag_open(int n, const complex<double> a[], 
									 const complex<double> b[], 
									 const complex<double> c[], 
									 const complex<double> d[], 
									 complex<double>x[]);
	friend void tridiag_closed(int n, const complex<double> a[], 
									  const complex<double> b[], 
									  const complex<double> c[], 
									  const complex<double> d[], 
									  complex<double>x[]);
	friend complex<double> getCtrlPt(const PlanarPH &pph, int i, int j);
	friend int getIter(const PlanarPH &pph);
	friend double getParaSpeed(const PlanarPH &pph, int i, double t);
	friend double getArcLength(const PlanarPH &pph, int i, double t);
	friend double getEnergy(const PlanarPH &pph, int i);
	double getEnergy();
	friend complex<double> * getOffset(const PlanarPH &pph, int i, double d);
};
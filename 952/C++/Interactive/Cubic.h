#pragma once
#include <complex>
#include "PointPH.h"
using namespace std;

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CCubic
{
protected:
	void cubicspline(BOOL isclose);

public:
	complex<double> p[Max+1], dp[Max+1];
	double t[Max+1];
	int num;

	CCubic(void);
	CCubic(const CPointPH m_pt[], int index);
	~CCubic(void);
	friend complex<double> getSpline(const CCubic &cubic, int i, double u);
};
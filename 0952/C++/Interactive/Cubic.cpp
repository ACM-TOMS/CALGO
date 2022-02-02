#include "stdafx.h"
#include "Cubic.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

CCubic::CCubic(void)
{
}

CCubic::CCubic(const CPointPH m_pt[], int index)
{
	num = index;
	for(int i = 0; i <= num; i++)
		p[i] = complex<double>(m_pt[i].xPos, m_pt[i].yPos);
	cubicspline((p[num] == p[0]));
}

CCubic::~CCubic(void)
{
}

complex<double> getSpline(const CCubic &cubic, int i, double u)
{
	if(cubic.p[cubic.num] == cubic.p[0])
	{
		complex<double> b[4][Max+3];
		int j, k;
		double tau;

		for(j = i-3; j <= i; j++)
		{
			if(j > cubic.num-1)
				b[0][j] = cubic.dp[j-cubic.num];
			else
				b[0][j] = cubic.dp[j];
		}
		for(k = 1; k <= 3; k++)
			for(j = i-3+k; j <= i; j++)
			{
				tau = (u-(j-3.0))/(((j+4.0-k)-3.0)-(j-3.0));
				b[k][j] = (1.0-tau)*b[k-1][j-1]+tau*b[k-1][j];
			}

		return b[3][i];
	}
	else
	{
		double alpha0, alpha1, belta0, belta1;

		alpha0 = 1-3*u*u+2*u*u*u;
		alpha1 = 3*u*u-2*u*u*u;
		belta0 = u-2*u*u+u*u*u;
		belta1 = 0-u*u+u*u*u;

		return (cubic.p[i-1]*alpha0+cubic.p[i]*alpha1+cubic.t[i]*(cubic.dp[i-1]*belta0+cubic.dp[i]*belta1));
	}
}

void CCubic::cubicspline(BOOL isclose)
{
	int i;
	complex<double> d[Max+1], delta[Max+1];
	double a[Max+1], b[Max+1], c[Max+1], beta[Max+1], epsilon[Max+1], theta, m;

	if(isclose)
	{
		for(i = 0; i <= num-1; i++)
		{
			a[i] = 1.0/6.0;
			b[i] = 2.0/3.0;
			c[i] = 1.0/6.0;
			d[i] = p[i];
		}
		beta[0] = b[0];
		delta[0] = d[0];
		epsilon[0] = a[0];
		for(i = 0; i <= num-3; i++)
		{
			m = a[i+1]/beta[i];
			beta[i+1] = b[i+1]-m*c[i];
			delta[i+1] = d[i+1]-m*delta[i];
			epsilon[i+1] = -m*epsilon[i];
		}
		epsilon[num-2] = epsilon[num-2]+c[num-2];
		m = a[num-1]/beta[num-2];
		beta[num-1] = b[num-1]-m*epsilon[num-2];
		delta[num-1] = d[num-1]-m*delta[num-2];
		theta = c[num-1];
		for(i = 0; i <= num-2; i++)
		{
			m = theta/beta[i];
			beta[num-1] = beta[num-1]-m*epsilon[i];
			delta[num-1] = delta[num-1]-m*delta[i];
			theta = -m*c[i];
		}
		dp[num-1] = delta[num-1]/beta[num-1];
		dp[num-2] = (delta[num-2]-epsilon[num-2]*dp[num-1])/beta[num-2];
		for(i = num-3; i >= 0; i--)
			dp[i] = (delta[i]-c[i]*dp[i+1]-epsilon[i]*dp[num-1])/beta[i];
	}
	else
	{
		for(i = 1; i <= num; i++)
			t[i] = abs(p[i]-p[i-1]);
		for(i = 0; i <= num; i++)
		{
			if(i == 0)
			{
				b[i] = 1.0;
				c[i] = 1.0;
				d[i] = 2.0*(p[i+1]-p[i])/t[i+1];
			}
			else
				if(i == num)
				{
					b[i] = 1.0;
					a[i] = 1.0;
					d[i] = 2.0*(p[i]-p[i-1])/t[i];
				}
				else
				{
					a[i] = t[i+1];
					b[i] = 2.0*(t[i]+t[i+1]);
					c[i] = t[i];
					d[i] = 3.0*(t[i]*(p[i+1]-p[i])/t[i+1]+t[i+1]*(p[i]-p[i-1])/t[i]);
				}
		}
		beta[0] = b[0];
		delta[0] = d[0];
		for(i = 0; i < num; i++)
		{
			m = a[i+1]/beta[i];
			beta[i+1] = b[i+1]-m*c[i];
			delta[i+1] = d[i+1]-m*delta[i];
		}
		dp[num] = delta[num]/beta[num];
		for(i = num-1; i >= 0; i--)
			dp[i] = (delta[i]-c[i]*dp[i+1])/beta[i];
	}
}
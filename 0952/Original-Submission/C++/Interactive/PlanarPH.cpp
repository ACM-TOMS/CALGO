#include "stdafx.h"
#include "PlanarPH.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

PlanarPH::PlanarPH(void)
{
}

PlanarPH::PlanarPH(const CPointPH m_pt[], int index)
{
	num = index;
	for(int i = 0; i <= num; i++)
		q[i] = complex<double>(m_pt[i].xPos, m_pt[i].yPos);
	Spline((q[num] == q[0]));
}

PlanarPH::PlanarPH(const CPointPH m_pt[])
{
	num = 1;
	complex<double> p0, p1, p4, p5;
	p0 = complex<double>(m_pt[0].xPos, m_pt[0].yPos);
	p1 = complex<double>(m_pt[1].xPos, m_pt[1].yPos);
	p4 = complex<double>(m_pt[3].xPos, m_pt[3].yPos);
	p5 = complex<double>(m_pt[2].xPos, m_pt[2].yPos);
	Hermite(p0, p1, p4, p5);
}

PlanarPH::~PlanarPH(void)
{
}

//
/*
Computes value of a degree-n polynomial with Bernstein coefficients 
b[] at the value t of the independent variable
using the de Casteljau algorithm.
n <= 2*Dgr+1.
Maximum allowed degree is MAXDEGREE.
Return complex double.
*/
complex<double> beval(int n, const complex<double> b[], double t)
{
	int i, j;
	complex<double> blast[MAXDEGREE+1], bnext[MAXDEGREE+1];

	if(n == 0)return b[0];
	for(i = 0; i <= n; i++)
		blast[i] = b[i];
	for(i = 1; i <= n; i++)
	{
		for(j = i; j <= n; j++)
			bnext[j] = (1.0-t)*blast[j-1]+t*blast[j];
		for(j = i; j <= n; j++)
			blast[j] = bnext[j];
	}
	return bnext[n];
}

/* 
Computes value of a degree-n polynomial with Bernstein coefficients 
b[] at the value t of the independent variable
using the de Casteljau algorithm.
n <= 2*Dgr+1.
Maximum allowed degree is MAXDEGREE.
Return double.
*/
double beval(int n, const double b[], double t)
{
	int i, j;
	double blast[MAXDEGREE+1], bnext[MAXDEGREE+1];

	if(n == 0)return b[0];
	for(i = 0; i <= n; i++)
		blast[i] = b[i];
	for(i = 1; i <= n; i++)
	{
		for(j = i; j <= n; j++)
			bnext[j] = (1.0-t)*blast[j-1]+t*blast[j];
		for(j = i; j <= n; j++)
			blast[j] = bnext[j];
	}
	return bnext[n];
}

/*
Solve tridiagonal system appropriate to open PH quintic spline with cubic end conditions.
Array a[], b[], c[] define the lower, main, and upper diagonal matrix elements,
array d[] defines the right hand side values. 
The solutions are returned in array x[].
n <= Max.
*/
void tridiag_open(int n, const complex<double> a[],
				   const complex<double> b[],
				   const complex<double> c[],
				   const complex<double> d[],
				   complex<double> x[])
{
	int i;
	complex<double> beta[Max], rho[Max];

	beta[1] = b[1];
	rho[1] = d[1];
	for(i = 1; i <= n-1; i++)
	{
		beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i];
		rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i];
    }
	x[n] = rho[n]/beta[n];
	for(i = n-1; i >= 1; i--)x[i] = (rho[i]-c[i]*x[i+1])/beta[i];
}

/*
Solve tridiagonal system appropriate to closed PH quintic spline with periodic end conditions.
Array a[], b[], c[] define the lower, main, and upper diagonal matrix elements,
array d[] defines the right hand side values. 
The solutions are returned in array x[].
n <= Max.
*/
void tridiag_closed(int n, const complex<double> a[],
					const complex<double> b[],
					const complex<double> c[],
					const complex<double> d[],
					complex<double> x[])
{
	int i;
	complex<double> beta[Max], rho[Max], zeta[Max], theta;
  
	beta[1] = b[1]; 
	rho[1] = d[1]; 
	zeta[1] = a[1];
	for(i = 1; i <= n-2; i++)
    {
		beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i];
		rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i];
		zeta[i+1] = -a[i+1]*zeta[i]/beta[i];
    }
	zeta[n-1] += c[n-1];
	beta[n] = b[n]-a[n]*zeta[n-1]/beta[n-1];
	rho[n] = d[n]-a[n]*rho[n-1]/beta[n-1];
	theta = c[n];
	for(i = 1; i <= n-1; i++)
    {
		beta[n] -= theta*zeta[i]/beta[i];
		rho[n] -= theta*rho[i]/beta[i];
		theta *= -c[i]/beta[i];
    }
	x[n] = rho[n]/beta[n];
	x[n-1] = (rho[n-1]-zeta[n-1]*x[n])/beta[n-1];
	for(i = n-2; i >= 1; i--)
    {
		x[i] = (rho[i]-c[i]*x[i+1]-zeta[i]*x[n])/beta[i];
    }
}

complex<double> getCtrlPt(const PlanarPH &pph, int i, int j)
{
	return pph.spline[i].p[j];
}

int getIter(const PlanarPH &pph)
{
	return pph.iter;
}

double getParaSpeed(const PlanarPH &pph, int i, double t)
{
	return beval(Dgr-1, pph.spline[i].sigma, t);
}

double getArcLength(const PlanarPH &pph, int i, double t)
{
	if(i == 1)
		return beval(Dgr, pph.spline[i].s, t);
	else
		return beval(Dgr, pph.spline[i].s, t)+getArcLength(pph, i-1, 1.0);
}

/* 
Computes the bending energy of the planar PH quintic specified 
in struct PH quintic "curve". The energy is returned as the value of
the function. The parameter "sml" is used to determine when w(t) 
degenerates to a linear polynomial. 
*/
double getEnergy(const PlanarPH &pph, int i)
{
	complex<double> I(0.0, 1.0);
	complex<double> w0, w1, w2, a1, a2, a3, b1, b2, b3, a, b, k, ac, bc, kc;
	double alpha, beta, energy, f0, f1;
	double sml = 0.000000000001;

	w0 = pph.spline[i].w[0];
	w1 = pph.spline[i].w[1];
	w2 = pph.spline[i].w[2];

	k = w2-2.0*w1+w0;

	if(abs(k) < sml)
	{
		if(abs(w2-w0) < sml)
			energy = 0.0;
		else
		{
			k = w2-w0;
			a = w0/(w0-w2);
			alpha = imag(a);
			f0 = abs(a)*abs(a);
			f1 = abs(1.0-a)*abs(1.0-a);

			energy = real(a)/(f0*f0) + real(1.0-a)/(f1*f1)+
					1.5*(real(a)/f0+real(1.0-a)/f1)/(alpha*alpha)+
					1.5*(atan((1.0-real(a))/alpha)+atan(real(a)/alpha))/pow(alpha, 3);
			energy /= norm(k);
		}
	}
	else
	{
		a = (w0-w1+sqrt(w1*w1-w0*w2))/k;
		b = (w0-w1-sqrt(w1*w1-w0*w2))/k;
		alpha = imag(a); 
		beta = imag(b);
		ac = conj(a); 
		bc = conj(b); 
		a3 = I/(8.0*alpha*(a-b)*(a-bc));
		b3 = I/(8.0*beta*(a-b)*(ac-b));
		a2 = (1.5*I/alpha+1.0/(a-b)-3.0/(a-bc))*a3;
		b2 = (1.5*I/beta-1.0/(a-b)+3.0/(ac-b))*b3;
		a1 = 1.5*I*a2/alpha+(0.75/(alpha*alpha)-2.0/((a-b)*(a-b))
	     + 6.0/((a-bc)*(a-bc))-(1.0-2.0*beta/alpha)/((a-b)*(a-bc)))*a3;
		b1 = 1.5*I*b2/beta+(0.75/(beta*beta)-2.0/((a-b)*(a-b))
	      + 6.0/((ac-b)*(ac-b))-(1.0-2.0*alpha/beta)/((a-b)*(ac-b)) )*b3;

		energy = 2.0*real(a1)*log(abs((1.0-a)/a))-2.0*imag(a1)*(arg(1.0-a)-arg(-a))
			   + 2.0*real(b1)*log(abs((1.0-b)/b))-2.0*imag(b1)*(arg(1.0-b)-arg(-b))
			   - real(2.0*a2/(a*(1.0-a))+2.0*b2/(b*(1.0-b))+(2.0*a-1.0)*a3/(a*a*(1.0-a)*(1.0-a))
			   + (2.0*b-1.0)*b3/(b*b*(1.0-b)*(1.0-b))); 
		energy *= 4.0/norm(k);
	}

	return energy;
}

double PlanarPH::getEnergy()
{
	complex<double> I(0.0, 1.0);
	complex<double> w0, w1, w2, a1, a2, a3, b1, b2, b3, a, b, k, ac, bc, kc;
	double alpha, beta, energy[Max], energysum;
	int i, loopbegin, loopend;

	loopbegin = 1;
	loopend = num;
	
	if(q[0] != q[num])
	{
		loopbegin = 2;
		loopend = num;
	}
	
	energysum = 0.0;
	for(i = loopbegin; i <= loopend; i++)
	{
		w0 = spline[i].w[0];
		w1 = spline[i].w[1];
		w2 = spline[i].w[2];

		k = w2-2.0*w1+w0;
		a = (w0-w1+sqrt(w1*w1-w0*w2))/k;
		b = (w0-w1-sqrt(w1*w1-w0*w2))/k;
		alpha = imag(a); 
		beta = imag(b) ;
		ac = conj(a); 
		bc = conj(b); 
		kc = conj(k) ;
		a3 = I/(8.0*alpha*(a-b)*(a-bc));
		b3 = I/(8.0*beta*(a-b)*(ac-b));
		a2 = (1.5*I/alpha+1.0/(a-b)-3.0/(a-bc))*a3;
		b2 = (1.5*I/beta-1.0/(a-b)+3.0/(ac-b))*b3;
		a1 = 1.5*I*a2/alpha+(0.75/(alpha*alpha)-2.0/((a-b)*(a-b))
		   + 6.0/((a-bc)*(a-bc))-(1.0-2.0*beta/alpha)/((a-b)*(a-bc)))*a3;
		b1 = 1.5*I*b2/beta+(0.75/(beta*beta)-2.0/((a-b)*(a-b))
		   + 6.0/((ac-b)*(ac-b))-(1.0-2.0*alpha/beta)/((a-b)*(ac-b)) )*b3;

		energy[i] = (double)(2.0*real(a1)*log(abs((1.0-a)/a))-2.0*imag(a1)*(arg(1.0-a)-arg(-a))
				  + 2.0*real(b1)*log(abs((1.0-b)/b))-2.0*imag(b1)*(arg(1.0-b)-arg(-b))
				  - real(2.0*a2/(a*(1.0-a))+2.0*b2/(b*(1.0-b))+(2.0*a-1.0)*a3/(a*a*(1.0-a)*(1.0-a))
				  + (2.0*b-1.0)*b3/(b*b*(1.0-b)*(1.0-b)))); 
		energy[i] *= 4.0/abs(k*kc);
		energysum += energy[i];
	}

	return energysum;
}

/* 
Computes the offset curve at distance d from the PH quintic
defined in the struct "curve". The positive sense of the normal
vector is assumed to be to the right when traversing the curve 
with increasing parameter value.
The homogeneous coordinates for the weights and control points
of the degree 9 rational Bezier curve defining the offset are
stored in the arrays W[], X[], Y[].
Return coefficients of offset curve.
*/
complex<double> * getOffset(const PlanarPH &pph, int i, double d)
{
	int j, jmin, jmax, k, binom[10][10];
	double P[6][3], deltaP[6][3], X[10], Y[10], W[10], f;
	complex<double> * offsetCtrl = new complex<double>[10];

	binom[0][0] = 1;
	for(j = 1; j <= 9; j++)
	{
		binom[j][0] = 1;
		binom[j][j] = 1;
	}
	for(j = 2; j <= 9; j++)
		for(k = 1; k <= j-1; k++)
			binom[j][k] = binom[j-1][k-1]+binom[j-1][k];
	for(j = 0; j <= 5; j++)
    {
		P[j][0] = 1.0 ;
		P[j][1] = real(pph.spline[i].p[j]) ;
		P[j][2] = imag(pph.spline[i].p[j]) ;
    }
	for(j = 0; j <= 4; j++)
		for(k = 0; k <= 2; k++)
			deltaP[j][k] = P[j+1][k]-P[j][k];
	for(k = 0; k <= 9; k++)
    {
		W[k] = 0.0;
		X[k] = 0.0;
		Y[k] = 0.0;
		jmax = 4; 
		if(k < jmax)
			jmax = k;
		jmin = 0;
		if(k-5 > 0)
			jmin = k-5;
		for(j = jmin; j <= jmax; j++)
		{
			f = ((double)binom[k][j]*binom[9-k][4-j])/126.0;
			W[k] += f*pph.spline[i].sigma[j];
			X[k] += f*(pph.spline[i].sigma[j]*P[k-j][1]+5.0*d*deltaP[j][2]) ;
			Y[k] += f*(pph.spline[i].sigma[j]*P[k-j][2]-5.0*d*deltaP[j][1]) ;
		}
		offsetCtrl[k] = complex<double>(X[k]/W[k], Y[k]/W[k]);
    }

	return offsetCtrl;
}

/*
Computes a planar C^2 PH quintic spline interpolating n+1 points with
***cubic/periodic*** end conditions. 
The coordinates of the points to be interpolated are represented as 
array q[].
These points should have reasonably even spacing. 
The data for the resulting n PH quintic segments are returned as 
spline[] of the struct PHquintic spline array.
*/
void PlanarPH::Spline(BOOL closed)
{
	int i, j;
	complex<double> a[Max], b[Max], c[Max], d[Max], dz[Max], dq[Max], z[Max+1];
	complex<double> beta[Max], rho[Max], zeta[Max], w0, w1, w2;
	complex<double> theta;
	complex<double> complex0, complex1, complex4, complex6, complex8;
	double eta, znorm, dznorm, err;

	complex0 = complex<double>(0.0, 0.0);
	complex1 = complex<double>(1.0, 1.0);
	complex4 = complex<double>(4.0, 4.0);
	complex6 = complex<double>(6.0, 6.0);
	complex8 = complex<double>(8.0, 8.0);
/*
Construct the ordinary C^2 cubic spline.
The complex values dq[] are the computed nodal derivatives.
*/
	if(closed)
	{
		a[0] = complex1;
        b[0] = complex4;
        c[0] = complex1; 
        d[0] = 3.0*(q[1]-q[num-1]);
        zeta[0] = a[0];
	}
	else
	{
		b[0] = complex1;
		c[0] = complex1;
		a[num] = complex1;
		b[num] = complex1;
		d[0] = 2.0*(q[1]-q[0]);
		d[num] = 2.0*(q[num]-q[num-1]);
	}
	for(i = 1; i < num; i++) 
	{
		a[i] = complex1;
		b[i] = complex4;
		c[i] = complex1; 
		d[i] = 3.0*(q[i+1]-q[i-1]);
	}
	beta[0] = b[0]; 
	rho[0] = d[0]; 
	if(closed)
	{
		for(i = 0; i < num-2; i++)
		{
        	beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i];
        	rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i];
        	zeta[i+1] = complex0-a[i+1]*zeta[i]/beta[i];
    	}
    	zeta[num-2] += c[num-2]; 
    	rho[num-1] = d[num-1]-a[num-1]*rho[num-2]/beta[num-2];
    	beta[num-1] = b[num-1]-a[num-1]*zeta[num-2]/beta[num-2];
    	theta = c[num-1];
    	for(i = 0; i < num-1; i++)
    	{
        	beta[num-1] -= theta*zeta[i]/beta[i];
        	rho[num-1] -= theta*rho[i]/beta[i];
        	theta *= complex0-c[i]/beta[i];
    	}
    	dq[num-1] = rho[num-1]/beta[num-1];
    	dq[num-2] = (rho[num-2]-zeta[num-2]*dq[num-1])/beta[num-2];
    	for(i = num-3; i >= 0; i--)
    		dq[i] = (rho[i]-c[i]*dq[i+1]-zeta[i]*dq[num-1])/beta[i];
//    	q[num] = q[0]; 
    	dq[num] = dq[0];
	}
	else
	{
		for(i = 0; i < num; i++)
		{
			beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i];
			rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i];	
		}
		dq[num] = rho[num]/beta[num];
		for(i = num-1; i >= 0; i--)
			dq[i] = (rho[i]-c[i]*dq[i+1])/beta[i];
	}
/*
Compute the starting values for the iteration from the mid-point
derivatives of the ordinary C^2 cubic spline.
*/
	d[1] = 4.0*sqrt(6.0*(q[1]-q[0])-dq[0]-dq[1]);
	for(i = 2; i <= num; i++)
    {
		d[i] = 4.0*sqrt(6.0*(q[i]-q[i-1])-dq[i-1]-dq[i]);
		if(abs(d[i-1]+d[i]) < abs(d[i-1]-d[i]))d[i] = complex0-d[i];
    }
    for(i = 2; i < num; i++) 
    {
		a[i] = complex1; 
		b[i] = complex6; 
		c[i] = complex1; 
	}
    if(closed)
    {
    	if(abs(d[num]+d[1]) < abs(d[num]-d[1]))eta = -1.0;
    	else eta = 1.0;
    	a[1] = complex<double>(eta, eta); 
    	b[1] = complex6; 
    	c[1] = complex1;
    	a[num] = complex1;
    	b[num] = complex6;
    	c[num] = complex<double>(eta, eta);
    	tridiag_closed(num, a, b, c, d, z);
    }
    else
    {
    	b[1] = complex8; 
		c[1] = complex0;
		a[num] = complex0; 
		b[num] = complex8;
		tridiag_open(num, a, b, c, d, z);
    }
/* 
Begin Newton-Raphson iteration to compute the complex values 
z[] that defin the C^2 PH quintic spline.
*/
	err = 1.0; 
	for(iter = 1; err > Epsilon && iter <= 20; iter++)
    {
		for(i = 2; i < num; i++)
		{
			a[i] = 6.0*z[i-1]+13.0*z[i]+z[i+1];
			b[i] = 13.0*z[i-1]+54.0*z[i]+13.0*z[i+1];
			c[i] = z[i-1]+13.0*z[i]+6.0*z[i+1];
			d[i] = 60.0*(q[i]-q[i-1])-3.0*z[i-1]*z[i-1]-27.0*z[i]*z[i]-3.0*z[i+1]*z[i+1]
					-z[i-1]*z[i+1]-13.0*z[i-1]*z[i]-13.0*z[i]*z[i+1];
		}
		if(closed)
		{
			a[1] = 6.0*z[num]+eta*(13.0*z[1]+z[2]);
        	b[1] = eta*13.0*z[num]+54.0*z[1]+13.0*z[2];
        	c[1] = eta*z[num]+13.0*z[1]+6.0*z[2];
        	d[1] = 60.0*(q[1]-q[0])-3.0*z[num]*z[num]-27.0*z[1]*z[1]-3.0*z[2]*z[2]
                -eta*(z[num]*z[2]+13.0*z[num]*z[1])-13.0*z[1]*z[2];
            a[num] = 6.0*z[num-1]+13.0*z[num]+eta*z[1];
        	b[num] = 13.0*z[num-1]+54.0*z[num]+eta*13.0*z[1];
        	c[num] = eta*(z[num-1]+13.0*z[num])+6.0*z[1];
        	d[num] = 60.0*(q[num]-q[num-1])-3.0*z[num-1]*z[num-1]-27.0*z[num]*z[num]-3.0*z[1]*z[1]
                -eta*(z[num-1]*z[1]+13.0*z[num]*z[1])-13.0*z[num-1]*z[num];
            tridiag_closed(num, a, b, c, d, dz);
		}
		else
		{
			b[1] = 130.0*z[1]-10.0*z[2];
			c[1] = 10.0*z[2]-10.0*z[1];
			d[1] = 60.0*(q[1]-q[0])-65.0*z[1]*z[1]-5.0*z[2]*z[2]+10.0*z[1]*z[2];
			a[num] = 10.0*z[num-1]-10.0*z[num];
			b[num] = 130.0*z[num]-10.0*z[num-1];
			d[num] = 60.0*(q[num]-q[num-1])-65.0*z[num]*z[num]-5.0*z[num-1]*z[num-1]+10.0*z[num-1]*z[num];
			tridiag_open(num, a, b, c, d, dz);
		}
		for(i = 1; i <= num; i++)z[i] += dz[i];
	    znorm = 0.0; 
		dznorm = 0.0;
		for(i = 1; i <= num; i++)
		{
			znorm += abs(z[i])*abs(z[i]);
			dznorm += abs(dz[i])*abs(dz[i]);
		}
		err = sqrt(dznorm/znorm);
	}
	iter--;
	if(closed)
	{
		z[0] = eta*z[num]; 
    	z[num+1] = eta*z[1];
	}
	else
	{
		z[0] = 2.0*z[1]-z[2];
		z[num+1] = 2.0*z[num]-z[num-1];
	}
/*
Construct the Bezier control points for each segment of C^2 PH quintic spline
*/
	for(i = 1; i <= num; i++)
    {
		w0 = 0.5*(z[i-1]+z[i]);
		w1 = z[i];
		w2 = 0.5*(z[i]+z[i+1]);

		spline[i].w[0] = w0;
		spline[i].w[1] = w1;
		spline[i].w[2] = w2;
		
		spline[i].p[0] = q[i-1];
		spline[i].p[1] = spline[i].p[0]+0.2*w0*w0;
		spline[i].p[2] = spline[i].p[1]+0.2*w0*w1;
		spline[i].p[3] = spline[i].p[2]+0.2*(2.0*w1*w1+w0*w2)/3.0;
		spline[i].p[4] = spline[i].p[3]+0.2*w1*w2;
		spline[i].p[5] = spline[i].p[4]+0.2*w2*w2;

		spline[i].sigma[0] = abs(w0)*abs(w0);
		spline[i].sigma[1] = real(w0*conj(w1));
		spline[i].sigma[2] = (2*abs(w1)*abs(w1)+real(w0*conj(w2)))/3;
		spline[i].sigma[3] = real(w1*conj(w2));
		spline[i].sigma[4] = abs(w2)*abs(w2);

		spline[i].s[0] = 0.0;
		for(j = 1; j <= Dgr; j++)
			spline[i].s[j] = spline[i].s[j-1]+0.2*spline[i].sigma[j-1];
    }
}      

/* 
Given the first and last pairs of Bezier control points (p0,p1) 
(p4,p5) this function constructs the middle control points (p2,p3) 
such that the resulting curve is a planar PH quintic.
The complete set of Bezier control points are returned in struct
PHquintic "curve" together with the Bernstein coefficients of the 
w(t) polynomial, the parametric speed polynomial sigma(t), and the
arc length polynomial s(t).
Among the four formal solutions, the "good" solution is identified 
as that which minimizes the absolute rotation index.
The initial and final points, p0 and p5, are assumed to be distinct.
Otherwise, the function returns a degenerate curve: the single point 
p0=p5
The parameter "sml" is used to determine when w(t) degenerates to a 
linear polynomial, and when the numerator of the curvature degenerates 
to a linear polynomial.
//the expression |w0-2*w1+w2| in the coefficients of w(t). 
*/
void  PlanarPH::Hermite(complex<double> p0, complex<double> p1, complex<double> p4, complex<double> p5)
{
	int eta0, eta2, eta0_min, eta2_min, i, intervals;
	complex<double> q0, q1, q4, q5, scale, w0, w1, w2, k, a, b;
	double s1, s2, s3, c0, c1, c2, anglea, angleb;
	double t1, t2, temp, t[4], rindex, rindex_min;
	double sml = 0.000000000001;

	if(p5 == p0)
    {
		w0 = w1 = w2 = 0.0 ;
    }
	else
    {
		//scale input data (p0,p1,p4,p5) to canonical form (q0,q1,q4,q5)
		scale = p5-p0;
		q0 = complex<double>(0.0, 0.0);
		q1 = (p1-p0)/scale;
		q4 = (p4-p0)/scale;
		q5 = complex<double>(1.0, 0.0);

		//loop over all four combinations of eta0, eta2 to identify
		//the case that gives the smallest absolute rotation index
		rindex_min = 1.0e12 ;
		for(eta0 = -1; eta0 <= 1; eta0 += 2)
			for(eta2 = -1; eta2 <= 1; eta2 += 2)
			{
				w0 = (double)eta0*sqrt(5.0*(q1-q0));
				w2 = (double)eta2*sqrt(5.0*(q5-q4));
				w1 = -0.75*(w0+w2)+0.25*sqrt(120.0*(q5-q0)-15.0*(w0*w0+w2*w2)+10.0*w0*w2);
				k = w0-2.0*w1+w2;
				if(abs(k) < sml)
				{
					//special case where degree of w(t) < 2
					if(abs(w2-w0) < sml)
						rindex = 0.0;
					else
					{
						a = w0/(w0-w2);
						s1 = abs(a);
						s2 = abs(1.0-a);
						s3 = 1.0;
						anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2));
						rindex = anglea;
					}
				}
				else
				{
					//generic case where degree of w(t) = 2
					a = (w0-w1+sqrt(w1*w1-w0*w2))/k;
					b = (w0-w1-sqrt(w1*w1-w0*w2))/k;
					if(imag(a)*imag(b) > 0.0)
					{ 
						intervals = 1;
						s1 = abs(a); 
						s2 = abs(1.0-a); 
						s3 = 1.0;
						anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2));
						s1 = abs(b); 
						s2 = abs(1.0-b); 
						s3 = 1.0;
						angleb = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2));
						rindex = anglea+angleb;
					}
					else
					{
						intervals = 0; 
						t[0] = 0.0;
						c2 = imag(a+b) ;
						c1 = -2.0*imag(a*b);
						c0 = abs(a)*abs(a)*imag(b)+abs(b)*abs(b)*imag(a);
						if(abs(c2) > sml)
						{
							t1 = (-c1-sqrt(c1*c1-4.0*c2*c0))/(2.0*c2);
							t2 = (-c1+sqrt(c1*c1-4.0*c2*c0))/(2.0*c2);
							if(t2 < t1) 
							{ 
								temp = t2; 
								t2 = t1; 
								t1 = temp; 
							}
							if(t1 > 0.0 && t1 < 1.0) 
								t[++intervals] = t1;
							if(t2 > 0.0 && t2 < 1.0) 
								t[++intervals] = t2;
						}
						else if(abs(c1) > sml)
						{
							t1 = 0.0 - c0/c1;
							if(t1 > 0.0 && t1 < 1.0)
								t[++intervals] = t1;
						}
						t[++intervals] = 1.0;
						rindex = 0.0;
						for(i = 1; i <= intervals; ++i)
						{
							s1 = abs(a-t[i]); 
							s2 = abs(t[i+1]-a); 
							s3 = t[i+1]-t[i];
							anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2));
							s1 = abs(b-t[i]); 
							s2 = abs(t[i+1]-b); 
							s3 = t[i+1]-t[i];
							angleb = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2));
							rindex += fabs(anglea-angleb);
						}
					}
				}
				if(rindex < rindex_min)
				{ 
					eta0_min = eta0; 
					eta2_min = eta2; 
					rindex_min = rindex; 
				}
			}
		w0 = (double)eta0_min*sqrt(5.0*(q1-q0));
		w2 = (double)eta2_min*sqrt(5.0*(q5-q4));
		w1 = -0.75*(w0+w2)+0.25*sqrt(120.0*(q5-q0)-15.0*(w0*w0+w2*w2)+10.0*w0*w2);
		w0 *= sqrt(scale);
		w1 *= sqrt(scale);
		w2 *= sqrt(scale);
		}
	spline[1].w[0] = w0; 
	spline[1].w[1] = w1; 
	spline[1].w[2] = w2;
	spline[1].p[0] = p0 ;
	spline[1].p[1] = spline[1].p[0]+0.2*w0*w0;
	spline[1].p[2] = spline[1].p[1]+0.2*w0*w1;
	spline[1].p[3] = spline[1].p[2]+0.2*(2.0*w1*w1+w0*w2)/3.0;
	spline[1].p[4] = spline[1].p[3]+0.2*w1*w2;
	spline[1].p[5] = spline[1].p[4]+0.2*w2*w2;
	spline[1].sigma[0] = norm(w0);
	spline[1].sigma[1] = real(w0*conj(w1));
	spline[1].sigma[2] = (2.0*norm(w1)+real(w0*conj(w2)))/3.0;
	spline[1].sigma[3] = real(w1*conj(w2));
	spline[1].sigma[4] = norm(w2);
	spline[1].s[0] = 0.0;
	for(i = 1; i <= 5; i++)
		spline[1].s[i] = spline[1].s[i-1]+0.2*spline[1].sigma[i-1];
}

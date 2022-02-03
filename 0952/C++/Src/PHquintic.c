#include<stdio.h> 
#include<math.h> 
#include<complex.h> 

/* Author: Rida T. Farouki, University of California, Davis (2014)

This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. */
 
#define MAXPOINTS 1001   /* maximum number of points for PH spline */

#define EPSILON 0.000000000001 /* convergence tolerance for construction
				  of PH quintic splines */

#define MAXDEGREE 100 /* maximum polynomial degree for beval() function */

#define DEBUG 1 /* set to 1 to print run-time information, 0 otherwise */

struct PHquintic { 

  complex double p[6] ; /* Bezier control points of PH quintic */

  complex double w[3] ; /* Bernstein coefficients of w polynomial */ 

  double sigma[5] ; /* parametric speed Bernstein coefficients */ 

  double s[6] ; /* arc length Bernstein coefficients */ 

} ; 
 
main()

{
  struct PHquintic curve ;  /* struct to define single PH quintic curve */

  struct PHquintic spline[MAXPOINTS] ; /* struct array for PH quintic spline */

  void construct_PHquintic() ; /* constructs planar PH quintic from given
				  initial and final control point pairs */

  void open_PHquintic_spline() ; /* constructs open PH quintic spline */

  void closed_PHquintic_spline() ; /* constructs closed PH quintic spline */

  void PHquintic_offset() ; /* computes offset to a PH quintic curve */

  double PHquintic_energy() ; /* computes bending energy of a PH quintic */

  complex double p0 , p1 , p4 , p5 ; /* initial and final pairs of control 
					points for constructing PH quintic */

  complex double q[MAXPOINTS] ; /* array for spline interpolation points */

  int n ; /* n+1 is number of spline interpolation points q[0],...,q[n] */

  double d ; /* offset distance to a PH quintic curve */

  double W[10] , X[10] , Y[10] ; /* offset curve homogeneous coordinates */

  double energy ; /* elastic bending energy of PH quintic curve */

  int i , j ; /* general purpose loop counters */

  double beval() ; /* function to evaluate Bernstein-form polynomials */

  /* construct the planar PH quintic with the assigned control points */

  p0 = 0.8 + 1.4*I ;
  p1 = 2.4 + 0.8*I ;
  p4 = 2.6 + 4.2*I ; 
  p5 = 4.2 + 3.6*I ;

  construct_PHquintic(p0,p1,p4,p5,&curve) ;

  if ( DEBUG )
    { printf("control points of PH quintic: \n\n") ;
      for ( i = 0 ; i <= 5 ; ++i ) 
	printf("p[%d] = (%f,%f) \n",i,creal(curve.p[i]),cimag(curve.p[i])) ; }

  /* compute bending energy of the planar PH quintic */

  energy = PHquintic_energy(&curve) ;

  if ( DEBUG ) printf("\nenergy of PH quintic = %f \n",energy) ;

  /* compute offset curve to planar PH quintic */

  d = 1.25 ;

  PHquintic_offset(d,&curve,W,X,Y) ;

  if ( DEBUG )
    { printf("\ncontrol points of rational offset curve: \n\n") ;
      for ( i = 0 ; i <= 9 ; ++i ) 
	printf("P[%d] = (%f,%f,%f) \n",i,W[i],X[i],Y[i]) ; }

  /* compute an open C^2 PH quintic spline curve */

  n = 6 ;

  q[0] = -2.1 + 1.8*I ;
  q[1] = -3.1 + 0.0*I ;
  q[2] = -0.3 - 0.8*I ;
  q[3] =  0.7 + 2.2*I ;
  q[4] =  3.4 + 0.5*I ;
  q[5] =  1.1 - 0.6*I ;
  q[6] =  2.3 - 2.4*I ;

  open_PHquintic_spline(n,q,spline) ;

  if ( DEBUG )
    { 
      printf("\ncontrol points of open PH quintic spline: \n") ;
      for ( i = 1 ; i <= n ; ++i ) 
	{ printf("\nsegment %d \n",i) ;
	  for ( j = 0 ; j <= 5 ; ++j )
	    printf("p[%d] = (%f,%f) \n",
		   j,creal(spline[i].p[j]),cimag(spline[i].p[j])) ; }
    }

  /* compute a closed C^2 PH quintic spline curve */

  n = 9 ;

  q[0] = -4.1 - 0.8*I ;
  q[1] = -1.5 - 1.5*I ;
  q[2] = -0.6 - 3.6*I ;
  q[3] =  1.2 - 1.5*I ;
  q[4] =  4.1 + 0.4*I ;
  q[5] =  1.2 + 3.3*I ;
  q[6] =  0.9 + 0.4*I ;
  q[7] = -1.4 - 0.2*I ;
  q[8] = -2.3 + 1.7*I ;
  q[9] = -4.1 - 0.8*I ;

  closed_PHquintic_spline(n,q,spline) ;

  if ( DEBUG )
    { 
      printf("\ncontrol points of closed PH quintic spline: \n") ;
      for ( i = 1 ; i <= n ; ++i ) 
	{ printf("\nsegment %d \n",i) ;
	  for ( j = 0 ; j <= 5 ; ++j )
	    printf("p[%d] = (%f,%f) \n",
		   j,creal(spline[i].p[j]),cimag(spline[i].p[j])) ; }
    }

}

/********************************************************************/

void construct_PHquintic( complex double p0 , 
			  complex double p1 ,
 			  complex double p4 , 
			  complex double p5 ,
			  struct PHquintic *curve )

/* Given the first and last pairs of Bezier control points p0,p1 
and p4,p5 this function constructs the middle control points p2,p3 
such that the resulting curve is a planar PH quintic.

The complete set of Bezier control points is returned in struct
PHquintic "curve" together with the Bernstein coefficients of the 
w(t) polynomial, the parametric speed polynomial sigma(t), and the
arc length polynomial s(t).

Among the four formal solutions, the "good" solution is identified 
as that which minimizes the absolute rotation index. 

The initial and final points, p0 and p5, are assumed to be distinct.
Otherwise, the function returns a degenerate curve: the single point 
p0=p5.

The parameter "small" is used to determine when w(t) degenerates to a 
linear polynomial, and when the numerator of the curvature degenerates 
to a linear polynomial. */

{
  int eta0 , eta2 , eta0_min , eta2_min , i , intervals ;

  complex double q0 , q1 , q4 , q5 , scale , w0 , w1 , w2 , k , a , b ;

  double s1 , s2 , s3 , c0 , c1 , c2 , anglea , angleb ;

  double t1 , t2 , temp , t[4] , rindex , rindex_min ;

  double small = 0.000000000001 ;

  if ( p5 == p0 )
    {
      w0 = w1 = w2 = 0.0 ;
    }
  else
    {
      /* scale input data (p0,p1,p4,p5) to canonical form (q0,q1,q4,q5) */

      scale = p5-p0 ;

      q0 = 0.0+0.0*I ; 
      q1 = (p1-p0)/scale ;
      q4 = (p4-p0)/scale ; 
      q5 = 1.0+0.0*I ;

      /* loop over all four combinations of eta0, eta2 to identify
	 the case that gives the smallest absolute rotation index */

      rindex_min = 1.0e12 ;

      for ( eta0 = -1 ; eta0 <= 1 ; eta0 += 2 )
	for ( eta2 = -1 ; eta2 <= 1 ; eta2 += 2 )
	  {
	    w0 = eta0*csqrt(5.0*(q1-q0)) ;
	    w2 = eta2*csqrt(5.0*(q5-q4)) ;
	    
	    w1 = - 0.75*(w0+w2) 
	      + 0.25*csqrt(120.0*(q5-q0)-15.0*(w0*w0+w2*w2)+10.0*w0*w2) ;

	    k = w2 - 2.0*w1 + w0 ;

	    if ( cabs(k) < small )
	      {
		/* special case where degree of w(t) < 2 */

		if ( cabs(w2-w0) < small )
		  rindex = 0.0 ;
		else
		  {
		    a = w0/(w0-w2) ;
		    s1 = cabs(a) ; 
		    s2 = cabs(1.0-a) ; 
		    s3 = 1.0 ;
		    anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2)) ;
		    rindex = anglea ;		    
		  }
	      }
	    else
	      {
		/* generic case where degree of w(t) = 2 */

		a = ( w0-w1 + csqrt(w1*w1-w0*w2) ) / k ;
		b = ( w0-w1 - csqrt(w1*w1-w0*w2) ) / k ;

		if ( cimag(a)*cimag(b) > 0.0 )
		  {
		    intervals = 1 ;
		    s1 = cabs(a) ; 
		    s2 = cabs(1.0-a) ; 
		    s3 = 1.0 ;
		    anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2)) ;
		    s1 = cabs(b) ; 
		    s2 = cabs(1.0-b) ; 
		    s3 = 1.0 ;
		    angleb = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2)) ;
		    rindex = anglea+angleb ;
		  }
		else
		  {
		    intervals = 0 ; t[0] = 0.0 ;
		    c2 = cimag(a+b) ;
		    c1 = -2.0*cimag(a*b) ;
		    c0 = cabs(a)*cabs(a)*cimag(b)+cabs(b)*cabs(b)*cimag(a) ;
		    if ( fabs(c2) > small )
		      {
			t1 = (-c1-sqrt(c1*c1-4.0*c2*c0))/(2.0*c2) ;
			t2 = (-c1+sqrt(c1*c1-4.0*c2*c0))/(2.0*c2) ;
			if ( t2 < t1 ) { temp = t2 ; t2 = t1 ; t1 = temp ; }
			if ( t1 > 0.0 && t1 < 1.0 ) t[++intervals] = t1 ;
			if ( t2 > 0.0 && t2 < 1.0 ) t[++intervals] = t2 ;
		      }
		    else if ( fabs(c1) > small )
		      {
			t1 = -c0/c1 ;
			if ( t1 > 0.0 && t1 < 1.0 ) t[++intervals] = t1 ;
		      }
		    t[++intervals] = 1.0 ;
		    rindex = 0.0 ;
		    for ( i = 0 ; i < intervals ; ++i )
		      {
			s1 = cabs(a-t[i]) ; 
			s2 = cabs(t[i+1]-a) ; 
			s3 = t[i+1]-t[i] ;
			anglea = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2)) ;
			s1 = cabs(b-t[i]) ; 
			s2 = cabs(t[i+1]-b) ; 
			s3 = t[i+1]-t[i] ;
			angleb = acos((s1*s1+s2*s2-s3*s3)/(2.0*s1*s2)) ;
			rindex += fabs(anglea-angleb) ;
		      }
		  }
	      }
	    
	    if ( rindex < rindex_min )
	      { eta0_min = eta0 ; eta2_min = eta2 ; rindex_min = rindex ; }
	  }
      
      /* compute the good solution, restore original scaling, and
	 store definition of good solution in the struct "curve" */
      
      w0 = eta0_min*csqrt(5.0*(q1-q0)) ;
      w2 = eta2_min*csqrt(5.0*(q5-q4)) ;
      
      w1 = - 0.75*(w0+w2) 
	+ 0.25*csqrt(120.0*(q5-q0)-15.0*(w0*w0+w2*w2)+10.0*w0*w2) ;
      
      w0 *= csqrt(scale) ; w1 *= csqrt(scale) ; w2 *= csqrt(scale) ;
    }

  curve->w[0] = w0 ; 
  curve->w[1] = w1 ; 
  curve->w[2] = w2 ;

  curve->p[0] = p0 ;
  curve->p[1] = curve->p[0] + 0.2*w0*w0 ;
  curve->p[2] = curve->p[1] + 0.2*w0*w1 ;
  curve->p[3] = curve->p[2] + 0.2*(2.0*w1*w1+w0*w2)/3.0 ;
  curve->p[4] = curve->p[3] + 0.2*w1*w2 ;
  curve->p[5] = curve->p[4] + 0.2*w2*w2 ;

  curve->sigma[0] = w0*conj(w0) ;
  curve->sigma[1] = creal(w0*conj(w1)) ;
  curve->sigma[2] = (2.0*w1*conj(w1)+creal(w0*conj(w2)))/3.0 ;
  curve->sigma[3] = creal(w1*conj(w2)) ;
  curve->sigma[4] = w2*conj(w2) ;

  curve->s[0] = 0.0 ;
  for ( i = 1 ; i <= 5 ; ++i )
    curve->s[i] = curve->s[i-1] + 0.2*curve->sigma[i-1] ;
}      

/********************************************************************/

void open_PHquintic_spline( int n , 
			    complex double q[] , 
			    struct PHquintic spline[] )

/* Computes a planar C^2 PH quintic spline interpolating n+1 points 
with cubic end spans as end conditions. The coordinates of the points 
to be interpolated are represented as the complex-number array elements
q[0],...,q[n]. These points should have reasonably even spacing. The
data for the resulting n PH quintic segments are returned as elements 
spline[1] ... spline[n] of the struct PHquintic spline array. */

{
  int i , k , iter ;

  complex double a[MAXPOINTS] , b[MAXPOINTS] , c[MAXPOINTS] , d[MAXPOINTS] ; 

  complex double dz[MAXPOINTS] , dq[MAXPOINTS] , z[MAXPOINTS+1] ;

  complex double beta[MAXPOINTS] , rho[MAXPOINTS] , w0 , w1 , w2 ;

  double znorm , dznorm , err , cx[6] , cy[6] , beval() ;

  void tridiag_open() ;

  /* first construct the ordinary C^2 cubic spline: the complex values 
     dq[0],..,.,dq[n] are the computed nodal derivatives */

  b[0] = c[0] = 1.0 ;
  
  for ( i = 1 ; i <= n-1 ; ++i ) 
    { a[i] = 1.0 ; b[i] = 4.0 ; c[i] = 1.0 ; }
  
  a[n] = b[n] = 1.0 ;
  
  d[0] = 2.0*(q[1]-q[0]) ;

  for ( i = 1 ; i <= n-1 ; ++i ) 
    d[i] = 3.0*(q[i+1]-q[i-1]) ;
  
  d[n] = 2.0*(q[n]-q[n-1]) ;

  beta[0] = b[0] ; rho[0] = d[0] ; 

  for ( i = 0 ; i <= n-1 ; ++i )
    {
      rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i] ;
      beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i] ;
    }

  dq[n] = rho[n]/beta[n] ;

  for ( i = n-1 ; i >= 0 ; --i )
    dq[i] = (rho[i]-c[i]*dq[i+1])/beta[i] ;

  for ( k = 1 ; k <= n ; ++k )
    {
      cx[0] = creal(q[k-1]) ; 
      cy[0] = cimag(q[k-1]) ;
      cx[1] = creal(q[k-1]+dq[k-1]/3.0) ; 
      cy[1] = cimag(q[k-1]+dq[k-1]/3.0) ;
      cx[2] = creal(q[k]-dq[k]/3.0) ; 
      cy[2] = cimag(q[k]-dq[k]/3.0) ;
      cx[3] = creal(q[k]) ; 
      cy[3] = cimag(q[k]) ;
    }

  /* compute the starting values for the iteration from the mid-point
     derivatives of the ordinary C^2 cubic spline */

  d[1] = 4.0*csqrt(6.0*(q[1]-q[0])-dq[0]-dq[1]) ;
  for ( i = 2 ; i <= n ; ++i )
    {
      d[i] = 4.0*csqrt(6.0*(q[i]-q[i-1])-dq[i-1]-dq[i]) ;
      if ( cabs(d[i-1]+d[i]) < cabs(d[i-1]-d[i]) ) d[i] = -d[i] ;
    }

  b[1] = 8.0 ; c[1] = 0.0 ;
  
  for ( i = 2 ; i <= n-1 ; ++i ) 
    { a[i] = 1.0 ; b[i] = 6.0 ; c[i] = 1.0 ; }
  
  a[n] = 0.0 ; b[n] = 8.0 ;
  
  tridiag_open(n,a,b,c,d,z) ;

  /* now begin the Newton-Raphson iteration to compute the complex 
     values z[0],...,z[n] that defin the C^2 PH quintic spline */

  if ( DEBUG ) printf("\nconvergence of NR iterations \n\n") ;

  iter = 0 ; err = 1.0 ; 
  
  while ( err > EPSILON && iter <= 20 )
    {
      iter += 1 ;

      b[1] = 130.0*z[1]-10.0*z[2] ;
      c[1] = 10.0*z[2]-10.0*z[1] ;
      d[1] = 60.0*(q[1]-q[0])
	-65.0*z[1]*z[1]-5.0*z[2]*z[2]+10.0*z[1]*z[2] ;

      for ( i = 2 ; i <= n-1 ; ++i )
	{
	  a[i] = 6.0*z[i-1]+13.0*z[i]+z[i+1] ;
	  b[i] = 13.0*z[i-1]+54.0*z[i]+13.0*z[i+1] ;
	  c[i] = z[i-1]+13.0*z[i]+6.0*z[i+1] ;
	  d[i] = 60.0*(q[i]-q[i-1]) 
	    -3.0*z[i-1]*z[i-1]-27.0*z[i]*z[i]-3.0*z[i+1]*z[i+1]
	    -z[i-1]*z[i+1]-13.0*z[i-1]*z[i]-13.0*z[i]*z[i+1] ;
	}
      
      a[n] = 10.0*z[n-1]-10.0*z[n] ;
      b[n] = 130.0*z[n]-10.0*z[n-1] ;
      d[n] = 60.0*(q[n]-q[n-1])
	-65.0*z[n]*z[n]-5.0*z[n-1]*z[n-1]+10.0*z[n-1]*z[n] ;

      tridiag_open(n,a,b,c,d,dz) ;

      for ( i = 1 ; i <= n ; ++i ) z[i] += dz[i] ;

      znorm = 0.0 ; dznorm = 0.0 ;

      for ( i = 1 ; i <= n ; ++i )
	{
	  znorm += cabs(z[i])*cabs(z[i]) ;
	  dznorm += cabs(dz[i])*cabs(dz[i]) ;
	}

      err = sqrt(dznorm/znorm) ;

      if ( DEBUG) printf("iter = %d, err = %20.15f \n",iter,err) ;
    }

  z[0] = 2.0*z[1]-z[2] ;

  z[n+1] = 2.0*z[n]-z[n-1] ;

  /* construct the Bezier control points for each segment of the
     C^2 PH quintic spline */

  for ( i = 1 ; i <= n ; ++i )
    {
      w0 = 0.5*(z[i-1]+z[i]) ;
      w1 = z[i] ;
      w2 = 0.5*(z[i]+z[i+1]) ;

      spline[i].w[0] = w0 ; spline[i].w[1] = w1 ; spline[i].w[2] = w2 ;

      spline[i].p[0] = q[i-1] ;
      spline[i].p[1] = spline[i].p[0] + 0.2*w0*w0 ;
      spline[i].p[2] = spline[i].p[1] + 0.2*w0*w1 ;
      spline[i].p[3] = spline[i].p[2] + 0.2*(2.0*w1*w1+w0*w2)/3.0 ;
      spline[i].p[4] = spline[i].p[3] + 0.2*w1*w2 ;
      spline[i].p[5] = spline[i].p[4] + 0.2*w2*w2 ;

      spline[i].sigma[0] = w0*conj(w0) ;
      spline[1].sigma[1] = creal(w0*conj(w1)) ;
      spline[i].sigma[2] = (2.0*w1*conj(w1)+creal(w0*conj(w2)))/3.0 ;
      spline[i].sigma[3] = creal(w1*conj(w2)) ;
      spline[i].sigma[4] = w2*conj(w2) ;

      spline[i].s[0] = 0.0 ;
      for ( k = 1 ; k <= 5 ; ++k )
	spline[i].s[k] = spline[i].s[k-1] + 0.2*spline[i].sigma[k-1] ;
    }
}      

/********************************************************************/

void closed_PHquintic_spline( int n , 
			      complex double q[] , 
			      struct PHquintic spline[] )

/* Computes a planar C^2 PH quintic spline interpolating n+1 points with
periodic end conditions. The coordinates of the points to be interpolated 
are represented as the complex-number array elements q[0],...,q[n] where 
q[n]=q[0]. These points should have reasonably even spacing. The data for 
the resulting n PH quintic segments are returned as elements spline[1] 
... spline[n] of the struct PHquintic spline array. */

{
  int i , k , eta , iter ;

  complex double a[MAXPOINTS] , b[MAXPOINTS] , c[MAXPOINTS] , d[MAXPOINTS] ; 

  complex double z[MAXPOINTS+1] , dz[MAXPOINTS] , dq[MAXPOINTS] ;

  complex double beta[MAXPOINTS] , rho[MAXPOINTS] , zeta[MAXPOINTS] ; 

  complex double theta , w0 , w1 , w2 ;

  double znorm , dznorm , err , cx[6] , cy[6] , beval() ;

  void tridiag_closed() ;

  /* first construct the ordinary C^2 cubic spline: the complex values 
     dq[0],..,.,dq[n] are the computed nodal derivatives */

  for ( i = 0 ; i <= n-1 ; ++i ) 
    { a[i] = 1.0 ; b[i] = 4.0 ; c[i] = 1.0 ; }

  d[0] = 3.0*(q[1]-q[n-1]) ;

  for ( i = 1 ; i <= n-1 ; ++i ) 
    d[i] = 3.0*(q[i+1]-q[i-1]) ;
  
  beta[0] = b[0] ; rho[0] = d[0] ; zeta[0] = a[0] ;

  for ( i = 0 ; i <= n-3 ; ++i )
    {
      beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i] ;
      rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i] ;
      zeta[i+1] = -a[i+1]*zeta[i]/beta[i] ;
    }

  zeta[n-2] += c[n-2] ; 
  rho[n-1] = d[n-1]-a[n-1]*rho[n-2]/beta[n-2] ;
  beta[n-1] = b[n-1]-a[n-1]*zeta[n-2]/beta[n-2] ;

  theta = c[n-1] ;

  for ( i = 0 ; i <= n-2 ; ++i )
    {
      beta[n-1] -= theta*zeta[i]/beta[i] ;
      rho[n-1] -= theta*rho[i]/beta[i] ;
      theta *= -c[i]/beta[i] ;
    }

  dq[n-1] = rho[n-1]/beta[n-1] ;
  dq[n-2] = (rho[n-2]-zeta[n-2]*dq[n-1])/beta[n-2] ;

  for ( i = n-3 ; i >= 0 ; --i )
    dq[i] = (rho[i]-c[i]*dq[i+1]-zeta[i]*dq[n-1])/beta[i] ;

  q[n] = q[0] ; dq[n] = dq[0] ;

  for ( k = 1 ; k <= n ; ++k )
    {
      cx[0] = creal(q[k-1]) ; 
      cy[0] = cimag(q[k-1]) ;
      cx[1] = creal(q[k-1]+dq[k-1]/3.0) ; 
      cy[1] = cimag(q[k-1]+dq[k-1]/3.0) ;
      cx[2] = creal(q[k]-dq[k]/3.0) ; 
      cy[2] = cimag(q[k]-dq[k]/3.0) ;
      cx[3] = creal(q[k]) ; 
      cy[3] = cimag(q[k]) ;
    }

  /* compute the starting values for the iteration from the mid-point
     derivatives of the ordinary C^2 cubic spline */

  d[1] = 4.0*csqrt(6.0*(q[1]-q[0])-dq[0]-dq[1]) ;
  for ( i = 2 ; i <= n ; ++i )
    {
      d[i] = 4.0*csqrt(6.0*(q[i]-q[i-1])-dq[i-1]-dq[i]) ;
      if ( creal(d[i-1]*conj(d[i])) < 0.0 ) d[i] = -d[i] ;
    }

  eta = 1 ;

  if ( creal(d[n]*conj(d[1])) < 0.0 ) eta = -1 ;

  a[1] = eta ; b[1] = 6.0 ; c[1] = 1.0 ;
  
  for ( i = 2 ; i <= n-1 ; ++i ) 
    { a[i] = 1.0 ; b[i] = 6.0 ; c[i] = 1.0 ; }
  
  a[n] = 1.0 ; b[n] = 6.0 ; c[n] = eta ;
  
  tridiag_closed(n,a,b,c,d,z) ;

  /* now begin the Newton-Raphson iteration to compute the complex 
     values z[0],...,z[n] that defin the C^2 PH quintic spline */

  if ( DEBUG ) printf("\nconvergence of NR iterations \n\n") ;

  iter = 0 ; err = 1.0 ;
  
  while ( err > EPSILON && iter <= 6 )
    {
      iter += 1 ;

      a[1] = 6.0*z[n]+eta*(13.0*z[1]+z[2]) ;
      b[1] = eta*13.0*z[n]+54.0*z[1]+13.0*z[2] ;
      c[1] = eta*z[n]+13.0*z[1]+6.0*z[2] ;
      d[1] = 60.0*(q[1]-q[0])
	-3.0*z[n]*z[n]-27.0*z[1]*z[1]-3.0*z[2]*z[2]
	-eta*(z[n]*z[2]+13.0*z[n]*z[1])-13.0*z[1]*z[2] ;

      for ( i = 2 ; i <= n-1 ; ++i )
	{
	  a[i] = 6.0*z[i-1]+13.0*z[i]+z[i+1] ;
	  b[i] = 13.0*z[i-1]+54.0*z[i]+13.0*z[i+1] ;
	  c[i] = z[i-1]+13.0*z[i]+6.0*z[i+1] ;
	  d[i] = 60.0*(q[i]-q[i-1]) 
	    -3.0*z[i-1]*z[i-1]-27.0*z[i]*z[i]-3.0*z[i+1]*z[i+1]
	    -z[i-1]*z[i+1]-13.0*z[i-1]*z[i]-13.0*z[i]*z[i+1] ;
	}
      
      a[n] = 6.0*z[n-1]+13.0*z[n]+eta*z[1] ;
      b[n] = 13.0*z[n-1]+54.0*z[n]+eta*13.0*z[1] ;
      c[n] = eta*(z[n-1]+13.0*z[n])+6.0*z[1] ;
      d[n] = 60.0*(q[n]-q[n-1])
	-3.0*z[n-1]*z[n-1]-27.0*z[n]*z[n]-3.0*z[1]*z[1]
	-eta*(z[n-1]*z[1]+13.0*z[n]*z[1])-13.0*z[n-1]*z[n] ;

      tridiag_closed(n,a,b,c,d,dz) ;

      for ( i = 1 ; i <= n ; ++i ) z[i] += dz[i] ;

      znorm = 0.0 ; dznorm = 0.0 ;

      for ( i = 1 ; i <= n ; ++i )
	{
	  znorm += cabs(z[i])*cabs(z[i]) ;
	  dznorm += cabs(dz[i])*cabs(dz[i]) ;
	}

      err = sqrt(dznorm/znorm) ;

      if ( DEBUG ) printf("iter = %d, err = %20.15f \n",iter,err) ;
    }

  z[0] = eta*z[n] ; z[n+1] = eta*z[1] ;

  /* construct the Bezier control points for each segment of the
     C^2 PH quintic spline */

  for ( i = 1 ; i <= n ; ++i )
    {
      w0 = 0.5*(z[i-1]+z[i]) ;
      w1 = z[i] ;
      w2 = 0.5*(z[i]+z[i+1]) ;

      spline[i].w[0] = w0 ; spline[i].w[1] = w1 ; spline[i].w[2] = w2 ;

      spline[i].p[0] = q[i-1] ;
      spline[i].p[1] = spline[i].p[0] + 0.2*w0*w0 ;
      spline[i].p[2] = spline[i].p[1] + 0.2*w0*w1 ;
      spline[i].p[3] = spline[i].p[2] + 0.2*(2.0*w1*w1+w0*w2)/3.0 ;
      spline[i].p[4] = spline[i].p[3] + 0.2*w1*w2 ;
      spline[i].p[5] = spline[i].p[4] + 0.2*w2*w2 ;

      spline[i].sigma[0] = w0*conj(w0) ;
      spline[1].sigma[1] = creal(w0*conj(w1)) ;
      spline[i].sigma[2] = (2.0*w1*conj(w1)+creal(w0*conj(w2)))/3.0 ;
      spline[i].sigma[3] = creal(w1*conj(w2)) ;
      spline[i].sigma[4] = w2*conj(w2) ;

      spline[i].s[0] = 0.0 ;
      for ( k = 1 ; k <= 5 ; ++k )
	spline[i].s[k] = spline[i].s[k-1] + 0.2*spline[i].sigma[k-1] ;
    }
}      

/********************************************************************/

void PHquintic_offset( double d , 
		       struct PHquintic *curve ,
		       double W[] , double X[] , double Y[] )

/* Computes the offset curve at distance d from the PH quintic
defined in the struct "curve". The positive sense of the normal
vector is assumed to be to the right when traversing the curve 
with increasing parameter value.

The homogeneous coordinates for the weights and control points
of the degree 9 rational Bezier curve defining the offset are
returned in the arrays W[], X[], Y[]. */

{
  int i , j , jmin , jmax , k , binom[10][10] ;

  double P[6][3] , deltaP[6][3] , f ;

  binom[0][0] = 1 ;
  for ( k = 1 ; k <= 9 ; ++k ) 
    binom[k][0] = binom[k][k] = 1 ;
  for ( k = 2 ; k <= 9 ; ++k )
    for ( j = 1 ; j <= k-1 ; ++j ) 
      binom[k][j] = binom[k-1][j-1] + binom[k-1][j] ;

  for ( k = 0 ; k <= 5 ; ++k )
    {
      P[k][0] = 1.0 ;
      P[k][1] = creal(curve->p[k]) ;
      P[k][2] = cimag(curve->p[k]) ;
    }

  for ( k = 0 ; k <= 4 ; ++k )
    for ( j = 0 ; j <= 2 ; ++j )
      deltaP[k][j] = P[k+1][j]-P[k][j] ;

  for ( k = 0 ; k <= 9 ; ++k )
    {
      W[k] = X[k] = Y[k] = 0.0 ;

      jmax = 4 ; if ( k < jmax ) jmax = k ;
      jmin = 0 ; if ( k-5 > 0 ) jmin = k-5 ;

      for ( j = jmin ; j <= jmax ; ++j )
	{
	  f = ((double) binom[k][j]*binom[9-k][4-j]) / 126.0 ;
	  W[k] += f*curve->sigma[j] ;
	  X[k] += f*(curve->sigma[j]*P[k-j][1]+5.0*d*deltaP[j][2]) ;
	  Y[k] += f*(curve->sigma[j]*P[k-j][2]-5.0*d*deltaP[j][1]) ;
	}
    }
}

/********************************************************************/

double PHquintic_energy( struct PHquintic *curve )

/* Computes the bending energy of the planar PH quintic specified 
in struct PH quintic "curve". The energy is returned as the value of
the function. The parameter "small" is used to determine when w(t) 
degenerates to a linear polynomial. */

{
  complex double w0 , w1 , w2 , k , a , b , kc , ac , bc ;

  complex double a1 , a2 , a3 , b1 , b2 , b3 ;

  double alpha , beta , energy ;

  double f0 , f1 ;

  double small = 0.000000000001 ;

  w0 = curve->w[0] ; w1 = curve->w[1] ; w2 = curve->w[2] ;

  k = w2 - 2.0*w1 + w0 ;

  if ( fabs(k) < small )
    {
      if ( cabs(w2-w0) < small )
	energy = 0.0 ;
      else
	{
	  k = w2-w0 ;
	  a = w0/(w0-w2) ;
	  alpha = cimag(a) ;
	  kc = conj(k) ;
	  f0 = cabs(a)*cabs(a) ;
	  f1 = cabs(1.0-a)*cabs(1.0-a) ;

	  energy = creal(a)/(f0*f0) + creal(1.0-a)/(f1*f1)
	    + 1.5*(creal(a)/f0+creal(1.0-a)/f1)/(alpha*alpha)
	    + 1.5*(atan((1.0-creal(a))/alpha)+atan(creal(a)/alpha))
	    /(alpha*alpha*alpha) ;

	  energy /= k*kc ;
	}
    }
  else
    {
      a = ( w0-w1 + csqrt(w1*w1-w0*w2) ) / k ;
      b = ( w0-w1 - csqrt(w1*w1-w0*w2) ) / k ;

      alpha = cimag(a) ; beta = cimag(b) ;

      ac = conj(a) ; bc = conj(b) ; kc = conj(k) ;

      a3 = I / (8.0*alpha*(a-b)*(a-bc)) ;
      b3 = I / (8.0*beta*(a-b)*(ac-b)) ;

      a2 = ( 1.5*I/alpha + 1.0/(a-b) - 3.0/(a-bc) )*a3 ;
      b2 = ( 1.5*I/beta - 1.0/(a-b) + 3.0/(ac-b) )*b3 ;

      a1 = 1.5*I*a2/alpha + ( 0.75/(alpha*alpha) - 2.0/((a-b)*(a-b))
	+ 6.0/((a-bc)*(a-bc)) - (1.0-2.0*beta/alpha)/((a-b)*(a-bc)) ) * a3 ;
      b1 = 1.5*I*b2/beta + ( 0.75/(beta*beta) - 2.0/((a-b)*(a-b))
	+ 6.0/((ac-b)*(ac-b)) - (1.0-2.0*alpha/beta)/((a-b)*(ac-b)) ) * b3 ;

      energy = 2.0*creal(a1)*log(cabs((1.0-a)/a))
	     - 2.0*cimag(a1)*(carg(1.0-a)-carg(-a))
	     + 2.0*creal(b1)*log(cabs((1.0-b)/b))
	     - 2.0*cimag(b1)*(carg(1.0-b)-carg(-b))
	     - creal( 2.0*a2/(a*(1.0-a)) + 2.0*b2/(b*(1.0-b))
		   + (2.0*a-1.0)*a3/(a*a*(1.0-a)*(1.0-a))
		   + (2.0*b-1.0)*b3/(b*b*(1.0-b)*(1.0-b)) ) ; 

      energy *= 4.0/(k*kc) ;
    }

  return(energy) ;
}

/********************************************************************/

void tridiag_open( int n , 
		   complex double a[] , 
		   complex double b[] , 
		   complex double c[] , 
		   complex double d[] , 
		   complex double x[] )

/* Solves tridiagonal system appropriate to open PH quintic spline, with
cubic end spans. The arrays a[], b[], c[] define the lower, main, and upper
diagonal matrix elements, and the array d[] defines the right hand side
values. The solutions are returned in the array x[].  */

{
  int i ;

  complex double beta[MAXPOINTS] , rho[MAXPOINTS] ;

  beta[1] = b[1] ; rho[1] = d[1] ;

  for ( i = 1 ; i <= n-1 ; ++i )
    {
      beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i] ;
      rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i] ;
    }

  x[n] = rho[n]/beta[n] ;

  for ( i = n-1 ; i >= 1 ; --i )
    {
      x[i] = (rho[i]-c[i]*x[i+1])/beta[i] ;
    }
}

/********************************************************************/

void tridiag_closed( int n , 
		     complex double a[] , 
		     complex double b[] , 
		     complex double c[] , 
		     complex double d[] , 
		     complex double x[] )

/* Solves tridiagonal system appropriate to closed PH quintic spline, with
periodic end conditions. The arrays a[], b[], c[] define the lower, main, 
and upper diagonal matrix elements, and the array d[] defines the right hand 
side values. The solutions are returned in the array x[].  */

{
  int i ;

  complex double beta[MAXPOINTS] , rho[MAXPOINTS] , zeta[MAXPOINTS] , theta ;

  beta[1] = b[1] ; rho[1] = d[1] ; zeta[1] = a[1] ;

  for ( i = 1 ; i <= n-2 ; ++i )
    {
      beta[i+1] = b[i+1]-a[i+1]*c[i]/beta[i] ;
      rho[i+1] = d[i+1]-a[i+1]*rho[i]/beta[i] ;
      zeta[i+1] = -a[i+1]*zeta[i]/beta[i] ;
    }

  zeta[n-1] += c[n-1] ;

  beta[n] = b[n]-a[n]*zeta[n-1]/beta[n-1] ;
  rho[n] = d[n]-a[n]*rho[n-1]/beta[n-1] ;
  theta = c[n] ;
  
  for ( i = 1 ; i <= n-1 ; ++i )
    {
      beta[n] -= theta*zeta[i]/beta[i] ;
      rho[n] -= theta*rho[i]/beta[i] ;
      theta *= -c[i]/beta[i] ;
    }

  x[n] = rho[n]/beta[n] ;
  x[n-1] = (rho[n-1]-zeta[n-1]*x[n])/beta[n-1] ;
  
  for ( i = n-2 ; i >= 1 ; --i )
    {
      x[i] = (rho[i]-c[i]*x[i+1]-zeta[i]*x[n])/beta[i] ;
    }
}

/********************************************************************/

double beval( int n , 
	      double b[] , 
	      double t )

/* Computes value of a degree-n polynomial with Bernstein coefficients 
b[0],...,b[n] at value t of the independent variable by the de Casteljau 
algorithm. Maximum allowed degree is n=100. */

{
  int j , k ; 
  double blast[MAXDEGREE+1] , bnext[MAXDEGREE+1] ;

  if ( n == 0 ) return(b[0]) ;

  for ( k = 0 ; k <= n ; ++k ) blast[k] = b[k] ;

  for ( j = 1 ; j <= n ; ++j )
    { 
      for ( k = j ; k <= n ; ++k )
	bnext[k] = (1.0-t)*blast[k-1] + t*blast[k] ;
      for ( k = j ; k <= n ; ++k ) blast[k] = bnext[k] ;
    }
  
  return(bnext[n]) ;
}

/********************************************************************/

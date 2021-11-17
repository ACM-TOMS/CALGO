#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>

// ----------------------------------------------------------------------------------
// MODULE Polynomial234RootSolvers
//
//    This module contains the three root solvers for quadratic,
//    cubic and quartic polynomials.
//
// ----------------------------------------------------------------------------------
void solve_cubic_analytic(double *coeff, complex double sol[3]);

void quadraticRoots (double q1, double q0, int *nReal, complex double root[2]);

double max2(double a, double b)
{
  if (a >= b)
    return a;
  else
    return b;
}
 double max3(double a, double b, double c)
{
  double t;
  t = max2(a,b);
  return max2(t,c);
}
 double max4(double a, double b, double c, double d)
{
  double t;
  t = max3(a,b,c);
  return max2(t,d);
}
 double min2(double a, double b)
{
  if (a <= b)
    return a;
  else
    return b;
}
 double min3(double a, double b, double c)
{
  double t;
  t = min2(a,b);
  return min2(t,c);
}
 double min4(double a, double b, double c, double d)
{
  double t;
  t = min3(a,b,c);
  return min2(t,d);
}
//
// CUBIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call cubicRoots (real,              intent (in)  :: c2,
//                   real,              intent (in)  :: c1,
//                   real,              intent (in)  :: c0,
//                   integer,           intent (out) :: nReal,
//                   real,              intent (out) :: root (1:3,1:2),
//                   logical, optional, intent (in)  :: printInfo)
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the cubic polynomial:
//
//                 x^3 + c2 * x^2 + c1 * x + c0
//
//  The first real root (which always exists) is obtained using an optimized
//  Newton-Raphson scheme. The other remaining roots are obtained through
//  composite deflation into a quadratic. An option for printing detailed info
//  about the intermediate stages in solving the cubic is available.
//
//  The cubic root solver can handle any size of cubic coefficients and there is
//  no danger of overflow due to proper rescaling of the cubic polynomial.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) Since there can be only one complex conjugate pair root, no order
//           is necessary.
//
//        3) All real roots preceede the complex ones.
//
// ARGUMENTS
//
//  c2         : coefficient of x^2 term
//  c1         : coefficient of x term
//  c0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//  printInfo  : if given and true, detailed info will be printed about intermediate stages
//
// NOTES
//
//***

void cubicRoots (double c2, double c1, double c0, int *nReal, complex double root[3])
{
  int bisection, converged;
  int cubicType, deflateCase, oscillate;

  enum costanti {allzero = 0, linear = 1, quadratic = 2, general   = 3};

  double a0, a1, a2, a=0, b=0, c=0, k, s=0, t=0, u=0, x, y, z, xShift=0.0;

  double macheps =2.2204460492503131E-16, one27th = 1.0 / 27.0;
  double two27th = 2.0 / 27.0, third   = 1.0 /  3.0;
  double p1 = 1.09574,q1 = 3.23900E-1,r1=3.23900e-1;      //
  double s1 = 9.57439E-2;      //

  double p3 = 1.14413;         //
  double q3 = 2.75509E-1;      // Newton-Raphson coeffs for class 3
  double r3 = 4.45578E-1;      //
  double s3 = 2.59342E-2;      //

  double q4 = 7.71845E-1;      // Newton-Raphson coeffs for class 4
  double s4 = 2.28155E-1;      //

  double p51 = 8.78558E-1;     //
  double p52 = 1.92823E-1;     //
  double p53 = 1.19748;        //
  double p54 = 3.45219E-1;     //
  double q51 = 5.71888E-1;     //
  double q52 = 5.66324E-1;     //
  double q53 = 2.83772E-1;     // Newton-Raphson coeffs for class 5 and 6
  double q54 = 4.01231E-1;     //
  double r51 = 7.11154E-1;     //
  double r52 = 5.05734E-1;     //
  double r53 = 8.37476E-1;     //
  double r54 = 2.07216E-1;     //
  double s51 = 3.22313E-1;     //
  double s52 = 2.64881E-1;     //
  double s53 = 3.56228E-1;     //
  double s54 = 4.45532E-3;     //
  //
  //
  //     ...Start.
  //
  //
  //
  //
  //     ...Handle special cases.
  //
  //            1) all terms zero
  //            2) only quadratic term is nonzero -> linear equation.
  //            3) only independent term is zero -> quadratic equation.
  //
  //
  if (c0 == 0.0 && c1 == 0.0 && c2 == 0.0) {

    cubicType = allzero;

  } else if (c0 == 0.0 && c1 == 0.0) {

    k  = 1.0;
    a2 = c2;

    cubicType = linear;

  } else if (c0 == 0.0) {

    k  = 1.0;
    a2 = c2;
    a1 = c1;

    cubicType = quadratic;

  } else {
    //
    //
    //     ...The general case. Rescale cubic polynomial, such that largest fabsolute coefficient
    //        is (exactly!) equal to 1. Honor the presence of a special cubic case that might have
    //        been obtained during the rescaling process (due to underflow in the coefficients).
    //
    //
    x = fabs (c2);
    y = sqrt (fabs (c1));
    z = pow(fabs (c0), third);
    u = max3(x,y,z);

    if (u == x) {

      k  = 1.0 / x;
      a2 = copysign (1.0 , c2);
      a1 = (c1 * k) * k;
      a0 = ((c0 * k) * k) * k;

    } else if (u == y) {

      k  = 1.0 / y;
      a2 = c2 * k;
      a1 = copysign (1.0 , c1);
      a0 = ((c0 * k) * k) * k;

    } else {

      k  = 1.0 / z;
      a2 = c2 * k;
      a1 = (c1 * k) * k;
      a0 = copysign (1.0 , c0);

    }

    k = 1.0 / k;

    if (a0 == 0.0 && a1 == 0.0 && a2 == 0.0) {
      cubicType = allzero;
    } else if (a0 == 0.0 && a1 == 0.0) {
      cubicType = linear;
    } else if (a0 == 0.0) {
      cubicType = quadratic;
    } else {
      cubicType = general;
    }

  }
  //
  //
  //     ...Select the case.
  //
  //        1) Only zero roots.
  //
  //
  switch (cubicType)
    {
    case allzero:

      *nReal = 3;

      root[0]=root[1]=root[2] = 0.0+I*0.0;
      //root (:,Im) = 0.0;
      break;
      //
      //
      //     ...2) The linear equation case -> additional 2 zeros.
      //
      //
    case linear:

      x = - a2 * k;

	*nReal = 3;

	root [0] = max2 (0.0, x)+I*0.0;
	root [1] = 0.0+I*0.0;
	root [2] = min2 (0.0, x)+I*0.0;
	//
	//
	//     ...3) The quadratic equation case -> additional 1 zero.
	//
	//
	break;
    case quadratic:

      quadraticRoots (a2, a1, nReal, root);

      if (*nReal == 2) {

	x = creal(root[0]) * k;         // real roots of quadratic are ordered x >= y
	y = creal(root[1]) * k;

	*nReal = 3;

	root[0] = max2 (x, 0.0)+I*0.0;
	root[1] = max2 (y, min2 (x, 0.0))+I*0.0;
	root[2] = min2 (y, 0.0)+I*0.0;
	//root (:,Im) = 0.0;

      } else {

	*nReal = 1;

	root [2] = root [1] * k;
	root [1] = root [0] * k;
	root [0] = 0.0+I*0.0;
	//root (3,Im) = root (2,Im) * k;
	//root (2,Im) = root (1,Im) * k;
	//root (1,Im) = 0.0;

      }
      break;
      //
      //
      //     ...3) The general cubic case. Set the best Newton-Raphson root estimates for the cubic.
      //           The easiest and most robust conditions are checked first. The most complicated
      //           ones are last and only done when fabsolutely necessary.
      //
      //
    case general:

      if (a0 == 1.0) {

	x = - p1 + q1 * a1 - a2 * (r1 - s1 * a1);

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a0 == - 1.0) {

	x = p1 - q1 * a1 - a2 * (r1 - s1 * a1);

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a1 == 1.0) {

	if (a0 > 0.0) {
	  x = a0 * (- q4 - s4 * a2);
	} else {
	  x = a0 * (- q4 + s4 * a2);
	}

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a1 == - 1.0) {

	y = - two27th;
	y = y * a2;
	y = y * a2 - third;
	y = y * a2;

	if (a0 < y) {
	  x = + p3 - q3 * a0 - a2 * (r3 + s3 * a0);               // + guess
	} else {
	  x = - p3 - q3 * a0 - a2 * (r3 - s3 * a0);               // - guess
	}

	a = a2;
	b = a1;
	c = a0;
	xShift = 0.0;

      } else if (a2 == 1.0) {

	b = a1 - third;
	c = a0 - one27th;

	if (fabs (b) < macheps && fabs (c) < macheps) {        // triple -1/3 root

	  x = - third * k;

	  *nReal = 3;

	  root[0]=root[1]=root[2] = x + I*0.0;
	  return;

	} else {

	  y = third * a1 - two27th;

	  if (a1 <= third) {
	    if (a0 > y) {
	      x = - p51 - q51 * a0 + a1 * (r51 - s51 * a0);   // - guess
	    } else {
	      x = + p52 - q52 * a0 - a1 * (r52 + s52 * a0);   // + guess
	    }
	  } else {
	    if (a0 > y) {
	      x = - p53 - q53 * a0 + a1 * (r53 - s53 * a0);   // <-1/3 guess
	    } else {
	      x = + p54 - q54 * a0 - a1 * (r54 + s54 * a0);   // >-1/3 guess
	    }
	  }

	  if (fabs (b) < 1.e-2 && fabs (c) < 1.e-2) {  // use shifted root
	    c = - third * b + c;
	    if (fabs (c) < macheps) c = 0.0;                  // prevent random noise
	    a = 0.0;
	    xShift = third;
	    x = x + xShift;
	  } else {
	    a = a2;
	    b = a1;
	    c = a0;
	    xShift = 0.0;
	  }

	}

      } else if (a2 == - 1.0) {

	b = a1 - third;
	c = a0 + one27th;

	if (fabs (b) < macheps && fabs (c) < macheps) {        // triple 1/3 root

	  x = third * k;

	  *nReal = 3;

	  root [0]=root[1]=root[2] = x + I*0.0;
	  //root (:,Im) = 0.0;

	  return;

	} else {

	  y = two27th - third * a1;

	  if (a1 <= third) {
	    if (a0 < y) {
	      x = + p51 - q51 * a0 - a1 * (r51 + s51 * a0);   // +1 guess
	    } else {
	      x = - p52 - q52 * a0 + a1 * (r52 - s52 * a0);   // -1 guess
	    }
	  } else {
	    if (a0 < y) {
	      x = + p53 - q53 * a0 - a1 * (r53 + s53 * a0);   // >1/3 guess
	    } else {
	      x = - p54 - q54 * a0 + a1 * (r54 - s54 * a0);   // <1/3 guess
	    }
	  }

	  if (fabs (b) < 1.e-2 && fabs (c) < 1.e-2) {  // use shifted root
	    c = third * b + c;
	    if (fabs (c) < macheps) c = 0.0;                  // prevent random noise
	    a = 0.0;
	    xShift = - third;
	    x = x + xShift;
	  } else {
	    a = a2;
	    b = a1;
	    c = a0;
	    xShift = 0.0;
	  }

	}

      }
      //
      //
      //     ...Perform Newton/Bisection iterations on x^3 + ax^2 + bx + c.
      //
      //
      z = x + a;
      y = x + z;
      z = z * x + b;
      y = y * x + z;       // C'(x)
      z = z * x + c;       // C(x)
      t = z;               // save C(x) for sign comparison
      x = x - z / y;       // 1st improved root

      oscillate = 0;
      bisection = 0;
      converged = 0;

      while (!converged && !bisection)    // Newton-Raphson iterates
	{
	  z = x + a;
	  y = x + z;
	  z = z * x + b;
	  y = y * x + z;
	  z = z * x + c;

	  if (z * t < 0.0) {                       // does Newton start oscillating ?
	    if (z < 0.0) {
	      oscillate = oscillate + 1;              // increment oscillation counter
	      s = x;                                  // save lower bisection bound
	    } else {
	      u = x;                                  // save upper bisection bound
	    }
	    t = z;                                      // save current C(x)
	  }

	  y = z / y;                                      // Newton correction
	  x = x - y;                                      // new Newton root

	  bisection = (oscillate > 2)?1:0;                      // activate bisection
	  converged = (fabs (y) <= fabs (x) * macheps)?1:0;       // Newton convergence indicator


	}

      if (bisection) {

	t = u - s;                                     // initial bisection interval
	while (fabs (t) > fabs (x) * macheps)        // bisection iterates
	  {
	    z = x + a;                                  //
	    z = z * x + b;                              // C (x)
	    z = z * x + c;                              //

	    if (z < 0.0) {                       //
	      s = x;                                  //
	    } else {                                       // keep bracket on root
	      u = x;                                  //
	    }                                     //

	    t = 0.5 * (u - s);                       // new bisection interval
	    x = s + t;                                  // new bisection root


	  }
      }


      x = x - xShift;                                   // unshift root
      //
      //
      //     ...Forward / backward deflate rescaled cubic (if needed) to check for other real roots.
      //        The deflation analysis is performed on the rescaled cubic. The actual deflation must
      //        be performed on the original cubic, not the rescaled one. Otherwise deflation errors
      //        will be enhanced when undoing the rescaling on the extra roots.
      //
      //
      z = fabs (x);
      s = fabs (a2);
      t = fabs (a1);
      u = fabs (a0);

      y = z * max2 (s,z);           // take maximum between |x^2|,|a2 * x|

      deflateCase = 1;             // up to now, the maximum is |x^3| or |a2 * x^2|

      if (y < t) {             // check maximum between |x^2|,|a2 * x|,|a1|
	y = t * z;               // the maximum is |a1 * x|
	deflateCase = 2;         // up to now, the maximum is |a1 * x|
      } else {
	y = y * z;               // the maximum is |x^3| or |a2 * x^2|
      }

      if (y < u) {             // check maximum between |x^3|,|a2 * x^2|,|a1 * x|,|a0|
	deflateCase = 3;         // the maximum is |a0|
      }

      y = x * k;                   // real root of original cubic

      switch (deflateCase)
	{
	case 1:
	  x = 1.0 / y;
	  t = - c0 * x;              // t -> backward deflation on unscaled cubic
	  s = (t - c1) * x;          // s -> backward deflation on unscaled cubic
	  break;	
	case 2:
	  s = c2 + y;                // s ->  forward deflation on unscaled cubic
	  t = - c0 / y;              // t -> backward deflation on unscaled cubic
	  break;
	case 3:
	  s = c2 + y;                // s ->  forward deflation on unscaled cubic
	  t = c1 + s * y;            // t ->  forward deflation on unscaled cubic
	}

      quadraticRoots (s, t, nReal, root);

      if (*nReal == 2) {

	x = creal(root[0]);         // real roots of quadratic are ordered x >= z
	z = creal(root[1]);         // use 'z', because 'y' is original cubic real root


	*nReal = 3;

	root [0] = max2 (x, y)+I*0.0;
	root [1] = max2 (z, min2 (x, y))+I*0.0;
	root [2] = min2 (z, y)+I*0.0;
	//root (:,Im) = 0.0;
      } else {

	*nReal = 1;

	root [2] = root[1];
	root [1] = root[0];
	root [0] = y + I*0.0;
	//root (3,Im) = root (2,Im);
	//root (2,Im) = root (1,Im)
	  //root (1,Im) = 0.0;

	//printf("QUI root[0]=%.15G %.15G root[1]= %.15G %.15G\n", y, creal(root[0]), cimag(root[0]),creal(root[1]), cimag(root[1]));
      }
      break;
    }
  // end select
  //
  //
  //     ...Ready!
  //
  //
  return;
}
//end subroutine cubicRoots



//-----------------------------------------------------------------------------------
//
// QUADRATIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call quadraticRoots (real,    intent (in)  :: q1,
//                       real,    intent (in)  :: q0,
//                       integer, intent (out) :: nReal,
//                       real,    intent (out) :: root (1:2,1:2))
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the quadratic polynomial:
//
//                 x^2 + q1 * x + q0
//
//  The code checks internally, if rescaling of the coefficients is needed to
//  avoid overflow.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) Since there can be only one complex conjugate pair root, no order
//           is necessary.
//
// ARGUMENTS
//
//  q1         : coefficient of x term
//  q0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//
// NOTES
//
//***

void quadraticRoots (double q1, double q0, int *nReal, complex double root[2])
{

  int rescale;

  double a0, a1;
  double k, x, y, z;

  const double LPN = 1.7976931348623157E+308;   // the (L)argest (P)ositive (N)umber
  const double sqrtLPN = sqrt (LPN);      // and the square root of it
//
//
//     ...Handle special cases.
//
//
  if (q0 == 0.0 && q1 == 0.0) {

      *nReal = 2;

      root [0] = root[1] = 0.0+I*0.0;
      //root (:,Im) = 0.0;

  } else if (q0 == 0.0) {

      *nReal = 2;

      root [0] = max2 (0.0, - q1)+I*0.0;
      root [1] = min2 (0.0, - q1)+I*0.0;

  } else if (q1 == 0.0) {

      x = sqrt (fabs (q0));

      if (q0 < 0.0) {

          *nReal = 2;

          root [0] = x+I*0.0;
          root [1] = - x+I*0.0;
	  // root (:,Im) = 0.0;

      } else {

          *nReal = 0;

          root [0] = 0.0+I*x;
          root [1] = 0.0 - x*I;

      }

  } else {
//
//
//     ...The general case. Do rescaling, if either squaring of q1/2 or evaluation of
//        (q1/2)^2 - q0 will lead to overflow. This is better than to have the solver
//        crashed. Note, that rescaling might lead to loss of accuracy, so we only
//        invoke it when fabsolutely necessary.
//
//
      rescale = (q1 > (sqrtLPN + sqrtLPN))?1:0;     // this detects overflow of (q1/2)^2

      if (!rescale) {
           x = q1 * 0.5;                      // we are sure here that x*x will not overflow
           rescale = (q0 < (x * x - LPN));      // this detects overflow of (q1/2)^2 - q0
      }

      if (rescale) {

          x = fabs (q1);
          y = sqrt (fabs (q0));

          if (x > y) {
              k  = x;
              z  = 1.0 / x;
              a1 = copysign (1.0 , q1);
              a0 = (q0 * z) * z;
          } else {
              k  = y;
              a1 = q1 / y;
              a0 = copysign (1.0 , q0);
          }

      } else {
          a1 = q1;
          a0 = q0;
      }
//
//
//     ...Determine the roots of the quadratic. Note, that either a1 or a0 might
//        have become equal to zero due to underflow. But both cannot be zero.
//
//
      x = a1 * 0.5;
      y = x * x - a0;

      if (y >= 0.0) {

          y = sqrt (y);

          if (x > 0.0) {
              y = - x - y;
          } else {
              y = - x + y;
          }

          if (rescale) {
              y = y * k;                     // very important to convert to original
              z = q0 / y;                    // root first, otherwise complete loss of
          } else {                              // root due to possible a0 = 0 underflow
              z = a0 / y;
          }

          *nReal = 2;

          root [0] = max2 (y,z)+I*0.0;           // 1st real root of x^2 + a1 * x + a0
          root [1] = min2 (y,z)+I*0.0;           // 2nd real root of x^2 + a1 * x + a0

      } else {

          y = sqrt (- y);

          *nReal = 0;

          root [0] = - x+I*y;
          root [1] = - x - y*I;

          if (rescale) {
              root[0] = root[0] * k;
	      root[1] = root[1] * k;
          }

      }

  }
//
//
//     ...Ready!
//
//
  return;
}

//-----------------------------------------------------------------------------------
//
// QUARTIC POLYNOMIAL ROOT SOLVER
//
// SYNOPSIS
//
//  call quarticRoots (real,              intent (in)  :: q3,
//                     real,              intent (in)  :: q2,
//                     real,              intent (in)  :: q1,
//                     real,              intent (in)  :: q0,
//                     integer,           intent (out) :: nReal,
//                     real,              intent (out) :: root (1:4,1:2),
//                     logical, optional, intent (in)  :: printInfo)
//
// DESCRIPTION
//
//  Calculates all real + complex roots of the quartic polynomial:
//
//                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
//
//  An option for printing detailed info about the intermediate stages in solving
//  the quartic is available. This enables a detailed check in case something went
//  wrong and the roots obtained are not proper.
//
//  The quartic root solver can handle any size of quartic coefficients and there is
//  no danger of overflow, due to proper rescaling of the quartic polynomial.
//
//  The order of the roots is as follows:
//
//        1) For real roots, the order is according to their algebraic value
//           on the number scale (largest positive first, largest negative last).
//
//        2) For complex conjugate pair roots, the order is according to the
//           algebraic value of their real parts (largest positive first). If
//           the real parts are equal, the order is according to the algebraic
//           value of their imaginary parts (largest first).
//
//        3) All real roots preceede the complex ones.
//
// ARGUMENTS
//
//  q3         : coefficient of x^3 term
//  q2         : coefficient of x^2 term
//  q1         : coefficient of x term
//  q0         : independent coefficient
//  nReal      : number of real roots found
//  root (n,1) : real part of n-th root
//  root (n,2) : imaginary part of n-th root
//  printInfo  : if given and true, detailed info will be printed about intermediate stages
//
// NOTES
//
//***

void CquarticRoots (double cc[5], int *nReal, complex double root[4])
{
  int bisection;
  int converged;
  int iterate;
  int minimum;
  int overshoot;
  int dsignflip;
  int notZero;

  int deflateCase;
  int quarticType;

  double q1, q2, q3, q0;
  //const int biquadratic = 2, cubic = 3, general = 4; 
  enum casi {biquadratic=2, cubic=3, general=4 };
  double a0, a1, a2, a3;
  double a, b, c, d, k, s, t, u, x, y, z;

  const double macheps = 2.2204460492503131E-16; 
  const double third   = 1.0 / 3.0;
  //
  //
  //     ...Start.
  //
  //
#if 0
  if (present (printInfo)) {
    doPrint = printInfo
  } else {
    doPrint = 0
  }
#endif

#if 0
  if (doPrint) {
    write (*,wpformat) ' initial quartic q3    = ',q3
      write (*,wpformat) ' initial quartic q2    = ',q2
      write (*,wpformat) ' initial quartic q1    = ',q1
      write (*,wpformat) ' initial quartic q0    = ',q0
      write (*,wpformat) ' ------------------------------------------------'
  }
#endif
  //
  //
  //     ...Handle special cases. Since the cubic solver handles all its
  //        special cases by itself, we need to check only for two cases:
  //
  //            1) independent term is zero -> solve cubic and include
  //               the zero root
  //
  //            2) the biquadratic case.
  //
  //
  q3=cc[3]/cc[4];
  q2=cc[2]/cc[4];
  q1=cc[1]/cc[4];
  q0=cc[0]/cc[4];
  if (q0 == 0.0) {

    k  = 1.0;
    a3 = q3;
    a2 = q2;
    a1 = q1;

    quarticType = cubic;

  } else if (q3 == 0.0 && q1 == 0.0) {

    k  = 1.0;
    a2 = q2;
    a0 = q0;

    quarticType = biquadratic;

  } else {
    //
    //
    //     ...The general case. Rescale quartic polynomial, such that largest fabsolute coefficient
    //        is (exactly!) equal to 1. Honor the presence of a special quartic case that might have
    //        been obtained during the rescaling process (due to underflow in the coefficients).
    //
    //
    s = fabs (q3);
    t = sqrt (fabs (q2));
    u = pow(fabs (q1), third);
    x = sqrt (sqrt (fabs (q0)));
    y = max4 (s,t,u,x);

    if (y == s) {

      k  = 1.0 / s;
      a3 = copysign (1.0 , q3);
      a2 = (q2 * k) * k;
      a1 = ((q1 * k) * k) * k;
      a0 = (((q0 * k) * k) * k) * k;

    } else if (y == t) {

      k  = 1.0 / t;
      a3 = q3 * k;
      a2 = copysign (1.0 , q2);
      a1 = ((q1 * k) * k) * k;
      a0 = (((q0 * k) * k) * k) * k;

    } else if (y == u) {

      k  = 1.0 / u;
      a3 = q3 * k;
      a2 = (q2 * k) * k;
      a1 = copysign (1.0 , q1);
      a0 = (((q0 * k) * k) * k) * k;

    } else {

      k  = 1.0 / x;
      a3 = q3 * k;
      a2 = (q2 * k) * k;
      a1 = ((q1 * k) * k) * k;
      a0 = copysign (1.0 , q0);

    }

    k = 1.0 / k;

#if 0
    if (doPrint) {
      write (*,wpformat) ' rescaling factor      = ',k
	write (*,wpformat) ' ------------------------------------------------'
	write (*,wpformat) ' rescaled quartic q3   = ',a3
	write (*,wpformat) ' rescaled quartic q2   = ',a2
	write (*,wpformat) ' rescaled quartic q1   = ',a1
	write (*,wpformat) ' rescaled quartic q0   = ',a0
	write (*,wpformat) ' ------------------------------------------------'
    }
#endif
    if (a0 == 0.0) {
      quarticType = cubic;
    } else if (a3 == 0.0 && a1 == 0.0) {
      quarticType = biquadratic;
    } else {
      quarticType = general;
    }

  }
  //
  //
  //     ...Select the case.
  //
  //        1) The quartic with independent term = 0 -> solve cubic and add a zero root.
  //
  //
  switch (quarticType)
    {
    case cubic:

      cubicRoots (a3, a2, a1, nReal, root);

      if (*nReal == 3) {

	x = creal(root[0]) * k;       // real roots of cubic are ordered x >= y >= z
	y = creal(root[1])* k;
	z = creal(root[2]) * k;

	*nReal = 4;

	root [0] = max2 (x, 0.0) + I*0.0;
	root [1] = max2 (y, min2 (x, 0.0))+I*0.0;
	root [2] = max2 (z, min2 (y, 0.0))+I*0.0;
	root [3] = min2 (z, 0.0)+I*0.0;
	//root (:,Im) = 0.0;

      } else {                          // there is only one real cubic root here

	x = creal(root[0]) * k;

	*nReal = 2;

	root [3] = root [2] * k;
	root [2] = root [1] * k;
	root [1] = min2 (x, 0.0)+I*0.0;
	root [0] = max2 (x, 0.0)+I*0.0;

      }
      break;
      //
      //
      //     ...2) The quartic with x^3 and x terms = 0 -> solve biquadratic.
      //
      //
      //
    case biquadratic:

      quadraticRoots (q2, q0, nReal, root);

	if (*nReal == 2) {

	  x = creal(root [0]);         // real roots of quadratic are ordered x >= y
	  y = creal(root [1]);

	    if (y >= 0.0) {

	      x = sqrt (x) * k;
      	      y = sqrt (y) * k;

      	      *nReal = 4;

      	      root [0] = x+I*0.0;
	      root [1] = y+I*0.0;
	      root [2] = - y+I*0.0;
	      root [3] = - x+I*0.0;

	    } else if (x >= 0.0 && y < 0.0) {

	      x = sqrt (x)       * k;
      	      y = sqrt (fabs (y)) * k;

      	      *nReal = 2;

      	      root [0] = x+I*0.0;
	      root [1] = - x+I*0.0;
	      root [2] = 0.0+I*y;
	      root [3] = 0.0-I*y;

	    } else if (x < 0.0) {

	      x = sqrt (fabs (x)) * k;
      	      y = sqrt (fabs (y)) * k;

      	      *nReal = 0;

      	      root [0] = 0+I*y;
	      root [1] = 0+I*x;
	      root [2] = 0- I*x;
	      root [3] = 0- I*y;

	    }

	} else {          // complex conjugate pair biquadratic roots x +/- iy.

	  x = creal(root [0]) * 0.5;
	  y = cimag(root[0]) * 0.5;
	  z = sqrt (x * x + y * y);
	  y = sqrt (z - x) * k;
	  x = sqrt (z + x) * k;

	  *nReal = 0;

	  root [0] = x+I*y;
	  root [1] = x-I*y;
	  root [2] = - x+I*y;
	  root [3] = - x-I*y;

	}
      break;
      //
      //
      //     ...3) The general quartic case. Search for stationary points. Set the first
      //           derivative polynomial (cubic) equal to zero and find its roots.
      //           Check, if any minimum point of Q(x) is below zero, in which case we
      //           must have real roots for Q(x). Hunt down only the real root, which
      //           will potentially converge fastest during Newton iterates. The remaining
      //           roots will be determined by deflation Q(x) -> cubic.
      //
      //           The best roots for the Newton iterations are the two on the opposite
      //           ends, i.e. those closest to the +2 and -2. Which of these two roots
      //           to take, depends on the location of the Q(x) minima x = s and x = u,
      //           with s > u. There are three cases:
      //
      //              1) both Q(s) and Q(u) < 0
      //                 ----------------------
      //
      //                 The best root is the one that corresponds to the lowest of
      //                 these minima. If Q(s) is lowest -> start Newton from +2
      //                 downwards (or zero, if s < 0 and a0 > 0). If Q(u) is lowest
      //                 -> start Newton from -2 upwards (or zero, if u > 0 and a0 > 0).
      //
      //              2) only Q(s) < 0
      //                 -------------
      //
      //                 With both sides +2 and -2 possible as a Newton starting point,
      //                 we have to avoid the area in the Q(x) graph, where inflection
      //                 points are present. Solving Q''(x) = 0, leads to solutions
      //                 x = -a3/4 +/- discriminant, i.e. they are centered around -a3/4.
      //                 Since both inflection points must be either on the r.h.s or l.h.s.
      //                 from x = s, a simple test where s is in relation to -a3/4 allows
      //                 us to avoid the inflection point area.
      //
      //              3) only Q(u) < 0
      //                 -------------
      //
      //                 Same of what has been said under 2) but with x = u.
      //
      //
    case general:

      x = 0.75 * a3;
      y = 0.50 * a2;
      z = 0.25 * a1;

      cubicRoots (x, y, z, nReal, root);

      s = creal(root[0]);        // Q'(x) root s (real for sure)
      x = s + a3;
      x = x * s + a2;
      x = x * s + a1;
      x = x * s + a0;         // Q(s)

      y = 1.0;             // dual info: Q'(x) has more real roots, and if so, is Q(u) < 0 ? 

      if (*nReal > 1) {
	u = creal(root [2]);   // Q'(x) root u
	y = u + a3;
	y = y * u + a2;
	y = y * u + a1;
	y = y * u + a0;     // Q(u)
      }

#if 0
      if (doPrint) {
	write (*,wpformat) ' dQ(x)/dx root s       = ',s
	  write (*,wpformat) ' Q(s)                  = ',x
	  write (*,wpformat) ' dQ(x)/dx root u       = ',u
	  write (*,wpformat) ' Q(u)                  = ',y
	  write (*,wpformat) ' ------------------------------------------------'
      }
#endif
      if (x == 0.0 && y == 0.0) 
	{
	  if (fabs(s) > fabs(u))
	    x=s;
	  else
	    x=u;
	  *nReal = 1;
	  iterate = 0;
	} 
      else if (x==0)
	{
	  x=s;
	  *nReal = 1;
	  iterate = 0;
	}
      else if (y==0)
	{
	  x = u;
	  *nReal = 1;
	  iterate = 0;
	}
      else if (x < 0.0 && y < 0.0)
	{
	  if (s < 0.0 && a0 > 0.0) 
	    x = 0.0;
          else
	    x = 2.0;
          
          if (u > 0.0 && a0 > 0.0) 
              y = 0.0;
          else
              y = -2.0;

          a = x + a3;
          b = x + a;
          a = a * x + a2;
          b = b * x + a;
          a = a * x + a1;
          b = b * x + a;     // b = Q'(x)

          c = y + a3;
          d = y + c;
          c = c * y + a2;
          d = d * y + c;
          c = c * y + a1;
          d = d * y + c;     //! d = Q'(y)

	  if (fabs (b) > fabs (d))   
	    {                           // if Q'(y) < Q'(x),
              x = y;                     // take root u for Newton iterations
              s = u;                     // save for lower bisecion bound just in case
	    }

          *nReal = 1;
          iterate = 1;
	}
      else if (x < 0.0) {

	if (s < - a3 * 0.25) {
	  if (s > 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = - 2.0;
	  }
	} else {
	  if (s < 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = 2.0;
	  }
	}

	*nReal = 1;
	iterate = 1;

      } else if (y < 0.0) {

	if (u < - a3 * 0.25) {
	  if (u > 0.0 && a0 > 0.0) {
	    x =  0.0;
	  } else {
	    x = -2.0;
	  }
	} else {
	  if (u < 0.0 && a0 > 0.0) {
	    x = 0.0;
	  } else {
	    x = 2.0;
	  }
	}
	s = u;
	*nReal = 1;
	iterate = 1;
      } else {
	*nReal = 0;
      }
      //
      //
      //     ...Do all necessary Newton iterations. In case we have more than 2 oscillations,
      //        exit the Newton iterations and switch to bisection. Note, that from the
      //        definition of the Newton starting point, we always have Q(x) > 0 and Q'(x)
      //        starts (-ve/+ve) for the (-2/+2) starting points and (increase/decrease) smoothly
      //        and staying (< 0 / > 0). In practice, for extremely shallow Q(x) curves near the
      //        root, the Newton procedure can overshoot slightly due to rounding errors when
      //        approaching the root. The result are tiny oscillations around the root. If such
      //        a situation happens, the Newton iterations are abandoned after 3 oscillations
      //        and further location of the root is done using bisection starting with the
      //        oscillation brackets.
      //
      //
      if (*nReal > 0) {
	if (iterate)
	  {
	    y = x + a3;                                     //
	    z = x + y;                                     //
	    y = y * x + a2;                                 // y = Q(x)
	    z = z * x + y;                                  //
	    y = y * x + a1;                                 // z = Q'(x)
	    z = z * x + y;                                  //
	    y = y * x + a0;                                 //
	    t = z;
	    y = y / z;
	    x = x - y;

	    dsignflip = 0;
	    overshoot = 0;
	    bisection = 0;
	    converged = 0;

	    while (!converged && !bisection)    // Newton-Raphson iterates
	      {
		y = x + a3;                                      //
		z = x + y;                                      //
		y = y * x + a2;                                 // y = Q(x)
		z = z * x + y;                                  //
		y = y * x + a1;                                 // z = Q'(x)
		z = z * x + y;                                  //
		y = y * x + a0;                                 //

		if (y < 0.0)    
		  {                                            // does Newton start overshooting ?
		    overshoot = overshoot + 1;                 // increment overshoot counter
		    s = x;
		  }                                            // save lower bisection bound
		else
		  u = x;                                      // save upper bisection bound

		if (z * t < 0.0)  
		  {  // does Q'(x) have a sign flip ?
		    dsignflip = dsignflip + 1;                  // increment sign flip counter
		    t = z;
		  }  // save Q'(x) for next sign check

		if (z == 0.0)                           // safeguard against accidental
		  {
		    bisection = 1;                        // Q'(x) = 0 due to roundoff
		    break;
		  } // errors -> activate bisection
		// with current bracket [s,u]

		y = y / z;                                      // Newton correction
		x = x - y;                                      // new Newton root

		bisection = (overshoot > 2 || dsignflip > 2)?1:0;   // activate bisection
		converged = (fabs (y) <= fabs (x) * macheps)?1:0;   // Newton convergence indicator
	      }

	    if (bisection) 
	      {
		t = u - s;                                     // initial bisection interval
		while (fabs (t) > fabs (x) * macheps)        // bisection iterates
		  {
		    y = x + a3;                                 //
		    y = y * x + a2;                             // y = Q(x)
		    y = y * x + a1;                             //
		    y = y * x + a0;                             //

		    if (y < 0.0)                       //
		      s = x;                                  //
		    else                                       // keep bracket on root
		      u = x;                                  //

		    t = 0.5 * (u - s);                       // new bisection interval
		    x = s + t;                                  // new bisection root
		  }
	      }
	  }

	//
	//
	//     ...Find remaining roots -> reduce to cubic. The reduction to a cubic polynomial
	//        is done using composite deflation to minimize rounding errors. Also, while
	//        the composite deflation analysis is done on the reduced quartic, the actual
	//        deflation is being performed on the original quartic again to avoid enhanced
	//        propagation of root errors.
	//
	//
	z = fabs (x);            //
	a = fabs (a3);           //
	b = fabs (a2);           // prepare for composite deflation
	c = fabs (a1);           //
	d = fabs (a0);           //

	y = z * max2 (a,z);      // take maximum between |x^2|,|a3 * x|

	deflateCase = 1;        // up to now, the maximum is |x^4| or |a3 * x^3|

	if (y < b) {        // check maximum between |x^2|,|a3 * x|,|a2|
	  y = b * z;          // the maximum is |a2| -> form |a2 * x|
	  deflateCase = 2;    // up to now, the maximum is |a2 * x^2|
	} else {
	  y = y * z;          // the maximum is |x^3| or |a3 * x^2|
	}

	if (y < c) {        // check maximum between |x^3|,|a3 * x^2|,|a2 * x|,|a1|
	  y = c * z;          // the maximum is |a1| -> form |a1 * x|
	  deflateCase = 3;    // up to now, the maximum is |a1 * x|
	} else {
	  y = y * z;          // the maximum is |x^4|,|a3 * x^3| or |a2 * x^2|
	}

	if (y < d) {        // check maximum between |x^4|,|a3 * x^3|,|a2 * x^2|,|a1 * x|,|a0|
	  deflateCase = 4;    // the maximum is |a0|
	}

	x = x * k;              // 1st real root of original Q(x)

	switch (deflateCase)
	  {
	  case 1:
	    z = 1.0 / x;
	    u = - q0 * z;         // u -> backward deflation on original Q(x)
	    t = (u - q1) * z;     // t -> backward deflation on original Q(x)
	    s = (t - q2) * z;     // s -> backward deflation on original Q(x)
	    break;
	  case 2:
	    z = 1.0 / x;
	    u = - q0 * z;         // u -> backward deflation on original Q(x)
	    t = (u - q1) * z;     // t -> backward deflation on original Q(x)
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    break;
	  case 3:
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    t = q2 + s * x;       // t ->  forward deflation on original Q(x)
	    u = - q0 / x;         // u -> backward deflation on original Q(x)
	    break;
	  case 4:
	    s = q3 + x;           // s ->  forward deflation on original Q(x)
	    t = q2 + s * x;       // t ->  forward deflation on original Q(x)
	    u = q1 + t * x;       // u ->  forward deflation on original Q(x)
	  }

	cubicRoots (s, t, u, nReal, root);

	  if (*nReal == 3) {

	    s = creal(root [0]);    //
	    t = creal(root [1]);    // real roots of cubic are ordered s >= t >= u
	    u = creal(root [2]);    //

	    root [0] = max2 (s, x)+I*0.0;
	    root [1] = max2 (t, min2 (s, x))+I*0.0;
	    root [2] = max2 (u, min2 (t, x))+I*0.0;
	    root [3] = min2 (u, x)+I*0.0;

	    *nReal = 4;

	  } else {                   // there is only one real cubic root here

	    s = creal(root [0]);

	    root [3] = root [2];
	    root [2] = root [1];
	    root [1] = min2 (s, x)+I*0.0;
	    root [0] = max2 (s, x)+I*0.0;

	    *nReal = 2;

	  }

      } else {
	//
	//
	//     ...If no real roots have been found by now, only complex roots are possible.
	//        Find real parts of roots first, followed by imaginary components.
	//
	//
	s = a3 * 0.5;
	t =  s * s - a2;
	u =  s * t + a1;                   // value of Q'(-a3/4) at stationary point -a3/4


	notZero = (fabs (u) >= macheps)?1:0;    // H(-a3/4) is considered > 0 at stationary point

	if (a3 != 0.0) {
	  s = a1 / a3;
	  minimum = (a0 > (s * s))?1:0;                            // H''(-a3/4) > 0 -> minimum
	} else {
	  minimum = ((4 * a0) > (a2 * a2))?1:0;                      // H''(-a3/4) > 0 -> minimum
	}

	iterate = (notZero || (!notZero && minimum))?1:0;

	if (iterate) {

	  x = copysign (2.0,a3);                              // initial root -> target = smaller mag root

	  overshoot = 0;
	  bisection = 0;
	  converged = 0;

	  while (!converged && !bisection)    // Newton-Raphson iterates
	    {
	      a = x + a3;                                     //
	      b = x + a;                                      // a = Q(x)
	      c = x + b;                                      //
	      d = x + c;                                      // b = Q'(x)
	      a = a * x + a2;                                 //
	      b = b * x + a;                                  // c = Q''(x) / 2
	      c = c * x + b;                                  //
	      a = a * x + a1;                                 // d = Q'''(x) / 6
	      b = b * x + a;                                  //
	      a = a * x + a0;                                 //
	      y = a * d * d - b * c * d + b * b;              // y = H(x), usually < 0
	      z = 2 * d * (4 * a - b * d - c * c);            // z = H'(x)
	     
	      if (y > 0.0)                            // does Newton start oscillating ?
		{
		  overshoot = overshoot + 1;                  // increment oscillation counter
		  s = x;
		}
	      else
		{
		  u = x;	  
		}
	      if (z == 0.0)
		{
		  converged = 1;
		  bisection = 0;
		  break;
		}
	      y = y / z;                                      // Newton correction
	      x = x - y;                                      // new Newton root

	      bisection = (overshoot > 2)?1:0;                      // activate bisection
	      converged = (fabs (y) <= fabs (x) * macheps)?1:0;       // Newton convergence criterion

	      // if (doPrint) write (*,wpformat) ' Newton H(x) root      = ',x

	    }

	  if (bisection) {

	    t = u - s;                                     // initial bisection interval
	    while (fabs (t) > fabs (x) * macheps)        // bisection iterates
	      {
		a = x + a3;                                 //
		b = x + a;                                  // a = Q(x)
		c = x + b;                                  //
		d = x + c;                                  // b = Q'(x)
		a = a * x + a2;                             //
		b = b * x + a;                              // c = Q''(x) / 2
		c = c * x + b;                              //
		a = a * x + a1;                             // d = Q'''(x) / 6
		b = b * x + a;                              //
		a = a * x + a0;                             //
		y = a * d * d - b * c * d + b * b;          // y = H(x)

		if (y > 0.0) {                       //
		  s = x;                                  //
		} else {                                       // keep bracket on root
		  u = x;                                  //
		}                                     //

		t = 0.5 * (u - s);                       // new bisection interval
		x = s + t;                                  // new bisection root

	      }

	  }


	  a = x * k;                                         // 1st real component -> a
	  b = - 0.5 * q3 - a;                             // 2nd real component -> b
	  c = a * a;
	  d = b * b;

	  x = 4 * a + q3;                                    // Q'''(a)
	  y = x + q3 + q3;                                   //
	  y = y * a + q2 + q2;                               // Q'(a)
	  y = y * a + q1;                                    //
	  y = y / x;                                         // Q'(a) / Q'''(a)
	  s = c + y;                                         // magnitude^2 of (a + iy) root
	  x = 4 * b + q3;                                    // Q'''(b)
	  z = x + q3 + q3;                                   //
	  z = z * b + q2 + q2;                               // Q'(b)
	  z = z * b + q1;                                    //
	  z = z / x;                                         // Q'(b) / Q'''(b)
	  t = d + z;                                         // magnitude^2 of (b + iz) root

	  if (s > t) 
	    {                                   // minimize imaginary error
	      y = max2(y, 0.0);                           // ensure >= 0 for sqrt
      	      d = max2(q0 / s - d, 0.0);                  // small component using Vieta
	      c = sqrt(y);                                  // 1st imaginary component -> c
	      d = sqrt(d);
	    }// 2nd imaginary component -> d
	  else
	    {
	      c = max2(q0 / t - c, 0.0);                  // small component using Vieta
   	      z = max2(z, 0.0);                           // ensure >= 0 for sqrt
	      c = sqrt(c);                                  // 1st imaginary component -> c
	      d = sqrt(z);                                  // 2nd imaginary component -> d
	    }

	} else {                                                  // no bisection -> real components equal

	  a = - 0.25 * q3;                                // 1st real component -> a
	  b = a;                                             // 2nd real component -> b = a

	  x = a + q3;                                        //
	  x = x * a + q2;                                    // Q(a)
	  x = x * a + q1;                                    //
	  x = x * a + q0;                                    //
	  y = - 0.1875 * q3 * q3 + 0.5 * q2;           // Q''(a) / 2
	  z = max2 (y * y - x, 0.0);                       // force discriminant to be >= 0
	  z = sqrt (z);                                      // square root of discriminant
	  y = y + copysign (z,y);                                // larger magnitude root
	  if (y == 0.0)                                    // guard against larger magnitude root = 0
	    x = 0.0;                                       // in this case smaller magnitude root must be 0
	  else
	    x = x / y;                                     // smaller magnitude root
	  c = max2 (y, 0.0);                               // ensure root of biquadratic > 0
	  d = max2 (x, 0.0);                               // ensure root of biquadratic > 0
	  c = sqrt (c);                                      // large magnitude imaginary component
	  d = sqrt (d);                                      // small magnitude imaginary component

	}

	if (a > b) {

	  root [0] = a+I*c;
	  root [1] = a-I*c;
	  root [2] = b+I*d;
	  root [3] = b-I*d;

	} else if (a < b) {

	  root [0] = b+I*d;
	  root [1] = b-I*d;
	  root [2] = a+I*c;
	  root [3] = a-I*c;

	} else {
	  root [0] = a+I*c;
	  root [1] = a-I*c;
	  root [2] = a+I*d;
	  root [3] = a-I*d;

	}

      }    // # of real roots 'if'
    }
  //
  //
  //     ...Ready!
  //
  //
  return;
}
//end subroutine quarticRoots
//


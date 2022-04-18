#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/tests.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "tests.h"

/* These are the test functions from table 4.1 of the QUADPACK book */

/* f1(x) = x^alpha * log(1/x) */
/* integ(f1,x,0,1) = 1/(alpha + 1)^2 */

MpIeee f1(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(x,alpha) * log(MpIeee( "1" )/x) ;
}

/* f2(x) = 4^-alpha / ((x-pi/4)^2 + 16^-alpha) */
/* integ(f2,x,0,1) = arctan((4-pi)4^(alpha-1)) + arctan(pi 4^(alpha-1)) */

MpIeee f2(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(MpIeee( "4.0" ),-alpha) / (pow((x-M_PI/MpIeee( "4.0" )),MpIeee( "2.0" )) + pow(MpIeee( "16.0" ),-alpha)) ;
}

/* f3(x) = cos(2^alpha * sin(x)) */
/* integ(f3,x,0,pi) = pi J_0(2^alpha) */

MpIeee f3(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return cos(pow(MpIeee( "2.0" ),alpha) * sin(x)) ;
}

/* Functions 4, 5 and 6 are duplicates of functions  1, 2 and 3 */
/* ....                                                         */

/* f7(x) = |x - 1/3|^alpha */
/* integ(f7,x,0,1) = ((2/3)^(alpha+1) + (1/3)^(alpha+1))/(alpha + 1) */

MpIeee f7(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(fabs(x - (MpIeee( "1.0" )/MpIeee( "3.0" ))),alpha) ;
}

/* f8(x) = |x - pi/4|^alpha */
/* integ(f8,x,0,1) = 
   ((1 - pi/4)^(alpha+1) + (pi/4)^(alpha+1))/(alpha + 1) */

MpIeee f8(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(fabs(x - (M_PI/MpIeee( "4.0" ))),alpha) ;
}

/* f9(x) = sqrt(1 - x^2) / (x + 1 + 2^-alpha) */
/* integ(f9,x,-1,1) = pi/sqrt((1+2^-alpha)^2-1) */

MpIeee f9(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return MpIeee( "1" ) / ((x + MpIeee( "1" ) + pow(MpIeee( "2.0" ),-alpha)) * sqrt(MpIeee( "1" )-x*x)) ;
}

/* f10(x) = sin(x)^(alpha - 1) */
/* integ(f10,x,0,pi/2) = 2^(alpha-2) ((Gamma(alpha/2))^2)/Gamma(alpha) */

MpIeee f10(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(sin(x), alpha-MpIeee( "1" )) ;
}

/* f11(x) = log(1/x)^(alpha - 1) */
/* integ(f11,x,0,1) = Gamma(alpha) */

MpIeee f11(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(log(MpIeee( "1" )/x), alpha-MpIeee( "1" )) ;
}

/* f12(x) = exp(20*(x-1)) * sin(2^alpha * x) */
/* integ(f12,x,0,1) = 
   (20 sin(2^alpha) - 2^alpha cos(2^alpha) + 2^alpha exp(-20))
   /(400 + 4^alpha) */

MpIeee f12(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return exp(MpIeee( "20" )*(x-MpIeee( "1" ))) * sin(pow(MpIeee( "2.0" ),alpha) * x) ;
}

/* f13(x) = cos(2^alpha * x)/sqrt(x(1 - x)) */
/* integ(f13,x,0,1) = pi cos(2^(alpha-1)) J_0(2^(alpha-1))  */

MpIeee f13(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return cos(pow(MpIeee( "2.0" ),alpha)*x)/sqrt(x*(MpIeee( "1" )-x)) ;
}

MpIeee f14(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return exp(-pow(MpIeee( "2.0" ),-alpha)*x)*cos(x)/sqrt(x) ;
}

MpIeee f15(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return x*x * exp(-pow(MpIeee( "2.0" ),-alpha)*x) ;
}

MpIeee f16(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  if (x==MpIeee( "0" ) && alpha == MpIeee( "1" )) return MpIeee( "1" ) ;  /* make the function continuous in x */
  if (x==MpIeee( "0" ) && alpha > MpIeee( "1" )) return MpIeee( "0" ) ;   /* avoid problems with pow(0,1) */
  return pow(x,alpha-MpIeee( "1" ))/pow((MpIeee( "1" )+MpIeee( "10" )*x),MpIeee( "2.0" )) ;
}

MpIeee f17(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return pow(MpIeee( "2.0" ),-alpha)/(((x-MpIeee( "1" ))*(x-MpIeee( "1" ))+pow(MpIeee( "4.0" ),-alpha))*(x-MpIeee( "2" ))) ;
}

/* f454(x) = x^3 log|(x^2-1)(x^2-2)| */
/* integ(f454,x,0,inf) = 61 log(2) + (77/4) log(7) - 27 */

MpIeee f454(MpIeee x, void * params) {
  MpIeee x2=  x * x;
  MpIeee x3=  x * x2;
  params = 0 ;
  return x3 * log(fabs((x2 - MpIeee( "1.0" )) * (x2 - MpIeee( "2.0" )))) ;
}

/* f455(x) = log(x)/(1+100*x^2) */
/* integ(f455,x,0,inf) = -log(10)/20 */

MpIeee f455(MpIeee x, void * params) {
  params = 0 ;
  return log(x) / (MpIeee( "1.0" ) + MpIeee( "100.0" ) * x * x) ;
}

/* f456(x) = log(x) */
/* integ(f456*sin(10 pi x),x,0,1) = -(gamma + log(10pi) - Ci(10pi))/(10pi) */

MpIeee f456(MpIeee x, void * params) {
  params = 0 ;
  if (x == MpIeee( "0.0" ))
    {
      return MpIeee( "0" );
    }
  return log(x) ;
}

/* f457(x) = 1/sqrt(x) */
/* integ(f457*cos(pi x / 2),x,0,+inf) = 1 */

MpIeee f457(MpIeee x, void * params) {
  params = 0 ;
  if (x == MpIeee( "0.0" ))
    {
      return MpIeee( "0" );
    }
  return MpIeee( "1" )/sqrt(x) ;
}

/* f458(x) = 1/(1 + log(x)^2)^2 */
/* integ(log(x) f458(x),x,0,1) = (Ci(1) sin(1) + (pi/2 - Si(1)) cos(1))/pi 
                               = -0.1892752 */

MpIeee f458(MpIeee x, void * params) {
  params = 0 ;

  if (x == MpIeee( "0.0" )) 
    {
      return MpIeee( "0" );
    }
  else 
    {
      MpIeee u=  log(x);
      MpIeee v=  MpIeee( "1" ) + u * u;
      
      return MpIeee( "1.0" ) / (v * v) ;
    }
}

/* f459(x) = 1/(5 x^3 + 6) */
/* integ(f459/(x-0),x,-1,5) = log(125/631)/18 */

MpIeee f459(MpIeee x, void * params) {
  params = 0 ;
  return MpIeee( "1.0" ) / (MpIeee( "5.0" ) * x * x * x + MpIeee( "6.0" )) ;
}

/* myfn1(x) = exp(-x - x^2) */
/* integ(myfn1,x,-inf,inf) = sqrt(pi) exp(-1/4) */

MpIeee myfn1(MpIeee x, void * params) {
  params = 0;
  return exp(-x - x*x) ;
}

/* myfn2(x) = exp(alpha*x) */
/* integ(myfn2,x,-inf,b) = exp(alpha*b)/alpha */

MpIeee myfn2(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params ;
  return exp(alpha*x) ;
}

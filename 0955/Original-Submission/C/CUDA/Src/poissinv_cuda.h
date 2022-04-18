/*

This software was written by Mike Giles, copyright University of Oxford, 
and is provided under the terms of the GNU GPLv3 license:
http://www.gnu.org/licenses/gpl.html

Commercial users who would like to use the software under a more
permissive license, such as BSD, should contact the author:
mike.giles@maths.ox.ac.uk

*/

// include CUDA header which defines Inf and NaN constants

#include <math_constants.h>


//
// This single precision function computes the inverse
// of the Poisson CDF, using at most 10 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// For lam < 1e7,  max |error| no more than 1
//                 ave |error| < 2e-7 * max(1, sqrt(lam))
//
// For lam > 1e7, the errors will be about 1 ulp.
//

__device__ inline float poissinvf(float u, float lam) {

  float s, t, x=0.0f;

// handle exceptions

  if (u <0.0f) return CUDART_NAN_F;
  if (u==0.0f) return 0.0f;
  if (u==1.0f) return CUDART_INF_F;
  if (u >1.0f) return CUDART_NAN_F;

// large lam

  if (lam > 4.0f) {
    s = normcdfinvf(u)*rsqrtf(lam);

// use polynomial approximations in central region

    if ( (s>-0.6833501f) && (s<1.777993f) ) {;
      float rm;   
                                             
//  polynomial approximation to f^{-1}(s) - 1
                                             
      rm =  2.82298751e-07f;                 
      rm = -2.58136133e-06f + rm*s;          
      rm =  1.02118025e-05f + rm*s;          
      rm = -2.37996199e-05f + rm*s;          
      rm =  4.05347462e-05f + rm*s;          
      rm = -6.63730967e-05f + rm*s;          
      rm =  0.000124762566f + rm*s;          
      rm = -0.000256970731f + rm*s;          
      rm =  0.000558953132f + rm*s;          
      rm =  -0.00133129194f + rm*s;          
      rm =   0.00370367937f + rm*s;          
      rm =   -0.0138888706f + rm*s;          
      rm =     0.166666667f + rm*s;          
      rm =             s + s*(rm*s);         
                                                
//  polynomial approximation to correction c0(r)
                                                
      t  =   1.86386867e-05f;                   
      t  =  -0.000207319499f + t*rm;            
      t  =     0.0009689451f + t*rm;            
      t  =   -0.00247340054f + t*rm;            
      t  =    0.00379952985f + t*rm;            
      t  =   -0.00386717047f + t*rm;            
      t  =    0.00346960934f + t*rm;            
      t  =   -0.00414125511f + t*rm;            
      t  =    0.00586752093f + t*rm;            
      t  =   -0.00838583787f + t*rm;            
      t  =     0.0132793933f + t*rm;            
      t  =     -0.027775536f + t*rm;            
      t  =      0.333333333f + t*rm;            
                                    
//  O(1/lam) correction             
                                    
      x  =   -0.00014585224f;       
      x  =    0.00146121529f + x*rm;
      x  =   -0.00610328845f + x*rm;
      x  =     0.0138117964f + x*rm;
      x  =    -0.0186988746f + x*rm;
      x  =     0.0168155118f + x*rm;
      x  =     -0.013394797f + x*rm;
      x  =     0.0135698573f + x*rm;
      x  =    -0.0155377333f + x*rm;
      x  =     0.0174065334f + x*rm;
      x  =    -0.0198011178f + x*rm;
      x  = __fdividef(x,lam);

//    sum from smallest to largest to minimise rounding error;
//    use of __fadd_rd to round down final sum is important for
//    very large values of lambda to ensure correct rounding

      x = floorf( __fadd_rd(lam, (x+t)+lam*rm) );
    }

// otherwise use Newton iteration

    else if (s > -sqrtf(2.0f)) {
      float r, r2, s2;

      r = 1.0f + s;
      if (r<0.1f) r = 0.1f;

      do {
        t  = __logf(r);
        r2 = r;
        s2 = sqrtf(2.0f*((1.0f-r) + r*t));
        if (r<1.0f) s2 = -s2;
        r = r2 - (s2-s)*s2/t;
        if (r<0.1f*r2) r = 0.1f*r2;
      } while (fabsf(r-r2)>1e-5f);

      t = __logf(r);
      x = lam*r + __logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;
      x = floorf( x - 0.0218f/(x+0.065f*lam) );
    }
  }

// bottom-up summation

  if (x<10.0f) {
    float lami, del;

    lami = 1.0f/lam;
    x    = 0.0f;
    t    = expf(0.5f*lam);
    del  = 0.0f;
    if (u>0.5f)
      del  = t*(1e-6f*t);
    s    = 1.0f - t*(u*t) + del;

    while (s<0.0f) {
      x  += 1.0f;
      t   = x*lami;
      del = t*del;
      s   = t*s + 1.0f;
    }

// top-down summation if needed

    if (s < 2.0f*del) {
      del = 1e6f*del;
      t   = 1e7f*del;
      del = (1.0f-u)*del;

      while (del<t) {
        x   += 1.0f;
        del *= x*lami;
      }

      s = del;
      t = 1.0f;
      while (s>0.0f) {
        t *= x*lami;
        s -= t;
        x -= 1.0f;
      }
    }
  }

  return x;
}


//
// This double precision function computes the inverse
// of the Poisson CDF, using about 30 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// For lam < 1e15,  max |error| no more than 1
//  ave |error| < 1e-16*max(4,lam) for lam < 1e9
//              < 1e-6             for lam < 1e15
//
// For lam > 1e15, the errors will be about 1 ulp.
//

//
// naming convention: double precision variables are capitalised
//

__device__ inline double poissinv(double U, double Lam) {

  double X=0.0, Xi, S, T;

// handle exceptions

  if (U <0.0) return CUDART_NAN;
  if (U==0.0) return 0.0;
  if (U==1.0) return CUDART_INF;
  if (U >1.0) return CUDART_NAN;

// large lam

  if (Lam > 4.0) {
    float s, t, del;

    S   = normcdfinv(U)*rsqrt(Lam);
    s   = (float) S;
    del = 2.0e-6f;

// use polynomial approximations in central region
 
    if ( (s>-0.6833501f) && (s<1.777993f) ) {;
      float rm, x;
                                             
//  polynomial approximation to f^{-1}(s) - 1
                                             
      rm =  2.82298751e-07f;                 
      rm = -2.58136133e-06f + rm*s;          
      rm =  1.02118025e-05f + rm*s;          
      rm = -2.37996199e-05f + rm*s;          
      rm =  4.05347462e-05f + rm*s;          
      rm = -6.63730967e-05f + rm*s;          
      rm =  0.000124762566f + rm*s;          
      rm = -0.000256970731f + rm*s;          
      rm =  0.000558953132f + rm*s;          
      rm =  -0.00133129194f + rm*s;          
      rm =   0.00370367937f + rm*s;          
      rm =   -0.0138888706f + rm*s;          
      rm =     0.166666667f + rm*s;          

      S +=                 s*(rm*s);
      rm = (float) S;
                                                
//  polynomial approximation to correction c0(r)
                                                
      t  =   1.86386867e-05f;                   
      t  =  -0.000207319499f + t*rm;            
      t  =     0.0009689451f + t*rm;            
      t  =   -0.00247340054f + t*rm;            
      t  =    0.00379952985f + t*rm;            
      t  =   -0.00386717047f + t*rm;            
      t  =    0.00346960934f + t*rm;            
      t  =   -0.00414125511f + t*rm;            
      t  =    0.00586752093f + t*rm;            
      t  =   -0.00838583787f + t*rm;            
      t  =     0.0132793933f + t*rm;            
      t  =     -0.027775536f + t*rm;            
      t  =      0.333333333f + t*rm;            
                                    
//  O(1/lam) correction             
                                    
      x  =   -0.00014585224f;       
      x  =    0.00146121529f + x*rm;
      x  =   -0.00610328845f + x*rm;
      x  =     0.0138117964f + x*rm;
      x  =    -0.0186988746f + x*rm;
      x  =     0.0168155118f + x*rm;
      x  =     -0.013394797f + x*rm;
      x  =     0.0135698573f + x*rm;
      x  =    -0.0155377333f + x*rm;
      x  =     0.0174065334f + x*rm;
      x  =    -0.0198011178f + x*rm;
      x  = __fdividef(x,(float) Lam);

//    sum from smallest to largest to minimise rounding error;
//    use of __dadd_rd to round down final sum is important
//    for large values of lambda to ensure correct rounding

//    this corresponds to  (x+delta) in the paper
      S = __dadd_rd(Lam, ((x+del)+t)+Lam*S);
    }

// otherwise use Newton iteration

    else if (s > -sqrtf(2.0f)) {
      float r, r2, s2;

      r = 1.0f + s;
      if (r<0.1f) r = 0.1f;

      do {
        t  = __logf(r);
        r2 = r;
        s2 = sqrtf(2.0f*((1.0f-r) + r*t));
        if (r<1.0f) s2 = -s2;
        r = r2 - (s2-s)*s2/t;
        if (r<0.1f*r2) r = 0.1f*r2;
      } while (fabsf(r-r2)>1e-5f);

      t   = __logf(r);
      S   = Lam*r + __logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;
      //      S   = S - (8.2/405.0)/(S+0.025*Lam);
      S   = S - 0.0218/(S+0.065*Lam);
      del = 0.01/S;
      S   = S + del;
    }

// if x>10, round down to nearest integer, and check accuracy

    X = floor(S);

    if (S>10.0 && S<X+2.0*del) {

// correction procedure based on Temme approximation (double precision)

      if (X>0.5*Lam && X<2.0*Lam) {
        double Eta, B0, B1;

        Xi = 1.0 / X;
        Eta = X / Lam;
        Eta = sqrt(2.0*(1.0-Eta+Eta*log(Eta))/Eta);
        if (X>Lam) Eta = -Eta;

        B1 =  8.0995211567045583e-16;              S = B1;      
        B0 = -1.9752288294349411e-15;              S = B0 + S*Eta;
        B1 = -5.1391118342426808e-16 + 25.0*B1*Xi; S = B1 + S*Eta;
        B0 =  2.8534893807047458e-14 + 24.0*B0*Xi; S = B0 + S*Eta;
        B1 = -1.3923887224181616e-13 + 23.0*B1*Xi; S = B1 + S*Eta;
        B0 =  3.3717632624009806e-13 + 22.0*B0*Xi; S = B0 + S*Eta;
        B1 =  1.1004392031956284e-13 + 21.0*B1*Xi; S = B1 + S*Eta;
        B0 = -5.0276692801141763e-12 + 20.0*B0*Xi; S = B0 + S*Eta;
        B1 =  2.4361948020667402e-11 + 19.0*B1*Xi; S = B1 + S*Eta;
        B0 = -5.8307721325504166e-11 + 18.0*B0*Xi; S = B0 + S*Eta;
        B1 = -2.5514193994946487e-11 + 17.0*B1*Xi; S = B1 + S*Eta;
        B0 =  9.1476995822367933e-10 + 16.0*B0*Xi; S = B0 + S*Eta;
        B1 = -4.3820360184533521e-09 + 15.0*B1*Xi; S = B1 + S*Eta;
        B0 =  1.0261809784240299e-08 + 14.0*B0*Xi; S = B0 + S*Eta;
        B1 =  6.7078535434015332e-09 + 13.0*B1*Xi; S = B1 + S*Eta;
        B0 = -1.7665952736826086e-07 + 12.0*B0*Xi; S = B0 + S*Eta;
        B1 =  8.2967113409530833e-07 + 11.0*B1*Xi; S = B1 + S*Eta;
        B0 = -1.8540622107151585e-06 + 10.0*B0*Xi; S = B0 + S*Eta;
        B1 = -2.1854485106799979e-06 +  9.0*B1*Xi; S = B1 + S*Eta;
        B0 =  3.9192631785224383e-05 +  8.0*B0*Xi; S = B0 + S*Eta;
        B1 = -0.00017875514403292177 +  7.0*B1*Xi; S = B1 + S*Eta;
        B0 =  0.00035273368606701921 +  6.0*B0*Xi; S = B0 + S*Eta;
        B1 =   0.0011574074074074078 +  5.0*B1*Xi; S = B1 + S*Eta;
        B0 =   -0.014814814814814815 +  4.0*B0*Xi; S = B0 + S*Eta;
        B1 =    0.083333333333333329 +  3.0*B1*Xi; S = B1 + S*Eta;
        B0 =    -0.33333333333333331 +  2.0*B0*Xi; S = B0 + S*Eta;
        S  = S / (1.0 + B1*Xi);

        S = S*exp(-0.5*X*Eta*Eta)*rsqrt(2.0*3.141592653589793*X);
        if (X<Lam) {
          S += 0.5*erfc(Eta*sqrt(0.5*X));
          if (S > U) X -= 1.0;
        }
        else {
          S -= 0.5*erfc(-Eta*sqrt(0.5*X));
          if (S > U-1.0) X -= 1.0;
        }
      }

// sum downwards or upwards

      else {
        Xi = 1.0 / X;
        S = - (691.0/360360.0);
        S =   (1.0/1188.0) + S*Xi*Xi;
        S = - (1.0/1680.0) + S*Xi*Xi;
        S =   (1.0/1260.0) + S*Xi*Xi;
        S = - (1.0/360.0)  + S*Xi*Xi;
        S =   (1.0/12.0)   + S*Xi*Xi;
        S =                  S*Xi;
        S = (X - Lam) - X*log(X/Lam) - S;

        if (X<Lam) {
          T  = exp(-0.5*S);
          S  = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
          T  = 1.0;
          Xi = X;
          for (int i=1; i<50; i++) {
            Xi -= 1.0;
            T  *= Xi/Lam;
            S  += T;
          }
          if (S > 0.0) X -= 1.0;
        }

        else {
          T  = exp(-0.5*S);
          S  = 1.0 - T*((1.0-U)*T) * sqrt(2.0*3.141592653589793*X);
          Xi = X;
          for (int i=0; i<50; i++) {
            Xi += 1.0;
            S   = S*Xi/Lam + 1.0;
          }
          if (S < 0.0) X -= 1.0;
        }
      }
    }
  }

// bottom-up summation

  if (X<10.0) {
    double Del;

    X   = 0.0;
    T   = exp(0.5*Lam);
    Del = 0.0;
    if (U>0.5) Del = T*(1e-13*T);
    S   = 1.0 - T*(U*T) + Del;

    while (S<0.0) {
      X  += 1.0;
      T   = X/Lam;
      Del = T*Del;
      S   = T*S + 1.0;
    }

// top-down summation if needed

    if (S < 2.0*Del) {
      Del = 1e13*Del;
      T   = 1e17*Del;
      Del = (1.0-U)*Del;

      while (Del<T) {
        X   += 1.0;
        Del *= X/Lam;
      }

      S = Del;
      T = 1.0;
      while (S>0.0) {
        T *= X/Lam;
        S -= T;
        X -= 1.0;
      }
    }
  }

  return X;
}


/////////////////////////////////////////////////////////////////
//                                                             //
// This software was written by Mike Giles, 2014               //
//                                                             //
// It is copyright University of Oxford, and provided under    //
// the terms of the GNU GPLv3 license:                         //
// http://www.gnu.org/licenses/gpl.html                        //
//                                                             //
// Commercial users wanting to use the software under a more   //
// permissive license, such as BSD, should contact the author: //
// mike.giles@maths.ox.ac.uk                                   //
//                                                             //
/////////////////////////////////////////////////////////////////

// standard math header file

#include <math.h>

// declare prototype for inverse Normal CDF function
// defined at the bottom of this header file

double normcdfinv_as241(double);

//
// This double precision function computes the inverse
// of the Poisson CDF
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


// As described in the TOMS paper, there are two versions;
// the first is optimised for MIMD execution, whereas
// the second is designed for vector execution

double poissinv_core(double, double, double);

double poissinv(double U, double Lam) {
  return poissinv_core(U,1.0-U,Lam);
}

double poisscinv(double V, double Lam) {
  return poissinv_core(1.0-V,V,Lam);
}

inline double poissinv_core(double U, double V, double Lam) {
  int    i;
  double X=0.0, Xi, W, T, Del, R, R2, S, S2, Eta, B0, B1, Lami=1.0/Lam;

// handle exceptions -- constants defined in <math.h>

  if (U==0.0) return 0.0;
  if (V==0.0) return FP_INFINITE;
  if (!(U>0.0 && V>0.0)) return FP_NAN;   // handles NAN inputs as well

  if (Lam > 4.0) {
    W = normcdfinv_as241(fmin(U,V));
    if (U>V) W = -W;

// use polynomial approximations in central region
 
    if ( fabs(W)<3.0 ) {;
      double Lam_root = sqrt(Lam);

      S = Lam_root*W + (1.0/3.0 + (1.0/6.0)*W*W)*(1.0 - W/(12.0*Lam_root));

      Del = (1.0 /160.0);
      Del = (1.0 / 80.0) + Del*(W*W);
      Del = (1.0 / 40.0) + Del*(W*W);
      Del = Del * Lami;

      S = Lam + (S + Del);
    }

// otherwise use Newton iteration

    else {
      S = W / sqrt(Lam);
      R = 1.0 + S;
      if (R<0.1) R = 0.1;

      do {
        T  = log(R);
        R2 = R;
        S2 = sqrt(2.0*((1.0-R) + R*T));
        if (R<1.0) S2 = -S2;
        R = R2 - (S2-S)*S2/T;
        if (R<0.1*R2) R = 0.1*R2;
      } while (fabs(R-R2)>1e-8);

      T   = log(R);
      S   = Lam*R + log(sqrt(2.0*R*((1.0-R)+R*T))/fabs(R-1.0)) / T;
      S   = S - 0.0218/(S+0.065*Lam);
      Del = 0.01/S;
      S   = S + Del;
    }

// if x>10, round down to nearest integer, and check accuracy

    X = floor(S);

    if (S>10.0 && S<X+2.0*Del) {

// correction procedure based on Temme approximation

      if (X>0.5*Lam && X<2.0*Lam) {

        Xi = 1.0 / X;
        Eta = X * Lami;
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

        S = S*exp(-0.5*X*Eta*Eta)/sqrt(2.0*3.141592653589793*X);
        if (X<Lam) {
          S += 0.5*erfc(Eta*sqrt(0.5*X));
          if (S > U) X -= 1.0;
        }
        else {
          S -= 0.5*erfc(-Eta*sqrt(0.5*X));
          if (S > -V) X -= 1.0;
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
        S = (X - Lam) - X*log(X*Lami) - S;

        if (X<Lam) {
          T  = exp(-0.5*S);
          S  = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
          T  = 1.0;
          Xi = X;
          for (i=1; i<50; i++) {
            Xi -= 1.0;
            T  *= Xi*Lami;
            S  += T;
          }
          if (S > 0.0) X -= 1.0;
        }

        else {
          T  = exp(-0.5*S);
          S  = 1.0 - T*(V*T) * sqrt(2.0*3.141592653589793*X);
          Xi = X;
          for (i=0; i<50; i++) {
            Xi += 1.0;
            S   = S*Xi*Lami + 1.0;
          }
          if (S < 0.0) X -= 1.0;
        }
      }
    }
  }

// bottom-up summation

  if (X<10.0) {
    X   = 0.0;
    T   = exp(0.5*Lam);
    Del = 0.0;
    if (U>0.5) Del = T*(1e-13*T);
    S   = 1.0 - T*(U*T) + Del;

    while (S<0.0) {
      X  += 1.0;
      T   = X*Lami;
      Del = T*Del;
      S   = T*S + 1.0;
    }

// top-down summation if needed

    if (S < 2.0*Del) {
      Del = 1e13*Del;
      T   = 1e17*Del;
      Del = V*Del;

      while (Del<T) {
        X   += 1.0;
        Del *= X*Lami;
      }

      S = Del;
      T = 1.0;
      while (S>0.0) {
        T *= X*Lami;
        S -= T;
        X -= 1.0;
      }
    }
  }

  return X;
}



double poissinv_v(double U, double Lam) {
  int    i;
  double X=0.0, Xi, T, Del, Rm, R, R2, S, S2, Eta, B0, B1;

// handle exceptions -- constants defined in <math.h>

  if (U <0.0) return FP_NAN;
  if (U==0.0) return 0.0;
  if (U==1.0) return FP_INFINITE;
  if (U >1.0) return FP_NAN;

// large lam

  if (Lam > 4.0) {
    S   = normcdfinv_as241(U)/sqrt(Lam);
    Del = 2.0e-6;

// use polynomial approximations in central region
 
    if ( (S>-0.6833501) && (S<1.777993) ) {;
                                             
//  polynomial approximation to f^{-1}(s) - 1
                                             
      Rm =  2.82298751e-07;                 
      Rm = -2.58136133e-06 + Rm*S;          
      Rm =  1.02118025e-05 + Rm*S;          
      Rm = -2.37996199e-05 + Rm*S;          
      Rm =  4.05347462e-05 + Rm*S;          
      Rm = -6.63730967e-05 + Rm*S;          
      Rm =  0.000124762566 + Rm*S;          
      Rm = -0.000256970731 + Rm*S;          
      Rm =  0.000558953132 + Rm*S;          
      Rm =  -0.00133129194 + Rm*S;          
      Rm =   0.00370367937 + Rm*S;          
      Rm =   -0.0138888706 + Rm*S;          
      Rm =     0.166666667 + Rm*S;          
      S +=                S*(Rm*S);
      Rm = S;
                                                
//  polynomial approximation to correction c0(r)
                                                
      T  =   1.86386867e-05;                   
      T  =  -0.000207319499 + T*Rm;            
      T  =     0.0009689451 + T*Rm;            
      T  =   -0.00247340054 + T*Rm;            
      T  =    0.00379952985 + T*Rm;            
      T  =   -0.00386717047 + T*Rm;            
      T  =    0.00346960934 + T*Rm;            
      T  =   -0.00414125511 + T*Rm;            
      T  =    0.00586752093 + T*Rm;            
      T  =   -0.00838583787 + T*Rm;            
      T  =     0.0132793933 + T*Rm;            
      T  =     -0.027775536 + T*Rm;            
      T  =      0.333333333 + T*Rm;            
                                   
//  O(1/lam) correction             
                                    
      X  =   -0.00014585224;       
      X  =    0.00146121529 + X*Rm;
      X  =   -0.00610328845 + X*Rm;
      X  =     0.0138117964 + X*Rm;
      X  =    -0.0186988746 + X*Rm;
      X  =     0.0168155118 + X*Rm;
      X  =     -0.013394797 + X*Rm;
      X  =     0.0135698573 + X*Rm;
      X  =    -0.0155377333 + X*Rm;
      X  =     0.0174065334 + X*Rm;
      X  =    -0.0198011178 + X*Rm;
      X  =  X / Lam;

//    sum from smallest to largest to minimise rounding error

      S = Lam + (((X+Del)+T)+Lam*S);
    }

// otherwise use Newton iteration

    else if (S > -sqrt(2.0)) {

      R = 1.0 + S;
      if (R<0.1) R = 0.1;

      do {
        T  = log(R);
        R2 = R;
        S2 = sqrt(2.0*(1.0 - R + R*T));
        if (R<1.0) S2 = -S2;
        R = R2 - (S2-S)*S2/T;
        if (R<0.1*R2) R = 0.1*R2;
      } while (fabs(R-R2)>1e-5);

      T   = log(R);
      S   = Lam*R + log(sqrt(2.0*R*(1.0-R+R*T))/fabs(R-1.0)) / T;
      S   = S - (8.2/405.0)/(S+0.025*Lam);
      Del = 0.01/S;
      S   = S + Del;
    }

// if x>10, round down to nearest integer, and check accuracy

    X = floor(S);

    if (S>10.0 && S<X+2.0*Del) {

// correction procedure based on Temme approximation

      if (X>0.5*Lam && X<2.0*Lam) {

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

        S = S*exp(-0.5*X*Eta*Eta)/sqrt(2.0*3.141592653589793*X);
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
          for (i=1; i<50; i++) {
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
          for (i=0; i<50; i++) {
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

//////////////////////////////////////////////////////////////////////
//                                                                  //
// The routine below is a C version of the code in                  //
//                                                                  //
// ALGORITHM AS241: APPLIED STATS (1988) VOL. 37, NO. 3, 477-44.    //
// http://lib.stat.cmu.edu/apstat/241                               //
//                                                                  //
// The relative error is less than 1e-15, and the accuracy is       //
// verified in the accompanying MATLAB code as241.m                 //
//                                                                  //
//////////////////////////////////////////////////////////////////////

double normcdfinv_as241(double p) {
  
  double q, r, num, den, res;

  q = p - 0.5;
  if (fabs(q) <= 0.425) {
    r = 0.180625 - q*q;

    num =         2.5090809287301226727e+3;
    num = r*num + 3.3430575583588128105e+4;
    num = r*num + 6.7265770927008700853e+4;
    num = r*num + 4.5921953931549871457e+4;
    num = r*num + 1.3731693765509461125e+4;
    num = r*num + 1.9715909503065514427e+3;
    num = r*num + 1.3314166789178437745e+2;
    num = r*num + 3.3871328727963666080e0;

    den =         5.2264952788528545610e+3;
    den = r*den + 2.8729085735721942674e+4;
    den = r*den + 3.9307895800092710610e+4;
    den = r*den + 2.1213794301586595867e+4;
    den = r*den + 5.3941960214247511077e+3;
    den = r*den + 6.8718700749205790830e+2;
    den = r*den + 4.2313330701600911252e+1;
    den = r*den + 1.0000000000e+00;

    res = q * num / den; 

    return res;
  }

  else {

    if (q < 0.0)
      r = p;
    else
      r = 1.0 - p;

    r = sqrt(-log(r));

    if (r <= 5.0) {
      r = r - 1.6;

      num =         7.74545014278341407640e-4;
      num = r*num + 2.27238449892691845833e-2;
      num = r*num + 2.41780725177450611770e-1;
      num = r*num + 1.27045825245236838258e0;
      num = r*num + 3.64784832476320460504e0;
      num = r*num + 5.76949722146069140550e0;
      num = r*num + 4.63033784615654529590e0;
      num = r*num + 1.42343711074968357734e0;

      den =         1.05075007164441684324e-9;
      den = r*den + 5.47593808499534494600e-4;
      den = r*den + 1.51986665636164571966e-2;
      den = r*den + 1.48103976427480074590e-1;
      den = r*den + 6.89767334985100004550e-1;
      den = r*den + 1.67638483018380384940e0;
      den = r*den + 2.05319162663775882187e0;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    }

    else {
      r = r - 5.0;

      num =         2.01033439929228813265e-7;
      num = r*num + 2.71155556874348757815e-5;
      num = r*num + 1.24266094738807843860e-3;
      num = r*num + 2.65321895265761230930e-2;
      num = r*num + 2.96560571828504891230e-1;
      num = r*num + 1.78482653991729133580e0;
      num = r*num + 5.46378491116411436990e0;
      num = r*num + 6.65790464350110377720e0;

      den =         2.04426310338993978564e-15;
      den = r*den + 1.42151175831644588870e-7;
      den = r*den + 1.84631831751005468180e-5;
      den = r*den + 7.86869131145613259100e-4;
      den = r*den + 1.48753612908506148525e-2;
      den = r*den + 1.36929880922735805310e-1;
      den = r*den + 5.99832206555887937690e-1;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    }

    if (q < 0.0)
      res = - res;

    return res;
  }
}


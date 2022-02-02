//-*- C++ -*-
/*
  Random Number Generator
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include "cmplxtype.h"
#include "mat2.h"
#include "mat3.h"

class Random
{
private:
  // static const int RNDGODO; //=48828125;
  // static const int RNDMASK; //=2147483647; // pow(2,31) - 1
  enum { RNDGODO =48828125, RNDMASK =2147483647};

  unsigned long n;

public:
  Random(int seed = 1) : n(seed) {} // seed must be nonzero integer.
  void init(int seed) { n = seed; }
  unsigned long gen() {
    n = n * RNDGODO;
    n = n & RNDMASK;
    return(n);
  }   
  double gen(double m) {
    n = n * RNDGODO;
    n = n & RNDMASK;
    return( double(n)/(double(RNDMASK)+1.0e+0) * m );
  }
  int gen(int m) { return( (int)( gen((double)m) ) ); }
  Complex gen(Complex cm) {
    return( Complex(gen(real(cm)),gen(imag(cm))) );
  }
  
  double operator()(double m) {
    n = n * RNDGODO;
    n = n & RNDMASK;
    return( double(n)/(double(RNDMASK)+1.0e+0) * m );
  }
  int operator()(int m) { return( (int)( gen((double)m) ) ); }
  Complex operator()(Complex cm) {
    return( Complex(gen(real(cm)),gen(imag(cm))) );
  }
  Vec2 operator()(Vec2 m) {
    return ( Vec2(gen(m.x), gen(m.y)) );
  }
  Vec3 operator()(Vec3 m) {
    return ( Vec3(gen(m.x), gen(m.y), gen(m.z)) );
  }

  // Normal Distribution ( average = 0.0 , variance = m*m );
  // ref. "Computer Simulation of Liquids", p347
  double normal(double m=1.0) {
    double zeta = 0.0;
    // each variance is 1/12*m*m and each average is 1/2*m .
    for (int i=0; i<12; i++) zeta += gen(m); 
    // Thus total variance is m*m and total average is 6*m .

    return ( zeta - 6.0*m );
  }
};

#endif // RANDOM_H_



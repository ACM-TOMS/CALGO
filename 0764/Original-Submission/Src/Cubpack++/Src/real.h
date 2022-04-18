/////////////////////////////////////////////////////////
//                                                     //
//    Cubpack++                                        //
//                                                     //
//        A Package For Automatic Cubature             //
//                                                     //
//        Authors : Ronald Cools                       //
//                  Dirk Laurie                        //
//                  Luc Pluym                          //
//                                                     //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// File : real.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#ifndef REAL_H
#define REAL_H
/////////////////////////////////////////////////////////
#include <float.h>

#ifdef FLOAT
  typedef float real;

#define REAL_MAX  FLT_MAX
#define REAL_MIN  FLT_MIN
#define REAL_EPSILON  FLT_EPSILON
#define REAL_MAX_EXP  FLT_MAX_EPS
#define DEFAULT_REL_ERR_REQ  (1.0e-4)

#else
  typedef double real;

#define REAL_MAX  DBL_MAX
#define REAL_MIN  DBL_MIN
#define REAL_EPSILON  DBL_EPSILON
#define REAL_MAX_EXP  DBL_MAX_EXP
#define DEFAULT_REL_ERR_REQ  (1.0e-6)

#endif
/////////////////////////////////////////////////////////
#endif

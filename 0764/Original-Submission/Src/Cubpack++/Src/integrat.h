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
// File : integrat.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   28 Mar 1996     V0.1h(DefaultMaxEval = unsigned long)
/////////////////////////////////////////////////////////
// PURPOSE:
//   Integrate is the algorithm controller of Cubpack++.
//   4 different versions are supplied:
//
//   1)real Integrate(
//      Function f,
//      COMPOUND_REGION& TheCollection,
//      real RequestedAbsolute,
//      real RequestedRelative,
//      unsigned long  MaxEval,
//      );
//   -----------------------------------------------
//   input parameters:
//     RequestedAbsolute: absolute error request
//     RequestedRelative: relative error request
//     MaxEval: bound on the number of function
//           evaluations.
//     f:the integrand
//   return value:
//     the approximation of the integral is returned
//   in-out parameters:
//     TheCollection: collection of regions over which the
//                 integration is to take place.
//                 On return, TheCollection.Integral()
//                 will contain an approximation of the
//                 integral.
//
//
//   2)real Integrate(
//      COMPOUND_REGION& TheCollection,
//      real RequestedAbsolute,
//      real RequestedRelative,
//      unsigned long  MaxEval,
//      );
//   -----------------------------------------------
//   input parameters:
//     RequestedAbsolute: absolute error request
//     RequestedRelative: relative error request
//     MaxEval: bound on the number of function
//           evaluations.
//   return value:
//     the approximation of the integral is returned
//   in-out parameters:
//     TheCollection: collection of regions over which the
//                 integration is to take place.
//                 Before the call of Integrate()
//                 a local integrand must have been
//                 specified for all members of TheCollection
//                 On return, TheCollection.Integral()
//                 will contain an approximation of the
//                 integral.
//
//
//   3)void Integrate(
//      Function f,
//      COMPOUND_REGION& TheCollection,
//      real& Integral,
//      real& AbsErr,
//      real RequestedAbsolute,
//      real RequestedRelative,
//      unsigned long  MaxEval,
//      );
//   -----------------------------------------------
//   input parameters:
//     RequestedAbsolute: absolute error request
//     RequestedRelative: relative error request
//     MaxEval: bound on the number of function
//           evaluations.
//     f:the integrand
//   output parameters:
//     Integral: approximation of the integral
//     AbsError: estimate of the error
//   in-out parameters:
//     TheCollection: collection of regions over which the
//                 integration is to take place.
//                 On return, TheCollection.Integral()
//                 will contain an approximation of the
//                 integral.
//
//
//   4)void Integrate(
//      COMPOUND_REGION& TheCollection,
//      real& Integral,
//      real& AbsErr,
//      real RequestedAbsolute,
//      real RequestedRelative,
//      unsigned long  MaxEval,
//      );
//   -----------------------------------------------
//   input parameters:
//     RequestedAbsolute: absolute error request
//     RequestedRelative: relative error request
//     MaxEval: bound on the number of function
//           evaluations.
//   output parameters:
//     Integral: approximation of the integral
//     AbsError: estimate of the error
//   in-out parameters:
//     TheCollection: collection of regions over which the
//                 integration is to take place.
//                 Before the call of Integrate()
//                 a local integrand must have been
//                 specified for all members of TheCollection
//                 On return, TheCollection.Integral()
//                 will contain an approximation of the
//                 integral.
///////////////////////////////////////////////////////


#include <boolean.h>

#ifndef INTEGRAT_H
#define INTEGRAT_H

const unsigned long DefaultMaxEval =  100000;
const real DefaultAbsErrReq = 0.0;
const real DefaultRelErrReq = DEFAULT_REL_ERR_REQ;

extern
real Integrate( Function f,
                COMPOUND_REGION& R,
                real AbsErrReq = DefaultAbsErrReq,
                real RelErrReq = DefaultRelErrReq,
                unsigned long  MaxEval = DefaultMaxEval);

extern
real Integrate( COMPOUND_REGION& R,
                real AbsErrReq = DefaultAbsErrReq,
                real RelErrReq = DefaultRelErrReq,
                unsigned long  MaxEval = DefaultMaxEval);

extern
void Integrate( Function f,
                COMPOUND_REGION& R,
                real& Integral,
                real& AbsError,
                Boolean& Success,
                real AbsErrReq ,
                real RelErrReq ,
                unsigned long  MaxEval );

extern
void Integrate( COMPOUND_REGION& R,
                real& Integral,
                real& AbsError,
                Boolean& Success,
                real AbsErrReq ,
                real RelErrReq ,
                unsigned long  MaxEval );

#endif

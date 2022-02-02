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
/////////////////////////////////////////////////
//File integrat.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////

#include <eval_ctr.h>
#include <function.h>
#include <real.h>
#include <compreg.h>
#include <iostream.h>
#include <math.h>
////////////////////////////////////////////////


void
Integrate ( COMPOUND_REGION& R,
            real& Integral,
            real& AbsError,
            Boolean& Success,
            real AbsErrReq,
            real RelErrReq,
            unsigned long MaxEval)
        
  {
  EvaluationCounter E;
  E.Start();
  do
    {
    if (R.Hopeless()) break;
    R.Process();
    }
  while (R.AbsoluteError() > AbsErrReq
         && R.AbsoluteError() > fabs(R.Integral())*RelErrReq
         && E.Read() < MaxEval);
  Success = (Boolean)(R.AbsoluteError() <= AbsErrReq ||
                     R.AbsoluteError() <= fabs(R.Integral()*RelErrReq));
  Integral = R.Integral();
  AbsError = R.AbsoluteError();
  }
/////////////////////////////////////////////////////////
real
Integrate ( COMPOUND_REGION& R,
            real AbsErrReq,
            real RelErrReq,
            unsigned long MaxEval)
        
  {
  real Integral;
  real AbsError;
  Boolean Success;
  Integrate(R,Integral,AbsError,Success,AbsErrReq,RelErrReq,MaxEval);
  return Integral;
  }
/////////////////////////////////////////////////////////
void
Integrate ( Function f,
            COMPOUND_REGION& R,
            real& Integral,
            real& AbsError,
            Boolean& Success,
            real AbsErrReq,
            real RelErrReq,
            unsigned long MaxEval)
        
  {
  R.LocalIntegrand(new Integrand(f));
  Integrate(R,Integral,AbsError,Success,AbsErrReq,RelErrReq,MaxEval);
  }
/////////////////////////////////////////////////////////
real
Integrate ( Function f,
            COMPOUND_REGION& R,
            real AbsErrReq,
            real RelErrReq,
            unsigned long MaxEval)
        
  {
  R.LocalIntegrand(new Integrand(f));
  real Integral;
  real AbsError;
  Boolean Success;
  Integrate(R,Integral,AbsError,Success,AbsErrReq,RelErrReq,MaxEval);
  return Integral;
  }
/////////////////////////////////////////////////////////

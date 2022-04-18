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
///////////////////////////////////////////////////
//File integran.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//    9 Sep 1994     V0.1a(operator added)
//   10 Sep 1994     V0.1b(constructors of ReferenceCounting added)
//   25 Jan 1996     V0.1f(unused parameter removed)
//   30 Jan 1996     V0.1g(no longer compare signed and unsigned)
//   28 Mar 1996     V0.1h(long instead of int)
//////////////////////////////////////////////////

#include <integran.h>
#include <trnsfrm.h>
#include <error.h>
#include <math.h>
#include <iostream.h>
//////////////////////////////////////////////////
long Integrand::Number = 0;

//////////////////////////////////////////////////////
Integrand::Integrand (real (*f)(const Point&))
  :TheFunction(f),AppliedTransformations()
   ,ReferenceCounting()
  {
  }
////////////////////////////////////////////////
Integrand::Integrand(const Integrand& I, Transformation* T)
  : TheFunction(I.TheFunction),
    AppliedTransformations(I.AppliedTransformations)
   ,ReferenceCounting(I)
  {
  Pointer<Transformation> p=T;
  AppliedTransformations += p;
  }
/////////////////////////////////////////////////////////
Integrand::Integrand(const Integrand& I)
  : TheFunction(I.TheFunction),
    AppliedTransformations(I.AppliedTransformations)
   ,ReferenceCounting(I)
  {
  }
////////////////////////////////////////////////
long
Integrand::NumberOfEvaluations()
  {
  return Number;
  }
//////////////////////////////////////////////////////
real
Integrand::operator() (const Point& StartingPoint)
  {
  real JacProd =1;
  Point p =StartingPoint;
  for (unsigned int i=0; JacProd!=0 && i<AppliedTransformations.Size();
       i++)
    {
      AppliedTransformations[i]->Transform(JacProd,p);
    };
  if (JacProd ==  0) return  0.0;
  real R =  TheFunction(p)*fabs(JacProd);
  Number++;
  return R;

  }
////////////////////////////////////////////////
real ErrorMessage(const Point&)
  {
  Error(True,"Error : Attempt to integrate without knowing integrand");
  return 0;
  }

Integrand::Integrand()
  :TheFunction(ErrorMessage),AppliedTransformations()
  ,ReferenceCounting()
  {
  }

///////////////////////////////////////////////
Boolean
Integrand::operator==(const Integrand& i) const
  {
  if ((AppliedTransformations.Size() > 0 ) || (i.AppliedTransformations.Size() > 0 ))
     return False;
  return (Boolean)(TheFunction == i.TheFunction);
  }
///////////////////////////////////////////////


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
// File : T2tops.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <T2tops.h>
/////////////////////////////////////////////////////////
T2toPS::T2toPS(ParabolicSegment* p)
  : Transformation(),
    P(p->P()),
    M((p->A()+p->B())/2),
    H((p->B()-p->A())/2)
  {
  a= -0.5* (P-M).Length();
  k= H.X()*(P.Y()-M.Y())-H.Y()*(P.X()-M.X());
  k*=a;
  }
////////////////////////////////////////////////////////
void
T2toPS::Transform(real & w, Point& p)
  {
  if (p.X()<0)
    {
    real X = p.X();
    real Y = a*(X-1)*p.Y();
    p= M+X*H+Y*(P-M);
    w *= k*(X-1);
    }
  else
    {
    real X = p.X();
    real Y = -a*(X+1)*p.Y();
    p= M+X*H+Y*(P-M);
    w *= -k*(X+1);
    };
  }
/////////////////////////////////////////////////////////

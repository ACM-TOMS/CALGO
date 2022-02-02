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
//File E2tostrp.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////
#include <E2tostrp.h>
#include <point.h>
#include <math.h>
////////////////////////////////////////////
void
E2toIS::Transform(real& w, Point& p)
      {
      InfiniteStrip& s = *IS_ptr;
      Point D  = s.B()- s.A();
      Point C(-D.Y(),D.X());
      C = C/C.Length();
      w *= D.Length()*.5*(1-tanh(p.X())*tanh(p.X()));
      p = s.A() + p.Y()*C + (.5*(1+tanh(p.X())))*D;
      }
///////////////////////////////////////////////////
E2toIS::E2toIS(InfiniteStrip* g)
  : Transformation(),IS_ptr(g)
  {
  }
//////////////////////////////////////////////////

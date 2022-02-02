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
//File gr.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////
#include <C2togr.h>
#include <point.h>
////////////////////////////////////////////


void
C2toGR::Transform(real& w, Point& p)
      {
      GeneralizedRectangle& G = *GR_ptr;
      Point M = G.B()-G.A();
      Point C = G.A() + p.X()*M;
      real ml = M.Length();
      real dist = G.Boundary(C), ratio=dist/ml;
      Point P(-ratio*M.Y(),ratio*M.X());
      w *= dist*ml;
      p = C+p.Y()*P;
      }
///////////////////////////////////////////////////
C2toGR::C2toGR(GeneralizedRectangle* g)
  : Transformation(),GR_ptr(g)
  {
  }
//////////////////////////////////////////////////

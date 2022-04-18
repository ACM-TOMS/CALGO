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
//////////////////////////////////////////
//File C2toS2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////

#include <trnsfrm.h>
#include <C2toS2.h>
#include <math.h>

//////////////////////////////////////////
void
PolarToRectangular::Transform(real& w, Point& p)
  { Point P(p.R()*cos(p.Theta()),p.R()*sin(p.Theta()));
    w *= p.R();
    p = P;
  }
///////////////////////////////////////////
PolarToRectangular::PolarToRectangular()
  : Transformation()
  {
  }
///////////////////////////////////////////

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
////////////////////////////////////////////
//File semistrp.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////
#include <semistrp.h>
//////////////////////////////////////////////
SemiInfiniteStrip::SemiInfiniteStrip(const Point& a,const Point& b)
  :Geometry(2),TheA(a),TheB(b)
  {
  }
////////////////////////////////////////////
const Point&
SemiInfiniteStrip::A()
const
  {
  return TheA;
  }
///////////////////////////////////////////
const Point&
SemiInfiniteStrip::B()
const
  {
  return TheB;
  }
///////////////////////////////////////////

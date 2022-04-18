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
//File strip.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////
#include <strip.h>
#include <passbuck.h>
#include <C2.h>
//////////////////////////////////////////////
InfiniteStrip::InfiniteStrip(const Point& a,const Point& b)
  :Geometry(2),TheA(a),TheB(b)
  {
  }
////////////////////////////////////////////
const Point&
InfiniteStrip::A()
const
  {
  return TheA;
  }
///////////////////////////////////////////
const Point&
InfiniteStrip::B()
const
  {
  return TheB;
  }
///////////////////////////////////////////

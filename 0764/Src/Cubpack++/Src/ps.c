
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
// File : ps.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <ps.h>
/////////////////////////////////////////////////////////
ParabolicSegment::ParabolicSegment(const Point& a,
                                   const Point& b,
                                   const Point& p)
  :Geometry(2),TheA(a),TheB(b),TheP(p)
  {
  }
/////////////////////////////////////////////////////////
const Point&
ParabolicSegment::A()
const
  {
  return TheA;
  }
//////////////////////////////////////////////////////////
const Point&
ParabolicSegment::B()
const
  {
  return TheB;
  }
//////////////////////////////////////////////////////////
const Point&
ParabolicSegment::P()
const
  {
  return TheP;
  }
//////////////////////////////////////////////////////////

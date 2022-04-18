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
//////////////////////////////////////////////
//File translat.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(unused parameter removed)
////////////////////////////////////////////////

#include <translat.h>

///////////////////////////////////////////////
Translation::Translation( const Point& Offset)
  :Transformation(),TheOffset(Offset)
  {
  }
/////////////////////////////////////////////
void
Translation::Transform (real &, Point& p)
  { p += TheOffset; }
/////////////////////////////////////////////

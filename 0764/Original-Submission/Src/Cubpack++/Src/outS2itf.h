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
// File : outS2itf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS OUT_CIRCLE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<OutCircle>
//
// PURPOSE:
//   end-user interface to OutCircle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) OUT_CIRCLE(const Point& Center, const Point& B)
//     ----------------------------------------------
//     constructs a OutCircle with Center as its origin
//     and B a point on the boundary.
//
//     2) OUT_CIRCLE(const Point& Center, real Radius)
//     -------------------------------------------
//
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef OUTS2ITF_H
#define OUTS2ITF_H
////////////////////////////////////////////////
#include <userint.h>
#include <S2.h>
#include <outS2.h>
////////////////////////////////////////////////
class OUT_CIRCLE : public USERINTERFACE<OutCircle>
  {
  public:

  OUT_CIRCLE(const Point&,const Point&);
  OUT_CIRCLE(const Point&, real);
  };
////////////////////////////////////////////////
#endif

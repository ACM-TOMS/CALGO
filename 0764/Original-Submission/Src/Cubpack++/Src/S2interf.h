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
// File : S2interf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS CIRCLE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<Circle>
//
// PURPOSE:
//   end-user interface to circles
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) CIRCLE(const Point& Center, const Point& B)
//     ----------------------------------------------
//     constructs a circle with Center as its origin
//     and B a point on the boundary.
//
//     2) CIRCLE(const Point& Center, real Radius)
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
#ifndef S2INTERF_H
#define S2INTERF_H
////////////////////////////////////////////////
#include <userint.h>
#include <S2.h>
////////////////////////////////////////////////
class CIRCLE : public USERINTERFACE<Circle>
  {
  public:

  CIRCLE(const Point&,const Point&);
  CIRCLE(const Point&, real);
  };
////////////////////////////////////////////////
#endif

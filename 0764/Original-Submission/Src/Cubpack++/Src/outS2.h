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
// File : outS2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS OutCircle
// -------------------------------------------------
//
// BASECLASSES:
//   Circle
//
// PURPOSE:
//   the complement of a circle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) OutCircle(const Point& C, const Point& B)
//     --------------------------------------------
//     the complement of a circle centred at C and having
//     B as a point on its boundary.
//
//     2) OutCircle(const Point& C,real Radius)
//     ----------------------------------------
//     the complement of a circle centred at C and having
//     Radius as its radius.
//
//   SELECTORS:
//     1) Processor<OutCircle>* DefaultProcessor() const
//     -------------------------------------------------
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
#ifndef OUTS2_H
#define OUTS2_H
////////////////////////////////////////////////
#include <S2.h>
///////////////////////////////////////////////
class OutCircle : public Circle
  {
  public:

  OutCircle(const Point&,const Point&);
  OutCircle(const Point&,real);
  };
////////////////////////////////////////////////
#endif

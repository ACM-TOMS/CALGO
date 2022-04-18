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
// File : strip.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS InfiniteStrip Rectangle
// -----------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements an infinite strip
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) GeneralizedRectangle(
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     the line segment AB is the diameter of the strip
//
//   SELECTORS:
//     1) const Point& A() const
//     -------------------------
//     returns  the point A (see constructor)
//
//     2) const Point& B() const
//     -------------------------
//     returns the point B (see constructor)
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
//
///////////////////////////////////////////////////////////
#ifndef STRIP_H
#define STRIP_H
//////////////////////////////////////////////////
#include <point.h>
#include <geometry.h>
//////////////////////////////////////////////////
class InfiniteStrip : public  Geometry
  {
  public:

  InfiniteStrip ( const Point&,const Point&);
  const Point& A() const;
  const Point& B() const;

  private:

  Point TheA;
  Point TheB;
  };
//////////////////////////////////////////////////
#endif

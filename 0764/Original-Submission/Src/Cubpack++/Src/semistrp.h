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
// File : semistrp.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS SemiInfiniteStrip
// -----------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements a semi-infinite strip
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) SemiInfiniteStrip(
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     the region lies to the left of AB
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
#ifndef SEMISTRP_H
#define SEMISTRP_H
//////////////////////////////////////////////////
#include <point.h>
#include <geometry.h>
//////////////////////////////////////////////////
class SemiInfiniteStrip : public  Geometry
  {
  public:

  SemiInfiniteStrip ( const Point&,const Point&);
  const Point& A() const;
  const Point& B() const;

  private:

  Point TheA;
  Point TheB;
  };
//////////////////////////////////////////////////
#endif

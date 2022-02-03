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
// File : ps.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS ParabolicSegment
// -------------------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   the geometry of a parabolic segment
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) ParabolicSegment(const Point& A, const Point& B,
//                          const Point& P)
//     ---------------------------------------------------
//     the region is bounded by the straight line AB and the
//     parabolic arc through A and B such that AP and BP are
//     tangents
//
//
//   SELECTORS:
//     1) const Point& A() const
//     --------------------------
//
//     1) const Point& B() const
//     --------------------------
//
//     1) const Point& P() const
//     --------------------------
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
#ifndef PS_H
#define PS_H
/////////////////////////////////////////////////////////
#include <geometry.h>
#include <point.h>
/////////////////////////////////////////////////////////
class ParabolicSegment : public Geometry
  {
  public:

  ParabolicSegment(const Point&,const Point&,const Point&);
  const Point& A() const;
  const Point& B() const;
  const Point& P() const;

  private:

  Point TheA;
  Point TheB;
  Point TheP;
  };
/////////////////////////////////////////////////////////
#endif

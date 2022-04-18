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
// File : polC2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PolarRectangle
// -------------------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   the geometry of a polar rectangle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PolarRectangle(const Point& A, const Point& B,
//                       const Point& C)
//     -------------------------------------------------
//     constructs a Polar rectangle. AB should not be
//     perpendicular to BC.
//
//     2) PolarRectangle(const Point& O, real R1, real R2,
//                       real Theta1, real Theta2)
//     ---------------------------------------------------
//
//   SELECTORS:
//     1) const Point& Center() const
//     ------------------------------
//
//     2) real InnerRadius() const
//     ---------------------------
//
//     3) real OuterRadius() const
//     ---------------------------
//
//     4) real SmallAngle() const
//     --------------------------
//
//     5) real BigAngle() const
//     ------------------------
//
//     6)Processor<PolarRectangle*> DefaultProcessor() const
//     -----------------------------------------------------
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
#ifndef POLC2_H
#define POLC2_H
/////////////////////////////////////////
#include <geometry.h>
#include <point.h>
#include <real.h>
#include <regproc.h>
/////////////////////////////////////////

class PolarRectangle : public Geometry
  {
  public:

  PolarRectangle(
    const Point& A, const Point& B, const Point& C);
  PolarRectangle( const Point& O, real r1, real r2, real theta1,real
     theta2);
  const Point& Center() const;
  real InnerRadius() const;
  real OuterRadius() const;
  real SmallAngle () const;
  real BigAngle   () const;



  private:

  Point TheCenter;
  real TheInnerRadius, TheOuterRadius, TheSmallAngle, TheBigAngle;

  };
//////////////////////////////////////////////
#endif

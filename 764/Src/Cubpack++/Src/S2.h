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
// File : S2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Circle
// --------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//     Circle is a 2D circle with arbitrary center and
//     radius.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Circle(const Point& Center, real Radius)
//     --------------------------------
//     defines the region.
//
//     2) Circle( const Point& Center, Point Boundary)
//     --------------------------------
//     defines the region using its center and a point on the boundary.

//   SELECTORS:
//     1) const Point& Center()const
//     ------------------------------
//
//     2) real Radius() const
//     -----------------------
//
//     3) real Volume() const
//     ----------------------
//
//     4) Processor<Circle>* DefaultProcessor()const
//     ---------------------------------------------
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     None
//
///////////////////////////////////////////////////////////

#ifndef S2_H
#define S2_H
/////////////////////////////////////////
#include <geometry.h>
#include <regproc.h>
#include <point.h>
#include <real.h>
/////////////////////////////////////////

class Circle : public Geometry
  {
  public:

  Circle(const Point& Center,real Radius);
  Circle(const Point& Center,const Point& Boundary);
  real Volume()const ;
  const Point& Center()const;
  real Radius() const;

  private:
  Point TheCenter;
  real TheRadius;

  };
//////////////////////////////////////////////
#endif

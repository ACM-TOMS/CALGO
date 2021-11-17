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
// File : E2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Plane
// ----------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements the geometry of a plane
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Plane()
//     ----------------------------
//     constructs a plane with center=(0,0) and the unit
//     on both coordinate axis = 1.
//
//     2) Plane(const Point& center)
//     ---------------------------------------------------
//     Constructs a plane given an integrand an a center point.
//     the unit on both axis = 1.
//
//     3) Plane(const Point& center,
//              const real& xscale, const real& yscale)
//     ----------------------------------------------------
//     same as above, but now the scales can be set.
//
//
//   SELECTORS:
//     1) real ScaleX() const
//     ---------------------------
//     returns a scale factor for the X coordinate.
//
//     2) real ScaleX() const
//     ---------------------------
//     returns a scale factor for the X coordinate.
//
//     3) const Point& Center() const
//     ------------------------------
//
//     4) Processor<Plane>* DefaultProcessor() const
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

#ifndef E2_H
#define E2_H
//////////////////////////////////////////////////
#include <geometry.h>
#include <real.h>
#include <point.h>
#include <regproc.h>
//////////////////////////////////////////////////
class Plane : public  Geometry
  {
  public:

  Plane ();
  Plane (const Point&);
  Plane (const Point&, real, real);
  real ScaleX() const ;
  real ScaleY() const ;
  const Point& Center()const;

  private:
  real xscale,yscale;
  Point TheCenter;
  };
//////////////////////////////////////////////////

#endif

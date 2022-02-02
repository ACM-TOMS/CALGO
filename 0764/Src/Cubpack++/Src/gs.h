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
// File : gs.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Generalized Sector
// -----------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements a generalized sector
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) GeneralizedSector(Function f,
//        real alpha, real beta, const Point& Center)
//     ---------------------------------------------------
//      the region is bounded by 2 straight lines going
//      through Center and having an angle of alpha and
//      beta respectively with the positive X-axis. a
//      third boundary line is provided by f, which denotes
//      the distance from the Center to the boundary for
//      all angles between alpha and beta.
//
//   SELECTORS:
//     1) real Alpha() const
//     -------------------------
//     returns  alpha (see constructor)
//
//     2) real Beta() const
//     -------------------------
//     returns beta (see constructor)
//
//     3) real Boundary(real P)const
//     -------------------------------------
//     evaluates the boundary function at P(see constructor)
//
//    4) const Point& Center() const
//    -------------------------------
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
//
///////////////////////////////////////////////////////////
#ifndef GS_H
#define GS_H
//////////////////////////////////////////////////

#include <point.h>
#include <geometry.h>
//////////////////////////////////////////////////
typedef real (*RealFunction)(real);
//////////////////////////////////////////////////
class GeneralizedSector : public  Geometry
  {
  public:

  GeneralizedSector (real(*)(real),
       real,real,const Point&);
  real Alpha() const;
  real Beta() const;
  real Boundary(real) const;
  RealFunction Boundary() const;
  const Point& Center() const;

  private:

  Point TheCenter;
  real TheAlpha,TheBeta;
  real(*TheBoundary)(real);
  };
//////////////////////////////////////////////////
#endif

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
// File : polC2itf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS POLAR_RECTANGLE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<PolarRectangle>
//
// PURPOSE:
//   End-user interface to PolarRectangles
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) POLAR_RECTANGLE(const Point& A,
//                       const Point& B,
//                       const Point& C)
//     ---------------------------------
//     constructs a polar rectangle with vertices A,B,C
//     AB shouldn't be perpendicular to BC
//
//     2) POLAR_RECTANGLE( const& Point& O,
//         real r1, real r2, real theta1,real theta2)
//     ----------------------------------------------
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
#ifndef POLC2ITF_H
#define POLC2ITF_H
/////////////////////////////////////////
#include <userint.h>
#include <polC2.h>
/////////////////////////////////////////
class POLAR_RECTANGLE : public USERINTERFACE<PolarRectangle>
  {
  public:

  POLAR_RECTANGLE(const Point&,const Point&, const Point&);
  POLAR_RECTANGLE(const Point&,real,real,real,real);
  };
/////////////////////////////////////////
#endif

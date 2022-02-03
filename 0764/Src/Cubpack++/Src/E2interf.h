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
// File : E2interf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PLANE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<Plane>
//
// PURPOSE:
//   End-user interface to planes
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1)PLANE()
//     --------
//
//     2)PLANE(const Point& Center)
//     ----------------------------
//     Center indicates where in the plane the
//     activity is the greatest.
//
//     3)PLANE(const Point& Center,
//                real ScaleX, real ScaleY)
//     ------------------------------------
//     same as 2) but with an indication of scales.
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
#ifndef E2INTERF_H
#define E2INTERF_H
////////////////////////////////////////////////
#include <userint.h>
#include <E2.h>
#include <point.h>
#include <real.h>
////////////////////////////////////////////////
class PLANE : public USERINTERFACE<Plane>
  {
  public:

  PLANE();
  PLANE(const Point&);
  PLANE(const Point&, real, real);
  };
////////////////////////////////////////////////
#endif

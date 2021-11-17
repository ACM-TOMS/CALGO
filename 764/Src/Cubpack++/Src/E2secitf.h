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
// File : E2secitf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PLANE_SECTOR
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<PlaneSector>
//
// PURPOSE:
//   End-user interface to PlaneSectors
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PLANE_SECTOR(const Point& A,
//                       const Point& B,
//                       const Point& C)
//     ---------------------------------
//     constructs a plane sector with vertices A,B,C
//     AB shouldn't be perpendicular to BC
//
//     2) PLANE_SECTOR(const Point& O,real r,real theta1,real theta2)
//      ---------------------------------------------------------------
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
#ifndef E2SECITF_H
#define E2SECITF_H
/////////////////////////////////////////
#include <userint.h>
#include <E2sec.h>
/////////////////////////////////////////
class PLANE_SECTOR : public USERINTERFACE<PlaneSector>
  {
  public:

  PLANE_SECTOR(const Point&,const Point&, const Point&);
  PLANE_SECTOR(const Point&,real,real,real);
  };
/////////////////////////////////////////
#endif

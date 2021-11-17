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
// File : C2toS2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// DEFINITION OF CLASS PolarToRectangular
// -------------------------------
//
// BASECLASSES:
//   Transformation
//
//
// PURPOSE:
//   implements a transformation from polar to
//   rectangular coordinates
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PolarToRectangular()
//     ------------------------
//
//   SELECTORS:
//     None
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) void Transform (real& w, Point& p)
//     ---------------------------------------
//     Multiplies w by the Jacobian at p and then replaces p by
//     the transformed point.
//
///////////////////////////////////////////////////////////

#ifndef C2TOS2_H
#define C2TOS2_H
////////////////////////////////////////////////////////

#include <trnsfrm.h>
#include <point.h>

///////////////////////////////////////////////////////
class PolarToRectangular :public Transformation
  {
  public:

  PolarToRectangular();
  void Transform (real &w, Point& p);
  };
////////////////////////////////////////////////////////////
#endif

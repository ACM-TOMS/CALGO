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
// File : invert.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Invert
// --------------------------
//
// BASECLASSES:
//   Transformation
//
//
// PURPOSE:
//   Inverts with respect to a given (arbitrary) circle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Invert(Circle* )
//     --------------------
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

#ifndef INVERT_H
#define INVERT_H
////////////////////////////////////////////////////////

#include <trnsfrm.h>
#include <point.h>
#include <S2.h>

///////////////////////////////////////////////////////
class Invert :public Transformation
  {
  public:
  Invert(Circle*);
  void Transform (real &w, Point& p);
  private:

  Pointer <Circle>  C_ptr;
  real RadiusSq;
  };
////////////////////////////////////////////////////////////
#endif

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
// File : trnsfrm.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Transformation
// ----------------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
//
// PURPOSE:
//   a pure virtual base class providing a common interface
//   to all transformations.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Transformation()
//     --------------------
//
//   SELECTORS:
//     None
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) virtual void Transform (real& w, Point& p)
//     -----------------------------------------------
//     Multiplies w by the Jacobian at p and then replaces p by
//     the transformed point.
//
///////////////////////////////////////////////////////////
#ifndef TRNSFRM_H
#define TRNSFRM_H
#include <point.h>
#include <refcount.h>

class Transformation: public ReferenceCounting
  {
  public:

  Transformation();
  virtual void Transform (real& w, Point& p) =0;
  virtual ~Transformation();
  };

#endif

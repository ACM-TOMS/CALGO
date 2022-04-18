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
// File : E2tostrp.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS E2toIS
// -------------------------------------------------
//
// BASECLASSES:
//   Transformation
//
// PURPOSE:
//   transformation from plane to  infinite strip
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) E2toIS(InfiniteStrip*)
//     ------------
//     The target region has to be supplied
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
//     1) void Transform(real & w, Point & p)
//     --------------------------------------
//     multiplies w by the Jacobian at p and then
//     replaces p by the transformed  point
//
/////////////////////////////////////////////////////////
#ifndef E2TOSTRP_H
#define  E2TOSTRP_H
/////////////////////////////////////////////////////////
#include <trnsfrm.h>
#include <strip.h>
#include <E2.h>
#include <pointer.h>
///////////////////////////////////////////////////
class E2toIS : public Transformation
  {
  public:

  E2toIS( InfiniteStrip*);
  void Transform(real& w, Point& p);

  private:

  Pointer<InfiniteStrip> IS_ptr;
  };
//////////////////////////////////////////////////////
#endif

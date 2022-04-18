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
// File : sttosmst.h
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
//   transformation from strip to semistrip
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) IStoSIS(SemiInfiniteStrip*)
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
#ifndef STTOSMST_H
#define  STTOSMST_H
/////////////////////////////////////////////////////////
#include <trnsfrm.h>
#include <strip.h>
#include <semistrp.h>
#include <pointer.h>
///////////////////////////////////////////////////
class IStoSIS : public Transformation
  {
  public:

  IStoSIS( SemiInfiniteStrip*);
  void Transform(real& w, Point& p);

  private:

  Pointer<SemiInfiniteStrip> SIS_ptr;
  };
//////////////////////////////////////////////////////
#endif

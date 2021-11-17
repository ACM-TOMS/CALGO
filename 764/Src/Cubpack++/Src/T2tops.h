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
// File : T2tops.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member declarations)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS C2toGR
// -------------------------------------------------
//
// BASECLASSES:
//   Transformation
//
// PURPOSE:
//   transformation from triangle to parabolic segment
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) T2toPS(ParabolicSegment* g)
//     ------------
//     the target region has to be supplied
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
#ifndef T2TOPS_H
#define  T2TOPS_H
/////////////////////////////////////////////////////////
#include <trnsfrm.h>
#include <ps.h>
#include <T2.h>
#include <pointer.h>
///////////////////////////////////////////////////
class T2toPS : public Transformation
  {
  public:

  T2toPS( ParabolicSegment*);
  void Transform(real& w, Point& p);

  private:

  Point P,M,H;
  real a,k;
  };
//////////////////////////////////////////////////////
#endif

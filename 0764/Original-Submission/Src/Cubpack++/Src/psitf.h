
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
// File : psitf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PARABOLIC_SEGMENT
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<ParabolicSegment>
//
// PURPOSE:
//   end-user interface to a parabolic segment.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PARABOLIC_SEGMENT(const Point& A,
//                          const Point& B,
//                          const Point& P)
//     -------------------------------------
//     the region  is bounded by the straight line AB
//     and the parabolic arc through A and B such that
//     AP and BP are tangents.
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
#ifndef PSITF_H
#define PSITF_H
/////////////////////////////////////////////////////////
#include <userint.h>
#include <ps.h>
/////////////////////////////////////////////////////////
class PARABOLIC_SEGMENT :  public USERINTERFACE<ParabolicSegment>
  {
  public:

  PARABOLIC_SEGMENT(const Point& A,const Point& B,const Point& p);
  };
////////////////////////////////////////////////////////////
#endif

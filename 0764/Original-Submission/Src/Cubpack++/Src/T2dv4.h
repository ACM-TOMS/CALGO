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
// File : T2dv4.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Triangle_Divide4
// ------------------------------------------
//
// BASECLASSES:
//   SameShapeDivisor<Triangle>
//
//
// PURPOSE:
//   provides a means to cut Triangles in 4.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Triangle_Divide4()
//     ---------------------
//
//   SELECTORS:
//     1) NumberOfParts() const
//     ------------------------
//     returns 4
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) void Apply(Triangle& H,
//     Stack<Triangle>& S,const Vector<unsigned int>& D)
//     --------------------------------
//     divides H into parts and returns them in S
//
///////////////////////////////////////////////////////////

#ifndef T2DV4_H
#define T2DV4_H

//////////////////////////////////////////

#include <samediv.h>
#include <T2.h>

////////////////////////////////////////

class Triangle_Divide4 :public SameShapeDivisor<Triangle>
  {

  public:

  Triangle_Divide4();
  void Apply(const Triangle&, Stack<Triangle>&, const Vector<unsigned int>&);
  int NumberOfParts() const {return 4;};

  };
//////////////////////////////////////////
#endif

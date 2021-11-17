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
// File : C2dv2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Parallelogram_Divide2
// -------------------------------------------------
//
// BASECLASSES:
//   SameShapeDivisor<Parallelogram>
//
// PURPOSE:
//   to cut a Parallelogram in 2 new ones
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Parallelogram_Divide4()
//     ------------------------
//
//   SELECTORS:
//     1) int NumberOfParts() const
//     ----------------------------
//     returns 2
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     1) void Apply(const Parallelogram& P, Stack<Parallelogram>& S
//                   , const Vector<unsigned int>& DiffOrder)
//     ----------------------------------------------------------
//     divides P by halving its sides. the new parts are
//     returned in S. If DiffOrder[0] > DiffOrder[1] then P
//     is divided by halving the edge {Vertex(0),Vertex(1)}
//     otherwise the edge {Vertex(0),Vertex(2)} is halved.
//
/////////////////////////////////////////////////////////
#ifndef C2DV2_H
#define C2DV2_H

//////////////////////////////////////////

#include <samediv.h>
#include <C2.h>

////////////////////////////////////////

class Parallelogram_Divide2 :public SameShapeDivisor<Parallelogram>
  {

  public:

  Parallelogram_Divide2();
  void Apply(const Parallelogram&, Stack<Parallelogram>&, 
    const Vector<unsigned int>&);
  int NumberOfParts() const {return 2;};

  };
//////////////////////////////////////////
#endif

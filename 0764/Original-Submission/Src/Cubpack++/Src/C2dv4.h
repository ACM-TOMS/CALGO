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
// File : C2dv4.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Parallelogram_Divide4
// -------------------------------------------------
//
// BASECLASSES:
//   SameShapeDivisor<Parallelogram>
//
// PURPOSE:
//   Divides a Parallelogram in 4 new ones.
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
//     returns 4
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
//     returned in S,DiffOrder isn't used.
//
/////////////////////////////////////////////////////////
#ifndef C2DV4_H
#define C2DV4_H

//////////////////////////////////////////

#include <samediv.h>
#include <C2.h>

////////////////////////////////////////

class Parallelogram_Divide4 :public SameShapeDivisor<Parallelogram>
  {

  public:

  Parallelogram_Divide4();
  void Apply(const Parallelogram&, Stack<Parallelogram>&, 
    const Vector<unsigned int>&);
  int NumberOfParts() const {return 4;};

  };
//////////////////////////////////////////
#endif

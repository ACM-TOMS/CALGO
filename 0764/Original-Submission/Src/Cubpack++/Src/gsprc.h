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
// File : gsprc.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS GeneralizedSector_Processor
// -------------------------------------------------
//
// BASECLASSES:
//   Processor<GeneralizedSector>
//
// PURPOSE:
//   transforms the generalized sector into a rectangle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1)GeneralizedSector_Processor()
//     ----------------------------
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
//     1) void Process(Stack<AtomicRegion>&)
//     -------------------------------------------------
//     see Processor<>
//
//     2) virtual Processor<GEOMETRY>* NewCopy() const=0
//     -------------------------------------------------
//
//     makes a new copy (using the copy constructor) and
//     returns a pointer to it.
/////////////////////////////////////////////////////////
#ifndef GSPRC_H
#define GSPRC_H
//////////////////////////////////////////////
#include <gs.h>
#include <regproc.h>
/////////////////////////////////////////////
class GeneralizedSector_Processor :
     public Processor<GeneralizedSector>
  {
  public:

  GeneralizedSector_Processor();
  void Process(Stack<AtomicRegion>&);
  Processor<GeneralizedSector>* NewCopy() const ;
  };
/////////////////////////////////////////////
#endif

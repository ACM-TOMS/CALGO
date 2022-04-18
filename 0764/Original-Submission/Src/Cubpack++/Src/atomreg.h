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
// File : atomreg.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS AtomicRegion
// --------------------------------
//
// BASECLASSES:
//   Region
//
// PURPOSE:
//   common base class to all Atomic<> regions, merely
//   to be able to group them in a single set.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) AtomicRegion()
//     -----------------
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//     1) void LocalIntegrand(Integrand*)
//     ----------------------------------
//     to change the integrand
//
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     1) void Process(Stack<AtomicRegion>& Offspring)
//     -----------------------------------------------
//     attempts to compute integral and error over
//     the region. possible Offspring is returned.
//
/////////////////////////////////////////////////////////
#ifndef ATOMREG_H
#define ATOMREG_H
///////////////////////////////////////////////
#include <region.h>
#include <integran.h>
#include <stack.h>
///////////////////////////////////////////////
class AtomicRegion : public Region
  {
  public:

  AtomicRegion();
  virtual void LocalIntegrand(Integrand*)=0;
  virtual void Process(Stack<AtomicRegion>& Offspring)=0;
  virtual ~AtomicRegion();
  };
/////////////////////////////////////////////////
#endif

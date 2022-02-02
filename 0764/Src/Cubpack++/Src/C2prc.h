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
// File : C2prc.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(typedef introduced)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Parallelogram_Processor
// -------------------------------------------------
//
// BASECLASSES:
//   Processor<Parallelogram>
//
// PURPOSE:
//   A processor for two-dimensional parallelograms.
//   first a degree 13 rule is applied to the region
//   Fourth order divided differences along the axes are
//   computed. If they're very much different, the region
//   is cut in two, otherwise it is cut in four.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Parallelogram_Processor();
//     -----------------------------
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
//     1) void Process(Stack<AtomicRegion>&);
//     --------------------------------------------------
//     see class Processor
//
//     2) virtual Processor<GEOMETRY>* NewCopy() const=0
//     -------------------------------------------------
//
//     makes a new copy (using the copy constructor) and
//     returns a pointer to it.
/////////////////////////////////////////////////////////
#ifndef C2PRC_H
#define C2PRC_H
//////////////////////////////////////////////
#include <pointer.h>
#include <rule.h>
#include <samediv.h>
#include <C2.h>
#include <regproc.h>
/////////////////////////////////////////////

class Parallelogram_Processor :
     public Processor<Parallelogram>
  {
  public:
  typedef Rule<Parallelogram> RuleParallelogram;
  typedef SameShapeDivisor<Parallelogram> SameShapeDivisorParallelogram;

  Parallelogram_Processor();
  void Process(Stack<AtomicRegion>&);
  Processor<Parallelogram>* NewCopy() const;

  protected:

  static Pointer<RuleParallelogram> TheRule;
  static Pointer<SameShapeDivisorParallelogram> TheDivisor4;
  static Pointer<SameShapeDivisorParallelogram> TheDivisor2;
  unsigned int TimesCalled;
  Parallelogram_Processor* Descendant() const;
  Vector<real> Diffs;


  };
/////////////////////////////////////////////
#endif

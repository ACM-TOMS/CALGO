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
// File : s_adapt.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member declarations)
//   25 Jan 1996     V0.1f(typedef introduced, code from .c to .h)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS SimpleAdaptive
// -------------------------------------------------
//
// BASECLASSES:
//   Processor<T>
//
// PURPOSE:
//   Simple Adaptive means applying a rule once
//   (during the first call of Process())
//   and then cutting the region into subregions with
//   the same shape as the original (during the second
//   call of Process()).
//
// TEMPLATES:
//   the Geometry T to work with
//
// METHODS:
//   CONSTRUCTORS:
//     1) SimpleAdaptive( Rule<T>* ,
//                        SameShapeDivisor<T>* )
//     -------------------------------------------------
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
//     see Processor<T>
//
//     2) virtual Processor<GEOMETRY>* NewCopy() const=0
//     -------------------------------------------------
//
//     makes a new copy (using the copy constructor) and
//     returns a pointer to it.
/////////////////////////////////////////////////////////
#ifndef S_ADAPT_H
#define S_ADAPT_H
///////////////////////////////////////////
#include <rule.h>
#include <samediv.h>
#include <regproc.h>
///////////////////////////////////////////

template <class GEOMETRY>
class SimpleAdaptive : public Processor<GEOMETRY>
  {
  public:
  typedef Rule<GEOMETRY> RuleGEOMETRY;
  typedef SameShapeDivisor<GEOMETRY> SameShapeDivisorGEOMETRY;

  SimpleAdaptive(const Pointer<RuleGEOMETRY> &R,
                 const Pointer<SameShapeDivisorGEOMETRY>& D)
    :Processor<GEOMETRY>(),
    TimesCalled(0),
    Diffs(2),
    TheRule(R),
    TheDivisor(D)
  {
  }
  SimpleAdaptive(Rule<GEOMETRY>*,
                 SameShapeDivisor<GEOMETRY>*);
  void Process( Stack<AtomicRegion>& Offspring);
  Processor<GEOMETRY>* NewCopy() const;

  protected:

  SimpleAdaptive<GEOMETRY>* Descendant() const;
  unsigned int TimesCalled;
  Vector<real> Diffs;
  Pointer<RuleGEOMETRY> TheRule;
  Pointer<SameShapeDivisorGEOMETRY > TheDivisor;

  };
/////////////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <s_adapt.c>
#endif
/////////////////////////////////////////////////
#endif

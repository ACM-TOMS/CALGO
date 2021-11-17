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
// File : S2adapt.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member declarations)
//   25 Jan 1996     V0.1f(typedef introduced)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS CircleAdaptive
// -------------------------------------------------
//
// BASECLASSES:
//   Processor<Circle>
//
// PURPOSE:
//   A circle processor. First a rule is
//   applied to the circle. Afterwards the circle is
//   divided into a central circle and some
//   polar rectangles.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1)CircleAdaptive(Rule<Circle>*)
//     -------------------------------
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
//     1) void Process( Stack<AtomicRegion&>)
//       -------------------------------------------------
//     see class Processor
//
//     2) virtual Processor<GEOMETRY>* NewCopy() const=0
//     -------------------------------------------------
//
//     makes a new copy (using the copy constructor) and
//     returns a pointer to it.
/////////////////////////////////////////////////////////
#ifndef S2ADAPT_H
#define S2ADAPT_H
/////////////////////////////////////////////////
#include <regproc.h>
#include <S2.h>
#include <pointer.h>
#include <rule.h>
/////////////////////////////////////////////////
class CircleAdaptive : public Processor<Circle>
  {
  public:
  typedef Rule<Circle> RuleCircle;

  CircleAdaptive( Rule <Circle>* );
  void Process(Stack<AtomicRegion>& );
  Processor<Circle>* NewCopy() const;

  private:

  Pointer< RuleCircle > TheRule;
  CircleAdaptive* Descendant()const;
  unsigned int TimesCalled;
  unsigned int GenerationNumber;
  };
///////////////////////////////////////////////
#endif

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
// File : E2adapt.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member declarations)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PlaneAdaptive
// -------------------------------------------------
//
// BASECLASSES:
//   Processor<Plane>
//
// PURPOSE:
//   during the first call of Process() a degree 13 rule
//   is applied to the plane. during the second call 
//   the plane is divided into a CIRCLE and an OUT_CIRCLE
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PlaneAdaptive()
//     ----------------
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
//     1) void Process( Stack<AtomicRegion>&)
//     -------------------------------------------------------
//     see Processor<>
//
//     2) virtual Processor<GEOMETRY>* NewCopy() const=0
//     -------------------------------------------------
//
//     makes a new copy (using the copy constructor) and
//     returns a pointer to it.
/////////////////////////////////////////////////////////
#ifndef E2ADAPT_H
#define  E2ADAPT_H
/////////////////////////////////////////////////////
#include <regproc.h>
#include <E2.h>
#include <integran.h>
////////////////////////////////////////////////////
class PlaneAdaptive: public Processor<Plane>
  {
  public:

  PlaneAdaptive();
  void Process(Stack<AtomicRegion>&);
  Processor<Plane>* NewCopy() const;

  private:

  static void Rule(Integrand&,Plane&,real&,real&,real&);
  unsigned int TimesCalled;
  real HalfValueRadius;
  }
  ;
/////////////////////////////////////////////////
#endif

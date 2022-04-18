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
//File passbuck.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <passbuck.h>
#include <trnsfrm.h>
#include <integran.h>
#include <atomic.h>
//////////////////////////////////////////////////
template<class FROM,class TO,class VIA>
PassTheBuck<FROM,TO,VIA>::PassTheBuck(AtomicRegion* a)
  :Processor<TO>(), AR_ptr(a)
  {
  }
//////////////////////////////////////////////////
template<class FROM,class TO,class VIA>
void
PassTheBuck<FROM,TO,VIA>::Process( Stack<AtomicRegion>& Offspring)
  {
  Transformation* T = new VIA(&Geometry());
  Integrand* I = new Integrand(LocalIntegrand(),T);
  AR_ptr->LocalIntegrand(I);
  Offspring.Push(AR_ptr);
  }
/////////////////////////////////////////////////
template <class FROM,class TO, class VIA>
Processor<TO>*
PassTheBuck<FROM,TO,VIA>::NewCopy()
const
  {
  return new PassTheBuck<FROM,TO,VIA>(*this);
  }
////////////////////////////////////////////////////

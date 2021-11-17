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
//////////////////////////////////////////////
//File atomic.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////
#include <atomic.h>
/////////////////////////////////////////////
template <class GEOMETRY>
Atomic<GEOMETRY>::Atomic(GEOMETRY* g, Processor<GEOMETRY>* p)
  :AtomicRegion(),G_ptr(g),RP_ptr(p),I_ptr( new Integrand)
  {
  p->LocalAtomic(this);
  }
//////////////////////////////////////////////
template <class GEOMETRY>
void
Atomic<GEOMETRY>::Use(Processor<GEOMETRY>* p)
  {
  RP_ptr =p;
  p->LocalAtomic(this);
  }
////////////////////////////////////////////
template <class GEOMETRY>
void
Atomic<GEOMETRY>::Process(Stack<AtomicRegion>& Offspring)
  {
  RP_ptr->Process(Offspring);
  }
////////////////////////////////////////////
template <class GEOMETRY>
void
Atomic<GEOMETRY>::LocalIntegrand(Integrand* i)
  {
  I_ptr = i;
  }
//////////////////////////////////////////////
template <class GEOMETRY>
RegionInfo*
Atomic<GEOMETRY>::LocalRegionInfo()
  {
  return &(LocalInfo());
  }
///////////////////////////////////////////////

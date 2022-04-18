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
///////////////////////////////////////////////
//File regproc.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
///////////////////////////////////////////////
#include <regproc.h>
#include <atomic.h>
///////////////////////////////////////////////
template <class GEOMETRY>
Processor<GEOMETRY>::Processor()
  :ReferenceCounting()
   {
   }
///////////////////////////////////////////////
template <class GEOMETRY>
Processor<GEOMETRY>::~Processor()
  {
  }
////////////////////////////////////////////
template <class GEOMETRY>
Integrand&
Processor<GEOMETRY>::LocalIntegrand()
const
  {
  return (*(A_ptr->I_ptr));
  }
///////////////////////////////////////////
template <class GEOMETRY>
GEOMETRY&
Processor<GEOMETRY>::Geometry()
const
  {
  return (*(A_ptr->G_ptr));
  }
///////////////////////////////////////////
template <class GEOMETRY>
void
Processor<GEOMETRY>::LocalAtomic(Atomic<GEOMETRY>* a)
  {
  A_ptr = a;
  }
////////////////////////////////////////////
template <class GEOMETRY>
Atomic<GEOMETRY>&
Processor<GEOMETRY>::LocalAtomic()
const
  {
  return *A_ptr;
  }
//////////////////////////////////////////////
template <class GEOMETRY>
real&
Processor<GEOMETRY>::Integral()
  {
  return A_ptr->LocalRegionInfo()->Integral();
  }
//////////////////////////////////////////////////
template <class GEOMETRY>
real&
Processor<GEOMETRY>::AbsoluteError()
  {
  return A_ptr->LocalRegionInfo()->AbsoluteError();
  }
/////////////////////////////////////////////////
template <class GEOMETRY>
RegionInfo&
Processor<GEOMETRY>::LocalRegionInfo()
const
  {
  return *(A_ptr->LocalRegionInfo());
  }
////////////////////////////////////////////////

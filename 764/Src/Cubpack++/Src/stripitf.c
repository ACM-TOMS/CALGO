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
////////////////////////////////////////////////
//File stripitf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////
#include <stripitf.h>
#include <passbuck.h>
#include <E2interf.h>
#include <E2tostrp.h>
//////////////////////////////////////////////
INFINITE_STRIP::INFINITE_STRIP(const Point& a,const Point& b)
  :USERINTERFACE<InfiniteStrip>()
  {
  PLANE P;
  StoreAtomic(new InfiniteStrip(a,b),
      new PassTheBuck<Plane ,InfiniteStrip,
               E2toIS>((AtomicRegion*)P));
  }
///////////////////////////////////////////////

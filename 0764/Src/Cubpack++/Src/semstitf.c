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
//////////////////////////////////////////////
#include <semstitf.h>
#include <passbuck.h>
#include <stripitf.h>
#include <sttosmst.h>
//////////////////////////////////////////////
SEMI_INFINITE_STRIP::SEMI_INFINITE_STRIP(const Point& a,const Point& b)
  :USERINTERFACE<SemiInfiniteStrip>()
  {
  Point or(0,0),one(1,0);
  INFINITE_STRIP I(or,one);
  StoreAtomic(new SemiInfiniteStrip(a,b),
      new PassTheBuck<InfiniteStrip,SemiInfiniteStrip,
               IStoSIS>((AtomicRegion*)I));
  }
///////////////////////////////////////////////

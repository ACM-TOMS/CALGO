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
//File gritf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////
#include <gritf.h>
#include <passbuck.h>
#include <C2.h>
#include <C2interf.h>
#include <C2togr.h>
#include <point.h>
//////////////////////////////////////////////
GENERALIZED_RECTANGLE::GENERALIZED_RECTANGLE(Function f,const Point& a,const Point& b)
  :USERINTERFACE<GeneralizedRectangle>()
  {
  Point A(0,0),B(1,0),C(0,1);
  PARALLELOGRAM R(A,B,C);
  StoreAtomic(new GeneralizedRectangle(f,a,b),
      new PassTheBuck<Parallelogram,GeneralizedRectangle,
               C2toGR>((AtomicRegion*)R));
  }
///////////////////////////////////////////////

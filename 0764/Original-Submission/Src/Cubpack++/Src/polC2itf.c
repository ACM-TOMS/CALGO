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
//File polC2itf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////
#include <polC2itf.h>
#include <polC2prc.h>
//////////////////////////////////////////////
POLAR_RECTANGLE::POLAR_RECTANGLE
  (const Point& A, const Point& B,const Point& C)
  {
  StoreAtomic(new PolarRectangle(A,B,C),new PolarRectangle_Processor);
  }
//////////////////////////////////////////////
POLAR_RECTANGLE::POLAR_RECTANGLE
  (const Point& O, real r1, real r2, real t1, real t2)
  {
  StoreAtomic(new PolarRectangle(O,r1,r2,t1,t2),new PolarRectangle_Processor);
  }
//////////////////////////////////////////////

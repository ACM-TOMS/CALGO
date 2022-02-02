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
///////////
//File E2interf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////
#include <E2interf.h>
#include <E2adapt.h>
///////////////////////
PLANE::PLANE()
  {
  StoreAtomic(new Plane,new PlaneAdaptive);
  }
//////////////////////////////////////////////
PLANE::PLANE(const  Point& Center)
  {
  StoreAtomic(new Plane(Center),new PlaneAdaptive);
  }
//////////////////////////////////////////////
PLANE::PLANE(const Point& Center, real ScaleX, real ScaleY)
  {
  StoreAtomic(new Plane(Center,ScaleX, ScaleY),new PlaneAdaptive);
  }
///////////////////////////////////////////////

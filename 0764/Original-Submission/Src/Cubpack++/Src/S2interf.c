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
//File S2interf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <S2interf.h>
#include <S2adapt.h>
#include <S2rule13.h>
/////////////////////////////////////////////////////////
CIRCLE::CIRCLE(const Point& c,const Point& b)
  {
        StoreAtomic(new Circle(c,b),new CircleAdaptive(new Circle_Rule13));
  }
/////////////////////////////////////////////////////////
CIRCLE::CIRCLE(const Point& c, real Radius)
  {
        StoreAtomic(new Circle(c,Radius),new CircleAdaptive(new Circle_Rule13));
  }
/////////////////////////////////////////////////////////

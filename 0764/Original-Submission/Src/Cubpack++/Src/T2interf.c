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
//File T2interf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(typedef introduced)
////////////////////////////////////////////////
#include <T2interf.h>
#include <s_adapt.h>
#include <T2rule13.h>
#include <T2dv4.h>
///////////////////////////////////////////////
typedef Rule<Triangle> RuleTriangle;
typedef SameShapeDivisor<Triangle> SameShapeDivisorTriangle;
TRIANGLE::TRIANGLE(const Point& p1,const  Point&  p2,
                   const Point& p3)
  { Pointer<RuleTriangle> R13 (new Triangle_Rule13);
    Pointer<SameShapeDivisorTriangle> D4 (new Triangle_Divide4);
    StoreAtomic(new Triangle(p1,p2,p3),
       new SimpleAdaptive<Triangle>(R13,D4) );
  }
///////////////////////////////////////////////

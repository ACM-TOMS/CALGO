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
//File E2secitf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////
#include <E2secitf.h>
#include <E2secprc.h>
//////////////////////////////////////////////
PLANE_SECTOR::PLANE_SECTOR
  (const Point& A, const Point& B,const Point& C)
  {
  StoreAtomic(new PlaneSector(A,B,C),new PlaneSector_Processor);
  }
//////////////////////////////////////////////
PLANE_SECTOR::PLANE_SECTOR
  (const Point& O, real r  , real theta1, real theta2)
  {
  StoreAtomic(new PlaneSector(O,r,theta1,theta2),
              new PlaneSector_Processor);
  }
///////////////////////////////////////////////////

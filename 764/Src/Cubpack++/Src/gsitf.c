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
//File gsitf.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////
#include <gsitf.h>
#include <point.h>
#include <gsprc.h>
//////////////////////////////////////////////
GENERALIZED_SECTOR::GENERALIZED_SECTOR(real (*F)(real),
     real a, real b, const Point& Center)
  :USERINTERFACE<GeneralizedSector>()
  {
  StoreAtomic(new GeneralizedSector(F,a,b,Center),
  new GeneralizedSector_Processor);
  }
///////////////////////////////////////////////

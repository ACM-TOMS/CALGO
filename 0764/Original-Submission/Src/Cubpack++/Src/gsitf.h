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
// File : gritf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS GENERALIZED_SECTOR
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<GeneralizedSector>
//
// PURPOSE:
//   End-user interface to Generalized Rectangles
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) GENERALIZED_SECTOR(real (*F) real,
//        real Alpha, real Beta, const Point& Center)
//     ---------------------------------------------------
//      the region is bounded by 2 straight lines going
//      through Center and having an angle of alpha and
//      beta respectively with the positive X-axis. a
//      third boundary line is provided by f, which denotes
//      the distance from the Center to the boundary for
//      all angles between alpha and beta.
//
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef GSITF_H
#define GSITF_H
////////////////////////////////////////////
#include <userint.h>
#include <gs.h>
////////////////////////////////////////////
class GENERALIZED_SECTOR : public
   USERINTERFACE<GeneralizedSector>
   {
   public:

   GENERALIZED_SECTOR( real(*)(real), real,real,const Point&);
   };
///////////////////////////////////////

#endif

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
// File : semstitf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS SEMI_INFINITE_STRIP
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<SemiInfiniteStrip>
//
// PURPOSE:
//   End-user interface to semi-infinite strips
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) SEMI_INFINITE_STRIP(
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     the region will be to the left of AB.
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
#ifndef SEMSTITF_H
#define SEMSTITF_H
////////////////////////////////////////////
#include <userint.h>
#include <semistrp.h>
////////////////////////////////////////////
class SEMI_INFINITE_STRIP : public
   USERINTERFACE<SemiInfiniteStrip>
   {
   public:

   SEMI_INFINITE_STRIP(const Point&,const Point&);
   };
////////////////////////////////////////////
#endif

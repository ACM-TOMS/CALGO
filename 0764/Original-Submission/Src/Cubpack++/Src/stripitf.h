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
// File : stripitf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS INFINITE_STRIP
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<InfiniteStrip>
//
// PURPOSE:
//   End user interface to  infinite strips
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) INFINITE_STRIP(
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     AB will be the diameter of the strip
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
#ifndef STRIPITF_H
#define STRIPITF_H
////////////////////////////////////////////
#include <userint.h>
#include <strip.h>
////////////////////////////////////////////
class INFINITE_STRIP : public
   USERINTERFACE<InfiniteStrip>
   {
   public:

   INFINITE_STRIP(const Point&,const Point&);
   };
////////////////////////////////////////////
#endif

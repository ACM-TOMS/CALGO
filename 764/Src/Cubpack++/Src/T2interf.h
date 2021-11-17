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
// File :  T2interf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS TRIANGLE
// ----------------------------
//
// BASECLASSES:
//   USERINTERFACE<Triangle>
//
// PURPOSE:
//   an end-user interface to  Triangle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1)TRIANGLE(const Point& , const Point&, const Point&)
//     ---------------------------------------------------
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
#ifndef T2INTERF_H
#define T2INTERF_H
///////////////////////////////////////////////////
#include <userint.h>
#include <T2.h>
///////////////////////////////////////////////////
class TRIANGLE : public USERINTERFACE<Triangle>
  {
  public:

  TRIANGLE(const Point&,const Point&,const Point&);
  };

//////////////////////////////////////////////////
#endif

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
// DEFINITION OF CLASS GENERALIZED_RECTANGLE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<GeneralizedRectangle>
//
// PURPOSE:
//   End-user interface to Generalized Rectangles
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) GENERALIZED_RECTANGLE(Function f,
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     constructs a generalized rectangle with two points
//     and a function;  f(P) is the length of the
//     perpendicular to
//     the boundary from a point P on the line AB.
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
#ifndef GRITF_H
#define GRITF_H
////////////////////////////////////////////
#include <userint.h>
#include <gr.h>
////////////////////////////////////////////
class GENERALIZED_RECTANGLE : public
   USERINTERFACE<GeneralizedRectangle>
   {
   public:

   GENERALIZED_RECTANGLE(Function, const Point&,const Point&);
   };
////////////////////////////////////////////
#endif

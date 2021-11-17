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
// File : C2interf.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS PARALLELOGRAM
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<Parallelogram>
//
// PURPOSE:
//   End-user interface to a parallelogram.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) PARALLELOGRAM(const Point& A,
//                      const Point& B,
//                      const Point& C)
//     ------------------------------
//     constructs a parallelogram with 4 vertices
//     A, B, C, and D=B+C-A
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
// DEFINITION OF CLASS RECTANGLE
// -------------------------------------------------
//
// BASECLASSES:
//   USERINTERFACE<Parallelogram>
//
// PURPOSE:
//   End-user interface to a parallelogram.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) RECTANGLE(const Point& A,
//                      const Point& B,
//                      const Point& C)
//     ------------------------------
//     constructs a parallelogram with 4 vertices
//     A, B, C, and D=B+C-A
//     AB should be orthogonal to BC
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
#ifndef C2INTERF_H
#define C2INTERF_H
///////////////////////////////////////////////////
#include <userint.h>
#include <C2.h>
///////////////////////////////////////////////////
class PARALLELOGRAM : public USERINTERFACE<Parallelogram>
  {
  public:

  PARALLELOGRAM(const Point&,const Point&,const Point&);
  };

//////////////////////////////////////////////////
class RECTANGLE : public USERINTERFACE<Parallelogram>
  {
  public:

  RECTANGLE(const Point&,const Point&,const Point&);
  };

//////////////////////////////////////////////////
#endif

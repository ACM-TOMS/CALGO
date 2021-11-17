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
// File : translat.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Translation
// -------------------------------
//
// BASECLASSES:
//   Transformation
//
//
// PURPOSE:
//   implements translations.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Translation(const Point& Offset)
//     -----------------------------------
//     constructs the translation by specifying the
//     Offset that has to be added by T()
//
//   SELECTORS:
//     None
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) void Transform (real& w, Point& p)
//     ---------------------------------------
//     Leaves w unchanged and replaces p by p+Offset.
//
///////////////////////////////////////////////////////////

#ifndef TRANSLAT_H
#define TRANSLAT_H
////////////////////////////////////////////////////////

#include <trnsfrm.h>
#include <point.h>

///////////////////////////////////////////////////////
class Translation :public Transformation
  {
  public:

  Translation( const Point& Offset);
  //Translation::T(p) ADDS Offset to p
  void Transform(real& w, Point& p);

  private:

  Point TheOffset;
  };
////////////////////////////////////////////////////////////
#endif

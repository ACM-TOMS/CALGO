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
// File : div.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Divisor
// ---------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
//
// PURPOSE:
//   Divisor is a pure virtual base class for implementing
//   dividing-procedures. Information on the number of
//   parts the divisor produces, can be retrieved.
//
// TEMPLATES:
//   The type of region T which has to be divided is
//   templated. Currently it is assumed that the parts
//   after division also are of type T. T should be
//   derived from class Region.
//
// METHODS:
//   CONSTRUCTORS:
//     1) Divisor()
//     ------------
//
//   SELECTORS:
//     1) virtual int NumberOfParts() const=0
//     ---------------------------------------------------
//     returns the number of parts the divisor produces.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////

#ifndef DIV_H
#define DIV_H
/////////////////////////////////////////
#include <refcount.h>

//////////////////////////////////////////
template <class GEOMETRY>
class Divisor: public ReferenceCounting
  {

  public:

  Divisor();
  virtual int NumberOfParts() const =0;
  virtual ~Divisor();


  };
///////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <div.c>
#endif
//////////////////////////////////////////

#endif

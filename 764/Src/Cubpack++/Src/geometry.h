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
// File : geometry.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Geometry
// -------------------------------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
// PURPOSE:
//   common base for all geometries
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Geometry( unsigned int Dim)
//     ------------------------------
//     the Dimension must be provided.
//
//   SELECTORS:
//     1) unsigned int Dimension() const
//     ---------------------------------
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
#ifndef GEOMETRY_H
#define GEOMETRY_H
/////////////////////////////////////////////////////////

#include <refcount.h>
//////////////////////////////////////////////////

class Geometry : public ReferenceCounting
  {
  public:

  Geometry (unsigned int Dim);
  unsigned int Dimension()const;

  private:

  unsigned int TheDimension;
  };
//////////////////////////////////////////////////
#endif

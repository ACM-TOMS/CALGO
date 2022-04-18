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
// File : T2.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Triangle
// -------------------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements the geometry of a simple triangle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Triangle (const Point&,const Point&,const Point&)
//     ----------------------------------------------------
//     constructs a 2D Triangle given its 3 vertices
//
//   SELECTORS:
//     1) const Point& Vertex(int i) const
//     -----------------------------------
//     returns the i-th vertex, 0<=i<3
//
//     2) real Volume() const
//     -----------------------
//
//     3) Processor<Triangle>* DefaultProcessor() const
//     ------------------------------------------------
//
//   MODIFIERS:
//     1) void Volume(real)
//     --------------------
//     informs the triangle about its volume, to avoid
//     future computations of it
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef T2_H
#define T2_H
//////////////////////////////////////////////////
#include <geometry.h>
#include <real.h>
#include <point.h>
#include <vector.h>
#include <boolean.h>
#include <regproc.h>

//////////////////////////////////////////////////

class Triangle : public Geometry
  {
  public :      

  Triangle(const Point&,const Point&,const Point&);
  const Point& Vertex (int) const;
  real Volume() const;
  void Volume(real);

  private:

  Vector<Point> Vertices;
  real TheVolume;
  Boolean TheVolumeKnown;
  void ComputeVolume();
  };

////////////////////////////////////////////////////

#endif

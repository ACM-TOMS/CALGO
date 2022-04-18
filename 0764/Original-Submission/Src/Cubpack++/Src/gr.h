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
// File : gr.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member declarations)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Generalized Rectangle
// -----------------------------------------
//
// BASECLASSES:
//   Geometry
//
// PURPOSE:
//   implements a generalized rectangle
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) GeneralizedRectangle(Function f,
//        const Point& A, const Point& B)
//     ---------------------------------------------------
//     constructs a generalized rectangle with two points
//     and a function;  f(P) is the length of the
//     perpendicular to
//     the boundary from a point P on the line AB.
//
//   SELECTORS:
//     1) const Point& A() const
//     -------------------------
//     returns  the point A (see constructor)
//
//     2) const Point& B() const
//     -------------------------
//     returns the point B (see constructor)
//
//     3) real Boundary(const Point& P)const
//     -------------------------------------
//     evaluates the boundary function at P(see constructor)
//
//     4) Processor<GeneralizedRectangle>*
//             DefaultProcessor() const
//     -----------------------------------
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
//
///////////////////////////////////////////////////////////
#ifndef gr2_H
#define gr2_H
//////////////////////////////////////////////////

#include <point.h>
#include <geometry.h>
#include <function.h>
#include <regproc.h>
//////////////////////////////////////////////////
class GeneralizedRectangle : public  Geometry
  {
  public:

  GeneralizedRectangle (Function,
       const Point&,const Point&);
  const Point& A() const;
  const Point& B() const;
  real Boundary(const Point&) const;

  private:

  const Function TheBoundary;
  Point TheA;
  Point TheB;
  };
//////////////////////////////////////////////////
#endif

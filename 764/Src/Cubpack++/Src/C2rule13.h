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
// File : C2rule13.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Parallelogram_Rule13
// ------------------------------------
//
// BASECLASSES:
//   Rule<Parallelogram>
//
//
// PURPOSE:
//   A degree 13 rule for a rectangle.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     None
//   SELECTORS:
//     1) int Degree() const
//     ---------------------
//     returns 13
//
//     2) int NumberOfPoints(int ) const
//     --------------------------------------
//     returns 37.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) ApplyWithDiffs(Integrand & I,
//              Parallelogram& H,real& Result,real&  Error
//              Vector<real>& D)
//     -------------------------------------------------
//     Input parameters:
//      I: Integrand
//      H: rectangle to be integrated. H.Dimension()
//      should be 2.
//     Output parameters:
//      Result: approximation to the integral
//      Error: absolute error estimation
//      D: Fourth order divided differences in each direction
//         D[0] refers to {Vertex(0),Vertex(1)};
//         D[1] refers to {Vertex(0),Vertex(2)};
//
//
///////////////////////////////////////////////////////////

#ifndef C2RULE13_H
#define C2RULE13_H
/////////////////////////////////////////
#include <C2.h>
#include <rule.h>
//////////////////////////////////////////

class Parallelogram_Rule13 : public Rule<Parallelogram>
  {
  public:

  Parallelogram_Rule13();
  void ApplyWithDiffs(Integrand&,Parallelogram&,real& Result,
    real& Error,Vector<real>&);
  int NumberOfPoints() const {return 37;};
  int Degree () const {return 13;};
  };
////////////////////////////////////////////
#endif

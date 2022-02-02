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
// File : T2rule13.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Triangle_Rule13
// ------------------------------------
//
// BASECLASSES:
//   Rule<Triangle>
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
//     1)Triangle_Rule13()
//     -------------------
//
//   SELECTORS:
//     1) int Degree() const
//     ---------------------
//     returns 13
//
//     2) int NumberOfPoints() const
//     --------------------------------------
//     returns 37.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) Apply(Integrand& I,Triangle& H,real& Result,real&
//     Error)
//     -------------------------------------------------
//     Input parameters:
//      I: the integrand
//      H: Triangle to be integrated.
//     Output parameters:
//      Result: approximation to the integral
//      Error: absolute error estimation
//
///////////////////////////////////////////////////////////

#ifndef T2_RULE13
#define T2_RULE13
///////////////////////////////////////////

#include <rule.h>
#include <T2.h>

///////////////////////////////////////////
class Triangle_Rule13 :public Rule<Triangle>
  {
  public:

  Triangle_Rule13();
  void Apply(Integrand&,Triangle& ,real&,real& );
  int Degree() const {return 13;};
  int NumberOfPoints () const {return 37;};

  };
///////////////////////////////////////////

#endif

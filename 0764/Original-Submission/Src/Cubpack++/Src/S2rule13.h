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
// File : S2rule13.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Circle_Rule13
// ------------------------------------
//
// BASECLASSES:
//   Rule<Circle>
//
//
// PURPOSE:
//   A degree 13 rule for a circle.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1)Circle_Rule13()
//     ------------------
//
//   SELECTORS:
//     1) int Degree() const
//     ---------------------
//     returns 13
//
//     2) int NumberOfPoints() const
//     --------------------------------------
//     returns 36.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) Apply(Integrand& I,Circle& H,
//              real& Result,real& Error)
//     -------------------------------------------------
//     Inputparameters:
//      I: the integrand
//      H: circle to be integrated. H.Dimension()
//      should be 2.
//     Outputparameters:
//      Result: approximation to the integral
//      Error: absolute error estimation
//
///////////////////////////////////////////////////////////

#ifndef S2RULE13_H
#define S2RULE13_H
/////////////////////////////////////////
#include <S2.h>
#include <rule.h>
//////////////////////////////////////////

class Circle_Rule13 : public Rule<Circle>
  {
  public:

  Circle_Rule13();
  void Apply(Integrand&,Circle&,real& Result ,real& Error);
  int NumberOfPoints() const {return 36;};
  int Degree () const {return 13;};
  };
////////////////////////////////////////////
#endif

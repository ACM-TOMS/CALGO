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
// File : rule.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Rule
// ------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
//
// PURPOSE:
//   provides a pure virtual base class as a common
//   interface to all rules.
//
// TEMPLATES:
//   T: type of region on which the rule has to be
//     applied. T should be derived from
//     class Region .
//
// METHODS:
//   CONSTRUCTORS:
//     None
//   SELECTORS:
//     1) int Degree() const
//     ---------------------
//     returns the degree of the formula
//
//     2) int NumberOfPoints() const
//     ------------------------------------------
//     returns the number of points the rule uses.
//     this might depend on Dimension
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     None
//   SPECIAL:
//     1) Apply(Integrand& I,T& H,real& Result,real&
//     Error)
//     -------------------------------------------------
//     Input parameters:
//      I: the integrand
//      H: Region to be integrated. H.Dimension()
//      should be 2.
//     Output parameters:
//      Result: approximation to the integral
//      Error: absolute error estimation
//
//     2) ApplyWithDiffs(Integrand& I,T& H,real& Result,real&
//     Error,Vector<unsigned int>& D)
//     -------------------------------------------------
//     Input parameters:
//      I: the integrand
//      H: Region to be integrated. H.Dimension()
//      should be 2.
//     Output parameters:
//      Result: approximation to the integral
//      Error: absolute error estimation
//      D: order of indices regarding a higher order difference.
///////////////////////////////////////////////////////////

#ifndef RULE_H
#define RULE_H

/////////////////////////////////////////
#include <refcount.h>
#include <integran.h>
#include <real.h>
//////////////////////////////////////////

template <class GEOMETRY>
class Rule : public ReferenceCounting
  {
  public:

  Rule();
  virtual void ApplyWithDiffs(Integrand&,GEOMETRY&,real&,real&,Vector<real>&) ;
  virtual void Apply(Integrand&,GEOMETRY&,real&,real&);
  virtual int Degree() const =0;
  virtual int NumberOfPoints () const =0;
  virtual ~Rule();

  };
//////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <rule.c>
#endif
//////////////////////////////////////////

#endif

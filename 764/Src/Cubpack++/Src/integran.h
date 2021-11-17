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
// File : integran.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//    8 Sep 1994     V0.1a(operator added)
//   25 Jan 1996     V0.1f(typedef introduced)
//   28 Mar 1996     V0.1h(long instead of int)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Integrand
// -----------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
// PURPOSE:
//   implements an integrand class. Apart from evaluating
//   integrands, it takes care of counting the number of
//   evaluations and it permits the integrand to be
//   transformed.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Integrand( Function f)
//     ---------------------------------
//     constructs the Integrand. The function to be
//     integrated is f.
//
//     2) Integrand( const Integrand& I)
//     ---------------------------------
//     copy constructor
//
//     3) Integrand( const Integrand& I,
//                const Transformation& T)
//     --------------------------------------
//     constructs a new Integrand by applying T to I
//
//     4) Integrand()
//     --------------
//     constructs an Integrand without the function
//     to be integrated being known. an evaluation of
//     this Integrand results in an error.
//
//   SELECTORS:
//     1)static long NumberOfEvaluations()
//     ----------------------------------
//     returns the total number of function evaluations
//     since the program was started.
//     See also class EvaluationCounter.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     1) real operator() (const Point& p)
//     -------------------------------------
//     evaluates the Integrand at Point p.
//
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////


#ifndef INTEGRAN_H
#define INTEGRAN_H
#include <point.h>
#include <trnsfrm.h>
#include <vstack.h>
#include <real.h>
#include <function.h>
#include <refcount.h>
#include <pointer.h>
///////////////////////////////////////////////////////
class Integrand : public ReferenceCounting
  {
  public:
  typedef Pointer< Transformation > PointerTransformation;

  Integrand();
  Integrand(Function);
  Integrand(const Integrand&,Transformation* );
  Integrand(const Integrand&);
  real operator()(const Point&);
  Boolean operator==(const Integrand&) const;
  static long    NumberOfEvaluations();


  private:

  static long    Number;
  Function TheFunction;
  VectorStack< PointerTransformation >  AppliedTransformations;
  };

#endif

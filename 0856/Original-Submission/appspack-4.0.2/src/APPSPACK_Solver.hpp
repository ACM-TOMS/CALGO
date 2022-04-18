// $Id: APPSPACK_Solver.hpp,v 1.14 2004/11/24 19:26:39 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Solver.hpp,v $ 

//@HEADER
// ************************************************************************
// 
//          APPSPACK: Asynchronous Parallel Pattern Search
//                 Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.                                                                           .
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*!
  \file APPSPACK_Solver.hpp
  \brief Class definition for APPSPACK::Solver
*/

#ifndef APPSPACK_SOLVER_HPP
#define APPSPACK_SOLVER_HPP

#include "APPSPACK_Point.hpp"
#include "APPSPACK_Constraints_Interface.hpp"
#include "APPSPACK_Conveyor.hpp"
#include "APPSPACK_Directions.hpp"
#include "APPSPACK_List.hpp"
#include "APPSPACK_Parameter_List.hpp"
#include "APPSPACK_Print.hpp"

namespace APPSPACK
{

//! The solver itself

class Solver
{
  
public:
  
  //! State of the solver
  enum State {
    Continue,			//!< Default state
    StepConverged,		//!< Step length converged
    FunctionConverged,		//!< Function tolerance converged
    EvaluationsExhausted        //!< Number of function evaluations exhausted
  };

  //! Constructor 
  Solver(const Parameter::List& params_in, 
	 Executor::Interface& executor_in, 
	 const Constraints::Interface& constraints_in);

  /*! \brief Called by the Solver constructor to parse the parameter list and
    create the initial best point. */
  static Point* initializeBestPointPtr(Parameter::List& params, 
				       const Constraints::Interface& constraints);


  //! Destructor 
  ~Solver();


  // --------------------------------------------------------------------------------

  //@{ \name Accessors (for after solve() is called)
  
  //! Return the x-vector corresponding to the best point
  const Vector& getBestX() const;

  //! Retun true is there is a finite function value associated with the best point
  bool isBestF() const;

  //! Return the finite function value (if any) associated with the best point
  double getBestF() const;

  //@}

  // --------------------------------------------------------------------------------

  //@{ \name Manipulators

  //! Find the minimum of the function that was set up in the constructor
  APPSPACK::Solver::State solve();

  //! Do a single APPS iteration 
  APPSPACK::Solver::State iterate();

  //@}

private:

  /*!  

    Process a new best point - Delete the old point (if any) and
    replace it with the new.

    If the function tolerance convergence test is being employed, we
    check for convergence at this point, and return if convergence is
    detected.

    Otherwise, we generate new search directions, re-initialize the
    #step array, and reset various bookkeeping variables (#trueStep
    and #tag).

    All the steps start the same, and they are calculated as follows:
    \f[
    {\rm step}_{\rm new} = \max \{ {\rm step}_{\rm best}, 2 * {\rm stepTolerance} \}
    \f]

   */
  void processNewBestPoint(Point* newBestPointPtr = NULL);

  /*!
    
    Generate trial points for any directions that are not converged
    and do not already have an associated trial point.

    Note that a direction is considered converged if its corresponding
    step length is strictly smaller than #stepTolerance.

  */
  void generateTrialPoints();


  /*!  

    Create a new trial point corresponding to the direction with the
    given index, idx. Let \f$d_i\f$ (see #directions) denote the given
    direction and \f$\Delta_i\f$ (see #steps) the corresponding step
    length. Then we calculate the new trial point as:
    \f[
    x = x_{\rm best} + \Delta_i d_i
    \f]

    If this point is infeasible, then we reset \f$\Delta_i\f$ to the
    longest possible step according to the following formula. Let
    \f$a^T x = b\f$ denote a linear constraint that is violated. 
    Then
    \f[
    \tilde \Delta_i = \frac{b - a^T x_{\rm best}}{a^T d_i}
    \f]
    
    The value of \f$\tilde \Delta_i\f$ is stored in #trueSteps, and
    the trial point is then calculated as:
    \f[
    x = x_{\rm best} + \tilde \Delta_i d_i
    \f]

   */
  void createTrialPoint(int idx);

  /*!

    Process a list of trial points that has been returned from the
    UberEvaluator.

    First, check to see if the best of all these points is better than
    the current best point. If so, replace that best point (see
    processNewBestPoint()) and throw away the rest of the list.

    Otherwise, process each point and delete it.

    Finally, check if all the directions have converged and check if
    we have exhausted the maximum number of function evaluations.

   */
  void processEvaluatedTrialPoints();

  //@{ \name Print Functions

  void printInitializationInformation() const;
  void printBestPoint(const string label = "") const;
  
  //@}

private:

  //! Constraints
  const Constraints::Interface& constraints;

  //! Parameters
  Parameter::List params;

  //! Print
  Print print;

  //! Pointer to the best trial point thus far
  Point* bestPointPtr;

  //! The search directions 
  Directions directions;

  //! Trial Point Evaluator
  Conveyor conveyor;

  //! List of trial points to be processed in some way
  List exchangeList;

  //! The state of the solver
  State state;

  //@{ \name Stopping Criteria

  //! Enforce function value tolerance
  bool isFunctionTolerance;

  //! Function value
  Value functionTolerance;

  //! Enforce function evaluation budget
  bool isMaxEvaluations;

  //! Function evaluation budget
  int maxEvaluations;

  //@}

  //! Tolerance for saying whether or not we are on a boundary
  double boundsTolerance;

  //! A temporary vector
  Vector tmpVector;

  //! Machine epsilon
  double epsMach;

};

}

//! Printing an APPSPACK::Solver::State
ostream& operator<<(ostream& stream, APPSPACK::Solver::State state);

#endif

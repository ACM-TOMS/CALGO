// $Id: APPSPACK_Directions.hpp,v 1.11 2004/11/23 22:26:01 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Directions.hpp,v $ 

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
  \file APPSPACK_Directions.hpp
  \brief Class description for APPSPACK::Directions
*/

#ifndef APPSPACK_DIRECTIONS_HPP
#define APPSPACK_DIRECTIONS_HPP

#include "APPSPACK_Point.hpp"
#include "APPSPACK_Constraints_Interface.hpp"
#include "APPSPACK_Parameter_List.hpp"

namespace APPSPACK
{

//! The search directions and associated information
/*!  
  The job of this object is to generate an appropriate set of search
  directions for given point, and then to track and update the
  associated step lengths.
  
  For each direction, we track the desired step length, the true step
  length (if the desired step length takes us outside the feasible
  region we instead substitute a step length that keeps us feasible),
  and the associated trial point tag. 

  The true step and tag are actually calculated by
  Solver::generateTrialPoints and just stored in this object for
  information purposes.

  <b>Parameters</b>

  These parameters are are stored in the Parameter::List that is
  passed to the Solver constructor. See \ref pageParameters for full details on
  these parameters.

  <ul>
  <li>"Step Tolerance"
  <li>"Minimum Step"
  <li>"Contraction Factor"
  </ul>

*/
class Directions
{

public:

  //! Constructor 
  /*!  
    Reads the "Step Tolerance", "Minimum Step", and "Contraction
    Factor" from the parameter list.
  */
  Directions(Parameter::List& params, const Constraints::Interface& constraints_in);
  
  //! Destructor 
  ~Directions();
  
  //@{ \name Accessors

  //! Returns the indices of directions that are ready for new trial points.
  /*!  
    A direction is ready for a new trial point if its associated
    step length is greater than or equal to the convergence tolerance
    and it doesn't already have an active trial point.
  */
  const vector<int>& getDirectionIndices() const;

  //! Returns the i-th search direction
  const Vector& getDirection(int i) const;

  //! Returns the i-th step
  double getStep(int i) const;

  //! Returns true if \e every step length is less than the step convergence tolerance 
  /*!
    The steps are stored in the #step vector.
    The convergence tolerance is stored in #stepTolerance.
  */
  bool isStepConverged() const;
  
  //@}


  //@{ \name Manipulators

  //! Computes a new set of search directions based on a new best point. 
  /*!
    Also resets the step information.
    The initial step (which is the same for all directions) is calculated as
    \f[
    \Delta_{\rm new} = \max \{ \Delta_{\rm old}, \Delta_{\min} \}.
    \f]
    Here, \f$\Delta_{\min}\f$ is #minStep.
  */
  void computeNewDirections(const Point& newPoint);

  //! Set the true step and tag for the trial point corresponding to direction i.
  void setTrueStepAndTag(int i, double trueStep_in, int tag_in);

  //! Reduce the step corresponding to direction i.
  /*! 
    The new step is calculated as 
    \f[
    \Delta_{\rm new} = 0.5 \Delta_{\rm old}.
    \f]
    Also resets the corresponding trueStep and tag. 
  */
  void reduceStep(int i);

  //@}

  //@{ \name Printing

  //! %Print the directions, preceeded by the given label
  void print(const string label) const;

  //@}

private:

  //! The constraints object
  const Constraints::Interface& constraints;

  //! The number of dimensions
  const int nDimensions;

  //! The vector of all zeros
  const Vector zero;

  //! The step tolerance used to test convergence, etc.
  const double stepTolerance;

  //! Minimum size of the step after a new minimum is found
  const double minStep;

  //! Constraction parameter
  const double theta;

  //! The current number of search directions
  int nDirections;

  //! The search directions
  vector< Vector > direction;

  //! The steps associated with each search direction
  Vector step;

  //! The actual step associated with each search direction
  /*! The step represents the longest <em>feasible</em> step that is
    less than or equal to the corresponding delta.
  */
  Vector trueStep;

  //! The trial point tags corresponding to the search directions
  /*! A value of -1 means that there is no trial point associated with
    the given direction.
  */
  vector<int> tag;

  //! A temporary vector
  Vector tmpVector;

  //! A temporary integer vector
  mutable vector<int> idxVector;

};

}

#endif

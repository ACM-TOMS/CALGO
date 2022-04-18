// $Id: APPSPACK_Directions.cpp,v 1.13 2004/11/23 22:26:01 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Directions.cpp,v $ 

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
  \file APPSPACK_Directions.cpp
  \brief Implemtation of APPSPACK::Directions
*/

#include "APPSPACK_Directions.hpp"
#include "APPSPACK_Print.hpp"
#include "APPSPACK_Utils.hpp"

APPSPACK::Directions::Directions(Parameter::List& params, const Constraints::Interface& constraints_in) :
  constraints(constraints_in),
  nDimensions(constraints_in.getScaling().size()),
  zero(APPSPACK::createZeroVector(nDimensions)),
  stepTolerance(params.getParameter("Step Tolerance", 0.01)),
  minStep(params.getParameter("Minimum Step", 2 * stepTolerance)),
  theta(params.getParameter("Contraction Factor", 0.5))
{
  // Check parameters
  if (stepTolerance <= 0)
  {
    cout << "APPSPACK::Directions::Directions - Error: \"Step Tolerance\" cannot be negative." << endl;
    throw "APPSACK Error";
  }

  if (minStep <= stepTolerance)
  {
    cout << "APPSPACK::Directions::Directions - Error: \"Minimum Step\" must be greater than \"Step Tolerance\"." << endl;
    throw "APPSACK Error";
  }

  if ((theta <= 0) || (theta >= 1))
  {
    cout << "APPSPACK::Directions::Directions - Error: \"Contraction Factor\" must be strictly between zero and one." << endl;
    throw "APPSACK Error";
  }

  nDirections = 0;
  direction.reserve(2*nDimensions);
  step.reserve(2*nDimensions);
  trueStep.reserve(2*nDimensions);
  tag.reserve(2*nDimensions);

}


APPSPACK::Directions::~Directions()
{
}

const APPSPACK::Vector& APPSPACK::Directions::getDirection(int i) const
{
  return direction[i];
}

double APPSPACK::Directions::getStep(int i) const
{
  return step[i];
}

const vector<int>& APPSPACK::Directions::getDirectionIndices() const
{
  idxVector.resize(0);
  for (int i = 0; i < nDirections; i ++)
    if ((step[i] >= stepTolerance) && (tag[i] == -1))
      idxVector.push_back(i);
  return idxVector;
}

void APPSPACK::Directions::computeNewDirections(const Point& newPoint)
{
  // Compute the appropriate directions (only works for bounds right now)

  const Vector& scaling = constraints.getScaling();
  const Vector& lower = constraints.getLower();
  const Vector& upper = constraints.getUpper();
  const vector<bool>& isLower = constraints.getIsLower();
  const vector<bool>& isUpper = constraints.getIsUpper();

  const Vector& x = newPoint.getX();
  direction.resize(0);
  for (int i = 0; i < nDimensions; i ++)
  {
    tmpVector = zero;

    if ( (!isUpper[i]) || (x[i] < upper[i]) )
    {
      // Add +e_i
      tmpVector[i] = scaling[i];
      direction.push_back(tmpVector);
    }
    
    if ( (!isLower[i]) || (x[i] > lower[i]) )
    {
      // Add -e_i
      tmpVector[i] = -1 * scaling[i];
      direction.push_back(tmpVector);
    }

  }

  // Update the step, trueStep, and tag information

  nDirections = direction.size();
  step.resize(nDirections);
  trueStep.resize(nDirections);
  tag.resize(nDirections);

  double newStep = max(newPoint.getStep(), minStep);

  for (int i = 0; i < nDirections; i ++)
  {
    step[i] = newStep;
    trueStep[i] = -1;
    tag[i] = -1;
  }

}
    
void APPSPACK::Directions::setTrueStepAndTag(int i, double trueStep_in, int tag_in)
{
  trueStep[i] = trueStep_in;
  tag[i] = tag_in;
}


void APPSPACK::Directions::print(const string label) const
{
  if (!label.empty())
    cout << "\n" << label << ":\n";

  for (int i = 0; i < nDirections; i ++)
  {
    cout << setw(4) << i << " : ";

    cout << "d = " << direction[i] << " ";

    cout<< "step = " << Print::formatPositiveDouble(step[i]) << " ";

    if (tag[i] != -1) 
    {
      cout << "tag = " << setw(6) << tag[i] << " ";
      cout << "trueStep = " <<  Print::formatPositiveDouble(trueStep[i]);
    }

    cout << "\n";
  }
}

bool APPSPACK::Directions::isStepConverged() const
{
  for (int i = 0; i < nDirections; i ++)
  {
    if (step[i] >= stepTolerance)
      return false;
  }

  return true;
}


void APPSPACK::Directions::reduceStep(int i)
{
  double tmpStep = theta * step[i];

  step[i] = tmpStep;
  trueStep[i] = -1;
  tag[i] = -1;
}

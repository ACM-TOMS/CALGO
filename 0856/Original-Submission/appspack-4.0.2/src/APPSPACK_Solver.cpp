// $Id: APPSPACK_Solver.cpp,v 1.30 2004/11/23 22:36:18 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Solver.cpp,v $ 

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
  \file APPSPACK_Solver.cpp
  \brief Implementation of APPSPACK::Solver
*/
#include "APPSPACK_Solver.hpp"
#include "APPSPACK_Print.hpp"


APPSPACK::Solver::Solver(const Parameter::List& params_in,
			 Executor::Interface& executor_in,
			 const Constraints::Interface& constraints_in) :
  constraints(constraints_in),
  params(params_in),
  print(params),		// modifies params and sets globals
  bestPointPtr(initializeBestPointPtr(params, constraints)), // modifies params
  directions(params, constraints), // modifies params
  conveyor(params, constraints.getScaling(), executor_in), // modifies params
  exchangeList(),
  state(Continue),
  isFunctionTolerance(false),
  isMaxEvaluations(false)
{

  cout << "\n"
       << "-----------------------------------------------------\n"
       << "APPSPACK: Asynchronous Parallel Pattern Search\n"
       << "Written by T. G. Kolda et al., Sandia National Labs\n"
       << "For more information visit \n"
       << "http://software.sandia.gov/appspack\n"
       << "-----------------------------------------------------\n"
       << "\n";

  // Stopping Criteria
  if ((params.isParameterDouble("Function Tolerance")) || (params.isParameterValue("Function Tolerance")))
  {
    isFunctionTolerance = true;
    if (params.isParameterDouble("Function Tolerance"))
    {
      functionTolerance = params.getDoubleParameter("Function Tolerance");
    }
    else if (params.isParameterValue("Function Tolerance"))
    {
      functionTolerance = params.getValueParameter("Function Tolerance");
    }
  }

  if (params.isParameter("Maximum Evaluations"))
  {
    isMaxEvaluations = true;
    maxEvaluations = params.getParameter("Maximum Evaluations", 0);
  }

  boundsTolerance = params.getParameter("Bounds Tolerance", params.getDoubleParameter("Step Tolerance") / 2.0);
  if (boundsTolerance <= 0)
  {
    cout << "APPSPACK::Solver::Solver - Invalid non-positive value for \"Bounds Tolerance\"." << endl;
    throw "APPSPACK Error";
  }

  // Machine Epsilon
  epsMach = params.getParameter("Machine Epsilon", 1.0e-14);

  //! Print
  if (Print::doPrint(Print::InitialData))
    printInitializationInformation();
  
  processNewBestPoint();
}

APPSPACK::Solver::~Solver()
{
  delete bestPointPtr;
}

const APPSPACK::Vector& APPSPACK::Solver::getBestX() const
{
  return bestPointPtr->getX();
}

bool APPSPACK::Solver::isBestF() const
{
  return bestPointPtr->getF().getIsValue(); 
}

double APPSPACK::Solver::getBestF() const
{
  return bestPointPtr->getF().getValue(); 
}


APPSPACK::Point* APPSPACK::Solver::initializeBestPointPtr(Parameter::List& params, const Constraints::Interface& constraints)
{
  // Note that this needs to be called before the Cache constructor because of the Step Tolerance!

  const Vector& lower = constraints.getLower();
  const Vector& upper = constraints.getUpper();
  const vector<bool>& isLower = constraints.getIsLower();
  const vector<bool>& isUpper = constraints.getIsUpper();

  // Initial point
  Vector nominalX;
  if (!params.isParameter("InitialX"))
  {
    nominalX.resize(constraints.getScaling().size());
    for (int i = 0; i < nominalX.size(); i ++)
    {
      if ((isUpper[i]) && (isLower[i]))
	nominalX[i] = ( upper[i] + lower[i] ) / 2.0 ;
      else if (isUpper[i])
	nominalX[i] = upper[i];
      else if (isLower[i])
	nominalX[i] = lower[i];
      else
	nominalX[i] = 0;
    }
  }
  Vector initialX = params.getParameter("Initial X", nominalX);

  // Check that the initial point is the right size
  if (initialX.size() != lower.size())
  {
    cerr << "Error: Size mismatch\n";
    cerr << "The size of the initial X is not the same as the lower bounds.\n";
    throw "APPSPACK Error";
  }

  // Check that initial point is feasible
  for (int i = 0; i < initialX.size(); i ++)
  {
    if (isUpper[i])
      if (initialX[i] > upper[i])
      {
	cerr << "ERROR: Infeasible initial point\n";
	cerr << "The " << i+1 << " component of the initial X (" << initialX[i] 
	     << ") is greater than the upper bound (" << upper[i] << ")" << endl;
	throw "APPSPACK Error";
      }
    if (isLower[i])
      if (initialX[i] < lower[i])
      {
	cerr << "ERROR: Infeasible initial point\n";
	cerr << "The " << i+1 << " component of the initial X (" << initialX[i] 
	     << ") is lower than the lower bound (" << lower[i] << ")" << endl;
	throw "APPSPACK Error";
      }
  }

  // Initial F
  Value initialF;
  if (params.isParameterDouble("Initial F"))
  {
    initialF = params.getDoubleParameter("Initial F");
  }
  else if (params.isParameterValue("Initial F"))
  {
    initialF = params.getValueParameter("Initial F");
  }
  
  // Initial step
  double initialStep = params.getParameter("Initial Step", 1.0);
  
  // Point
  double alpha = params.getParameter("Sufficient Decrease Factor", 0.01);
  
  // Create first best point
  return new Point(initialX, initialF, initialStep, alpha);
}

APPSPACK::Solver::State APPSPACK::Solver::solve()
{
  while (state == Continue)
  {
    iterate();
  }

  // Print the final solution
  if (Print::doPrint(Print::FinalSolution))
  {
    cout << "\nFinal State: " << state << "\n";
    printBestPoint("Final Min");
    directions.print("Final Directions");
    conveyor.getCounter().print();
  }

  return state;
}

APPSPACK::Solver::State APPSPACK::Solver::iterate()
{
    generateTrialPoints();

    conveyor.exchange(exchangeList);

    processEvaluatedTrialPoints();

    return state;
}

// PRIVATE
void APPSPACK::Solver::processNewBestPoint(Point* newBestPointPtr)
{
  // Update the best point
  if (newBestPointPtr != NULL)
  {
    delete bestPointPtr;
    bestPointPtr = newBestPointPtr;
  }

  // Print
  if (Print::doPrint(Print::NewBestPoint))
    printBestPoint("New Min");
  
  // Check for convergence based on the function value
  if ((isFunctionTolerance) && (bestPointPtr->getF() < functionTolerance))
  {
    state = FunctionConverged;
    return;
  }

  // Update the (scaled) search directions
  directions.computeNewDirections(*bestPointPtr);

  // Print
  if (Print::doPrint(Print::NewBestDirections))
    directions.print("New Directions");
}

// PRIVATE
void APPSPACK::Solver::generateTrialPoints()
{
  // Local references
  const Vector& parentX = bestPointPtr->getX();
  //const Vector& scaling = constraints.getScaling();
  const Vector& lower = constraints.getLower();
  const Vector& upper = constraints.getUpper();
  const vector<bool>& isLower = constraints.getIsLower();
  const vector<bool>& isUpper = constraints.getIsUpper();

  Vector& x = tmpVector;
  int n = parentX.size();
  x.resize(n);

  const vector<int> indices = directions.getDirectionIndices();

  for (int i = 0; i < indices.size(); i ++)
  {
    int idx = indices[i];
    const Vector& direction = directions.getDirection(idx);
    const double step = directions.getStep(idx);
    
    // Used in the calculation of the longest feasible step
    double tmpStep = step;
    
    // Calculate the new point
    for (int k = 0; k < n; k ++)
      x[k] = parentX[k] + tmpStep * direction[k];
    
    // Correct violations of the constraints
    for (int j = 0; j < n; j ++) 
    {
      if ((isLower[j]) && ((lower[j] - x[j]) > 0))
      {
	if (lower[j] - x[j] > epsMach)
	{
	  tmpStep = ( lower[j] - parentX[j] ) / direction[j];
	  for (int k = 0; k < n; k ++)
	    x[k] = parentX[k] + tmpStep * direction[k];
	}
	else
	  x[j] = lower[j];	// nudge onto exact boundary
      }
      else if ((isUpper[j]) && ((x[j] - upper[j]) > 0))
      {
	if ((x[j] - upper[j]) > 0)
	{
	  tmpStep = ( upper[j] - parentX[j] ) / direction[j];
	  for (int k = 0; k < n; k ++)
	    x[k] = parentX[k] + tmpStep * direction[k];
	}
	else
	  x[j] = upper[j];	// nudge onto exact boundary
      }
    }

    /*

    // Snap nearby points to the constraints
    for (int j = 0; j < n; j ++) 
    {
      if ( (isLower[j]) && 
	   (abs(lower[j] - x[j]) < (boundsTolerance * scaling[j]) ) )
      {
	x[j] = lower[j];
      }
      else if ( (isUpper[j]) && 
		(abs(upper[j] - x[j]) < (boundsTolerance * scaling[j]) ) )
      {
	x[j] = upper[j];
      }
    }

    */
    
    // Create a new trial point
    Point* newPointPtr = new Point(x, step, *bestPointPtr, idx);
    
    // Save off trial point information
    directions.setTrueStepAndTag(idx, tmpStep, newPointPtr->getTag());
    
    // Push this trial point onto the new trial point list.
    // Ownership of the pointer transfers to the list.
    exchangeList.push(newPointPtr);
    
  }

  if (Print::doPrint(Print::Directions))
    directions.print("Directions After Trial Point Generation");

  if (Print::doPrint(Print::UnevaluatedPoints))
    exchangeList.print("Trial Points Before Conveyor");

}

void APPSPACK::Solver::processEvaluatedTrialPoints()
{
  // Print
  if (Print::doPrint(Print::EvaluatedPoints))
    exchangeList.print("Trial Points After Conveyor");

  // Check for a new best point
  if (exchangeList.best() < *bestPointPtr)
  {
    processNewBestPoint(exchangeList.popBest());
    exchangeList.prune();
    conveyor.prune();
  }
  else 
  {
    // Otherwise, just process the list
    Point* ptr;
    while (!exchangeList.isEmpty())
    {
      ptr = exchangeList.pop();
      if (ptr->getParentTag() == bestPointPtr->getTag())
	directions.reduceStep(ptr->getIndex());
      delete ptr;
    }

    // Check for step length convergence
    if (directions.isStepConverged())
      state = StepConverged;
  }

  // Check for number of evaluations
  if ((state == Continue) && 
      (isMaxEvaluations) && 
      (conveyor.getCounter().getNumEvaluations() >= maxEvaluations))
    state = EvaluationsExhausted; 
}

void APPSPACK::Solver::printInitializationInformation() const
{
  cout << "\n";
  cout << "###########################################" << "\n";
  cout << "###   APPSPACK Initialization Results   ###" << endl;
  cout << "\n";
  
  // Paramters
  cout << "*** Parameter List ***" << "\n";
  params.print();

  // Constraints
  cout << "\n" << "*** Constraints ***" << "\n";
  constraints.print();
  
  // Conveyor
  cout << "\n" << "*** Conveyor ***" << "\n";
  conveyor.print();

  cout << "\n";
  cout << "### End APPSPACK Initialization Results ###" << endl;
  cout << "###########################################" << "\n";
}

void APPSPACK::Solver::printBestPoint(const string label) const
{
  cout << "\n" << label << ": " << *bestPointPtr << endl;
}

ostream& operator<<(ostream& stream, APPSPACK::Solver::State state) 
{
  switch (state)
  {
  case APPSPACK::Solver::Continue:
    stream << "Continue";
    break;
  case APPSPACK::Solver::StepConverged:
    stream << "Step Converged";
    break;
  case APPSPACK::Solver::FunctionConverged:
    stream << "Function Converged";
    break;
  case APPSPACK::Solver::EvaluationsExhausted:
    stream << "Evaluations Exhausted";
    break;
  }

  return stream;
}


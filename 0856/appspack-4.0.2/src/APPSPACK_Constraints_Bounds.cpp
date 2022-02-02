// $Id: APPSPACK_Constraints_Bounds.cpp,v 1.12 2004/10/21 23:06:53 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Constraints_Bounds.cpp,v $ 

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
  \file APPSPACK_Constraints_Bounds.hpp
  \brief Implementation of non-virtual functions in APPSPACK::Constraints::Bounds
*/

#include "APPSPACK_Vector.hpp"
#include "APPSPACK_Constraints_Bounds.hpp"
#include "APPSPACK_Print.hpp"

void APPSPACK::Constraints::Bounds::error(const string& fname, const string& msg)
{
  cout << "APPSPACK::Constraints::" << fname << " - " << msg << endl;
  throw "APPSPACK Error";
}

//static
APPSPACK::Vector APPSPACK::Constraints::Bounds::setup(Parameter::List& params)
{

  Vector s,l,u,isl,isu;
  
  if (params.isParameterVector("Scaling"))
    s = params.getVectorParameter("Scaling");

  if (params.isParameterVector("Upper"))
    u = params.getVectorParameter("Upper");

  if (params.isParameterVector("Lower"))
    l = params.getVectorParameter("Lower");

  if (params.isParameterVector("Is Lower"))
    isl = params.getVectorParameter("Is Lower");

  if (params.isParameterVector("Is Upper"))
    isu = params.getVectorParameter("Is Upper");


  // Determine the vector length
  int n;
  if (s.size() > 0)
    n = s.size();
  else if (l.size() > 0)
    n = l.size();
  else
    error("setup", "Neither \"Scaling\" nor \"Lower\" is defined in the paramater list");


  // Fill l and isl (if necessary)
  if ((l.empty()) && (isl.empty()))
  {
    for (int i = 0; i < n; i ++)
    {
      l.push_back(0);
      isl.push_back(0);
    }

    params.getParameter("Lower", l);
    params.getParameter("Is Lower", isl);
  }
  else if (l.empty())
  {
    for (int i = 0; i < n; i ++)
      if (isl[i] != 0)
	error("setup", "\"Is Lower\" has a true (nonzero) entry but \"Lower\" is not defined");
    
    for (int i = 0; i < n; i ++)
      l.push_back(0);
    
    params.getParameter("Lower", l);
  }
  else if (isl.empty())
  {
    for (int i = 0; i < n; i ++)
      isl.push_back(1);

    params.getParameter("Is Lower", isl);
  }

  // Fill u and isu (if necessary)
  if ((u.empty()) && (isu.empty()))
  {
    for (int i = 0; i < n; i ++)
    {
      u.push_back(0);
      isu.push_back(0);
    }

    params.getParameter("Upper", u);
    params.getParameter("Is Upper", isu);
  }
  else if (u.empty())
  {
    for (int i = 0; i < n; i ++)
      if (isl[i] != 0)
	error("setup", "\"Is Upper\" has a true (nonzero) entry but \"Upper\" is not defined");
    
    for (int i = 0; i < n; i ++)
      u.push_back(0);
    
    params.getParameter("Upper", u);
  }
  else if (isu.empty())
  {
    for (int i = 0; i < n; i ++)
      isu.push_back(1);

    params.getParameter("Is Upper", isu);
  }
   
  // Check sizes
  if (l.size() != n)
    error("setup", "Size mismatch for lower bounds");

  if (u.size() != n)
    error("setup", "Size mismatch for upper bounds");

  if (isu.size() != n)
    error("setup", "Size mismatch for isupper bounds");

  if (isl.size() != n)
    error("setup", "Size mismatch for islower bounds");


  // Fix the scaling (if necessary)
  if (s.empty())
  {  
    // Check that all the bounds are specified
    for (int i = 0; i < n; i ++)
      if ((isu[i] == 0) || (isl[i] == 0))
	error("setup", "\"Scaling\" must be specified if upper and lower bounds are not fully specified");

    for (int i = 0; i < n; i ++)
      s.push_back(u[i] - l[i]);

    params.getParameter("Scaling", s);
  }

  if (s.size() != n)
    error("setup", "Size mismatch for scaling");

  return s;
}

//static
vector<bool> APPSPACK::Constraints::Bounds::convertToBool(const Vector& v)
{
  vector<bool> b;
  b.resize(v.size());
  for (int i = 0; i < v.size(); i ++)
    b[i] = (v[i] == 0) ? false : true;
  return b;
}

//static
void APPSPACK::Constraints::Bounds::checkVector(const string& name, const Vector& v, int n)
{
  if (v.size() != n)
  {
    cout << "APPSPACK::Constraints::Bounds::checkVector - \"" << name << "\" is the wrong size." << endl;
    throw "APPSPACK Error";
  }
}


//static
void APPSPACK::Constraints::Bounds::checkVector(const string& name, const vector<bool>& v, int n)
{
  if (v.size() != n)
  {
    cout << "APPSPACK::Constraints::Bounds::checkVector - \"" << name << "\" is the wrong size." << endl;
    throw "APPSPACK Error";
  }
}


APPSPACK::Constraints::Bounds::Bounds(Parameter::List& params) :
  scaling(setup(params)),
  lower(params.getVectorParameter("Lower")),
  upper(params.getVectorParameter("Upper")),
  isLower(convertToBool(params.getVectorParameter("Is Lower"))),
  isUpper(convertToBool(params.getVectorParameter("Is Upper")))
{
  errorCheck();
}

void APPSPACK::Constraints::Bounds::errorCheck() const
{
  int n = scaling.size();

  if (n == 0)
  {
    cout << "APPSPACK::Constraints::Bounds::errorCheck() - No \"Scaling\" defined." << endl;
    throw "APPSPACK Error";
  }

  checkVector("Lower", lower, n);
  checkVector("Upper", upper, n);
  checkVector("Is Lower", isLower, n);
  checkVector("Is Upper", isUpper, n);

  for (int i = 0; i < n; i ++)
  {
    if (scaling[i] <= 0)
    {
      cout << "APPSPACK::Constraints::Bounds::errorCheck() - Negative \"Scaling\" value." << endl;
      throw "APPSPACK Error";
    }

    if ( (isLower[i]) && (isUpper[i]) && (lower[i] > upper[i]) )
    {
      cout << "APPSPACK::Constraints::Bounds::errorCheck() - Upper bound is less than lower bound." << endl;
      throw "APPSPACK Error";
    }
  }

}

const APPSPACK::Vector& APPSPACK::Constraints::Bounds::getScaling() const
{
  return scaling;
}

const APPSPACK::Vector& APPSPACK::Constraints::Bounds::getLower() const
{
  return lower;
}

const APPSPACK::Vector& APPSPACK::Constraints::Bounds::getUpper() const
{
  return upper;
}

const vector<bool>& APPSPACK::Constraints::Bounds::getIsLower() const
{
  return isLower;
}

const vector<bool>& APPSPACK::Constraints::Bounds::getIsUpper() const
{
  return isUpper;
}

void APPSPACK::Constraints::Bounds::print() const
{
  cout << "\n";
  cout << "Bound Constraints " << endl;

  cout << "lower = ";
  cout << "[";
  for (int i = 0; i < scaling.size(); i ++)
  {
    if (isLower[i])
    {
      cout  << APPSPACK::Print::formatDouble(lower[i]) << " ";
    }
    else 
    {
      cout << "<null>";
      for (int i = 0; i < APPSPACK::Print::precision + 1; i ++)
	cout << " ";
    }

  }
  cout << "]\n";

  cout << "upper = ";
  cout << "[";
  for (int i = 0; i < scaling.size(); i ++)
  {
    if (isUpper[i])
    {
      cout  << APPSPACK::Print::formatDouble(upper[i]) << " ";
    }
    else 
    {
      cout << "<null>";
      for (int i = 0; i < APPSPACK::Print::precision + 1; i ++)
	cout << " ";
    }

  }
  cout << "]\n";

  cout << "scaling = " << scaling << endl;
}


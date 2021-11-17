// $Id: APPSPACK_Print.hpp,v 1.9 2003/11/26 16:27:11 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Print.hpp,v $ 

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
  \file APPSPACK_Print.hpp
  \brief Defines global variables (though all in namespace APPSPACK)
*/

#ifndef APPSPACK_PRINT_HPP
#define APPSPACK_PRINT_HPP

#include "APPSPACK_Common.hpp"
#include "APPSPACK_Value.hpp"

namespace APPSPACK 
{

// forward declaration
namespace Parameter {
class List;
}

/*!
\brief Printing utilities

  <b>Parameters</b>

  These parameters are are stored in the Parameter::List that is
  passed to the constructor. See \ref pageParameters for full details on
  these parameters.

  <ul>
  <li>"Debug"
  <li>"Precision"
  </ul>

*/


class Print 
{

public:
  
  //! Set #precision and #debug from parameter list.
  Print(Parameter::List& params);

  //! Does nothing
  ~Print();

  //! What debugging information to print, see #debug
  enum PrintType 
    {
      //! Final Solution
      FinalSolution = 1,
      //! Every new minimum
      NewBestPoint = 2,
      //! Initial data
      InitialData = 3,
      //! All evaluated trial points
      EvaluatedPoints = 4,
      //! Directions corresponding to new minimum
      NewBestDirections = 5,
      //! All unevaluated trial points
      UnevaluatedPoints = 6,
      //! Directions at every iteration
      Directions = 7
    };
  
  //! Debug level 
  /*!  
    This controls how much information is printed to the
    screen. The higher the value, the more information that's
    printed. 
    See "Debug" in \ref pageParameters.
  */    
  static unsigned int debug;		

  //! Precision for output of real numbers 
  /*!
    See "Precision" in \ref pageParameters.
  */    
  static int precision;		

  //! Class used for printing formatted positive doubles
  class PrintablePositiveDouble
  {
  public:
    PrintablePositiveDouble(double d_in, int precision_in = -1) :  
      d(d_in), 
      precision(precision_in < 0 ? APPSPACK::Print::precision : precision_in)
    {};
    ~PrintablePositiveDouble() {};
    double d;
    int precision;
  };
  
  //! Class used for printing formatted doubles
  class PrintableDouble
  {
  public:
    PrintableDouble(double d_in, int precision_in = -1) :  
      d(d_in), 
      precision(precision_in < 0 ? APPSPACK::Print::precision : precision_in)
    {};
    
    ~PrintableDouble() {};
    double d;
    int precision;
  };
  
  //! Prints a positive double using the specified #precision
  static PrintablePositiveDouble formatPositiveDouble(double d, int precision_in = -1);
  
  //! Prints a double using the specified #precision
  static PrintableDouble formatDouble(double d, int precision_in = -1);
  
  //! Return true if the type of item specified should be printed given the current value of #debug.
  static bool doPrint(enum PrintType type);

};
}

//! Print a positive double with the correct spacing
ostream& operator<< (ostream& stream, const APPSPACK::Print::PrintablePositiveDouble value);
//! Print a double with the correct spacing
ostream& operator<< (ostream& stream, const APPSPACK::Print::PrintableDouble value);


#endif


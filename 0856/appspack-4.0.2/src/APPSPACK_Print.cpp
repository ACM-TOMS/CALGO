// $Id: APPSPACK_Print.cpp,v 1.7 2003/11/26 16:27:11 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Print.cpp,v $ 

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
  \file APPSPACK_Print.cpp
  \brief Actual definition of global variables
*/

#include "APPSPACK_Print.hpp"
#include "APPSPACK_Parameter_List.hpp"

int APPSPACK::Print::precision;
unsigned int APPSPACK::Print::debug;

APPSPACK::Print::Print(Parameter::List& params)
{
  // Statics
  debug = params.getParameter("Debug", 3);
  precision = params.getParameter("Precision", 3);
}

APPSPACK::Print::~Print()
{

}

APPSPACK::Print::PrintablePositiveDouble APPSPACK::Print::formatPositiveDouble(double d, int precision_in)
{
  return PrintablePositiveDouble(d, precision_in);
}

APPSPACK::Print::PrintableDouble APPSPACK::Print::formatDouble(double d, int precision_in)
{
  return PrintableDouble(d, precision_in);
}

/*
APPSPACK::Print::PrintableValue APPSPACK::Print::formatValue(Value v, int precision_in)
{
  return PrintableValue(v, precision_in);
}
*/

bool APPSPACK::Print::doPrint(enum APPSPACK::Print::PrintType type)
{
  return (debug >= type);
}

ostream& operator<< (ostream& stream, const APPSPACK::Print::PrintablePositiveDouble value)
{
  stream.setf(ios::scientific);
  stream.precision(value.precision);
  stream << setw(APPSPACK::Print::precision + 6) << value.d;
  stream.unsetf(ios::scientific);
  return stream;
}

ostream& operator<< (ostream& stream, const APPSPACK::Print::PrintableDouble value)
{
  stream.setf(ios::scientific);
  stream.precision(value.precision);
  stream << setw(APPSPACK::Print::precision + 7) << value.d;
  stream.unsetf(ios::scientific);
  return stream;
}

/*
ostream& operator<< (ostream& stream, const APPSPACK::Print::PrintableValue value)
{
  if (value.v.getIsValue())
    stream << APPSPACK::Print::formatDouble(value.v.getValue(), value.precision);
  else
  {
    stream << "<null>";
    for (int i = 0; i < value.precision + 1; i ++)
      stream << " ";
  }
  return stream;
}
*/

/*
ostream& operator<< (ostream& os, const APPSPACK::PrintablePositiveDouble value)
{
  ostringstream s;
  s.setf(ios::scientific);
  s.precision(APPSPACK::Print::precision);
  s << setw(APPSPACK::Print::precision + 6) << value.d;
  return os << s.str();
}

ostream& operator<< (ostream& os, const APPSPACK::PrintableDouble value)
{
  ostringstream s;
  s.setf(ios::scientific);
  s.precision(APPSPACK::Print::precision);
  os << setw(APPSPACK::Print::precision + 7) << value.d;
  return os << s.str();
}
*/



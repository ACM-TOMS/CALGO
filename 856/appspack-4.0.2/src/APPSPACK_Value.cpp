// $Id: APPSPACK_Value.cpp,v 1.8 2003/11/26 16:27:11 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Value.cpp,v $ 

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
  \file APPSPACK_Value.cpp
  \brief Implementation of APPSPACK::Value
*/
#include "APPSPACK_Value.hpp"
#include "APPSPACK_Print.hpp"

APPSPACK::Value::Value() :
  isValue(false),
  value(0)
{
}

APPSPACK::Value::Value(bool isValue_in, double value_in) :
  isValue(isValue_in),
  value(value_in)
{
}

APPSPACK::Value::Value(double value_in) :
  isValue(true),
  value(value_in)
{
}

APPSPACK::Value::Value(const Value& source) :
  isValue(source.isValue),
  value(source.value)
{
}

APPSPACK::Value::~Value()
{

}

bool APPSPACK::Value::getIsValue() const
{
  return isValue;
}

double APPSPACK::Value::getValue() const
{
  return value;
}

void APPSPACK::Value::operator=(const Value& source)
{
  isValue = source.isValue;
  value = source.value;
}

void APPSPACK::Value::operator=(double source)
{
  isValue = true;
  value = source;
}

void  APPSPACK::Value::setValueTo(bool isValue_in, double value_in)
{
  isValue = isValue_in;
  value = value_in;
}

void  APPSPACK::Value::setValueTo(double value_in)
{
  isValue = true;
  value = value_in;
}

void APPSPACK::Value::setValueToUnknown()
{
  isValue = false;
  value = 0;
}

bool APPSPACK::Value::operator<(const Value& source) const
{
  // Case I: f(x) = infty
  if (!isValue)
    return false;

  // Case II: f(y) = infty and f(x) is finite
  if (!source.isValue)
    return true;

  // Case III: f(x) and f(y) are finite
  return (value < source.value);
}

bool APPSPACK::Value::operator>(const Value& source) const
{
  // Case I: f(y) = infty
  if (!source.isValue)
    return false;

  // Case II: f(x) = infty and f(y) is finite
  if (!isValue)
    return true;

  // Case III: f(x) and f(y) are finite
  return (value > source.value);
}

bool APPSPACK::Value::operator==(const Value& source) const
{
  // Case I: f(x) = infty, so return true if f(y) is also infty and
  // false otherwise.
  if (!isValue)
    return (!source.isValue);

  // Case II: f(x) is finite and f(y) = infty.
  if (!source.isValue)
    return false;

  // Case III: Both f(x) and f(y) are finite.
  return (value == source.value);
}

bool APPSPACK::Value::isSufficientDecrease
(const Value& source, double rho) const
{
  // Case I: f(x) = infty
  if (!isValue)
    return false;

  // Case II: f(x) is finite and f(y) = infty
  if (!source.isValue)
    return true;

  // Case III: f(x) and f(y) are finite
  return (value < (source.value - rho));
}

  
ostream& APPSPACK::Value::leftshift(ostream& stream) const
{
  // Print point to the given stream. 

  if (isValue)
    stream << Print::formatDouble(value);
  else 
  {
    stream << "<null>";
    for (int i = 0; i < APPSPACK::Print::precision + 1; i ++)
      stream << " ";
  }

  return stream;
}

ostream& operator<<(ostream& stream, const APPSPACK::Value& value)
{
  return value.leftshift(stream);
}

// $Id: APPSPACK_Value.hpp,v 1.8 2004/03/29 05:25:34 wehart Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Value.hpp,v $ 

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
  \file APPSPACK_Value.hpp
  \brief Class definition of APPSPACK::Value
*/
#ifndef APPSPACK_VALUE_HPP
#define APPSPACK_VALUE_HPP

#include "APPSPACK_Common.hpp"

namespace APPSPACK
{
//! Stores and Manipulates a (function) value in the range \f$(-\infty,+\infty]\f$.
/*! 
  Stores a function value, v, which may be infinite (or unknown).
  
  Two function values \f$a\f$ and \f$b\f$ are compared as follows:
  \f[
  \begin{array}{ll}
  a < b & \mbox{ if $a$ and $b$ are finite and $a < b$, or $a$ is finite and $b = +\infty$}.\\
  a > b & \mbox{ if $a$ and $b$ are finite and $a > b$, or $a = +\infty$ and $b$ is finite}.\\
  a = b & \mbox{ if $a$ and $b$ are finite and $a = b$, or $a = +\infty$ and $b = +\infty$}.
  \end{array}
  \f]

*/
class Value
{
public:

  //@{ \name Constructors

  /*! Sets v = \f$+\infty\f$.*/
  Value();

  /*!
    \param isValue_in - True if a finite value exists. False otherwise.
    \param value_in - The function value (only has meaning if isKnown is true)
  */
  Value(bool isValue_in, double value_in);

  /*!
    \param value_in - The finite function value

    (Also sets #isValue to true.)
  */
  Value(double value_in);

  /*! Copy constructor */
  Value(const Value& source);

  //@}

  //@{ \name Destructor

  /*! Destructor */
  ~Value();

  //@}

  //@{ \name Accessors
  
  bool getIsValue() const;
  double getValue() const;
  
  //@

  //@}

  //@{ \name Manipulators

  /*! Copy operator */
  void operator=(const Value& source);

  /*! Copy from double */
  void operator=(double source);

  /*! See Value(bool isValue_in, double value_in) */
  void setValueTo(bool isValue_in, double value_in);

  /*! See Value(double value_in) */
  void setValueTo(double value_in);

  /*! Reset to \f$+\infty\f$. */
  void setValueToUnknown();

  //@}

  //@{ \name Comparisons

  /*! Return true if v < w; false otherwise. */
  bool operator<(const Value& w) const;

  /*! Return true if v > w; false otherwise. */
  bool operator>(const Value& w) const;

  /*! Return true if v = w; false otherwise. */
  bool operator==(const Value& w) const;

  /*! Return true if v < (w - rho); false otherwise. */
  bool isSufficientDecrease(const Value& w, double rho) const;

  //@}

  //! Print object to the given stream
  ostream& leftshift(ostream& stream) const;

private:

  //! If true, there is a finite function value. Otherwise the function value is \f$+\infty\f$.
  bool isValue;

  //! The function value (only has meaning if #isValue is true).
  double value;

};

}

//! Printing an APPSPACK::Value
ostream& operator<<(ostream& stream, const APPSPACK::Value& value);

#endif

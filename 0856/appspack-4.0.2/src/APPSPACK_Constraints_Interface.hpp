// $Id: APPSPACK_Constraints_Interface.hpp,v 1.10 2004/10/21 23:06:53 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Constraints_Interface.hpp,v $ 

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
  \file APPSPACK_Constraints_Interface.hpp
  \brief Class definition for APPSPACK::Constraints::Interface
*/

#ifndef APPSPACK_CONSTRAINTS_INTERFACE_HPP
#define APPSPACK_CONSTRAINTS_INTERFACE_HPP

#include "APPSPACK_Vector.hpp"

namespace APPSPACK
{

//! Namespace for constraint-related objects
namespace Constraints
{

//! Abstract interface for constraints
class Interface
{

public:
  /*! Constructor */
  Interface() {};
  
  /*! Destructor */
  virtual ~Interface() {} ;

  //! Return the scaling vector. 
  /*! The scaling vector is typically defined as
    \f[
    s_i = u_i - \ell_i
    \f]
    where \f$s_i\f$ represents the i-th entry of the scaling vector,
    \f$u_i\f$ represents the i-th upper bound, and
    \f$\ell_i\f$ represents the i-th lower bound.
  */
  virtual const Vector& getScaling() const = 0;
  
  //! Return vector of lower bounds
  virtual const Vector& getLower() const = 0;

  //! Return vector of upper bounds
  virtual const Vector& getUpper() const = 0;

  //! Return boolean vector where each vector is true if the corresponding lower bound is defined
  virtual const vector<bool>& getIsLower() const = 0;

  //! Return boolean vector where each vector is true if the corresponding upper bound is defined
  virtual const vector<bool>& getIsUpper() const = 0;

  //! Optional print function.
  /*! Defaults to doing nothing. Can optionally print out relevant
    information about the constraints. */
  virtual void print() const {};

};

}
}

#endif

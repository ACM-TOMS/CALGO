// $Id: APPSPACK_Executor_Serial.hpp,v 1.9 2004/11/23 22:26:01 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Executor_Serial.hpp,v $ 

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
  \file APPSPACK_Executor_Serial.hpp
  \brief Class definition for APPSPACK::Executor::Serial
*/
#ifndef APPSPACK_EXECUTOR_SERIAL_HPP
#define APPSPACK_EXECUTOR_SERIAL_HPP

#include "APPSPACK_Executor_Interface.hpp"
#include "APPSPACK_Evaluator_Interface.hpp"

namespace APPSPACK 
{

namespace Executor
{

//! Serial Executor
/*!  Serial implementation of the Executor::Interface. Mainly useful
  for testing.  There are no workers; instead, each evaluations is
  performed when the spawn() function is called and stored until the
  recv() function is called to retrieve the result.
*/
class Serial : public Interface
{
public:

  //! Constructor 
  Serial(Evaluator::Interface& evaluator_in);

  //! Destructor 
  virtual ~Serial() {};

  // derived
  virtual bool isWaiting() const;

  // derived
  virtual bool spawn(const Vector& x_in, int tag_in);

  // derived
  virtual int recv(int& tag_out, bool& isF_out, double& f_out, string& msg_out);

  // Prints out the Evaluator information.
  virtual void print() const;

private:

  //! Interface to object that computes the actual function evaluation
  Evaluator::Interface& evaluator;

  //! Is the executor currently busy?
  bool isFree;
  
  //! Tag of current point, if any
  int tag;

  //! True if a function value exists
  bool isF;

  //! The function value
  double f;

  //! The message
  string msg;

};

}

}

#endif

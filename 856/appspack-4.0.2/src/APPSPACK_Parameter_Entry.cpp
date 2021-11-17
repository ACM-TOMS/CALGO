// $Id: APPSPACK_Parameter_Entry.cpp,v 1.5.2.1 2005/06/29 17:07:42 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Parameter_Entry.cpp,v $ 

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
  \file APPSPACK_Parameter_Entry.cpp
  \brief Implementation of APPSPACK::Parameter::Entry
*/

#include "APPSPACK_Parameter_Entry.hpp"
#include "APPSPACK_Parameter_List.hpp"
#include "APPSPACK_GCI.hpp"

using namespace APPSPACK::Parameter;

Entry::Entry() : 
  type(APPSPACK_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(false)
{
}

Entry::Entry(const Entry& source) :
  type(APPSPACK_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(false)
{
  operator=(source);
}

Entry& Entry::operator=(const Entry& source)
{
  if (&source == this)
    return *this;

  reset();

  type = source.type;
  bval = source.bval;
  ival = source.ival;
  dval = source.dval;
  sval = source.sval;
  
  if ((type == APPSPACK_LIST) && (source.lval != NULL)) 
  {
    lval = new List(*source.lval);
  }
  
  valueval = source.valueval;
  vectorval = source.vectorval;

  isGotten = source.isGotten;
  isSetByGet = source.isSetByGet;

  return *this;
}

Entry::Entry(bool value, bool isCreatedByGet) : 
  type(APPSPACK_BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(isCreatedByGet)
{
}

Entry::Entry(int value, bool isCreatedByGet) : 
  type(APPSPACK_INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(double value, bool isCreatedByGet) : 
  type(APPSPACK_DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const string& value, bool isCreatedByGet) : 
  type(APPSPACK_STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  lval(NULL),
  valueval(),
  vectorval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const Value& value, bool isCreatedByGet) : 
  type(APPSPACK_VALUE),
  bval(false),
  ival(0),
  dval(0),
  sval("" ),
  lval(NULL),
  valueval(value),
  vectorval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const Vector& value, bool isCreatedByGet) : 
  type(APPSPACK_VECTOR),
  bval(false),
  ival(0),
  dval(0),
  sval("" ),
  lval(NULL),
  valueval(),
  vectorval(value),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::~Entry() 
{
  reset();
}

void Entry::reset()
{
  type = APPSPACK_NONE;

  delete lval;
  lval = NULL;

  isGotten = false;
  isSetByGet = false;
}

void Entry::setValue(bool value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_BOOL;
  bval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(int value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_INT;
  ival = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(double value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_DOUBLE;
  dval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const char* value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const string& value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const Value& value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_VALUE;
  valueval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const Vector& value, bool isCreatedByGet)
{
  reset();
  type = APPSPACK_VECTOR;
  vectorval = value;
  isSetByGet = isCreatedByGet;
}

APPSPACK::Parameter::List& Entry::setList(bool isCreatedByGet)
{
  reset();
  type = APPSPACK_LIST;
  lval = new List();
  isSetByGet = isCreatedByGet;
  isGotten = true;
  return *lval;
}


bool Entry::isBool() const
{
  return (type == APPSPACK_BOOL);
}

bool Entry::isInt() const
{
  return (type == APPSPACK_INT);
}

bool Entry::isDouble() const
{
  return (type == APPSPACK_DOUBLE);
}

bool Entry::isString() const
{
  return (type == APPSPACK_STRING);
}

bool Entry::isList() const
{
  return (type == APPSPACK_LIST);
}

bool Entry::isValue() const
{
  return (type == APPSPACK_VALUE);
}

bool Entry::isVector() const
{
  return (type == APPSPACK_VECTOR);
}

bool Entry::getBoolValue() const
{
  isGotten = true;
  return bval;
}

int Entry::getIntValue() const
{
  isGotten = true;
  return ival;
}

double Entry::getDoubleValue() const
{
  isGotten = true;
  return dval;
}

const string& Entry::getStringValue() const
{
  isGotten = true;
  return sval;
}

APPSPACK::Parameter::List& Entry::getListValue() 
{
  isGotten = true;
  return *lval;
}

const APPSPACK::Parameter::List& Entry::getListValue() const
{
  isGotten = true;
  return *lval;
}

const APPSPACK::Value& Entry::getValueValue() const
{
  isGotten = true;
  return valueval;
}

const APPSPACK::Vector& Entry::getVectorValue() const
{
  isGotten = true;
  return vectorval;
}

ostream& Entry::leftshift(ostream& stream) const
{
  switch(type) {
  case APPSPACK_BOOL: 
    stream << (bval ? "true" : "false");
    break;
  case APPSPACK_INT:
    stream << ival;
    break;
  case APPSPACK_DOUBLE:
    stream << dval;
    break;
  case APPSPACK_STRING:
    stream << "\"" << sval << "\"";
    break;
  case APPSPACK_LIST:
    break;
  case APPSPACK_VALUE:
    stream << valueval;
    break;
  case APPSPACK_VECTOR:
    stream << vectorval;
    break;
  default:
    stream << "(empty non-typed parameter)";
    break;
  }

  if (isSetByGet)
    stream << "   [default]";
  else if (!isGotten)
    stream << "   [unused]";
  

  return stream;
}


void Entry::pack() const
{
  switch(type) {
  case APPSPACK_BOOL: 
    GCI::pack(APPSPACK_BOOL);
    GCI::pack(bval);
    break;
  case APPSPACK_INT:
    GCI::pack(APPSPACK_INT);
    GCI::pack(ival);
    break;
  case APPSPACK_DOUBLE:
    GCI::pack(APPSPACK_DOUBLE);
    GCI::pack(dval);
    break;
  case APPSPACK_STRING:
    GCI::pack(APPSPACK_STRING);
    GCI::pack(sval);
    break;
  case APPSPACK_LIST:
    GCI::pack(APPSPACK_LIST);
    lval->pack();
    break;
  case APPSPACK_VALUE:
    {
      GCI::pack(APPSPACK_VALUE);
      bool isF = valueval.getIsValue();
      double f = valueval.getValue();
      GCI::pack(isF);
      GCI::pack(f);
      break;
    }
  case APPSPACK_VECTOR:
    GCI::pack(APPSPACK_VECTOR);
    GCI::pack(vectorval);
    break;
  default:
    cerr << "APPSPACK::Parameter::Entry::pack - Empty non-typed parameter";
    throw "APPSPACK Error";
    break;
  }

  GCI::pack(isGotten);
  GCI::pack(isSetByGet);
}

void Entry::unpack()
{
  int itype;
  GCI::unpack(itype);

  switch(itype) {
  case APPSPACK_BOOL: 
    type = APPSPACK_BOOL;
    GCI::unpack(bval);
    break;
  case APPSPACK_INT:
    type = APPSPACK_INT;
    GCI::unpack(ival);
    break;
  case APPSPACK_DOUBLE:
    type = APPSPACK_DOUBLE;
    GCI::unpack(dval);
    break;
  case APPSPACK_STRING:
    type = APPSPACK_STRING;
    GCI::unpack(sval);
    break;
  case APPSPACK_LIST:
    type = APPSPACK_LIST;
    lval = new List();
    lval->unpack();
    break;
  case APPSPACK_VALUE:
    type = APPSPACK_VALUE;
    bool isF;
    double f;
    GCI::unpack(isF);
    GCI::unpack(f);
    valueval.setValueTo(isF,f);
    break;
  case APPSPACK_VECTOR:
    type = APPSPACK_VECTOR;
    GCI::unpack(vectorval);
    break;
  default:
    cerr << "APPSPACK::Parameter::Entry::pack - Empty non-typed parameter";
    throw "APPSPACK Error";
    break;
  }

  GCI::unpack(isGotten);
  GCI::unpack(isSetByGet);
}

ostream& operator<<(ostream& stream, const Entry& e)
{
  return e.leftshift(stream);
}



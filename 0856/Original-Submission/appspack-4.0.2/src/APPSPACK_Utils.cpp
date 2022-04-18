// $Id: APPSPACK_Utils.cpp,v 1.17.2.2 2005/06/29 17:49:35 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/APPSPACK_Utils.cpp,v $ 

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
  \file APPSPACK_Utils.cpp
  \brief Implementations of utility functions in the APPSPACK namespace
*/

#include "APPSPACK_Utils.hpp"
#include "APPSPACK_Common.hpp"


#ifndef SNL_TFLOPS_ENV

bool APPSPACK::getNextQuotedString(const string& line, string::size_type& pos, string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  string::size_type linelength = line.length();

  // Find the location of the first quote
  string::size_type pos1 = line.find('"', pos); 
  
  // Check that the operation succeeded and that we're not at the end of the line
  if ((pos1 == string::npos) || (pos1 == linelength - 1))
  {
    pos = string::npos;
    return false;
  }

  // Advance to first character after the first quote
  pos1 = pos1 + 1;

  // Find the location of the second quote
  string::size_type pos2 = line.find('"', pos1); 

  // Check that the operation was successful
  if (pos2 == string::npos)
  {
    pos = string::npos;
    return false;
  }

  // Compute the length of the expression
  string::size_type length = pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool APPSPACK::getNextString(const string& line, string::size_type& pos, string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  string::size_type linelength = line.length();

  // Find the location of the first non-space character
  string::size_type pos1 = line.find_first_not_of(' ', pos); 
  
  // Check that the operation succeeded 
  if (pos1 == string::npos) 
  {
    pos = string::npos;
    return false;
  }

  // Find the location of the next space, if any
  string::size_type pos2 = line.find(' ', pos1); 

  // Compute the length of the expression
  string::size_type length = (pos2 == string::npos) ? linelength - pos1 : pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool APPSPACK::getNextInt(const string& line, string::size_type& pos, int& value)
{
  string field;
  if ((!getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  return ((sscanf(field.c_str(), "%d", &value)) == 1);
}

bool APPSPACK::getNextDouble(const string& line, string::size_type& pos, double& value)
{
  string field;
  if ((!getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  return ((sscanf(field.c_str(), "%le", &value)) == 1);
}

bool APPSPACK::processTextInputFileLine(const string& line, Parameter::List& params, Parameter::List*& subPtr)
{
  
  string::size_type pos;	// current position inside line

  string field;			// name of parameter

  bool tmpbool;			// for reading in a boolean
  int tmpint;			// for reading in an int
  double tmpdouble;		// for reading in a double
  string tmpstring;		// for reading in a string
  Value tmpvalue;		// for reading in an APPSPACK::Value
  Vector tmpvector;		// for reading in an APPSPACK::Vector

  string type;			// parameter type (used in new style)

  static int n = -1;			// size of vector (used in old style)
  char c;			// used to read in '=' character in old style

  if (line.size() == 0) 	// empty line
  {
    return false;
  }
  else if (line[0] == '#')	// comment
  {
    return false;
  }
  else if (line[0] == '@')	// sublist command
  {
    subPtr = &params;
    
    pos = 0;
    
    while (pos != string::npos)
    {
      if ((getNextQuotedString(line, pos, field)) && (field.size() > 0))
	subPtr = &(subPtr->sublist(field));
    }

    return true;
  }
  else if (line[0] == '"')	// new style
  {

    // Get the name
    pos = 0;
    if ((!getNextQuotedString(line, pos, field)) || (field.empty()))
      return false;
    
    // Read in the type
    if (!getNextString(line, pos, type))
      return false;

    if (type == "int")
    {
      if (!getNextInt(line, pos, tmpint))
      {
	return false;
      }
      else
      {
	subPtr->setParameter(field, tmpint);
	return true;
      }
    }
    else if (type == "bool")
    {
      if (!getNextString(line, pos, tmpstring))
      {
	return false;
      }
      else
      {
	subPtr->setParameter(field, (tmpstring == "true") );
	return true;
      }
    }
    else if (type == "double")
    {
      if (!getNextDouble(line, pos, tmpdouble))
      {
	return false;
      }
      else
      {
	subPtr->setParameter(field, tmpdouble);
	return true;
      }
    }
    else if (type == "string")
    {
      if ((!getNextQuotedString(line, pos, tmpstring)) || (tmpstring.empty()))
      {
	return false;
      }
      else
      {
	subPtr->setParameter(field, tmpstring);
	return true;
      }
    }
    else if (type == "vector")
    {
      // get the size
      if (!getNextInt(line, pos, tmpint))
	return false;

      if (tmpint < 0)
	return false;

      tmpvector.resize(tmpint);
      for (int i = 0; i < tmpint; i ++)
      {
	if (!getNextDouble(line, pos, tmpvector[i]))
	  return false;
      }
      
      subPtr->setParameter(field, tmpvector);	
      return true;
    }
    else
      return false;

  } // end new style

  else				// old style
  {
    istringstream linein(line);
    linein >> field;
    
    if ((field == "n_parameters") || (field == "continuous_design"))
    {
      linein >> c;
      linein >> n;
    }    
    
    else if ((field == "initial_point") || (field == "cdv_initial_point"))
    {
      if (n == -1)
      {
	cout << "Error: initial_point specified before n_parameters " << endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Solver").setParameter("Initial X", tmpvector);
    }
    
    else if (field == "initial_f") 
    {
      linein >> c;
      linein >> tmpdouble;
      tmpvalue.setValueTo(true, tmpdouble);
      params.sublist("Solver").setParameter("Initial F", tmpvalue);
    }
    
    else if ((field == "lower_bounds") || (field == "cdv_lower_bounds"))
    {
      if (n == -1)
      {
	cout << "Error:lower_bounds specified before n_parameters " << endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Bounds").setParameter("Lower", tmpvector);
      
    }
    
    else if ((field == "upper_bounds") || (field == "cdv_upper_bounds"))
    {
      if (n == -1)
      {
	cout << "Error:upper_bounds specified before n_parameters " << endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Bounds").setParameter("Upper", tmpvector);
      
    }
    
    else if (field == "executable") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").setParameter("Executable Name", tmpstring);
    }
    
    else if (field == "params_prefix") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").setParameter("Input Prefix", tmpstring);
    }
    
    else if (field == "result_prefix") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").setParameter("Output Prefix", tmpstring);
    }
    
    else 
    {
      cout << "Ignoring unrecognized field." << field << "\n";
    }
    
  } /* old style */

  return true;
}

bool APPSPACK::parseTextInputFile(const string filename, Parameter::List& params)
{
  // Open the input file
  ifstream fin;
  fin.open(filename.c_str());

  if (!fin) 
  {
    cout << "Error: Cannot find input file " << filename << endl;
    return false;
  }

  string line;			// one line from input file
  Parameter::List* subPtr;	// sublist pointer

  subPtr = &params;

  while (!fin.eof())
  {
    getline(fin, line);
    if (!processTextInputFileLine(line, params, subPtr))
    {
      // cout << "Ignoring input file line: " << line << endl;
    }

  }

  fin.close();

  return true;
}


#endif

APPSPACK::Vector APPSPACK::createZeroVector(int n)
{
  Vector z;
  z.resize(n);
  for (int i = 0; i < n ; i ++)
    z[i] = 0;
  return z;
}


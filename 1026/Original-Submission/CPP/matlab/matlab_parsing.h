//@HEADER
// ******************************************************************************
//
//  CP-CALS: Software for computing the Canonical Polyadic Decomposition using
//  the Concurrent Alternating Least Squares Algorithm.
//
//  Copyright (c) 2020, Christos Psarras
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//  3. Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  ***************************************************************************************************
//
//  This file includes modified versions of code found in Genten,
//  which is covered by the following license:
//
//      Genten: Software for Generalized Tensor Decompositions
//      by Sandia National Laboratories
//
//  Sandia National Laboratories is a multimission laboratory managed
//  and operated by National Technology and Engineering Solutions of Sandia,
//  LLC, a wholly owned subsidiary of Honeywell International, Inc., for the
//  U.S. Department of Energy's National Nuclear Security Administration under
//  contract DE-NA0003525.
//
//  Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
//  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
//  Government retains certain rights in this software.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  1. Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ***************************************************************************************************
//@HEADER

#ifndef CALS_MATLAB_PARSING_H
#define CALS_MATLAB_PARSING_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "cals.h"

using std::vector;
using std::string;
using std::find;
using std::cerr;
using std::endl;

namespace cals::matlab::parsing
{
  template<typename T>
  T parse_enum(vector<string> &args, const string &cl_arg, T default_value, unsigned num_values, const T *values,
                   const char *const *names)
  {
    auto it = find(args.begin(), args.end(), cl_arg);
    // If not found, try removing the '--'
    if ((it == args.end()) && (cl_arg.size() > 2) &&
        (cl_arg[0] == '-') && (cl_arg[1] == '-'))
    {
      it = find(args.begin(), args.end(), cl_arg.substr(2));
    }
    if (it != args.end())
    {
      auto arg_it = it;
      // get next cl_arg
      ++it;
      if (it == args.end())
      {
        args.erase(arg_it);
        return default_value;
      }
      // convert to string
      std::string arg_val = *it;
      // Remove argument from list
      args.erase(arg_it, ++it);
      // find name in list of names
      for (unsigned i = 0; i < num_values; ++i)
      {
        if (arg_val == names[i])
          return values[i];
      }
      // if we got here, name wasn't found
      std::ostringstream error_string;
      error_string << "Bad input: " << cl_arg << " " << arg_val << ",  must be one of the values: ";
      for (unsigned i = 0; i < num_values; ++i)
      {
        error_string << names[i];
        if (i != num_values - 1)
          error_string << ", ";
      }
      error_string << "." << endl;
      cerr << error_string.str() << endl;
      exit(1);
    }
    // return default value if not specified on command line
    return default_value;
  }

  int parse_indx(vector<string> &args, const string &cl_arg, int default_value, int min, int max);

  double parse_real(vector<string> &args, const string &cl_arg, double default_value, double min, double max);

  bool parse_bool(vector<string>& args, const string& cl_arg_on, const string& cl_arg_off, bool default_value);

  /** Parse options from matlab call (defaults are the struct default values).
   *
   * This function is used to parse arguments from Matlab.
   *
   * @param params CalsParams object containing the parameters for a call to cals.
   * @param args Vector of strings containing CALS parameters obtained by the command line.
   */
  void parse(cals::CalsParams &params, vector<string> &args);
}
#endif //CALS_MATLAB_PARSING_H

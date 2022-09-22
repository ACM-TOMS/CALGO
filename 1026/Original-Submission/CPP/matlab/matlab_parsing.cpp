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

#include "matlab_parsing.h"

namespace cals::matlab::parsing
{
  int parse_indx(vector<string> &args, const string &cl_arg, int default_value, int min, int max)
  {
    int tmp = default_value;
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
        return tmp;
      }
      // convert to ttb_indx
      char *cend = 0;
      tmp = std::strtol(it->c_str(), &cend, 10);
      // check if cl_arg is actually a ttb_indx
      if (it->c_str() == cend)
      {
        std::ostringstream error_string;
        error_string << "Unparseable input: " << cl_arg << " " << *it << ", must be an integer" << endl;
        cerr << error_string.str() << endl;
        exit(1);
      }
      // Remove argument from list
      args.erase(arg_it, ++it);
    }
    // check if ttb_real is within bounds
    if (tmp < min || tmp > max)
    {
      std::ostringstream error_string;
      error_string << "Bad input: " << cl_arg << " " << tmp << ",  must be in the range (" << min << ", " << max
                   << ")" << endl;
      cerr << error_string.str() << endl;
      exit(1);
    }
    return tmp;
  }

  double parse_real(vector<string> &args, const string &cl_arg, double default_value, double min, double max)
  {
    double tmp = default_value;
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
        return tmp;
      }
      // convert to ttb_real
      char *cend = 0;
      tmp = std::strtod(it->c_str(), &cend);
      // check if cl_arg is actually a ttb_real
      if (it->c_str() == cend)
      {
        std::ostringstream error_string;
        error_string << "Unparseable input: " << cl_arg << " " << *it << ", must be a double" << endl;
        cerr << error_string.str() << endl;
        exit(1);
      }
      // Remove argument from list
      args.erase(arg_it, ++it);
    }
    // check if ttb_real is within bounds
    if (tmp < min || tmp > max)
    {
      std::ostringstream error_string;
      error_string << "Bad input: " << cl_arg << " " << tmp << ",  must be in the range (" << min << ", " << max
                   << ")" << endl;
      cerr << error_string.str() << endl;
      exit(1);
    }
    return tmp;
  }

  bool parse_bool(vector<string> &args, const string &cl_arg_on, const string &cl_arg_off, bool default_value)
  {
    // return true if arg_on is found
    auto it = find(args.begin(), args.end(), cl_arg_on);
    // If not found, try removing the '--'
    if ((it == args.end()) && (cl_arg_on.size() > 2) &&
        (cl_arg_on[0] == '-') && (cl_arg_on[1] == '-'))
    {
      it = find(args.begin(), args.end(), cl_arg_on.substr(2));
    }
    if (it != args.end())
    {
      args.erase(it);
      return true;
    }

    // return false if arg_off is found
    it = find(args.begin(), args.end(), cl_arg_off);
    // If not found, try removing the '--'
    if ((it == args.end()) && (cl_arg_off.size() > 2) &&
        (cl_arg_off[0] == '-') && (cl_arg_off[1] == '-'))
    {
      it = find(args.begin(), args.end(), cl_arg_off.substr(2));
    }
    if (it != args.end())
    {
      args.erase(it);
      return false;
    }

    // return default value if not specified on command line
    return default_value;
  }

  void parse(cals::CalsParams &params, vector<string> &args)
  {
    // Read factor matrix update method
    const cals::update::UPDATE_METHOD update_types[] = {cals::update::UPDATE_METHOD::UNCONSTRAINED,
                                                        cals::update::UPDATE_METHOD::NNLS};
    const char *update_names[] = {"unconstrained", "nnls"};
    params.update_method =
        cals::matlab::parsing::parse_enum(args, "update-method", params.update_method, 2, update_types, update_names);

    // Read MTTKRP method
    const cals::mttkrp::MTTKRP_METHOD mttkrp_types[] = {cals::mttkrp::MTTKRP_METHOD::MTTKRP,
                                                        cals::mttkrp::MTTKRP_METHOD::TWOSTEP0,
                                                        cals::mttkrp::MTTKRP_METHOD::TWOSTEP1,
                                                        cals::mttkrp::MTTKRP_METHOD::AUTO};
    const char *mttkrp_names[] = {"mttkrp", "twostep0", "twostep1", "auto"};
    params.mttkrp_method =
        cals::matlab::parsing::parse_enum(args, "mttkrp-method", params.mttkrp_method, 4, mttkrp_types, mttkrp_names);

    // Read rest of parameters
    params.max_iterations = cals::matlab::parsing::parse_indx(args, "maxiters", params.max_iterations, 1, INT32_MAX);
    params.buffer_size = cals::matlab::parsing::parse_indx(args, "buffer-size", params.buffer_size, 1, INT32_MAX);
    params.tol = cals::matlab::parsing::parse_real(args, "tol", params.tol, 0.0, FLT_MAX);

    params.cuda = cals::matlab::parsing::parse_bool(args, "cuda", "no-cuda", params.cuda);

    params.line_search = cals::matlab::parsing::parse_bool(args, "ls", "no-ls", params.line_search);
    params.line_search_interval = cals::matlab::parsing::parse_indx(args, "ls-interval", params.line_search_interval, 2, INT32_MAX);
    params.line_search_step = cals::matlab::parsing::parse_real(args, "ls-step", params.line_search_step, 0.001, FLT_MAX);
  }
}

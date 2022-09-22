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

#include "matlab.h"
#include "cals.h"
#include "timer.h"
#include "matlab_parsing.h"

using std::cout;
using std::endl;

extern "C" {

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  try
  {
    if (nrhs < 3)
    {
      cout << "Expected at least 3 command line arguments" << endl;
      return;
    }

    bool debug = false;
    DEBUG(debug = true;)

    // Parse inputs
    auto args = mxBuildArgList(nrhs, 3, prhs);
    cals::CalsParams cals_params;
    cals::matlab::parsing::parse(cals_params, args);
    cals_params.print();
    // TODO Implement the check_print_unused_args function
    //    if (cals_params.check_and_print_unused_args(args, cout))
    //    {
    //      cals_params.print_help(cout);
    //      throw std::string("Invalid command line arguments.");
    //    }

    // Get tensor
    cals::Tensor X = mxGetTensor(prhs[0], debug);

    // Receive a vector of ranks
    double *ranks = mxGetDoubles(prhs[1]);
    int n_ranks = mxGetNumberOfElements(prhs[1]);

    // Get initial Ktensors
    vector<cals::Ktensor> ktensor_init(n_ranks);
    const mxArray *arg = prhs[2];
    if (mxIsCell(arg))
    {
      auto i = 0;
      for (auto &ktensor : ktensor_init)
        ktensor = cals::Ktensor(mxGetKtensor(mxGetCell(arg, (mwIndex) i++), debug));
    } else if (mxIsChar(arg) && mxGetStdString(arg) == "random")
    {
      auto i = 0;
      for (auto &ktensor : ktensor_init)
      {
        ktensor = cals::Ktensor(static_cast<int>(ranks[i]), X.get_modes());
        ktensor.randomize();
      }
    } else
      throw std::string("Invalid type for initial guess specification.");

    auto cals_input(ktensor_init);

    cals::KtensorQueue cals_queue;
    for (auto &p : cals_input) cals_queue.emplace(p);

    // Call driver
    cout << "OpenMP threads: " << get_threads() << endl;
    auto report = cals::cp_cals(X, cals_queue, cals_params);

    // Return results
    if (nlhs >= 1)
    {
      mxArray *cell_array_ptr = mxCreateCellMatrix((mwSize) n_ranks, (mwSize) 1);
      for (int k = 0; k < n_ranks; k++)
      {
        mxArray *mat_ptr = mxSetKtensor(cals_input[k]);
        mxSetCell(cell_array_ptr, (mwIndex) k, mat_ptr);
      }
      plhs[0] = cell_array_ptr;
    }
    if (nlhs >= 2)
    {
      mxArray *cell_array_ptr = mxCreateCellMatrix((mwSize) n_ranks, (mwSize) 1);
      for (int k = 0; k < n_ranks; k++)
      {
        mxArray *mat_ptr = mxSetKtensor(ktensor_init[k]);
        mxSetCell(cell_array_ptr, (mwIndex) k, mat_ptr);
      }
      plhs[1] = cell_array_ptr;
    }
    cout << "Done." << endl;
  }
  catch (std::string &sExc)
  {
    cout << "Call to CALS threw an exception:" << endl;
    cout << "  " << sExc << endl;
  }
}

}



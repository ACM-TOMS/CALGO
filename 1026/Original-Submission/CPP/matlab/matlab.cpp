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

#include <iostream>

#ifndef NDEBUG
#include <numeric>
#endif

using std::cerr;
using std::endl;


cals::Tensor mxGetTensor(const mxArray *ptr, bool print)
{
  if (!mxIsClass(ptr, "tensor")) cerr << "Arg is not a tensor!" << endl;

  mxArray *data_field = mxGetField(ptr, 0, "data");
  mxArray *size_field = mxGetField(ptr, 0, "size");

  double *data = mxGetDoubles(data_field);
  DEBUG(int n_elements = mxGetNumberOfElements(data_field);)

  double *modes = mxGetDoubles(size_field);
  auto n_modes = mxGetNumberOfElements(size_field);

  std::vector<dim_t> modes_vector;
  modes_vector.reserve(n_modes);
  for (int i = 0; i < n_modes; i++)
    modes_vector.push_back(static_cast<dim_t>(modes[i]));

  assert(std::accumulate(modes_vector.cbegin(), modes_vector.cend(), 1lu, std::multiplies<>()) == n_elements);

  auto X = cals::Tensor(modes_vector, data);

  if (print) X.print();

  return X;
}

mxArray *mxSetKtensor(const cals::Ktensor &ktensor)
{
  const auto n_modes = ktensor.get_n_modes();

  // Create factor matrix array (cell array in matlab)
  mxArray *cell_array_ptr = mxCreateCellMatrix((mwSize) n_modes, (mwSize) 1);
  for (dim_t i = 0; i < n_modes; i++)
  {
    const auto rows = ktensor.get_factor(i).get_rows();
    const auto cols = ktensor.get_factor(i).get_cols();
    mxArray *mat_ptr = mxCreateDoubleMatrix((mwSize) rows, (mwSize) cols, mxREAL);
    double *mat = mxGetDoubles(mat_ptr);
    cals::Matrix(rows, cols, mat).copy(ktensor.get_factor(i));
    mxSetCell(cell_array_ptr, (mwIndex) i, mat_ptr);
  }

  // Create lambda array
  mxArray *lambda_ptr = mxCreateDoubleMatrix((mwSize) ktensor.get_components(), (mwSize) 1, mxREAL);
  double *lambda = mxGetDoubles(lambda_ptr);
  int i = 0;
  for (const auto &l : ktensor.get_lambda()) lambda[i++] = l;

  // Create Ktensor class by calling Ktensor constructor
  mxArray *lhs[1];
  mxArray *rhs[2] = {lambda_ptr, cell_array_ptr};
  mexCallMATLAB(1, lhs, 2, rhs, "ktensor");

  return lhs[0];
}

cals::Ktensor mxGetKtensor(const mxArray *ptr, bool print)
{
  if (!mxIsStruct(ptr)) cerr << "Arg is not a struct!" << endl;

  // Get lambda array
  mxArray *lambda_field = mxGetField(ptr, 0, "lambda");
  double *lambda = mxGetDoubles(lambda_field);
  auto rank = mxGetNumberOfElements(lambda_field);

  // Get Cell array of factor matrices
  mxArray *u_field = mxGetField(ptr, 0, "u");
  auto n_modes = mxGetNumberOfElements(u_field);

  // Create and initialize modes vector
  std::vector<dim_t> modes;
  modes.reserve(n_modes);
  for (int i = 0; i < n_modes; ++i)
  {
    mxArray *mat = mxGetCell(u_field, (mwIndex) i);
    modes.push_back(mxGetM(mat));
  }

  // Create Ktensor
  cals::Ktensor ktensor(static_cast<int>(rank), modes);

  // Copy lambda
  ktensor.set_lambda(lambda);

  // Copy factor matrices
  for (int i = 0; i < n_modes; ++i)
  {
    mxArray *mat = mxGetCell(u_field, (mwIndex) i);
    double *data_ptr = mxGetDoubles(mat);
    ktensor.set_factor(i, data_ptr);
  }

  if (print) ktensor.print();

  return ktensor;
}

std::string mxGetStdString(const mxArray *ptr)
{
  const mwSize str_len = mxGetNumberOfElements(ptr);
  char *c_str = new char[str_len + 1];
  int ret = mxGetString(ptr, c_str, str_len + 1);
  if (ret != 0) cerr << "mxGetString failed!" << endl;
  std::string str(c_str);
  delete[] c_str;
  return str;
}

std::vector<std::string> mxBuildArgList(int nargs, int offset, const mxArray *margs[])
{
  std::vector<std::string> args(static_cast<size_t>(nargs - offset));
  for (int i = 0; i < nargs - offset; i++)
  {
    const mxArray *arg = margs[i + offset];
    if (mxIsScalar(arg))
      args[i] = std::to_string(mxGetScalar(arg));
    else if (mxIsChar(arg))
      args[i] = mxGetStdString(arg);
    else
      cerr << "Unknown argument type for argument " << i + offset << endl;
  }
  return args;
}

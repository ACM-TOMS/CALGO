#ifndef CP_CALS_BENCH_MTTKRP_PLANC_H
#define CP_CALS_BENCH_MTTKRP_PLANC_H

#include <vector>
#include <fstream>
#include <iostream>

#include "tensor.h"
#include "matrix.h"
#include "timer.h"
#include "utils/utils.h"
#include "utils/mttkrp.h"

#include "experiments/bench_mttkrp.h"
#include "experiments/bench_utils.h"

#if WITH_PLANC

#include <armadillo>
#include "ddt.hpp"

#endif

namespace cals::bench
{
  BenchResult benchmark_mttkrp_planc(int rank, int mode, vector<dim_t> modes)
  {
#if WITH_PLANC
    arma::uvec dims;
    dims.resize(modes.size());
    for (auto i = 0; i < dims.size(); i++)
      dims[i] = static_cast<unsigned long long int>(modes[i]);

    planc::Tensor T(dims);
    T.rand();
    planc::NCPFactors ncp_factors(dims, rank, false);
    ncp_factors.randu(100);

//    cals::Tensor Tc(modes, T.m_data.data());
//    T.print();
//    Tc.print("Mine");
//    cals::Ktensor ktensor(rank, modes);
//    for (auto f = 0; f < dims.size(); f++)
//      for (auto i = 0; i < dims[f]; i++)
//        for (auto j = 0; j < rank; j++)
//          ktensor.get_factor(f)(i,j) = ncp_factors.factor(f).at(i, j);
//    for (auto f = 0; f < dims.size(); f++)
//      ncp_factors.factor(f).print("Factor " + std::to_string(f));
//    for (auto &f : ktensor.get_factors())
//      f.print("Mine");

//    std::vector<cals::Matrix> workspace;
//    workspace.emplace_back(Matrix(modes[0] * modes[1], rank));
//    cals::mttkrp::MttkrpParams params;
//    params.method = cals::mttkrp::MTTKRP_METHOD::AUTO;
//    auto &out_cals = cals::mttkrp::mttkrp(Tc, ktensor, workspace, tmode, params);

//    out.print();
//    out_cals.print();
//    exit(EXIT_SUCCESS);

    Matrix cache(2000, 2000);
    cache.randomize();

    cals::Matrix out(modes[mode], rank);

    vector<double> time;
    time.reserve(ITERATIONS);

    cals::Timer timer;

    int split;
    if (mode == 0)
      split = 0; // or 1 NOLINT(bugprone-branch-clone)
    else if (mode == 1)
      split = 0; // Only
    else if (mode == 2)
      split = 1; // Only
    else
      split = 0;

    auto ddt = DenseDimensionTree(T, ncp_factors, split);

    for (auto i = 0; i < ITERATIONS; i++)
    {
      double multittv_time, mttkrp_time;

      for (auto j = 0; j < cache.get_n_elements(); j++)
        cache[j] += 0.0001;

      timer.start();
      ddt.in_order_reuse_MTTKRP(mode, out.get_data(), true, multittv_time, mttkrp_time);
      timer.stop();
      time.push_back(timer.get_time());
    }

    return {*std::min_element(time.begin(), time.end()), 0, 0, 0, 0};
#endif

    return {0, 0, 0, 0, 0};
  }
}

#endif //CP_CALS_BENCH_MTTKRP_PLANC_H

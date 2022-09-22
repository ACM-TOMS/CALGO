#ifndef CP_CALS_BENCH_MTTKRP_CTF_H
#define CP_CALS_BENCH_MTTKRP_CTF_H

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

#if WITH_CTF
#include "mpi.h"
#include "ctf.hpp"
#endif

namespace cals::bench
{
  BenchResult benchmark_mttkrp_ctf(int rank, int mode, vector<dim_t> modes)
  {
#if WITH_CTF
    // CTF Initialization
    int lens[3] = {modes[0], modes[1], modes[2]};
    int lens_A[2] = {lens[0],rank};
    int lens_B[2] = {lens[1],rank};
    int lens_C[2] = {lens[2],rank};
    bool is_sparse = false;

    auto T = CTF::Tensor<double>(3, is_sparse, lens);

    auto A = CTF::Matrix<double>(lens_A[0], rank);
    auto B = CTF::Matrix<double>(lens_B[0], rank);
    auto C = CTF::Matrix<double>(lens_C[0], rank);


    Matrix cache(2000, 2000);
    cache.randomize();

    vector<double> time;
    time.reserve(ITERATIONS);

    cals::Timer timer;

    for (auto i = 0; i < ITERATIONS; i++)
    {
      T.fill_random(-1, 1);
      A.fill_random(-1, 1);
      B.fill_random(-1, 1);
      C.fill_random(-1, 1);
      for (auto j = 0; j < cache.get_n_elements(); j++)
        cache[j] += 0.0001;

      timer.start();
      if (mode == 0)
        A["ir"] = T["ijk"]*B["jr"]*C["kr"];
      else if (mode == 1)
        B["jr"] = T["ijk"]*A["ir"]*C["kr"];
      else if (mode == 2)
        C["kr"] = T["ijk"]*A["ir"]*B["jr"];
      timer.stop();
      time.push_back(timer.get_time());
    }

    return {*std::min_element(time.begin(), time.end()), 0, 0, 0, 0};
#endif

    return {0, 0, 0, 0, 0};
  }
}


#endif //CP_CALS_BENCH_MTTKRP_CTF_H

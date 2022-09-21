#ifndef CP_CALS_BENCH_MTTKRP_CALS_H
#define CP_CALS_BENCH_MTTKRP_CALS_H

#include <vector>
#include <fstream>
#include <iostream>

#include "tensor.h"
#include "matrix.h"
#include "timer.h"
#include "utils/utils.h"
#include "utils/mttkrp.h"

#include "experiments/bench_utils.h"

namespace cals::bench
{
  template<typename F>
  BenchResult benchmark_mttkrp_cals(int rank, int mode, vector<dim_t> modes, F &&f, cals::mttkrp::MttkrpParams &params)
  {
    auto T = Tensor(modes);
    auto cache = Matrix(2000, 2000);
    auto ktensor = Ktensor(rank, modes);

    std::sort(modes.begin(), modes.end(), std::greater<>());
    vector <Matrix> workspace;
    workspace.emplace_back(Matrix(modes[0] * modes[1], rank));

    if (params.cuda)
    {
#if CUDA_ENABLED
      T.allocate_cudata(T.get_n_elements());

      for (auto &f : ktensor.get_factors())
        f.allocate_cudata(f.get_n_elements());

      for (auto &w : workspace)
        w.allocate_cudata(w.get_n_elements());
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      abort();
#endif
    }

    vector<double> time;
    time.reserve(ITERATIONS);
    vector<double> f_time;
    f_time.reserve(ITERATIONS);
    vector<double> s_time;
    s_time.reserve(ITERATIONS);

    for (auto i = 0; i < ITERATIONS; i++)
    {
      T.randomize();
      workspace[0].randomize();
      ktensor.randomize();

      if (params.cuda)
      {
#if CUDA_ENABLED
      T.send_to_device();

      for (auto &f : ktensor.get_factors())
        f.send_to_device();

      for (auto &w : workspace)
        w.send_to_device();
#else
        std::cerr << "Not compiled with CUDA support" << std::endl;
        abort();
#endif
      }
      if (!params.cuda)
        for (dim_t j = 0; j < cache.get_n_elements(); j++)
          cache[j] += 0.0001;

      std::forward<decltype(f)>(f)(T, ktensor, workspace, mode, params);
      time.push_back(
          params.mttkrp_timers.timers[MttkrpTimers::TIMERS::MT_KRP].get_time() +
          params.mttkrp_timers.timers[MttkrpTimers::TIMERS::MT_GEMM].get_time() +
          params.mttkrp_timers.timers[MttkrpTimers::TIMERS::TS_GEMM].get_time() +
          params.mttkrp_timers.timers[MttkrpTimers::TIMERS::TS_GEMV].get_time());
      f_time.push_back(params.mttkrp_timers.timers[MttkrpTimers::TIMERS::MT_KRP].get_time() +
                       params.mttkrp_timers.timers[MttkrpTimers::TIMERS::TS_GEMM].get_time());
      s_time.push_back(params.mttkrp_timers.timers[MttkrpTimers::TIMERS::MT_GEMM].get_time() +
                       params.mttkrp_timers.timers[MttkrpTimers::TIMERS::TS_GEMV].get_time());
    }
    return {*std::min_element(time.begin(), time.end()),
            *std::min_element(f_time.begin(), f_time.end()),
            *std::min_element(s_time.begin(), s_time.end()),
            params.flops, params.memops};
  }
}


#endif //CP_CALS_BENCH_MTTKRP_CALS_H

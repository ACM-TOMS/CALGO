#ifndef CALS_BENCH_MTTKRP_H
#define CALS_BENCH_MTTKRP_H

#include <vector>
#include <fstream>
#include <iostream>

#include "tensor.h"
#include "matrix.h"
#include "timer.h"
#include "utils/utils.h"
#include "utils/mttkrp.h"

#include "experiments/bench_mttkrp_ctf.h"
#include "experiments/bench_mttkrp_cals.h"
#include "experiments/bench_mttkrp_planc.h"

namespace cals::bench
{
  void benchmark_mttkrp(const vector<dim_t> modes, int num_threads, int method, bool cuda_f=false)
  {
    set_threads(num_threads);
    std::cout << "Threads: " << get_threads() << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////
    //  Create CSV file name
    ////////////////////////////////////////////////////////////////////////////////////
    std::string folder = std::string(SOURCE_DIR) + "/data/" + std::string(CALS_BACKEND) + "/benchmark/";
    std::string is_cuda = (cuda_f) ? "CUDA_" : "";
    std::string modes_string = cals::utils::mode_string(modes);
    std::string file_name;
    file_name = folder
                + "benchmark_"
                + is_cuda
                + std::string(CALS_BACKEND) + "_"
                + modes_string + "_"
                + std::to_string(get_threads()) + "_"
                + std::to_string(method)
                + ".csv";

    initialize_csv(file_name);
    std::cout << "File name: " << file_name << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////
    //  Calculate target components
    ////////////////////////////////////////////////////////////////////////////////////
    vector<int> components;
    components.reserve(300);
    for (auto i = 1; i < 20; i += 1)
      components.push_back(i);
    for (auto i = 20; i < 100; i += 10)
      components.push_back(i);
    for (auto i = 100; i < 1000; i += 100)
      components.push_back(i);
    for (auto i = 1000; i <= 5000; i += 1000)
      components.push_back(i);

    int n_modes = static_cast<int>(modes.size());
    for (auto m = 0; m < n_modes; m++)
    {
      for (auto comp : components)
      {
        if (method <= cals::mttkrp::MTTKRP_METHOD::AUTO)
        {
          auto params = cals::mttkrp::MttkrpParams();
          params.method = static_cast<cals::mttkrp::MTTKRP_METHOD>(method);
          params.cuda = cuda_f;
          params.lut = mttkrp::read_lookup_table(modes, get_threads(), params.cuda);

          auto result = benchmark_mttkrp_cals(comp, m, modes, cals::mttkrp::mttkrp, params);
          append_csv(file_name, m, modes_string, comp, num_threads, result);
        } else if (method == EXTERNAL_MTTKRP_METHODS::CTF)
        {
          auto result = benchmark_mttkrp_ctf(comp, m, modes);
          append_csv(file_name, m, modes_string, comp, num_threads, result);
        } else if (method == EXTERNAL_MTTKRP_METHODS::PLANC)
        {
          auto result = benchmark_mttkrp_planc(comp, m, modes);
          append_csv(file_name, m, modes_string, comp, num_threads, result);
        } else
        {
          std::cerr << "MTTKRP method " << method << " not found. Exiting." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
}

#endif //CALS_BENCH_MTTKRP_H

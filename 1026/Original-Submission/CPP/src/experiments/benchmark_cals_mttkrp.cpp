#include <iostream>
#include <random>

#include "experiments/bench_mttkrp.h"

using cals::mttkrp::MTTKRP_METHOD;
using cals::mttkrp::MttkrpParams;
using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  cout << "=============================================================================" << endl;
  cout << "This executable performs experiments regarding MTTKRP for the paper." << endl;
  cout << "It benchmarks the different implementations of MTTKRP supported by CALS." << endl;
  cout << "The corresponding source file is located in `src/experiments/benchmark_cals_mttkrp.cpp`" << endl;
  cout << "The output is writen in CSV files in the `data` folder of the project." << endl;
  cout << "This executable accepts only one (mandatory) argument, the number of threads." << endl;
  cout << "=============================================================================" << endl;
  cout << endl;

  if (argc != 2) {
    cerr << "Not enough arguments. Give number of threads." << endl;
    cout << "USAGE: " << argv[0] << " <num_threads>" << endl;
    abort();
  } else {
    std::string arg = argv[1];
    if ((arg == "-h") || (arg == "--help")) {
      cout << "USAGE: " << argv[0] << " <num_threads>" << endl;
      return 0;
    }
  }

#if CUDA_ENABLED
  bool cuda_f = true;
#else
  bool cuda_f = false;
#endif

  int num_threads = std::atoi(argv[1]); // NOLINT(cert-err34-c)

  vector<MTTKRP_METHOD> methods_v;
  methods_v.push_back(MTTKRP_METHOD::MTTKRP);
  methods_v.push_back(MTTKRP_METHOD::TWOSTEP0);
  methods_v.push_back(MTTKRP_METHOD::TWOSTEP1);

  vector<vector<dim_t>> modes_v;
  modes_v.push_back({100, 100, 100});
  modes_v.push_back({200, 200, 200});
  modes_v.push_back({300, 300, 300});
  // modes_v.push_back({299, 301, 41});
  // modes_v.push_back({405, 136, 19});
  // modes_v.push_back({255, 281, 25});

  for (auto method : methods_v)
    for (auto &modes : modes_v)
      cals::bench::benchmark_mttkrp(modes, num_threads, method, cuda_f);
}

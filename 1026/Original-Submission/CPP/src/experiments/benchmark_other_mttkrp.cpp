#include <algorithm>
#include <iostream>

#include "experiments/bench_mttkrp.h"

#if WITH_CTF

#include "ctf.hpp"

#endif

using cals::bench::EXTERNAL_MTTKRP_METHODS;
using cals::mttkrp::MTTKRP_METHOD;
using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  cout << "=============================================================================" << endl;
  cout << "This executable performs experiments for the letter." << endl;
  cout << "It compares the implementation of MTTKRP in CALS with other implementations." << endl;
  cout << "The corresponding source file is located in `src/experiments/benchmark_other_mttkrp.cpp`" << endl;
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

#if WITH_CTF
  int mpi_rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
#endif

  int num_threads = std::atoi(argv[1]); // NOLINT(cert-err34-c)

  vector<int> methods_v;
  methods_v.push_back(MTTKRP_METHOD::AUTO);
  methods_v.push_back(EXTERNAL_MTTKRP_METHODS::CTF);
  methods_v.push_back(EXTERNAL_MTTKRP_METHODS::PLANC);

  vector<vector<dim_t>> modes_v;
  modes_v.push_back({100, 100, 100});
  modes_v.push_back({200, 200, 200});
  modes_v.push_back({300, 300, 300});

  for (auto method : methods_v)
    for (auto &modes : modes_v)
      cals::bench::benchmark_mttkrp(modes, num_threads, method);
}

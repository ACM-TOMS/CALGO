#include <cmath>
#include <iostream>
#include <numeric>

#include "matrix.h"
#include "timer.h"

const int iterations = 100;

using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {

  cout << "=============================================================================" << endl;
  cout << "This executable helps create the GEMM lines in figure 5." << endl;
  cout << "The corresponding source file is located in `src/experiments/peak_evaluator.cpp`" << endl;
  cout << "The output is writen in standard output." << endl;
  cout << "=============================================================================" << endl;

  cout << "USAGE: " << argv[0] << " 24 2000 32 1000" << endl
       << "for 24 threads, 2000 MHz AVX peak freq for 24 threads, 32 Flops/Cycle and matrix sizes 1000x1000." << endl;
  cout << endl;

  if (argc != 5) {
    cerr << "Not enough arguments!" << endl;
    cerr << "Give number of threads, AVX peak freq (in MHz) for those threads, Flops/Cycle and size of GEMM." << endl;
    cerr << "Example: " << argv[0] << " 24 2000 32 1000" << endl
         << "for 24 threads, 2000 MHz AVX peak freq for 24 threads, 32 Flops/Cycle and matrix sizes 1000x1000." << endl;
    abort();
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      cout << "Give number of threads, AVX peak freq (in MHz) for those threads, Flops/Cycle and size of GEMM." << endl;
      cout << "Example: " << argv[0] << " 24 2000 32 1000" << endl
           << "for 24 threads, 2000 MHz AVX peak freq for 24 threads, 32 Flops/Cycle and matrix sizes 1000x1000."
           << endl;
      return 0;
    }
  }

  typedef long long unsigned int lint;

  const lint num_threads = std::strtoull(argv[1], nullptr, 10);
  set_threads((int)num_threads);
  lint freq = std::strtoull(argv[2], nullptr, 10);
  freq *= 1000;
  freq *= 1000;
  const lint fpc = std::strtoull(argv[3], nullptr, 10);
  const lint tpp = num_threads * freq * fpc;

  lint size = std::strtoull(argv[4], nullptr, 10);
  cals::Matrix A(size, size), B(size, size), C(size, size);
  cals::Matrix cache(3000, 3000);
  vector<cals::Timer> vec_t(iterations, cals::Timer());

  A.randomize();
  B.randomize();
  C.randomize();
  const lint flops = 2llu * A.get_rows() * A.get_cols() * B.get_cols() + 2llu * C.get_rows() * C.get_cols();

  for (int i = 0; i < iterations; i++) {
    cache.randomize();
    auto &t = vec_t[i];
    t.start();
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A.get_rows(), B.get_cols(), A.get_cols(), 1.0, A.get_data(),
                A.get_col_stride(), B.get_data(), B.get_col_stride(), 0.0, C.get_data(), C.get_col_stride());
    t.stop();
  }

  vec_t.erase(vec_t.begin());

  cout << "Timings: " << endl;
  for (auto &i : vec_t)
    cout << i.get_time() << " ";
  cout << endl;

  std::vector<double> times(vec_t.size());
  for (auto i = 0lu; i < vec_t.size(); i++)
    times[i] = 1.0 * flops / vec_t[i].get_time() / tpp;

  auto sum = std::accumulate(times.begin(), times.end(), 0.0);
  auto mean = sum / times.size();
  double sq_sum = std::inner_product(times.begin(), times.end(), times.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / times.size() - mean * mean);
  auto max_t = *std::max_element(times.begin(), times.end());
  auto min_t = *std::min_element(times.begin(), times.end());

  std::nth_element(times.begin(), times.begin() + times.size() / 2, times.end());
  auto median = times[times.size() / 2];

  cout << endl;
  cout << "=============================================================" << endl;
  cout << "=============================================================" << endl;
  cout << "FLOPS:            " << flops << endl;
  cout << "-------------------------------------------------------------" << endl
       << "Min Performance:        " << min_t << endl
       << "Max Performance:        " << max_t << endl
       << "Mean Performance:       " << mean << endl
       << "Median Performance:     " << median << endl
       << "Stdev Performance:      " << stdev << endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Freq:    " << freq << endl;
  cout << "FPC:     " << fpc << endl;
  cout << "Size:    " << size << endl;
  cout << "Threads: " << num_threads << endl;
  cout << "{'mean': " << mean << ", 'median': " << median << ", 'std': " << stdev << "}" << endl;
  cout << "=============================================================" << endl;
  cout << "=============================================================" << endl;
  cout << endl;
}

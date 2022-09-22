#include "experiments/experiments_utils.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[]) {
  cout << "=============================================================================" << endl;
  cout << "This executable performs experiments related to the letter." << endl;
  cout << "The corresponding source file is located in `src/experiments/experiments_letter.cpp`" << endl;
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

  auto num_threads = static_cast<unsigned int>(std::strtol(argv[1], nullptr, 10));

  bool exp_defrag = true;

  if (exp_defrag) {
    vector<dim_t> modes = {200, 200, 200};
    cals::Tensor T(modes);
    T.randomize();

    auto ranks = generate_components(1, 20, 20);

    cals::CalsParams params;
    params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
    params.max_iterations = 1000;
    params.always_evict_first = true;
#if CUDA_ENABLED
    params.cuda = true;
#else
    params.cuda = false;
#endif

    run_cals(T, ranks, num_threads, params, true, "defrag");
  }
}
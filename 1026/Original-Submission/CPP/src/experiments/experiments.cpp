#include "experiments/experiments_utils.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[]) {
  cout << "=============================================================================" << endl;
  cout << "This executable performs the experiments in the paper." << endl;
  cout << "The corresponding source file is located in `src/experiments/experiments.cpp`" << endl;
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

  auto num_threads = static_cast<unsigned int>(std::strtoul(argv[1], nullptr, 10));

  bool experiment_6_1_1 = true;
  bool experiment_6_1_2_and_6_3 = true;
  bool experiment_6_2 = true;

  // Mock run to warmup everything.
  {
    auto components = generate_components(1, 2, 2);
    std::sort(components.begin(), components.end());

    cals::CalsParams params;
    params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
    params.max_iterations = 5;
    params.buffer_size = 10;
    params.tol = 1e-4;
#if CUDA_ENABLED
    params.cuda = true;
#else
    params.cuda = false;
#endif

    vector<dim_t> modes = {20, 20, 20};
    cals::Tensor T(modes);
    T.randomize();
    params.mttkrp_lut = cals::mttkrp::read_lookup_table({100, 100, 100}, num_threads, params.cuda);
    compare_als_cals(T, components, num_threads, params);
  }

  if (experiment_6_1_1) {
    // Set the parameters of the execution
    cals::CalsParams params;
    params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO; // Let the Lookup tables decide the best method.
    params.force_max_iter = true;                             // Force all models to reach max_iterations.
    params.max_iterations = 50;                               // Set maximum number of iterations.
#if CUDA_ENABLED
    params.cuda = true;
#else
    params.cuda = false;
#endif

    for (int comp = 1; comp <= 20; comp++) {
      auto components = generate_components(comp, comp, 20); // Generate 20 models for every number of components

      vector<vector<dim_t>> modes_v;
      modes_v.push_back({100, 100, 100});
      modes_v.push_back({200, 200, 200});
      modes_v.push_back({300, 300, 300});

      for (auto &modes : modes_v) {
        if (num_threads == 1 && modes[0] == 100 && modes[1] == 100 && modes[2] == 100)
          params.buffer_size = 90; // Selecting the optimal buffer size for MTTKRP (as explained in the paper)
        else
          params.buffer_size = static_cast<dim_t>(20) * comp; // Else, select the largest buffer size to fit all models
        int nt = (num_threads == 1) ? 1 : 24;
        params.mttkrp_lut = cals::mttkrp::read_lookup_table(modes, nt, params.cuda);
        cals::Tensor T(25, modes);
        T.randomize();
        compare_als_cals(T, components, num_threads, params, "speedup_" + std::to_string(comp));
      }
    }
  }

  if (experiment_6_1_2_and_6_3) {
    auto components = generate_components(1, 20, 20);
    std::sort(components.begin(), components.end());

    // Set the parameters of the execution
    cals::CalsParams params;
    params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO; // Let the Lookup tables decide the best method.
    params.force_max_iter = true;                             // Force all models to reach max_iterations.
    params.max_iterations = 50;                               // Set maximum number of iterations.
#if CUDA_ENABLED
    params.cuda = true;
#else
    params.cuda = false;
#endif

    vector<vector<dim_t>> modes_v;
    modes_v.push_back({100, 100, 100});
    modes_v.push_back({200, 200, 200});
    modes_v.push_back({300, 300, 300});

    for (auto &modes : modes_v) {
      if (num_threads == 1 && modes[0] == 100 && modes[1] == 100 && modes[2] == 100)
        params.buffer_size = 90; // Selecting the optimal buffer size for MTTKRP (as explained in the paper)
      else
        params.buffer_size = 4200; // Else, select the largest buffer size to fit all models
      int nt = (num_threads == 1) ? 1 : 24;
      params.mttkrp_lut = cals::mttkrp::read_lookup_table(modes, nt, params.cuda);
      cals::Tensor T(modes);
      T.randomize();
      compare_als_cals(T, components, num_threads, params);
    }
  }

  if (experiment_6_2) {
    auto components = generate_components(1, 20, 20);
    std::sort(components.begin(), components.end());

    cals::CalsParams params;
    params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
    params.max_iterations = 1000; // Set the maximum number of iterations.
    params.buffer_size = 4200;    // Set the buffer size to the sum of all number of components.
    params.tol = 1e-6;            // Set the tolerance.
#if CUDA_ENABLED
    params.cuda = true;
#else
    params.cuda = false;
#endif

    string file_name = string(SOURCE_DIR) + "/data/fluorescence_cancer_UD.txt";
    cals::Tensor T_cancer(file_name);
    int nt = (num_threads == 1) ? 1 : 24;
    params.mttkrp_lut = cals::mttkrp::read_lookup_table({299, 301, 41}, nt, params.cuda);
    compare_als_cals(T_cancer, components, num_threads, params, "", true);

    file_name = string(SOURCE_DIR) + "/data/eemdata.txt";
    cals::Tensor T_eemdata(file_name);
    params.mttkrp_lut = cals::mttkrp::read_lookup_table({405, 136, 19}, nt, params.cuda);
    compare_als_cals(T_eemdata, components, num_threads, params, "", true);
  }
}

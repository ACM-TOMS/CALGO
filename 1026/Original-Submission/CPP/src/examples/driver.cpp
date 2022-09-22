#include <iostream>
#include <random>
#include <sstream>

#include "als.h"
#include "cals.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

void show_usage(char *text) {
  cout << "USAGE: " << text << "[Options] " << endl;
  cout << endl;
  cout << "-n --nthreads THREADS" << endl;
  cout << "        Set the number of threads (defaults to OMP_NUM_THREADS, if set, or else to 10)" << endl;
  cout << "-c --components MIN:MAX:COPIES" << endl;
  cerr << "        Set the MIN and MAX number of components per model, as well as the randomly initialized COPIES of "
          "each model. (defaults to 1:10:10)"
       << endl;
  cout << "-t --tensor DIM0-DIM1-DIM2" << endl;
  cerr << "        Set the DIMensions of the target tensor." << endl;
  cout << endl;
}

void string_split(const std::string &str, vector<dim_t> &vect, char delimiter = ' ') {
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter))
    vect.push_back(std::strtoul(token.c_str(), nullptr, 10));
}

int main(int argc, char *argv[]) {
  cout << "GENERAL INFO" << endl;
  cout << "=======================================================================================" << endl;
  cout << "The `driver` executable demonstrates how CALS runs without specifying extra arguments. " << endl;
  cout << "The associated source file is located in `src/examples/driver.cpp`." << endl;
  cout << "=======================================================================================" << endl;
  cout << endl;

  /////////////////////////////////////////////////////////////////////////////////
  // Parse Arguments.
  /////////////////////////////////////////////////////////////////////////////////
  cout << "Parsing Arguments..." << endl;
  vector<dim_t> modes{210, 210, 210};
  int min_components{1};
  int max_components{10};
  int mul_components{10};
  unsigned int num_threads{10};
  bool f_threads = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } else if ((arg == "-n") || (arg == "--nthreads")) {
      if (i + 1 < argc) {
        num_threads = static_cast<unsigned int>(std::strtol(argv[++i], nullptr, 10));
        cout << "Threads specified: " << num_threads << endl;
        f_threads = true;
      } else {
        cerr << "--nthreads option requires one argument." << endl;
        return 1;
      }
    } else if ((arg == "-c") || (arg == "--components")) {
      if (i + 1 < argc) {
        vector<dim_t> split_components;
        string_split(string(argv[++i]), split_components, ':');
        if (split_components.size() != 3) {
          cerr << "--components option requires one argument of the form MIN:MAX:COPIES." << endl;
          return 1;
        } else {
          min_components = static_cast<int>(split_components[0]);
          max_components = static_cast<int>(split_components[1]);
          mul_components = static_cast<int>(split_components[2]);
        }
      } else {
        cerr << "--components option requires one argument of the form MIN:MAX:COPIES." << endl;
        return 1;
      }
    } else if ((arg == "-t") || (arg == "--tensor")) {
      if (i + 1 < argc) {
        vector<dim_t> split_tensor;
        string_split(string(argv[++i]), split_tensor, '-');

        if (split_tensor.size() < 3) {
          cerr << "--tensor option requires one argument of the form DIM0-DIM1-DIM2." << endl;
          return 1;
        } else {
          modes[0] = split_tensor[0];
          modes[1] = split_tensor[1];
          modes[2] = split_tensor[2];
          cout << "Tensor dimensions specified: ";
          for (auto &s : split_tensor)
            cout << s << " ";
          cout << endl;
        }
      } else {
        cerr << "--tensor option requires one argument of the form DIM0-DIM1-DIM2." << endl;
        return 1;
      }
    } else {
      cerr << "Unrecognized argument " << arg << endl;
      show_usage(argv[0]);
      return 1;
    }
  }
  cout << "Done Parsing Arguments..." << endl;
  cout << endl;

  /////////////////////////////////////////////////////////////////////////////////
  // Set the number of threads.
  /////////////////////////////////////////////////////////////////////////////////
  if (!f_threads) {
    const char *omp_num_threads = std::getenv("OMP_NUM_THREADS");
    if (omp_num_threads != nullptr) {
      num_threads = static_cast<unsigned int>(std::strtol(omp_num_threads, nullptr, 10));
      cout << "OMP_NUM_THREADS: " << num_threads << endl;
    } else {
      num_threads = 10;
      cout << "OMP_NUM_THREADS not set. Using default " << num_threads << " threads." << endl;
    }
  }

  set_threads(static_cast<int>(num_threads));

  /////////////////////////////////////////////////////////////////////////////////
  // Create target Tensor.
  /////////////////////////////////////////////////////////////////////////////////
  cals::Tensor X(modes);
  X.randomize();

  /////////////////////////////////////////////////////////////////////////////////
  // Create input Ktensors.
  /////////////////////////////////////////////////////////////////////////////////

  // Create a vector with all the numbers of components we want to run CP-ALS for.
  vector<int> components;
  // Components from min_components up to (and including) max_components
  for (auto comp = min_components; comp <= max_components; comp++)
    for (auto cp = 0; cp < mul_components; cp++) // mul_components copies of each number of components
      components.push_back(comp);

  // Create a vector containing the Ktensors (models) we want to fit to the target Tensor.
  vector<cals::Ktensor> cals_input(components.size());
  auto i = 0;
  for (auto &ktensor : cals_input) {
    ktensor = cals::Ktensor(components[i++], modes);
    ktensor.randomize();
  }
  cout << "Created " << components.size() << " Ktensors." << endl;
  cout << "Ktensor components: ";
  for (auto &r : components)
    cout << r << " ";
  cout << endl;

  // Copy the same Ktensors to be used as input to regular ALS to a new vector.
  auto als_input(cals_input);

  /////////////////////////////////////////////////////////////////////////////////
  // Run CALS.
  /////////////////////////////////////////////////////////////////////////////////

  // Create the input to CALS, as a Queue of references to the Ktensors created.
  cals::KtensorQueue cals_queue;
  for (auto &p : cals_input)
    cals_queue.emplace(p);

  // Set the CALS parameters (if necessary).
  cals::CalsParams cals_params;
  cals_params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
  cals_params.update_method = cals::update::UPDATE_METHOD::UNCONSTRAINED;
  // cals_params.force_max_iter = true;
  cals_params.max_iterations = 1000;
  cals_params.tol = 1e-5;

  cals_params.buffer_size = std::accumulate(components.cbegin(), components.cend(), static_cast<dim_t>(0));
  cals_params.print();

  cout << "Running CALS..." << endl;
  cals::Timer cals_timer;
  cals_timer.start();
  auto cals_report = cp_cals(X, cals_queue, cals_params);
  cals_timer.stop();
  cout << "Finished CALS." << endl;
  cout << endl;

  /////////////////////////////////////////////////////////////////////////////////
  // Run a for-loop of regular ALS for comparison.
  /////////////////////////////////////////////////////////////////////////////////

  // Set the ALS parameters (if necessary).
  cals::AlsParams als_params;
  als_params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
  als_params.update_method = cals::update::UPDATE_METHOD::UNCONSTRAINED;
  // als_params.force_max_iter = true;
  als_params.max_iterations = 1000;
  als_params.tol = 1e-5;

  // Suppress the warning that no Lookup table was found, to not clutter the standard output.
  als_params.suppress_lut_warning = true;
  als_params.print();

  cout << "Running ALS..." << endl;
  cals::Timer als_timer;
  als_timer.start();
  for (auto &kt : als_input)
    auto als_report = cp_als(X, kt, als_params);
  als_timer.stop();
  cout << "Finished ALS." << endl;
  cout << endl;

  /////////////////////////////////////////////////////////////////////////////////
  // Print timings.
  /////////////////////////////////////////////////////////////////////////////////

  cout << "======================================================================" << endl;
  cout << "ALS time: " << als_timer.get_time() << endl;
  cout << "CALS time: " << cals_timer.get_time() << endl;
  cout << "Speedup: " << als_timer.get_time() / cals_timer.get_time() << endl;
}

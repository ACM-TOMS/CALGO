#include "experiments/experiments_utils.h"

#include <cmath>
#include <iomanip>
#include <random>

#include "als.h"

using cals::Ktensor;
using cals::Tensor;
using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<int> generate_components(int min, int max, int copies) {
  vector<int> components;
  for (auto comp = min; comp <= max; comp++)
    for (auto cp = 0; cp < copies; cp++)
      components.push_back(comp);
  return components;
}

vector<cals::AlsReport>
regular_als(const Tensor &X, vector<Ktensor> &ktensor_v, cals::CalsParams &params, bool openmp_f = false) {
  vector<cals::AlsReport> reports(ktensor_v.size());

  cals::AlsParams als_params;
  als_params.update_method = params.update_method;
  als_params.mttkrp_method = params.mttkrp_method;
  als_params.mttkrp_lut = params.mttkrp_lut;
  als_params.max_iterations = params.max_iterations;
  als_params.force_max_iter = params.force_max_iter;
  als_params.tol = params.tol;
  als_params.line_search = params.line_search;
  als_params.line_search_interval = params.line_search_interval;
  als_params.line_search_step = params.line_search_step;
  als_params.cuda = params.cuda;

  cals::Timer als_time;
  als_time.start();
  if (openmp_f) {
    // cout << "MKL_DYNAMIC:     " << mkl_get_dynamic() << endl;
    // cout << "OMP_NUM_THREADS: " << omp_get_max_threads() << endl;
    reports = cals::cp_omp_als(X, ktensor_v, als_params);
  } else {
    for (auto i = 0lu; i < ktensor_v.size(); i++)
      reports[i] = cals::cp_als(X, ktensor_v[i], als_params);
  }
  als_time.stop();
  std::cout << ((openmp_f) ? "O" : " ") << "ALS Computation time: " << als_time.get_time() << std::endl;

  return reports;
}

cals::CalsReport concurrent_als(const Tensor &X, cals::KtensorQueue &kt_queue, cals::CalsParams &params) {
  cals::Timer cals_time;
  cals_time.start();

  // Fit the models using Concurrent ALS (CALS).
  auto report = cals::cp_cals(X, kt_queue, params);

  cals_time.stop();
  std::cout << "CALS Computation time: " << cals_time.get_time() << std::endl;

  return report;
}

void compare_als_cals(const Tensor &X,
                      vector<int> &components,
                      unsigned int num_threads,
                      cals::CalsParams &params,
                      const std::basic_string<char> &file_suffix,
                      bool reproducible) {
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "START ALS v OMP ALS v CALS" << endl;

  std::mt19937 reproducible_generator(0);
  ////////////////////////////////////////////////////////////////////////
  // Create Tensor, Ktensors and run experiments
  ////////////////////////////////////////////////////////////////////////
  set_threads(num_threads);
  cout << "Number of threads set: " << get_threads() << endl;

  auto const modes = X.get_modes();
  auto const n_ktensors = components.size();
  cout << "Modes: " << cals::utils::mode_string(modes) << endl;

  vector<Ktensor> ktensor_vector(n_ktensors);
  auto i = 0;
  for (auto &ktensor : ktensor_vector) {
    ktensor = Ktensor(components[i++], modes);
    if (reproducible) {
      std::uniform_real_distribution<double> dist(-1.0, 1.0);
      ktensor.fill(
          (function<double()> &&)[&dist, &reproducible_generator]()->double { return dist(reproducible_generator); });
    } else
      ktensor.randomize();
  }

  cout << "----------------------" << endl;
  cout << "Starting CALS" << endl;
  auto cals_input(ktensor_vector);
  cals::KtensorQueue cals_queue;
  for (auto p = 0lu; p < n_ktensors; p++)
    cals_queue.emplace(cals_input[p]);
  auto cals_report = concurrent_als(X, cals_queue, params);

  auto als_omp_report = vector<cals::AlsReport>();
  auto als_omp_input(ktensor_vector);
  if (num_threads != 1) {
    cout << "----------------------" << endl;
    cout << "Starting OMP ALS" << endl;
    auto enable_omp = true;
    als_omp_report = regular_als(X, als_omp_input, params, enable_omp);
  }

  cout << "----------------------" << endl;
  cout << "Starting ALS" << endl;
  auto als_input(ktensor_vector);
  auto enable_omp = false;
  auto als_report = regular_als(X, als_input, params, enable_omp);

  ////////////////////////////////////////////////////////////////////////
  // Check whether results are the same across all 3 methods
  ////////////////////////////////////////////////////////////////////////
  for (auto p = 0lu; p < n_ktensors; p++) {
    auto model_norm_diff = std::fabs(als_input[p].get_approximation_error() - cals_input[p].get_approximation_error());

    if (std::isnan(model_norm_diff) || std::isnan(cals_input[p].get_factor(0)[0]) ||
        std::isnan(cals_input[p].get_factor(1)[0]) || std::isnan(cals_input[p].get_factor(2)[0]))
      std::cerr << "NaN values in solutions" << std::endl;
    if (model_norm_diff > MODEL_DIFF_ACC)
      std::cerr << std::scientific << "Ktensor: " << p << " ALS v CALS: " << std::setw(9) << model_norm_diff
                << " ALS: " << std::setw(9) << als_input[p].get_approximation_error() << " CALS: " << std::setw(9)
                << cals_input[p].get_approximation_error() << endl;

    if (num_threads != 1) {
      auto model_norm_diff_1 =
          std::fabs(als_omp_input[p].get_approximation_error() - cals_input[p].get_approximation_error());
      if (std::isnan(model_norm_diff_1))
        std::cerr << "NaN values in solutions" << std::endl;
      if (model_norm_diff_1 > MODEL_DIFF_ACC)
        std::cerr << std::scientific << "Ktensor: " << p << " OALS v CALS: " << std::setw(9) << model_norm_diff_1
                  << " OALS: " << std::setw(9) << als_omp_input[p].get_approximation_error()
                  << " CALS: " << std::setw(9) << cals_input[p].get_approximation_error() << endl;
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Print results to file
  ////////////////////////////////////////////////////////////////////////
  string str_mode = cals::utils::mode_string(modes);
  string str_n_threads = std::to_string(cals_report.n_threads);
  string str_blas_vendor = CALS_BACKEND;
  string str_cuda = (params.cuda) ? "CUDA_" : "";
  string dir = string(SOURCE_DIR) + "/data/" + str_blas_vendor + "/";
  string suffix = ".csv";
  if (!file_suffix.empty())
    suffix = "_" + file_suffix + ".csv";

  string als_file_name = dir + "ALS_" + str_cuda + str_blas_vendor + "_" + str_mode + "_" + str_n_threads + suffix;
  als_report[0].print_header(als_file_name);
  for (auto &rm : als_report)
    rm.print_to_file(als_file_name);

  if (num_threads != 1) {
    string als_omp_file_name =
        dir + "ALS_" + "OMP_" + str_cuda + str_blas_vendor + "_" + str_mode + "_" + str_n_threads + suffix;
    als_report[0].print_header(als_omp_file_name);
    for (auto &rm : als_omp_report)
      rm.print_to_file(als_omp_file_name);
  }

  string cals_file_name = dir + "CALS_" + str_cuda + str_blas_vendor + "_" + str_mode + "_" + str_n_threads + suffix;
  cals_report.print_header(cals_file_name);
  cals_report.print_to_file(cals_file_name);

  // Print all Ktensor ids, components, errors and iterations
  string ktensor_file_name =
      dir + "Ktensors_" + str_cuda + str_blas_vendor + "_" + str_mode + "_" + str_n_threads + suffix;
  auto file = std::ofstream(ktensor_file_name, std::ios::out);
  file << "KTENSOR_ID;KTENSOR_RANK;ERROR;ITERS" << endl;
  for (auto &k : cals_input)
    file << k.get_id() << ";" << k.get_components() << ";" << k.get_approximation_error() << ";" << k.get_iters()
         << endl;

  cout << endl;
  cout << "STOP ALS v OMP ALS v CALS" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
}

void run_cals(const Tensor &X,
              vector<int> &components,
              unsigned int num_threads,
              cals::CalsParams &params,
              bool print_header,
              const std::basic_string<char> &file_suffix) {
  cout << "============================" << endl;
  cout << "START Single CALS call" << endl;

  set_threads(num_threads);
  cout << "Threads: " << get_threads() << endl;

  auto const modes = X.get_modes();
  cout << "Modes: " << cals::utils::mode_string(modes) << endl;

  vector<Ktensor> ktensor_vector(components.size());
  auto i = 0;
  for (auto &ktensor : ktensor_vector) {
    ktensor = Ktensor(components[i++], modes);
    ktensor.randomize();
  }
  auto cals_input(ktensor_vector);

  cals::KtensorQueue cals_queue;
  for (auto &p : cals_input)
    cals_queue.emplace(p);
  auto cals_report = concurrent_als(X, cals_queue, params);

  // Print results to file
  string str_mode = cals::utils::mode_string(modes);
  string str_n_threads = std::to_string(cals_report.n_threads);
  string str_blas_vendor = CALS_BACKEND;
  string dir = string(SOURCE_DIR) + "/data/" + str_blas_vendor + "/";
  string str_cuda = (params.cuda) ? "CUDA_" : "";
  string suffix = ".csv";
  if (!file_suffix.empty())
    suffix = "_" + file_suffix + ".csv";
  string cals_file_name = dir + "CALS_" + str_cuda + str_blas_vendor + "_" + str_mode + "_" + str_n_threads + suffix;

  cals_report.output_file_name = cals_file_name;
  if (print_header)
    cals_report.print_header(cals_file_name);
  cals_report.print_to_file(cals_file_name);

  cout << cals_report.output_file_name << endl;
  cout << "END Single CALS call" << endl;
  cout << "==========================" << endl;
}

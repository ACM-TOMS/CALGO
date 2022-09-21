#ifndef CALS_EXPERIMENTS_UTILS_H
#define CALS_EXPERIMENTS_UTILS_H

#include <cstdlib>

#include "cals.h"

#define MODEL_DIFF_ACC 1e-05

using std::vector;

vector<int> generate_components(int min, int max, int copies);

void compare_als_cals(const cals::Tensor &X,
                      vector<int> &components,
                      unsigned int num_threads,
                      cals::CalsParams &params,
                      const std::basic_string<char> &file_suffix = "",
                      bool reproducible = false);

void run_cals(const cals::Tensor &X,
              vector<int> &components,
              unsigned int num_threads,
              cals::CalsParams &params,
              bool print_header,
              const std::basic_string<char> &file_suffix = "");

#endif //CALS_EXPERIMENTS_UTILS_H

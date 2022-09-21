#include "gtest/gtest.h"

#include <cals.h>
#include <cmath>
#include <random>

#include "als.h"

#define MODEL_DIFF_ACC 1e-08

class CalsParametersTests : public ::testing::TestWithParam<bool> {};

TEST_P(CalsParametersTests, Correctness) {
  // parameter = GetParam()
  // You can pass a CalsParams object or a tuple of values and access them with std::get<1>(GetParam())

  vector<int> ranks;
  for (auto rank = 1; rank <= 10; rank++)
    for (auto copies = 0; copies < 30; copies++)
      ranks.push_back(rank);
  //  std::sort(ranks.begin(), ranks.end());
  std::default_random_engine generator;
  std::shuffle(ranks.begin(), ranks.end(), generator);

  cals::CalsParams cals_params;
  cals_params.mttkrp_method = cals::mttkrp::MTTKRP_METHOD::AUTO;
  cals_params.max_iterations = 100;
  cals_params.tol = 1e-4;
  cals_params.buffer_size = 30;
  cals_params.line_search = true;
  cals_params.line_search_interval = 5;
  cals_params.line_search_step = 0.9;
  cals_params.cuda = GetParam();

  cals::AlsParams als_params;
  als_params.mttkrp_method = cals_params.mttkrp_method;
  als_params.max_iterations = cals_params.max_iterations;
  als_params.line_search = cals_params.line_search;
  als_params.line_search_interval = cals_params.line_search_interval;
  als_params.line_search_step = cals_params.line_search_step;
  als_params.tol = cals_params.tol;
  als_params.cuda = cals_params.cuda;
  als_params.suppress_lut_warning = true;

  cals::Tensor T(5, {9, 4, 3});

  auto const modes = T.get_modes();

  int n_ktensors = static_cast<int>(ranks.size());
  vector<cals::Ktensor> ktensor_vector(n_ktensors);
  auto i = 0;
  for (auto &ktensor : ktensor_vector) {
    ktensor = cals::Ktensor(ranks[i++], modes);
    ktensor.randomize();
  }
  auto als_input(ktensor_vector);
  auto als_omp_input(ktensor_vector);
  auto cals_input(ktensor_vector);

  cals::KtensorQueue cals_queue;
  for (auto p = 0; p < n_ktensors; p++)
    cals_queue.emplace(cals_input[p]);
  cals::cp_cals(T, cals_queue, cals_params);

  for (auto p = 0; p < n_ktensors; p++)
    cals::cp_als(T, als_input[p], als_params);

  cals::cp_omp_als(T, als_omp_input, als_params);

  for (auto p = 0; p < n_ktensors; p++) {
    EXPECT_NEAR(als_input[p].get_approximation_error(), cals_input[p].get_approximation_error(), MODEL_DIFF_ACC);
    EXPECT_NEAR(als_omp_input[p].get_approximation_error(), cals_input[p].get_approximation_error(), MODEL_DIFF_ACC);
  }
};

INSTANTIATE_TEST_SUITE_P(CalsCorrectnessTests, CalsParametersTests, ::testing::Values(false));

#if CUDA_ENABLED
INSTANTIATE_TEST_SUITE_P(CalsCUDACorrectnessTests, CalsParametersTests, ::testing::Values(true));
#endif

int main(int argc, char **argv) {
  set_threads(4);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

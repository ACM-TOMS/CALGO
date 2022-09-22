#include "gtest/gtest.h"

#include <cmath>
#include <random>

#include "als.h"

#define MODEL_DIFF_ACC 1e-09

TEST(Als, ComputeCorrectResult3D) {
  vector<dim_t> modes = {9, 4, 2};
  cals::Tensor X(5, modes);
  cals::Ktensor k(5, modes);

  cals::mttkrp::MTTKRP_METHOD mttkrp_variants[] = {
      cals::mttkrp::MTTKRP_METHOD::MTTKRP, cals::mttkrp::MTTKRP_METHOD::TWOSTEP0, cals::mttkrp::MTTKRP_METHOD::TWOSTEP1,
      cals::mttkrp::MTTKRP_METHOD::AUTO};

  for (auto p = 0; p < 20; p++) {
    double errors[4] = {0.0};
    auto index = 0;
    k.randomize();

    for (auto &method : mttkrp_variants) {
      cals::AlsParams params;
      params.max_iterations = 100;
      params.mttkrp_method = method;
      params.line_search = false;
      params.line_search_interval = 5;
      params.line_search_step = 0.9;
      params.suppress_lut_warning = true;
#if CUDA_ENABLED
      std::random_device rd;
      std::default_random_engine eng(rd());
      std::uniform_int_distribution<int> distr(0, 1);
      if (distr(eng) == 0)
        params.cuda = true;
      else
        params.cuda = false;
#else
      params.cuda = false;
#endif
      cals::Ktensor k_copy(k);

      auto report = cp_als(X, k_copy, params);

      // Get the slow error calculated by reconstructing the tensor, elementwise subtracting and norming
      auto approx_X = k_copy.to_tensor();
      for (dim_t i = 0; i < X.get_n_elements(); i++)
        approx_X[i] = X[i] - approx_X[i];
      auto slow_error = approx_X.norm();
      errors[index++] = slow_error;

      EXPECT_FALSE(std::isnan(slow_error));
      EXPECT_LT(slow_error, 50);
    }
    for (auto e : errors)
      EXPECT_NEAR(e, errors[0], 1e-8);
  }
}

TEST(Als, ComputeCorrectResultConstrained3D) {
  vector<dim_t> modes = {18, 17, 16};
  cals::Tensor X(5, modes);
  cals::Ktensor k(5, modes);

  cals::mttkrp::MTTKRP_METHOD mttkrp_variants[] = {
      cals::mttkrp::MTTKRP_METHOD::MTTKRP, cals::mttkrp::MTTKRP_METHOD::TWOSTEP0, cals::mttkrp::MTTKRP_METHOD::TWOSTEP1,
      cals::mttkrp::MTTKRP_METHOD::AUTO};

  for (auto p = 0; p < 20; p++) {
    double errors[4] = {0.0};
    auto index = 0;
    k.randomize();

    for (auto &method : mttkrp_variants) {
      cals::AlsParams params;
      params.max_iterations = 100;
      params.mttkrp_method = method;
      params.update_method = cals::update::UPDATE_METHOD::NNLS;
      params.suppress_lut_warning = true;
      cals::Ktensor k_copy(k);

      auto report = cp_als(X, k_copy, params);

      for (auto &f : k_copy.get_factors())
        for (auto i = 0; i < f.get_n_elements(); i++)
          EXPECT_GE(f[i], 0.0);

      // Get the slow error calculated by reconstructing the tensor, elementwise subtracting and norming
      auto approx_X = k_copy.to_tensor();
      for (dim_t i = 0; i < X.get_n_elements(); i++)
        approx_X[i] = X[i] - approx_X[i];
      auto slow_error = approx_X.norm();
      errors[index++] = slow_error;

      EXPECT_FALSE(std::isnan(slow_error));
      EXPECT_LT(slow_error, 50);
    }
    for (auto e : errors)
      EXPECT_NEAR(e, errors[0], 1e-8);
  }
}

TEST(Als, ComputeCorrectResult4D) {
  cals::Tensor X(5, {3, 3, 3, 3});
  cals::Ktensor k(7, {3, 3, 3, 3});
  k.randomize();

  cals::AlsParams params;
  params.max_iterations = 100;
  params.suppress_lut_warning = true;
  auto report = cp_als(X, k, params);

  // Get the slow error calculated by reconstructing the tensor, elementwise subtracting and norming
  auto approx_X = k.to_tensor();
  for (dim_t i = 0; i < X.get_n_elements(); i++)
    approx_X[i] = X[i] - approx_X[i];
  auto slow_error = approx_X.norm();

  EXPECT_FALSE(std::isnan(slow_error));
  EXPECT_LT(slow_error, 1e-1);
}

TEST(Als, ComputeCorrectError) {
  cals::Tensor X(5, {9, 3, 2});
  cals::Ktensor k(5, {9, 3, 2});
  k.randomize();

  cals::AlsParams params;
  params.max_iterations = 3;
  params.suppress_lut_warning = true;
  auto report = cp_als(X, k, params);

  // Get the fast error reported by the algorithm
  auto fast_error = k.get_approximation_error();

  // Get the slow error calculated by reconstructing the tensor, elementwise subtracting and norming
  auto approx_X = k.to_tensor();
  for (dim_t i = 0; i < X.get_n_elements(); i++)
    approx_X[i] = X[i] - approx_X[i];
  auto slow_error = approx_X.norm();

  EXPECT_NEAR(fast_error, slow_error, 1e-10);
}

int main(int argc, char **argv) {
  set_threads(4);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
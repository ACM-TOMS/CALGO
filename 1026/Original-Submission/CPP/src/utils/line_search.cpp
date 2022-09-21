#include "utils/line_search.h"

#include <iostream>

#include "utils/error.h"
#include "utils/utils.h"

using cals::Ktensor;
using cals::Matrix;
using std::vector;

namespace cals::ls {
void line_search(Ktensor &ktensor,
                 Ktensor &ls_ktensor,
                 Matrix &krp_workspace,
                 Matrix &ten_workspace,
                 vector<Matrix> &gramians,
                 const Tensor &X,
                 double X_norm,
                 LineSearchParams &params) {
  for (dim_t n = 0; n < ls_ktensor.get_n_modes(); n++) {
    const auto &curr_factor = ktensor.get_factor(n);
    auto &old_factor = ls_ktensor.get_factor(n);

    krp_workspace.resize(curr_factor.get_rows(), curr_factor.get_cols());
    for (dim_t i = 0; i < old_factor.get_n_elements(); i++)
      krp_workspace[i] = curr_factor[i] - old_factor[i];

    old_factor.copy(curr_factor);
    cblas_daxpy(old_factor.get_n_elements(), params.step, krp_workspace.get_data(), 1, old_factor.get_data(), 1);
  }
  for (auto i = 0lu; i < ktensor.get_lambda().size(); i++)
    ls_ktensor.get_lambda()[i] = ktensor.get_lambda()[i];

  double error = cals::error::compute_error(X, ls_ktensor, krp_workspace, ten_workspace);

  const double old_error = ktensor.get_approximation_error();

  if (error < old_error) {
#if CUDA_ENABLED
    auto stream = cuda::create_stream();
#endif

    DEBUG(std::cout << "Fast forwarding... id: " << ktensor.get_id() << " Old: " << old_error << " New: " << error
                    << std::endl;);
    for (dim_t n = 0; n < ls_ktensor.get_n_modes(); n++) {
      ktensor.get_factor(n).copy(ls_ktensor.get_factor(n));

      if (params.cuda) {
#if CUDA_ENABLED
        ktensor.get_factor(n).send_to_device_async(stream);
#else
        std::cerr << "Not compiled with CUDA support" << std::endl;
        exit(EXIT_FAILURE);
#endif
      }
      cals::ops::update_gramian(ktensor.get_factor(n), gramians[n]);
    }

    ktensor.set_approximation_error(error);
    ktensor.calculate_new_fit(X_norm);

    if (params.cuda) {
#if CUDA_ENABLED
      cudaStreamSynchronize(stream);
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    }
#if CUDA_ENABLED
    cuda::destroy_stream(stream);
#endif
  }
}

void line_search_par(Ktensor &ktensor,
                     Ktensor &ls_ktensor,
                     Ktensor &ls_tr_ktensor,
                     vector<Matrix> &gramians,
                     const Tensor &X,
                     double X_norm,
                     LineSearchParams &params) {
#if CUDA_ENABLED
  auto stream = cuda::create_stream();
#endif
  for (dim_t n = 0; n < ls_ktensor.get_n_modes(); n++) {
    const auto &curr_factor = ktensor.get_factor(n);
    auto &old_factor = ls_ktensor.get_factor(n);

    for (dim_t i = 0; i < old_factor.get_n_elements(); i++)
      old_factor[i] = curr_factor[i] + params.step * (curr_factor[i] - old_factor[i]);
  }
  for (dim_t i = 0; i < ktensor.get_lambda().size(); i++)
    ls_ktensor.get_lambda()[i] = ktensor.get_lambda()[i];

  double error = cals::error::compute_error_par(X, ls_ktensor, ls_tr_ktensor);

  const double old_error = ktensor.get_approximation_error();

  if (error < old_error) {
    DEBUG(std::cout << "Fast forwarding... id: " << ktensor.get_id() << " Old: " << old_error << " New: " << error
                    << std::endl;);
    for (dim_t n = 0; n < ls_ktensor.get_n_modes(); n++) {
      ktensor.get_factor(n).copy(ls_ktensor.get_factor(n));

      if (params.cuda) {
#if CUDA_ENABLED
        ktensor.get_factor(n).send_to_device_async(stream);
#else
        std::cerr << "Not compiled with CUDA support" << std::endl;
        exit(EXIT_FAILURE);
#endif
      }
      cals::ops::update_gramian(ktensor.get_factor(n), gramians[n]);
    }

    ktensor.set_approximation_error(error);
    ktensor.calculate_new_fit(X_norm);

    if (params.cuda) {
#if CUDA_ENABLED
      cudaStreamSynchronize(stream);
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    }
#if CUDA_ENABLED
    cuda::destroy_stream(stream);
#endif
  }
}
} // namespace cals::ls

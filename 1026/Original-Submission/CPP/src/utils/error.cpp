#include "error.h"

#include "utils/mttkrp.h"

namespace cals::error
{

  double compute_error(const Tensor &X, Ktensor &ktensor, Matrix &krp_workspace, Matrix &ten_workspace)
  {
    cals::mttkrp::KrpParams krp_params;
    const auto &C = ktensor.get_factor(2);
    const auto &B = ktensor.get_factor(1);
    const auto &A = ktensor.get_factor(0);

    ktensor.denormalize();

    krp_workspace.resize(B.get_rows() * C.get_rows(), B.get_cols());
    cals::mttkrp::khatri_rao(C, B, krp_workspace, krp_params);

    ten_workspace.resize(A.get_rows(), krp_workspace.get_rows());
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, A.get_rows(), krp_workspace.get_rows(), krp_workspace.get_cols(),
                1.0, A.get_data(), A.get_col_stride(), krp_workspace.get_data(), krp_workspace.get_col_stride(),
                0.0, ten_workspace.get_data(), ten_workspace.get_col_stride());

#pragma omp parallel for default(none) shared(X, ten_workspace)
    for (dim_t i = 0; i < X.get_n_elements(); i++)
      ten_workspace[i] = X[i] - ten_workspace[i];

    ktensor.normalize();

    return ten_workspace.norm();
  }

  double compute_error_par(const Tensor &X, Ktensor &ktensor, Ktensor &tr_ktensor)
  {
    const auto &C = ktensor.get_factor(2);
    const auto &B = ktensor.get_factor(1);
    const auto &A = ktensor.get_factor(0);
    auto &C_tr = tr_ktensor.get_factor(2);
    auto &B_tr = tr_ktensor.get_factor(1);
    auto &A_tr = tr_ktensor.get_factor(0);
    const auto I = A.get_rows();
    const auto J = B.get_rows();
    const auto K = C.get_rows();
    const dim_t L = ktensor.get_components();
    const auto &lambda = ktensor.get_lambda();

    A_tr.transpose_copy(A);
    B_tr.transpose_copy(B);
    C_tr.transpose_copy(C);

    double nrm2 = 0;
#pragma omp parallel for shared(lambda, X, A_tr, B_tr, C_tr) reduction(+:nrm2)
    for (dim_t k=0; k < K; k++)
      for (dim_t j=0; j < J; j++)
        for (dim_t i=0; i < I; i++)
        {
          double tmp = 0.0;
          for (dim_t l=0; l < L; l++)
            tmp += lambda[l] * A_tr[l + i*L] * B_tr[l + j*L] * C_tr[l + k*L];
          tmp -= X[i + j*I + k*I*J];
          nrm2 += tmp * tmp;
        }

    return std::sqrt(nrm2);
  }

  double compute_fast_error(double X_norm,
                            const std::vector<double> &lambda,
                            const cals::Matrix &last_factor,
                            const cals::Matrix &last_mttkrp,
                            const cals::Matrix &gramian_hadamard)
  {
//    hadamard_all(gramians);

    // Sum the entries of the Hadamard product
    double term2 = 0.0;
    for (dim_t j = 0; j < gramian_hadamard.get_cols(); j++)
      for (dim_t i = 0; i < gramian_hadamard.get_rows(); i++)
        term2 += lambda[i] * lambda[j] * gramian_hadamard(i, j);

    // Compute the inner product of the last factor matrix and the last matrix G.
    double term3 = 0.0;
    for (dim_t j = 0; j < last_factor.get_cols(); j++)
      for (dim_t i = 0; i < last_factor.get_rows(); i++)
        term3 += lambda[j] * last_factor(i, j) * last_mttkrp(i, j);

    // Sum up the terms in the approximation error formula.
    double fast_error = 0.0;
    fast_error = std::fmax(X_norm * X_norm + term2 - 2 * term3, 0);
    fast_error = std::sqrt(fast_error);

    return fast_error;
  }
}

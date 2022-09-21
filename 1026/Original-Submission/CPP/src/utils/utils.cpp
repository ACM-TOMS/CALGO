#include "utils/utils.h"

namespace cals::utils {
std::string mode_string(vector<dim_t> const &modes) {
  std::string m_string;
  for (auto const &m : modes)
    m_string += std::to_string(m) + '-';
  m_string.pop_back();

  return m_string;
}
} // namespace cals::utils

namespace cals::ops {
void hadamard_all(vector<Matrix> &matrices) {
  for (auto i = 1lu; i < matrices.size(); i++) // skip first gramian
    matrices[0].hadamard(matrices[i]);
}

Matrix &hadamard_but_one(vector<Matrix> &matrices, int mode) {
  // Initialize target matrix
  matrices[mode].fill((function<double()> const &&)[]()->double { return static_cast<double>(1.0); });

  for (auto i = 0lu; i < matrices.size(); i++) {
    if (i == static_cast<unsigned long int>(mode)) // ...except mode...
      continue;
    else // ATA[n] := ATA[n] .* ATA[k]
      matrices[mode].hadamard(matrices[i]);
  }
  return matrices[mode];
}

void update_gramian(const Matrix &factor, Matrix &gramian) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, gramian.get_rows(), gramian.get_cols(), factor.get_rows(), 1.0,
              factor.get_data(), factor.get_col_stride(), factor.get_data(), factor.get_col_stride(), 0.0,
              gramian.get_data(), gramian.get_col_stride());
}
} // namespace cals::ops

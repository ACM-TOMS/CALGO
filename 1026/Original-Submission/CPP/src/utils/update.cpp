#include "utils/update.h"

#include <cmath>
#include <iostream>

using cals::Matrix;
using std::vector;

namespace cals::update {
struct ZeroPassiveSet : public std::exception {
  const char *what() const noexcept override { return "Passive set is zero."; }
};

struct CholFail : public std::exception {
  const char *what() const noexcept override { return "Cholesky failed due to Matrix not beind SPD."; }
};

Matrix &calculate_sp(const Matrix &y, Matrix &sp, const Matrix &G, Matrix &Gp, const vector<bool> &active) {
  int n_active = static_cast<int>(std::count(active.cbegin(), active.cend(), false));
  sp.resize(n_active, 1);
  Gp.resize(n_active, n_active);

  auto index = 0;
  for (int i = 0; i < static_cast<int>(active.size()); i++) {
    if (!active[i]) {
      sp[index] = y[i];
      auto index2 = 0;
      for (int j = 0; j < static_cast<int>(active.size()); j++)
        if (!active[j]) // TODO: Invert the index and index2++ for contiguous traversal
          Gp(index, index2++) = G(i, j);
      index++;
    }
  }

  int info = 0;
  int rhs = 1;
  int n = n_active;
  int lda = n_active;
  int ldb = n_active;
  dposv("L", &n, &rhs, Gp.get_data(), &lda, sp.get_data(), &ldb, &info);
  if (info < 0)
    std::cerr << "Parameter invalid. DPOSV returned info=" << info << std::endl;
  else if (info > 0)
    throw CholFail();

  return sp;
}

Matrix &calculate_lagrangian_multipliers(const Matrix &y, const Matrix &G, const Matrix &d, Matrix &w) {
  // TODO consider using DSYMV to potentially improve performance for larger (10?) matrices
  cblas_dgemv(CblasColMajor, CblasNoTrans, G.get_rows(), G.get_cols(), 1.0, G.get_data(), G.get_rows(), d.get_data(),
              d.get_cols(), 0.0, w.get_data(), w.get_cols());
  for (dim_t i = 0; i < w.get_n_elements(); i++)
    w[i] = y[i] - w[i];
  return w;
}

// Update the factor matrix by applying non-negativity constraints to the solution.
// factor is assumed to contain the result of the MTTKRP, while gramian contains the hadamard product of all gramians
// (except for the gramian of the target factor matrix)
Matrix &update_factor_non_negative_constrained(Matrix &factor, Matrix &gramian, vector<vector<bool>> &active_old) {
  Matrix &G = gramian;
  int n = factor.get_cols();

  const double eps = 2.2204e-16;
  const double tol = 10 * eps * gramian.one_norm() * n;

  // Create local buffers
  Matrix y(n, 1);
  Matrix d(n, 1);
  Matrix w(n, 1);
//  vector<bool> active(static_cast<unsigned long>(n), true);

  Matrix s(n, 1);
  Matrix sp(n, 1);
  Matrix Gp(n, n);

  // Iterative algorithm applied to each row of the factor matrix
  for (dim_t row = 0; row < factor.get_rows(); row++) {
    // Initialize local buffers
    d.zero();

    auto &active = active_old[row];
    // Get row and determine previous active set
    for (int i = 0; i < n; i++) {
      y[i] = factor[row + factor.get_rows() * i];
      if (y[i] > 0)
        active[i] = false;
    }

    // If all constraints are active, it's the first iteration, so no correction needed.
    if (std::find(active.begin(), active.end(), false) != active.end()) {
      try {
        calculate_sp(y, sp, G, Gp, active);

        auto index = 0;
        for (auto i = 0; i < n; i++)
          d[i] = active[i] ? 0.0 : sp[index++];

        // Modified inner loop (if a solution was negative, make it zero, update the active set, and retry)
        while (sp.min() <= tol) {
          for (int i = 0; i < n; i++)
            if (d[i] <= tol) {
              d[i] = 0.0;
              active[i] = (d[i] == 0.0);
            }

          if (std::find(active.begin(), active.end(), false) == active.end())
            throw ZeroPassiveSet();

          sp = calculate_sp(y, sp, G, Gp, active);

          index = 0;
          for (auto i = 0lu; i < active.size(); i++)
            d[i] = active[i] ? 0.0 : sp[index++];
        }
      } catch (std::exception &exception) {
        std::fill(active.begin(), active.end(), true);
        d.zero();
      }
    }

    w = calculate_lagrangian_multipliers(y, G, d, w);

    // Main loop
    while (std::find(active.begin(), active.end(), true) != active.end() && // R ≠ ∅
           w.max(active) > tol)                                             //  max(wn) > tol , n in R
    {
      auto m = w.max_id(active);
      active[m] = false;

      calculate_sp(y, sp, G, Gp, active);

      // Inner Loop
      while (sp.min() <= tol) {
        auto index = 0;
        for (auto i = 0lu; i < active.size(); i++)
          s[i] = active[i] ? 0.0 : sp[index++];

        auto a = std::numeric_limits<double>::max();
        for (int i = 0; i < n; i++)
          if (!active[i] && s[i] <= tol) {
            auto tmp = d[i] / (d[i] - s[i]);
            if (tmp < a)
              a = tmp;
          }
        if (a >= 1e20)
          std::cerr << "Bad a in NNLS" << std::endl;

        for (dim_t i = 0; i < d.get_n_elements(); i++) {
          d[i] = d[i] + a * (s[i] - d[i]);
          assert(d[i] >= 0 || std::fabs(d[i]) < tol);
          if (std::fabs(d[i]) < tol && !active[i]) {
            active[i] = true;
            d[i] = 0;
          }
        }
        sp = calculate_sp(y, sp, G, Gp, active);
      }

      if (sp.min() < 0)
        std::cerr << "CRITICAL ____________________________________" << std::endl;

      auto index = 0;
      for (auto i = 0lu; i < active.size(); i++)
        d[i] = active[i] ? 0.0 : sp[index++];

      w = calculate_lagrangian_multipliers(y, G, d, w);
    }

    for (int i = 0; i < n; i++)
      factor[row + factor.get_rows() * i] = d[i];
  }

  return factor;
}

Matrix &update_factor_unconstrained(Matrix &factor, Matrix &gramian) {
  // Solve U_n * H = G for U_n.
  int info = 0;
  int size = gramian.get_rows();
  int lda = gramian.get_col_stride();
  dpotrf("L", &size, gramian.get_data(), &lda, &info);
  if (info)
    std::cerr << "als_update_factor: DPORTF returned info=" << info << std::endl;

  cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, factor.get_rows(), factor.get_cols(),
              1.0, gramian.get_data(), gramian.get_col_stride(), factor.get_data(), factor.get_col_stride());
  cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, factor.get_rows(), factor.get_cols(),
              1.0, gramian.get_data(), gramian.get_col_stride(), factor.get_data(), factor.get_col_stride());
  return factor;
}

} // namespace cals::update

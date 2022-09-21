#ifndef CALS_KTENSOR_H
#define CALS_KTENSOR_H

#include <cfloat>
#include <cmath>
#include <vector>

#include "cals_blas.h"
#include "matrix.h"

using std::multiplies;
using std::vector;

static int universal_ktensor_id = 1;

namespace cals {
class Ktensor {
  int id{-1};         /*< Unique ID of the Ktensor */
  int components{-1}; /*< Number of components of the Ktensor */
  int iters{-1};      /*< Number of iterations spent in ALS */

  double fit{0.0};          /*< The fit to the target tensor */
  double old_fit{0.0};      /*< The fit in the previous iteration */
  double approx_error{0.0}; /*< The approximation error of the Ktensor to the target Tensor */

  bool normalized{false}; /*< Boolean indicating whether the Ktensor is normalized or not */

  vector<vector<vector<bool>>> active_set; /*< Active set for each factor matrix of last iteration (Used in NNLS) */

  vector<double> lambda{};  /*< Vector containing the lambda normalization factors per column of the factor matrices */
  vector<Matrix> factors{}; /*< Vector containing the factor matrices of the Ktensor */

  /** Recursive function that is used to compute the reconstruction of the Ktensor into a Tensor.
   */
  void rec_to_tensor(vector<dim_t> &modes, int &index, double *x_data, int level, dim_t *dim_ind);

public:
  Ktensor() = default;

  /** Allocates a Ktensor (no initialization).
   *
   * This constructor allocates a Ktensor of a specific \p rank and sizes of \p modes. The
   * contents of the Ktensor created are not initialized. One can use the randomize member function
   * to fill in the Ktensor with random values.
   *
   * @param rank The rank of the Ktensor.
   * @param vector A vector containing the sizes of each mode of the Ktensor.
   */
  Ktensor(int rank, const vector<dim_t> &modes)
      : id{universal_ktensor_id++}, components{rank},
        active_set{vector<vector<vector<bool>>>(modes.size(), vector<vector<bool>>())}, factors{vector<Matrix>(
                                                                                            modes.size())} {
    assert(rank > 0);
    for (dim_t n = 0; n < modes.size(); n++) {
      factors[n] = Matrix{modes[n], static_cast<dim_t>(rank)};
      active_set[n] = vector<vector<bool>>(modes[n], vector<bool>(rank, true));
    }

    // For each mode allocate a coefficients vector
    lambda.resize(get_components());
  }

  Ktensor &operator=(Ktensor &&rhs) = default;

  Ktensor(Ktensor &&rhs) = default;

  Ktensor(const Ktensor &rhs)
      : id{universal_ktensor_id++}, components{rhs.components}, active_set{vector<vector<vector<bool>>>(
                                                                    rhs.get_n_modes(), vector<vector<bool>>())},
        lambda{vector<double>(static_cast<size_t>(rhs.components))}, factors{vector<Matrix>(rhs.get_n_modes())} {
    for (dim_t n = 0; n < rhs.get_n_modes(); n++) {
      factors[n] = rhs.get_factor(n);
      active_set[n] = vector<vector<bool>>(get_factor(n).get_rows(), vector<bool>(components, true));
    }
    lambda = rhs.get_lambda();
  }

  Ktensor &operator=(const Ktensor &rhs) {
    if (this == &rhs) // Properly handle self assignment
      return *this;

    id = rhs.id;
    components = rhs.components;
    lambda = rhs.get_lambda();

    factors = vector<Matrix>(rhs.get_n_modes());
    for (dim_t n = 0; n < rhs.get_n_modes(); n++)
      factors[n] = rhs.get_factor(n);
    return *this;
  }

  ~Ktensor() = default;

  // Getters
  // TODO clean this. Components should either not exist or be updated properly when adjusting the dimensions of a
  // Ktensor.
  [[nodiscard]] inline int get_components() const noexcept { return factors[0].get_cols(); }

  [[nodiscard]] inline int get_iters() const noexcept { return iters; }

  [[nodiscard]] inline int get_id() const noexcept { return id; }

  inline double get_approximation_error() const noexcept { return approx_error; }

  inline vector<Matrix> &get_factors() noexcept { return factors; }

  inline vector<double> const &get_lambda() const noexcept { return lambda; }

  inline vector<double> &get_lambda() noexcept { return lambda; }

  // Setters
  inline void set_iters(int new_iters) noexcept { iters = new_iters; }

  inline void set_approximation_error(double new_error) noexcept { approx_error = new_error; }

  inline void set_factor(int index, double *data) noexcept {
    auto &target = get_factor(static_cast<dim_t>(index));
    cblas_dcopy(target.get_n_elements(), data, 1, target.get_data(), 1);
  }

  inline void set_lambda(double const *data) noexcept {
    auto i = 0;
    for (auto &l : lambda)
      l = data[i++];
  }

  /** Calculate the fit to the target Tensor.
   *
   * @param X_norm L2-Norm of the target Tensor.
   *
   * @return The fit to the target Tensor.
   */
  inline double calculate_new_fit(double X_norm) noexcept {
    assert(X_norm != 0);
    old_fit = fit;
    fit = 1 - std::fabs(approx_error) / X_norm;
    return fit;
  }

  /** Calculate the difference in fit between the last iterations.
   *
   * @return The difference in fit.
   */
  [[nodiscard]] inline double get_fit_diff() const noexcept { return std::fabs(old_fit - fit); }

  /** Get the number of modes of the Ktensor.
   */
  inline dim_t get_n_modes() const noexcept { return static_cast<dim_t>(factors.size()); }

  /** Get a reference to the last factor Matrix.
   */
  inline Matrix const &get_last_factor() const noexcept { return factors.back(); }

  /** Get a reference to a specific factor Matrix.
   *
   * @param mode The mode for which to get the factor Matrix.
   */
  inline Matrix const &get_factor(dim_t mode) const noexcept { return factors.at(mode); }

  inline Matrix &get_factor(dim_t mode) noexcept { return factors.at(mode); }

  inline vector<vector<bool>> &get_active_set(const int mode) noexcept { return active_set.at(mode); }

  [[nodiscard]] inline vector<Matrix> const &get_factors() const noexcept { return factors; }

  /** Print the contents of the Ktensor, together with some optional text.
   *
   * @param text (Optional) Text to display along with the contents of the Ktensor.
   */
  void print(const std::string &&text = "Ktensor") const;

  /** Attach every factor matrix to pointers (and copy their contents to the new locations).
   *
   * @param data_ptrs vector containing the pointers to which the factor matrices should point (must be the same length
   * as components)
   * @param multi_thread (Optional) Whether it should be performed using multiple threads or not (to avoid parallel copy
   * of overlaping memory regions).
   *
   * @return Reference to self.
   */
  Ktensor &attach(vector<double *> &data_ptrs, bool multi_thread = true);

  /** Reset the factor matrices to point to the data they own.
   */
  Ktensor &detach();

  /** Normalize the Ktensor.
   */
  Ktensor &normalize();

  /** Normalize a specific factor matrix.
   *
   * @param mode Specify the factor matrix to normalize.
   * @param iteration Depending on the iteration a different normalization function is used.
   *
   * @return Reference to self.
   */
  Ktensor &normalize(int mode, int iteration = 1);

  /** Remove normalization of the Ktensor.
   *
   * @return Reference to self.
   */
  Ktensor &denormalize();

  /** Randomize the Ktensor.
   *
   * @param r (Optional) Whether to use the global pre-seeded generator to create a reproducible set of values for
   * testing/experiments.
   *
   * @return Reference to self.
   */
  Ktensor &randomize();

  /** Fill the Ktensor with values.
   *
   * @param f function that returns a double every time it is invoked.
   *
   * @return Reference to self.
   */
  Ktensor &fill(function<double()> &&func);

  /** Convert the Ktensor to a Tensor.
   *
   * @return The reconstructed tensor.
   */
  Tensor to_tensor();

#if CUDA_ENABLED
  /** Attach every factor matrix to pointers (for GPU) (and copy their contents to the new locations).
   *
   * @param data_ptrs vector containing the pointers to which the factor matrices should point (must be the same length
   * as components)
   *
   * @return Reference to self.
   */
  Ktensor &cuattach(vector<double *> &cudata_ptrs, cudaStream_t &stream);

  /** Reset the factor matrices to point to the data they own.
   */
  Ktensor &cudetach();
#endif
};

} // namespace cals
#endif // CALS_KTENSOR_H

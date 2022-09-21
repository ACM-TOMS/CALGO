#ifndef CALS_MATRIX_H
#define CALS_MATRIX_H

#include "tensor.h"

namespace cals {
class Matrix : public Tensor {
  dim_t rows{};       /*< Number of rows of the Matrix */
  dim_t cols{};       /*< Number of columns of the Matrix */
  dim_t col_stride{}; /*< Column stride of the Matrix (all matrices stored in col major order) */

public:
  Matrix() = default;

  ~Matrix() = default;

  /** Generic constructor that allocates a Matrix (no initialization).
   *
   * This constructor allocates a Matrix with \p dim0 rows and \p dim 1 cols.
   * The contents of the Matrix are not initialized. One can use the randomize member function
   * to fill in the Matrix with random values.
   *
   * @param dim0 Desired number of rows.
   * @param dim1 Desired number of cols.
   */
  Matrix(dim_t dim0, dim_t dim1);

  /** Generic constructor that allocates a Matrix (no initialization).
   *
   * This constructor allocates a Matrix with \p dim0 rows and \p dim1 cols.
   * The contents of the Matrix are not initialized. One can use the randomize member function
   * to fill in the Matrix with random values.
   *
   * @param dim0 Desired number of rows.
   * @param dim1 Desired number of cols.
   * @param view_data (Optional) If set, the created Matrix object does not own data (RAII). Instead it points to
   * already existing data. In this case, the Matrix is a "view".
   */
  Matrix(dim_t dim0, dim_t dim1, double *view_data);

  // Move and Copy Constructors
  Matrix(Matrix &&rhs) = default;

  Matrix &operator=(Matrix &&rhs) = default;

  Matrix(const Matrix &rhs) = default;

  Matrix &operator=(const Matrix &rhs) = default;

  // Getters
  [[nodiscard]] inline dim_t get_rows() const noexcept { return rows; };

  [[nodiscard]] inline dim_t get_cols() const noexcept { return cols; };

  [[nodiscard]] inline dim_t get_col_stride() const noexcept { return col_stride; };

  /** Operator to access an element in the Matrix, giving its row and column. (all matrices are stored in col-major)
   *
   * @param row row index of the element.
   * @param col column index of the element.
   *
   * @return The element of the matrix
   */
  inline double &operator()(dim_t row, dim_t col) { return get_data()[row + col * col_stride]; };

  inline double operator()(dim_t row, dim_t col) const noexcept { return get_data()[row + col * col_stride]; };

  /** "Soft" resize, meaning the dimensions of the Matrix are only changed (not the underlying memory)
   *
   * The new requested sizes must fit inside the memory initially allocated for the Matrix (max_n_elements).
   *
   * @param new_rows new number of rows.
   * @param new_cols new number of columns.
   *
   * @return reference to self.
   */
  Matrix &resize(dim_t new_rows, dim_t new_cols) noexcept {
    vector<dim_t> modes = {new_rows, new_cols};
    Tensor::resize(new_rows * new_cols, modes);
    rows = new_rows;
    cols = new_cols;
    col_stride = new_rows;
    return *this;
  };

  /** Compute the hadamard product between the Matrix itself and another Matrix.
   *
   * The matrices must be of the same size.
   *
   * @param mat Matrix with which to compute the hadamard product.
   *
   * @return reference to self.
   */
  Matrix &hadamard(const Matrix &mat);

  /** Change the data member variable to point to a specific place in memory (other than the owned memory, if it
   * exists).
   *
   * @param data pointer to memory where the Matrix should point to.
   */
  inline void attach(double *data) { set_data(data); }

  /** Reset the data member variable to the memory "owned" by the Matrix (pointed to by the data_up variable).
   */
  inline void detach() { reset_data(); }

  /** Print the contents of the Matrix, together with some optional text.
   *
   * @param text (Optional) Text to display along with the contents of the Matrix.
   */
  void print(const std::string &&text = "Matrix") const;

  /** Print the size and other info regarding the Matrix.
   */
  void info() const;

  /** Compute the one-norm per column of the Matrix and return the max value.
   *
   * @return maximum value among the sum of magnitudes of elements per column of the Matrix.
   */
  [[nodiscard]] double one_norm() const {
    auto max = -DBL_MAX;
    for (dim_t col = 0; col < get_cols(); col++) {
      auto one_norm = cblas_dasum(get_rows(), get_data() + col * get_col_stride(), 1);
      if (one_norm > max)
        max = one_norm;
    }
    return max;
  };

  /** Read a transposed version of a Matrix.
   *
   * @param rhs Matrix, of which the transposed version to read.
   *
   * @return Reference to self.
   */
  Matrix &transpose_copy(const Matrix &rhs);

#if CUDA_ENABLED
  /** Change the data member variable (of the GPU) to point to a specific place in memory (other than the owned memory,
   * if it exists)
   *
   * @param cudata pointer to memory where the Matrix should point to.
   */
  inline void cuattach(double *cudata) { set_cudata(cudata); }

  /** Reset the data member variable (of the GPU) to the memory "owned" by the Matrix (pointed to by the data_up
   * variable)
   */
  inline void cudetach() { reset_cudata(); }
#endif
};
} // namespace cals

#endif // CALS_MATRIX_H

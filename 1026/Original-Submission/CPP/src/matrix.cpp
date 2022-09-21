#include "matrix.h"

#include <iomanip>
#include <iostream>

namespace cals {
Matrix::Matrix(dim_t dim0, dim_t dim1) : Tensor{dim0, dim1}, rows{dim0}, cols{dim1}, col_stride{dim0} {}

Matrix::Matrix(dim_t dim0, dim_t dim1, double *view_data)
    : Tensor{dim0, dim1, view_data}, rows{dim0}, cols{dim1}, col_stride{dim0} {}

Matrix &Matrix::hadamard(const Matrix &mat) {
  // TODO improve with element-wise multiplication in MKL
  for (dim_t i = 0; i < get_n_elements(); i++)
    get_data()[i] *= mat[i];
  return *this;
}

void Matrix::print(const std::string &&text) const {
  using std::cout;
  using std::endl;

  cout << "----------------------------------------" << endl;
  cout << text << endl;
  cout << "----------------------------------------" << endl;
  cout << "Rows: " << rows << ", Cols: " << cols << endl;
  for (dim_t row = 0; row < rows; row++) {
    for (dim_t col = 0; col < cols; col++)
      cout << "  " << std::setw(6) << (*this)(row, col) << "  ";
    cout << std::endl;
  }
  cout << "----------------------------------------" << endl;
}

void Matrix::info() const {
  std::cout << "nRows: " << get_rows() << ", nCols: " << get_cols() << ", nElements: " << get_n_elements()
            << ", maxNElements: " << get_max_n_elements() << std::endl;
}

Matrix &Matrix::transpose_copy(const Matrix &rhs) {
  auto const I = rows;
  auto const J = cols;
  for (dim_t i = 0; i < I; i++)
    for (dim_t j = 0; j < J; j++)
      get_data()[j + i * cols] = rhs.get_data()[i + j * rows];
  return *this;
}

} // namespace cals

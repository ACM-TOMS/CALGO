#include "ktensor.h"

#include <iostream>
#include <random>

#if CUDA_ENABLED
#include "cuda_utils.h"
#endif

namespace cals {
Ktensor &Ktensor::randomize() {
  for (auto &f : factors)
    f.randomize();
  normalize();
  return *this;
}

Ktensor &Ktensor::fill(function<double()> &&func) {
  for (auto &f : factors)
    f.fill(std::forward<decltype(func)>(func));
  normalize();
  return *this;
}

void Ktensor::rec_to_tensor(
    vector<dim_t> &modes, int &index, double *const x_data, int const level, dim_t *const dim_ind) {
  if (level >= 1)
    for (dim_t i = 0; i < modes[level - 1]; i++) {
      dim_ind[level - 1] = i;
      rec_to_tensor(modes, index, x_data, level - 1, dim_ind);
    }
  else if (level == 0) {
    double s = 0;
    for (auto r = 0; r < components; r++) {
      double m = 1.0;
      for (dim_t f = 0; f < get_n_modes(); f++)
        m *= factors[f](dim_ind[f], r);
      s += lambda[r] * m;
    }
    x_data[index++] = s;
  }
}

Tensor Ktensor::to_tensor() {
  vector<dim_t> modes(get_n_modes());
  for (dim_t i = 0; i < get_n_modes(); i++)
    modes[i] = factors[i].get_rows();

  Tensor X(modes);
  int index = 0;
  double *x_data = X.get_data();
  auto *dim_ind = new dim_t[get_n_modes()];
  int levels = static_cast<int>(get_n_modes());
  rec_to_tensor(modes, index, x_data, levels, dim_ind);
  delete[] dim_ind;
  return X;
}

Ktensor &Ktensor::normalize(int mode, int iteration) {
#pragma omp parallel for // NOLINT(openmp-use-default-none)
  for (dim_t col = 0; col < factors[mode].get_cols(); col++) {
    auto &factor = factors[mode];
    auto *data_ptr = factor.get_data() + col * factor.get_col_stride();

    if (iteration == 1)
      lambda[col] = cblas_dnrm2(factor.get_rows(), data_ptr, 1);
    else {
      auto index = cblas_idamax(factor.get_rows(), data_ptr, 1);
      assert(index >= 0 && index < factor.get_rows());
      lambda[col] = data_ptr[index];
    }
    if (lambda[col] != 0)
      cblas_dscal(factor.get_rows(), 1 / lambda[col], data_ptr, 1);
  }
  return *this;
}

Ktensor &Ktensor::normalize() {
  for (auto &l : lambda)
    l = 1.0;

  for (dim_t n = 0; n < get_n_modes(); n++) {
    auto &factor = factors[n];
    for (auto col = 0; col < get_components(); col++) {
      auto coeff = cblas_dnrm2(factor.get_rows(), factor.get_data() + col * factor.get_col_stride(), 1);
      cblas_dscal(factor.get_rows(), 1 / coeff, factor.get_data() + col * factor.get_col_stride(), 1);
      lambda[col] *= coeff;
    }
  }
  normalized = true;
  return *this;
}

Ktensor &Ktensor::denormalize() {
  auto &factor = factors[0];
  for (auto col = 0; col < get_components(); col++)
    cblas_dscal(factor.get_rows(), lambda[col], factor.get_data() + col * factor.get_col_stride(), 1);
  normalized = false;
  return *this;
}

Ktensor &Ktensor::attach(vector<double *> &data_ptrs, bool multi_thread) {
  assert(data_ptrs.size() == get_n_modes());

  auto th = get_threads();
  if (!multi_thread)
    set_threads(1);

  auto index = 0;
  for (auto &f : factors) {
    Matrix(f.get_rows(), f.get_cols(), data_ptrs[index]).copy(f);
    f.attach(data_ptrs[index++]);
  }

  if (!multi_thread)
    set_threads(th);
  return *this;
}

Ktensor &Ktensor::detach() {
  for (auto &f : factors) {
    auto *old_data = f.get_data();
    f.detach();
    f.copy(Matrix(f.get_rows(), f.get_cols(), old_data));
    Matrix(f.get_rows(), f.get_cols(), old_data).zero();
  }
  return *this;
}

void Ktensor::print(const std::string &&text) const {
  using std::cout;
  using std::endl;

  cout << "----------------------------------------" << endl;
  cout << "Ktensor:";
  cout << "----------------------------------------" << endl;

  cout << "Rank: " << get_components() << endl;
  cout << "Num Modes: " << get_n_modes() << endl;
  cout << "Modes: [ ";
  for (const auto &f : get_factors())
    cout << f.get_rows() << " ";
  cout << "]" << endl;

  cout << "Weights: [";
  for (const auto &l : lambda)
    cout << l << " ";
  cout << " ] " << endl;

  for (auto const &f : get_factors())
    f.print("factor");

  cout << "----------------------------------------" << endl;
}

#if CUDA_ENABLED
Ktensor &Ktensor::cuattach(vector<double *> &cudata_ptrs, cudaStream_t &stream) {
  assert(cudata_ptrs.size() == get_n_modes());

  auto index = 0;
  for (auto &f : factors) {
    f.cuattach(cudata_ptrs[index++]);
    f.send_to_device_async(stream);
  }

  return *this;
}

Ktensor &Ktensor::cudetach() {
  auto stream = cuda::create_stream();
  for (auto &f : factors) {
    cudaMemsetAsync(f.get_cudata(), 0, f.get_n_elements(), stream);
    f.cudetach();
  }
  cudaStreamSynchronize(stream);
  return *this;
}
#endif

} // namespace cals

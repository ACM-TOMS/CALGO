#include "tensor.h"
#include "ktensor.h"

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

using std::accumulate;
using std::cout;
using std::endl;
using std::make_unique;
using std::multiplies;

namespace cals {
Tensor::Tensor(const vector<dim_t> &modes)
    : n_elements{accumulate(modes.cbegin(), modes.cend(), static_cast<dim_t>(1), multiplies<>())},
      max_n_elements{n_elements}, modes{modes} {
  auto *host_ptr = (double *)operator new[](n_elements * sizeof(double), alignment);
  data_up = unique_ptr<double, decltype(&Dopnew)>(host_ptr, &Dopnew);
  data = data_up.get();
}

// Constructor for Matlab, to not copy the tensor twice.
Tensor::Tensor(const vector<dim_t> &modes, double *view_data)
    : n_elements{accumulate(modes.cbegin(), modes.cend(), static_cast<dim_t>(1), multiplies<>())},
      max_n_elements{n_elements}, modes{modes}, data_up{nullptr, &Dopnew}, data{view_data} {}

Tensor::Tensor(const std::string &file_name) {
  vector<dim_t> modes_read;
  cout << "Reading tensor: " << file_name << endl;

  auto file = std::ifstream(file_name, std::ios::in);

  std::string line;
  std::getline(file, line);
  std::stringstream ss(line);

  while (ss.good()) {
    std::string substr;
    getline(ss, substr, ' ');
    modes_read.push_back(static_cast<dim_t>(std::strtol(substr.c_str(), nullptr, 10)));
  }

  n_elements = accumulate(modes_read.cbegin(), modes_read.cend(), static_cast<dim_t>(1), multiplies<>());
  max_n_elements = n_elements;
  modes = modes_read;

  auto *host_ptr = (double *)operator new[](n_elements * sizeof(double), alignment);
  data_up = unique_ptr<double, decltype(&Dopnew)>(host_ptr, &Dopnew);
  data = data_up.get();

  double val = 0;
  dim_t index = 0;
  while (file >> val)
    data[index++] = val;

  assert(index == n_elements);
}

Tensor::Tensor(dim_t mode0, dim_t mode1, double *view_data)
    : n_elements{mode0 * mode1}, max_n_elements{n_elements}, modes{vector<dim_t>{mode0, mode1}} {
  if (view_data == nullptr) {
    auto *host_ptr = (double *)operator new[](n_elements * sizeof(double), alignment);
    data_up = unique_ptr<double, decltype(&Dopnew)>(host_ptr, &Dopnew);
    data = data_up.get();
  } else {
    data_up = nullptr;
    data = view_data;
  }
}

Tensor::Tensor(int rank, const vector<dim_t> &modes) {
  Ktensor P(rank, modes);
  P.randomize();

  *this = P.to_tensor();
  this->rank = rank;
}

Tensor::Tensor(const Tensor &rhs) : n_elements{rhs.n_elements}, modes{rhs.modes} {
  if (!rhs.is_view()) // If copying views, don't allocate new memory and don't copy over data
  {
    auto *host_ptr = (double *)operator new[](n_elements * sizeof(double), alignment);
    data_up = unique_ptr<double, decltype(&Dopnew)>(host_ptr, &Dopnew);
    data = data_up.get();
    cblas_dcopy(n_elements, rhs.get_data(), 1, get_data(), 1);
  } else {
    data_up = nullptr;
    data = rhs.get_data();
  }
  cout << "WARNING: Performed copy constructor =============================" << endl;
}

Tensor &Tensor::operator=(const Tensor &rhs) {
  if (this == &rhs) // Properly handle self assignment
    return *this;

  n_elements = rhs.n_elements;
  modes = rhs.modes;

  if (!rhs.is_view()) // If copying views, don't allocate new memory and don't copy over data
  {
    auto *host_ptr = (double *)operator new[](n_elements * sizeof(double), alignment);
    data_up = unique_ptr<double, decltype(&Dopnew)>(host_ptr, &Dopnew);
    data = data_up.get();
    cblas_dcopy(n_elements, rhs.get_data(), 1, get_data(), 1);
  } else {
    data_up = nullptr;
    data = rhs.get_data();
  }

  return *this;
}

Tensor &Tensor::randomize() {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  std::random_device device;
  std::mt19937 generator(device());
  fill((function<double()> const &&)[&dist, &generator ]()->double { return dist(generator); });

  return *this;
}

Tensor &Tensor::zero() {
  fill((function<double()> const &&)[]()->double { return static_cast<double>(0.0); });
  return *this;
}

Tensor &Tensor::fill(const std::function<double()> &&generator) {
  for (dim_t i = 0; i < n_elements; i++)
    data[i] = generator();
  return *this;
}

Unfolding Tensor::implicit_unfold(const dim_t mode) const {
  Unfolding unfolding{};

  auto prod_before = [this](dim_t mode) -> dim_t {
    auto prod = static_cast<dim_t>(1);
    for (dim_t n = 0; n < mode; n++)
      prod *= modes[n];
    return prod;
  };
  auto prod_after = [this](dim_t mode) -> dim_t {
    auto prod = static_cast<dim_t>(1);
    for (dim_t n = mode + 1; n < get_n_modes(); n++)
      prod *= modes[n];
    return prod;
  };

  if (mode == 0) {
    unfolding.n_blocks = 1;
    unfolding.block_offset = 0;
    unfolding.rows = modes[mode];
    unfolding.cols = prod_after(mode);
    unfolding.stride = modes[mode];
  } else if (mode == get_n_modes() - static_cast<dim_t>(1)) {
    unfolding.n_blocks = 1;
    unfolding.block_offset = 0;
    unfolding.rows = modes[mode];
    unfolding.cols = prod_before(mode);
    unfolding.stride = prod_before(mode);
  } else //(mode > 0 && mode < n_dims - 1)
  {
    unfolding.n_blocks = prod_after(mode);
    unfolding.block_offset = prod_before(mode) * modes[mode]; // include mode
    unfolding.rows = modes[mode];
    unfolding.cols = prod_before(mode);
    unfolding.stride = prod_before(mode);
  }
  return unfolding;
}

void Tensor::print(const std::string &&text) const {
  cout << "----------------------------------------" << endl;
  cout << "Modes: ";
  for (auto const &m : modes)
    cout << m << " ";
  cout << endl;
  cout << "data = [ ";
  for (dim_t i = 0; i < n_elements; i++)
    cout << data[i] << "  ";
  cout << "]" << endl;
  cout << "----------------------------------------" << endl;
}

} // namespace cals

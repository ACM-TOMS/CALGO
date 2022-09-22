#ifndef CALS_TENSOR_H
#define CALS_TENSOR_H

#include <cassert>
#include <cfloat>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "cals_blas.h"
#include "definitions.h"

#if CUDA_ENABLED

#include "cuda_utils.h"

#endif

using std::function;
using std::unique_ptr;
using std::vector;

namespace cals {

/* variable specifying the alignment with which to allocate performance critical memory segments */
constexpr std::align_val_t alignment = static_cast<std::align_val_t>(64);

#if CUDA_ENABLED

static void Ddev(double *dev_ptr) { cuda::deallocate(dev_ptr); }

#endif

static void Dopnew(double *host_ptr) { operator delete[](host_ptr, alignment); }

struct Unfolding {
  dim_t n_blocks;     /*< The number of matrix blocks. */
  dim_t block_offset; /*< The offset (# entries) between blocks in the tensor. */
  dim_t rows;         /*< The rows of one block. */
  dim_t cols;         /*< The cols of one block. */
  dim_t stride;       /*< The (row- or column-) stride of each block. */
};

// TODO Create separate class for "data_up" to support an iterator, if it does not inhibit performance.
//  or implement iterator for the tensor class to simplify certain C-like loops
class Tensor {
  int rank{0} /*< Rank of the tensor (if known) */;
  dim_t n_elements{} /*< Total number of elements in tensor */;
  dim_t max_n_elements{} /*< Total memory allocated size of the tensor (this cannot change after creation) */;
  vector<dim_t> modes; /*< Vector containing the dimensions of each mode of the tensor */

#if CUDA_ENABLED
  /* Unique pointer pointing to the data associated with the tensor (RAII) on the GPU */
  mutable unique_ptr<double, decltype(&Ddev)> cudata_up{nullptr, &Ddev};
  /* Pointer to the contents of the tensor on the GPU. This can be set to point elsewhere from the data owned by the
   * tensor (cudata_up) to facilitate "view" type functionality. */
  mutable double *cudata{};
#endif

  /* Unique pointer pointing to the data associated with the tensor (RAII) */
  unique_ptr<double, decltype(&Dopnew)> data_up{nullptr, &Dopnew};

  /* Pointer to the contents of the tensor. This can be set to point elsewhere from the data owned by the tensor
   * (data_up) to facilitate "view" type functionality. */
  double *data{};

public:
  Tensor() = default;

  ~Tensor() = default;

  /** Generic constructor that allocates a Tensor (no initialization).
   *
   * This constructor allocates a Tensor with specific sizes of \p modes.
   * The contents of the Tensor created are not initialized. One can use the randomize member function
   * to fill in the Tensor with random values.
   *
   * @param modes A vector containing the sizes of each mode of the Ktensor.
   */
  explicit Tensor(const vector<dim_t> &modes);

  /** Constructor that uses pre-existing Tensor.
   *
   * This constructor creates a Tensor object which points to data (does not own it) already existing in memory.
   * It is meant to be used, for example, with Matlab, to view a Tensor already initialized in Matlab.
   *
   * @param modes A vector containing the sizes of each mode of the Ktensor.
   * @param view_data Pointer to the location of the Tensor data.
   */
  explicit Tensor(const vector<dim_t> &modes, double *view_data);

  /** Constructor that reads a Tensor from a file.
   *
   * This constructor creates a Tensor, whose contents are read from a file.
   * The first line of the file must contain the dimensions, separated by spaces.
   * The rest of the file should contain the elements of the Tensor, separated by the newline character.
   *
   * @param file_name The name of the file to read the Tensor from.
   */
  explicit Tensor(const std::string &file_name);

  /** Constructor that simplifies the creation of Matrices (2D tensors).
   *
   * @param mode0 Number of rows of a Matrix.
   * @param mode1 Number of columns of a Matrix.
   * @param view_data (Optional) If set, the Matrix created points to data already existing in memory (does not
   * own that data). Otherwise, it allocates an (un-initialized) Matrix.
   */
  Tensor(dim_t mode0, dim_t mode1, double *view_data = nullptr); // Constructor to simplify creation of matrices

  /** Creates a Tensor of a specific rank.
   *
   * This constructor creates a Tensor of a specific rank by first creating a randomly initialized Ktensor
   * of a specific \p rank, and then transforms that Ktensor into a full Tensor.
   *
   * @param rank The desired rank of the Tensor.
   * @param modes A vector containing the sizes of each mode of the Ktensor.
   */
  Tensor(int rank, const vector<dim_t> &modes);

  // Move and Copy Constructors
  Tensor(Tensor &&rhs) = default;

  Tensor &operator=(Tensor &&rhs) = default;

  Tensor(const Tensor &);

  Tensor &operator=(const Tensor &rhs);

  // Getters
  [[nodiscard]] inline dim_t get_n_elements() const noexcept { return n_elements; };

  [[nodiscard]] inline dim_t get_max_n_elements() const noexcept { return max_n_elements; };

  [[nodiscard]] inline dim_t get_n_modes() const noexcept { return static_cast<dim_t>(modes.size()); };

  [[nodiscard]] inline vector<dim_t> get_modes() const noexcept { return modes; };

  [[nodiscard]] inline double *get_data() const noexcept { return data; };

  [[nodiscard]] inline int get_rank() const noexcept { return rank; };

  /** Change the data member variable to point to a specific place in memory (other than the owned memory, if it
   * exists).
   *
   * @param new_data pointer to memory where the Tensor should point to.
   */
  inline void set_data(double *new_data) noexcept { data = new_data; };

  /** Reset the data member variable to the memory "owned" by the Matrix (pointed to by the data_up variable).
   *
   * @return Reference to self.
   */
  inline Tensor &reset_data() noexcept {
    data = data_up.get();
    return *this;
  };

  /** Check whether the Tensor owns memory (RAII) or is it just a "view".
   *
   * @return True if the Tensor does not "own" data.
   */
  [[nodiscard]] inline bool is_view() const noexcept { return data_up == nullptr; };

  /** Operator to access an element in the Tensor, giving its index.
   *
   * @param index the index of the element.
   *
   * @return The element of the Tensor pointed to by the index.
   */
  inline double const &operator[](dim_t index) const noexcept { return data[index]; }

  inline double &operator[](dim_t index) noexcept { return data[index]; }

  /** "Soft" resize, meaning the dimensions of the Tensor are only changed (not the underlying memory)
   *
   * The new requested sizes must fit inside the memory initially allocated for the Tensor (max_n_elements).
   *
   * @param new_n_elements new number of total elements.
   * @param new_modes new dimensions of the Tensor.
   */
  void resize(dim_t new_n_elements, vector<dim_t> &new_modes) {
    assert(new_n_elements <= max_n_elements);
    assert(new_modes.size() == modes.size());

    n_elements = new_n_elements;
    modes = std::move(new_modes);
  }

  /** Compute the L2-norm of the Tensor.
   *
   * @return L2-norm of the Tensor.
   */
  [[nodiscard]] double norm() const { return cblas_dnrm2(n_elements, data, 1); };

  /** Fill the Tensor with values.
   *
   * @param f function that returns a double every time it is invoked.
   *
   * @return Reference to self.
   */
  Tensor &fill(const function<double()> &&f);

  /** Fill the Tensor with zeros.
   *
   * @return Reference to self.
   */
  Tensor &zero();

  /** Fill the Tensor with random values (uniform distribution between [-1.0--1.0]).
   *
   * @param r (Optional) Whether to generate a reproducible set of random numbers (for experiments and testing)
   *
   * @return Reference to self.
   */
  Tensor &randomize();

  /** Read (and copy) a Tensor.
   *
   * @param ten Tensor to read.
   */
  inline void copy(const Tensor &ten) noexcept { cblas_dcopy(ten.get_n_elements(), ten.get_data(), 1, data, 1); };

  /** Return the id of the maximum element (consulting a mask).
   *
   * @param mask vector of booleans specifying which ids to consider in finding the id of the maximum element.
   *
   * @return the id of the maximum element.
   */
  inline dim_t max_id(vector<bool> &mask) noexcept {
    dim_t max_id = 0;
    DEBUG(bool found = false;)
    auto max = -DBL_MAX;
    assert(mask.size() == n_elements);

    for (auto i = 0lu; i < mask.size(); i++)
      if (mask[i] && data[i] > max) {
        max = data[i];
        max_id = i;
        DEBUG(found = true;)
      }
    assert(found);
    return max_id;
  };

  /** Return the maximum element (consulting a mask).
   *
   * @param mask vector of booleans specifying which elements to consider in finding the maximum element.
   *
   * @return the maximum element of the Tensor.
   */
  inline double max(vector<bool> &mask) noexcept { return data[max_id(mask)]; };

  /** Return the minimum element.
   *
   * @return the minimum element of the Tensor.
   */
  inline double min() noexcept { return *std::min_element(data, data + n_elements); };

  /** Print the contents of the Tensor, together with some optional text.
   *
   * @param text (Optional) Text to display along with the contents of the Tensor.
   */
  void print(const std::string &&text = "Tensor") const;

  /** Calculate the variables associated with performing an MTTKRP without memory movement.
   *
   * An MTTKRP can be computing by explicitly computing the Khatri-Rao product and multiplying the unfolded Tensor
   * with this KRP. This function determines how to perform the multiplication with the Tensor in-place (without
   * transposing the Tensor).
   *
   * @param mode Calculate the implicit Unfolding to perform an MTTKRP for the \p mode.
   *
   * @return Unfolding object.
   */
  [[nodiscard]] Unfolding implicit_unfold(dim_t mode) const;

#if CUDA_ENABLED

  inline double *const &get_cudata() const noexcept { return cudata; };

  inline double *&get_cudata() noexcept { return cudata; };

  inline void allocate_cudata(dim_t size) const noexcept {
    cuda::allocate(cudata, size);
    cudata_up = unique_ptr<double, decltype(&Ddev)>(cudata, &Ddev);
  };

  inline void send_to_device() const noexcept { cuda::send_to_device(data, cudata, n_elements); };

  inline void receive_from_device() noexcept { cuda::receive_from_device(cudata, data, n_elements); };

  inline void send_to_device_async(cudaStream_t &stream) const noexcept {
    cuda::send_to_device_async(data, cudata, n_elements, stream);
  };

  inline void receive_from_device_async(cudaStream_t &stream) noexcept {
    cuda::receive_from_device_async(cudata, data, n_elements, stream);
  };

  inline void set_cudata(double *new_cudata) noexcept { cudata = new_cudata; };

  inline Tensor &reset_cudata() noexcept {
    cudata = cudata_up.get();
    return *this;
  };
#endif
};

} // namespace cals

#endif // CALS_TENSOR_H

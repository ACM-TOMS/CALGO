#ifndef CALS_MULTI_KTENSOR_H
#define CALS_MULTI_KTENSOR_H

#include "ktensor.h"
#include <map>

namespace cals {
/** Structure that defines all the relative information and workspaces to a model that has been entered to the
 * multi-Ktensor (buffer).
 * */
struct RegistryEntry {
  Ktensor &ktensor;        /*< Reference to the Ktensor. */
  vector<Matrix> gramians; /*< Gramians associated to the Ktensor. */
  int col;                 /*< The starting column in each factor matrix, where the Ktensor was placed. */
  Ktensor ls_ktensor{};    /*< Workspace related to line search (only allocated if line search is enabled). */
  Ktensor ls_tr_ktensor{}; /*< Workspace related to line search (only allocated if line search is enabled). */
};

/* A map associating every model currently in the buffer to data related to it */
typedef std::map<int, RegistryEntry> Registry;

class MultiKtensor : public Ktensor {
  int occupancy{0};          /*< Number of occupied cols.*/
  int start{0};              /*< starting column for every factor multi-matrix */
  int end{0};                /*< last column for every factor multi-matrix */
  vector<int> occupancy_vec; /*< Vector indicating the id of each ktensor (model) occupying every column */

  vector<dim_t> modes; /*< modes of the Ktensor (same as modes of the target tensor) */
  Registry registry;   /*< Keeps track of Ktensors currently in the buffer. */

  bool cuda{false};        /*< Whether CUDA is used (to also allocate GPU related memory) */
  bool line_search{false}; /*< Whether line search is used (to create workspace objects when inserting a new model) */

  /** Check whether a specific Ktensor fits into the buffer
   *
   * @throws BufferFull If the Ktensor does not fit
   * @return the index of the starting column, where the model can be placed in the buffer.
   * */
  int check_availability(Ktensor &ktensor);

  /** Adjust the \p end member variable to point to the end of the last Ktensor in the buffer.
   * */
  MultiKtensor &adjust_edges();

#if CUDA_ENABLED
  cudaStream_t stream;
#endif

public:
  MultiKtensor() = default;

  ~MultiKtensor() = default;

  MultiKtensor &operator=(MultiKtensor &&mk) = default;

  /** Create an empty MultiKtensor with a fixed maximum buffer_size (number of columns per factor matrix)
   *
   * @param modes The modes of the target Tensor.
   * @param buffer_size The maximum number of columns per factor matrix. Sets the upper limit of the sum of components
   * of the models that can concurrently reside on the buffer.
   * */
  explicit MultiKtensor(vector<dim_t> &modes, dim_t buffer_size);

  /** Add a Ktensor to the MultiKtensor.
   *
   * @throws BufferFull If the Ktensor does not fit.
   * @returns Reference to self.
   * */
  MultiKtensor &add(Ktensor &ktensor);

  /** Remove a Ktensor from the MultiKtensor.
   *
   * @param ktensor_id The unique global ID of the Ktensor.
   * @returns Reference to self.
   * */
  MultiKtensor &remove(int ktensor_id);

  // Getters & Setters
  Registry &get_registry() { return registry; }

  [[nodiscard]] inline int get_start() const noexcept { return start; }

  inline void set_cuda(bool value) { cuda = value; };

  inline void set_line_search(bool value) { line_search = value; };

  [[maybe_unused]] inline vector<dim_t> &get_modes() { return modes; };

  /** Return the id of the leftmost Ktensor in the buffer (usefull for experiments)
   *
   * */
  inline int get_leftmost_id() {
    if (!occupancy_vec.empty())
      return occupancy_vec[0];
    else
      return -1;
  };

  /** Compress the buffers to remove fragmentation.
   *
   * As models converge and get evicted from the \p occupancy_vec, fragmentation occurs in the factor matrices.
   * This function traverses the \p occupancy_vec, detects values set to 0 (meaning no model resides in the particular
   * column of every factor matrix). Then, it creates a set of operations that shift all the models to the start of
   * each factor matrix (buffer). Finally, the \p end is adjusted using the corresponding member function.
   * @returns Reference to self.
   * */
  MultiKtensor &compress();
};

/** Exception thrown when the buffer is full.
 * */
struct BufferFull : public std::exception {
  [[nodiscard]] const char *what() const noexcept override {
    return "Buffer is full, wait until some ktensors converge.";
  }
};
} // namespace cals

#endif // CALS_MULTI_KTENSOR_H

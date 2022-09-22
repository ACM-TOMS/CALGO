#ifndef CP_CALS_MTTKRP_H
#define CP_CALS_MTTKRP_H

#include <cstdint>
#include <map>
#include <vector>

#include "ktensor.h"
#include "matrix.h"
#include "tensor.h"
#include "timer.h"

namespace cals::mttkrp {
/* Define new type for a vector of lookup tables */
typedef std::vector<std::map<int, int>> LUT_v;

struct MttkrpLut {
  LUT_v lut_v{}; /*< Lookup Table containing the "best" method to compute MTTKRP for the size of the tensor, rank of
                    Ktensor and number of threads or GPU. */
  std::vector<int> keys_v{}; /*< Sorted keys of the Lookup Table, for searching. */
};

/* Enumeration of different MTTKRP methods.
 * */
enum MTTKRP_METHOD {
  MTTKRP = 0, /*< Compute KRP explicitly and multiply with unfolded tensor. (default choice for >3D tensors) */
  TWOSTEP0,   /*< Compute MTTKRP in two steps, involving one GEMM and multiple GEMVs. (only for 3D tensors) */
  TWOSTEP1, /*< Compute MTTKRP in two steps, involving multiple (parallelizable) GEMMs and multiple GEMVs. (only 3D) */
  AUTO,     /*< Automatically choose using Lookup Tables (if they exist) or our emprirical algorithm. */
  LENGTH
};

static const std::string mttkrp_method_names[MTTKRP_METHOD::LENGTH] = {"MTTKRP", "TWOSTEP0", "TWOSTEP1", "AUTO"};

struct KrpParams {
  uint64_t flops{0};  /*< Output: Number of FLOPs performed. */
  uint64_t memops{0}; /*< Output: Number of MemOps performed. */

  bool cuda{false}; /*< Input: Indicate whether cuda is used (to stream results back to the GPU when done). */

#if CUDA_ENABLED
  cuda::CudaParams cuda_params; /*< Parameters related to CUDA. */
#endif
};

struct MttkrpParams {
  MTTKRP_METHOD method{AUTO}; /*< Input: Method to use to compute MTTKRP. */
  KrpParams krp_params{};     /*< Input/Output:  Parameters to be used for the computation of the KRP (if computed). */

  cals::mttkrp::MttkrpLut
      lut{}; /*< Input: Lookup table with the best choice of computing the MTTKRP (used in AUTO with 3D tensors only) */

  bool cuda{false}; /*< Input: Indicate whether cuda is used (to stream results back to the GPU when done). */

  MttkrpTimers mttkrp_timers; /*< Output: Timer object containing a breakout of timings across the different phases. */

  u_int64_t flops{0};  /*< Output: Number of FLOPs performed. */
  u_int64_t memops{0}; /*< Output: Number of MemOps performed. */

#if CUDA_ENABLED
  cuda::CudaParams cuda_params{};
#endif
};

/** Function to compute the MTTKRP.
 *
 * @param X Target Tensor.
 * @param u Ktensor to use for MTTKRP.
 * @param workspace Vector of Matrix of appropriate size, to be used as workspace.
 * @param mode Mode for which to compute the MTTKRP.
 * @param params Specific parameters for MTTKRP.
 *
 * @return Reference to the factor matrix of the /p mode factor Matrix of /p u.
 */
cals::Matrix &mttkrp(const cals::Tensor &X,
                     cals::Ktensor &u,
                     std::vector<cals::Matrix> &workspace,
                     int mode,
                     cals::mttkrp::MttkrpParams &params);

/** Function to compute the Khatri-Rao product.
 *
 * @param[in] A left Matrix of the KRP.
 * @param[in] B right Matrix of the KRP.
 * @param[out] workspace workspace Matrix of size (at least) [A.rows*B.rows x B.cols].
 * @param[in, out] params KrpParams object specifying parameters for KRP.
 *
 * @return Reference to the workspace, where the result was stored.
 */
cals::Matrix &
khatri_rao(const cals::Matrix &A, const cals::Matrix &B, cals::Matrix &workspace, cals::mttkrp::KrpParams &params);

/** Read a LookUp Table (LUT) from a file.
 *
 * @param[in] modes The modes of the target tensor.
 * @param[in] threads Number of threads available.
 * @param[in] gpu (optional) Whether to look for GPU LUT.
 * @param[in] suppress_warning (optional) whether to suppress warning if LUT was not found.
 *
 * @return MttkrpLut object containing the LUT read from the file.
 */
MttkrpLut
read_lookup_table(std::vector<dim_t> const &modes, int threads, bool gpu = false, bool suppress_warning = false);
} // namespace cals::mttkrp

#endif // CP_CALS_MTTKRP_H

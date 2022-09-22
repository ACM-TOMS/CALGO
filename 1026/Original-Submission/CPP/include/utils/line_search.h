#ifndef CP_CALS_LINE_SEARCH_H
#define CP_CALS_LINE_SEARCH_H

#include "ktensor.h"

namespace cals::ls {
struct LineSearchParams {
  bool cuda{false}; /*< Indicate whether cuda is used (to stream results back to the GPU when done). */
  double step{0.0}; /*< Factor with which to extrapolate. */
};

/** Perform line search extrapolation. (only available for 3D tensors)
 *
 * The ktensor is extrapolated and the new approximation error is calculated (with respect to the target tensor). If the
 * new error decreases, substitute the ktensor with the extrapolation. Otherwise, leave the ktensor as is.
 *
 * @param[in, out] ktensor the Ktensor on which to perform line search. If the error decreases, the \p ktensor is
 * updated.
 * @param[in] ls_ktensor Ktensor object to temporarily store the extrapolated ktensor.
 * @param[in] krp_workspace Workspace of size (at least) [B.rows*C.rows x B.cols] to store the Khatri-rao product for
 * the reconstruction (necessary for error calculation).
 * @param[in] ten_workspace Workspace of size (at least) [I * J * K] to store the reconstructed tensor (necessary for
 * error calculation).
 * @param[out] gramians Vector containing the gramian Matrix objects, to be updated in case the \p ktensor is updated.
 * @param[in] X Target tensor, used for error calculation.
 * @param[in] X_norm L2 norm of the target tensor.
 * @param[in] params LineSearchParams object, containing parameters related to the line search extrapolation.
 */
void line_search(cals::Ktensor &ktensor,
                 cals::Ktensor &ls_ktensor,
                 cals::Matrix &krp_workspace,
                 cals::Matrix &ten_workspace,
                 std::vector<cals::Matrix> &gramians,
                 const cals::Tensor &X,
                 double X_norm,
                 cals::ls::LineSearchParams &params);

/** Perform line search extrapolation. This version favors concurrent calculation of line search, without high
 * memory usage. (only available for 3D tensors)
 *
 * The ktensor is extrapolated and the new approximation error is calculated (with respect to the target tensor). If the
 * new error decreases, substitute the ktensor with the extrapolation. Otherwise, leave the ktensor as is.
 *
 * @param[in, out] ktensor the Ktensor on which to perform line search. If the error decreases, the \p ktensor is
 * updated.
 * @param[in] ls_ktensor Ktensor object to temporarily store the extrapolated ktensor.
 * @param[in] ls_tr_ktensor Ktensor object to temporarily store the extrapolated ktensor.
 * @param[out] gramians Vector containing the gramian Matrix objects, to be updated in case the \p ktensor is updated.
 * @param[in] X Target tensor, used for error calculation.
 * @param[in] X_norm L2 norm of the target tensor.
 * @param[in] params LineSearchParams object, containing parameters related to the line search extrapolation.
 */
void line_search_par(cals::Ktensor &ktensor,
                     cals::Ktensor &ls_ktensor,
                     cals::Ktensor &ls_tr_ktensor,
                     std::vector<cals::Matrix> &gramians,
                     const cals::Tensor &X,
                     double X_norm,
                     cals::ls::LineSearchParams &params);

} // namespace cals::ls

#endif // CP_CALS_LINE_SEARCH_H

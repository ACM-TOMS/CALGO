#ifndef CP_CALS_ERROR_H
#define CP_CALS_ERROR_H

#include "ktensor.h"

namespace cals::error {

/** Compute the approximation error of the Ktensor based on the FastALS paper formula.
 *
 * Method implemented is described in the following paper:
 * A. Phan, P. Tichavský and A. Cichocki, "Fast Alternating LS Algorithms for High Order CANDECOMP/PARAFAC Tensor
 * Factorizations," in IEEE Transactions on Signal Processing, vol. 61, no. 19, pp. 4834-4846, Oct.1, 2013,
 * doi: 10.1109/TSP.2013.2269903.
 *
 * @param X_norm Norm of the target tensor.
 * @param lambda Vector of lambdas of Ktensor.
 * @param last_factor Last Factor matrix of Ktensor.
 * @param last_mttkrp The resulting Matrix of the last MTTKRP performed in the current iteration.
 * @param gramian_hadamard Matrix containing the hadamard product of the gramians of the Ktensor.
 *
 * @return The error of the Ktensor based on the FastALS paper formula.
 */
double compute_fast_error(double X_norm,
                          const std::vector<double> &lambda,
                          const cals::Matrix &last_factor,
                          const cals::Matrix &last_mttkrp,
                          const cals::Matrix &gramian_hadamard);

/** Compute the approximation error of the Ktensor by explicitly reconstructing a full tensor, elementwise subtracting
 * it from the real tensor \p X and returning the norm.
 *
 * (This function is used only in line search and thus only supports 3D tensors. The reconstruction is done by
 * first calculating the Khatri-Rao Product of factor matrices C and B, which is then multiplied with factor matrix A.)
 *
 * @param X The target tensor.
 * @param ktensor Ktensor for which to measure the approximation error.
 * @param krp_workspace Workspace capable of holding the Khatri-Rao product of factor matrices B and C.
 * @param ten_workspace Workspace capable of holding the reconstructed tensor.
 *
 * @return The approximation error of the Ktensor.
 */
double
compute_error(const cals::Tensor &X, cals::Ktensor &ktensor, cals::Matrix &krp_workspace, cals::Matrix &ten_workspace);

/** Compute the approximation error of the Ktensor by explicitly reconstructing a full tensor, elementwise subtracting
 * it from the real tensor \p X and returning the norm.
 *
 * (This function is used only in line search and thus only supports 3D tensors. The reconstruction is done on the fly,
 * without requiring extra workspace, using a custom quadruple loop implementation. It was designed to facilitate
 * concurrent line search invocations without running out of workspace memory).
 *
 * @param X The target tensor.
 * @param ktensor Ktensor for which to measure the approximation error.
 * @param tr_ktensor Ktensor to be used as workspace (of same size as \p ktensor).
 *
 * @return The approximation error of the Ktensor.
 */
double compute_error_par(const cals::Tensor &X, cals::Ktensor &ktensor, cals::Ktensor &tr_ktensor);

} // namespace cals::error
#endif // CP_CALS_ERROR_H

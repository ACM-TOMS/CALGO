#ifndef CALS_UTILS_H
#define CALS_UTILS_H

#include <string>

#include "matrix.h"

namespace cals::utils {
/** Create a string with the modes of a tensor.
 *
 * @param modes vector containing the modes of the tensor.
 *
 * @return String with the modes of a tensor.
 */
std::string mode_string(std::vector<dim_t> const &modes);
} // namespace cals::utils

namespace cals::ops {
/** Compute the gramian by performing an A^{T}A operation.
 *
 * @param[in] factor Factor Matrix for which to compute the gramian.
 * @param[out] gramian Matrix in which to store the gramian.
 */
void update_gramian(const cals::Matrix &factor, cals::Matrix &gramian);

/** Compute the hadamard product of a set of matrices except one (which will be the output).
 *
 * Hadamard of the all matrices in \p matrices, except the Matrix in position \p mode. (all matrices are assumed of
 * same size. Store the result in position \p mode.
 *
 * @param[in, out] matrices vector containing Matrix objects to be multiplied together.
 * @param[in] mode specify the index of the matrix in \p matrices, in which the output is going to be written.
 *
 * @return Reference to Matrix used as output.
 */
Matrix &hadamard_but_one(std::vector<cals::Matrix> &matrices, int mode);

/** Compute the hadamard product of all matrices and store the result in the first.
 *
 * @param[in, out] matrices vector containing Matrix objects to be multiplied together. The first matrix is used as
 * output.
 */
void hadamard_all(std::vector<cals::Matrix> &matrices);
} // namespace cals::ops

#endif // CALS_UTILS_H

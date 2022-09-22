#ifndef CP_CALS_UPDATE_H
#define CP_CALS_UPDATE_H

#include "matrix.h"

namespace cals::update {
enum UPDATE_METHOD { UNCONSTRAINED = 0, NNLS, LENGTH };

static const std::string update_method_names[UPDATE_METHOD::LENGTH] = {"unconstrained", "nnls"};

// Update the factor matrix by solving the SPD system factor * gramian = G
// hint: The factor matrix should initially hold the result of the MTTKRP
Matrix &update_factor_unconstrained(cals::Matrix &factor, cals::Matrix &gramian);

Matrix &update_factor_non_negative_constrained(cals::Matrix &factor,
                                               cals::Matrix &gramian,
                                               std::vector<std::vector<bool>> &active_old);
} // namespace cals::update

#endif // CP_CALS_UPDATE_H

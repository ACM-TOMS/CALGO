src/pcafast.m: Computes a nearly optimal low-rank approximation to a matrix.
src/diffsnorm.m: Computes the 2-norm accuracy of an approximation to a matrix.
src/eigen.m: Computes a low-rank eigendecomposition a self-adjoint matrix.
src/diffsnormschur.m: Computes 2-norm accuracy of a Schur decomp. of a matrix.
src/adaptivepca.m: Computes a low-rank approximation to a specified precision.
src/examples.m: Provides some example usages.

tests/alltests.m: Runs all tests listed below.
tests/pcafasttest.m: Tests pcafast for dense matrices.
tests/pcafasttestsparse.m: Tests pcafast for sparse matrices.
tests/eigenstest.m: Tests eigen on dense self-adjoint matrices.
tests/eigenstestsparse.m: Tests eigen on sparse self-adjoint matrices.
tests/eigenntest.m: Tests eigen on dense non-negative definite matrices.
tests/eigenntestsparse.m: Tests eigen on sparse non-negative definite matrices.
tests/diffsnormtest.m: Tests diffsnorm on dense matrices.
tests/diffsnormtestsparse.m: Tests diffsnorm on sparse matrices.
tests/diffsnormschurtest.m: Tests diffsnormschur on dense matrices.
tests/diffsnormschurtestsparse.m: Tests diffsnormschur on sparse matrices.
tests/diffsnormctest.m: Tests diffsnorm with centering on dense matrices.
tests/diffsnormctestsparse.m: Tests diffsnorm with centering on sparse matrices.
tests/adaptivepcatest.m: Tests adaptivepca for dense matrices.
tests/adaptivepcatestsparse.m: Tests adaptivepca for sparse matrices.

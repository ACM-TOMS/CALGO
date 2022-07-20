---
title: "randUTV SM-DM"
author: Basic installation and user guide
date: April 7, 2020
geometry: margin=1.5cm
fontsize: 9pt
output: pdf_document
---

`randUTV SM-DM` is a software package that accompanies the submission 
*Efficient algorithms for computing a rank-revealing UTV factorization
on parallel computing architectures* and includes an algorithm-by-blocks
implementation for shared-memory architectures (`SM` from now on) and
a distributed memory implementation (`DM` from now on) to compute 
a \textit{full} factorization of a given matrix that provides low-rank 
approximations with near-optimal error.

This guide complements the submitted paper and includes basic setup and 
execution steps for the package.  The provided software package includes 
all the algorithms described in the paper, namely:

* Directory `basic_randutv_sm` contains the *algorithm-by-blocks* for
shared-memory architectures that leverages the `libflame` SuperMatrix
infrastructure for the implementation.

* Directory `basic_randutv_dm` contains the distributed-memory implementation
based on ScaLAPACK.

# Basic installation steps for SM

## Prerequisites

The package has been tested with the Intel 2018 compiler suite, although
any other modern compiler would suffice.

The package requires an installation of the `libflame` library with
support for the SuperMatrix runtime task scheduler with OpenMP 
multithread support enabled.

To install the library, follow the next steps:

```sh
## Clone the library from the official Git repository.
https://github.com/figual/libflame.git

## Configure the library for appropriate SuperMatrix and OpenMP support. 
## Choose the desired library installation location.
./configure --prefix=/opt/libflame --enable-max-arg-list-hack \\
--enable-supermatrix --enable-multithreading=openmp

## Compile and install the library.
make && make install
```

## Package compilation

To compile the package, just tune the appropriate `libflame` installation
location in the Makefile (variable LIBFLAME) and execute `make`.

## Testing

A simple testing driver (basic_test.c) is provided to test correctness and
performance. The driver tests a small case (`m=n=8, nb=3`), although larger
matrices and different block sizes can be easily tested by modifying the code.
The macro `PRINT_DATA` can be activated to show input and output data prior and
after factorization.

```sh
# Compile software
make

# Run basic test
./basic_test.x
```
## Sample output for the basic test

The following is a sample output of the basic driver provided in the package
for `m=n=8, nb=3`.  The output includes the contents of the original matrix, 
and matrices U, T, V after factorization, together with residual computations
and report.

```sh
 Ai = [ 
5.658107e-01 2.533418e-02 2.689885e-01 2.234422e-01 1.470419e-01 9.781572e-01 3.210570e-01 4.331079e-01 
6.109299e-01 3.162596e-01 9.233823e-02 2.144171e-01 6.236489e-02 4.471605e-01 8.889882e-01 1.360038e-01 
5.057681e-01 6.127619e-02 8.809632e-01 1.198589e-02 7.724467e-02 3.548297e-01 2.567804e-01 9.622619e-01 
1.796469e-01 8.365902e-02 5.625732e-01 8.683906e-01 9.637282e-01 9.549481e-01 9.845709e-01 7.648950e-01 
8.166863e-01 9.767909e-01 6.635139e-01 3.317871e-01 2.458365e-01 4.252188e-01 8.704291e-01 6.721158e-01 
1.834717e-01 9.780583e-01 9.814409e-01 5.361120e-01 6.618976e-01 2.287188e-01 2.186908e-01 5.188587e-01 
5.846530e-01 8.738890e-01 9.619104e-01 5.565968e-01 3.858843e-01 8.024953e-03 1.240179e-01 6.624928e-01 
4.221561e-01 5.307681e-02 1.394470e-01 8.975978e-01 2.711707e-01 6.942072e-01 9.387133e-02 8.191577e-01 
 ];
% Just before computing factorization.
% Just after computing factorization.
 Af = [ 
4.028131e+00 0.000000e+00 0.000000e+00 1.572611e-03 -1.254069e-04 -3.974468e-05 2.133971e-06 -1.261850e-08 
0.000000e+00 1.560757e+00 0.000000e+00 9.447237e-02 8.679398e-03 -1.894195e-03 9.426621e-05 -4.129835e-07 
0.000000e+00 0.000000e+00 1.127758e+00 -8.138325e-02 8.320477e-02 -1.471131e-03 6.790575e-05 4.722842e-07 
0.000000e+00 0.000000e+00 0.000000e+00 1.025960e+00 0.000000e+00 0.000000e+00 2.974335e-04 -4.174366e-08 
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 7.202730e-01 0.000000e+00 -1.437927e-03 6.341683e-06 
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 4.495372e-01 -1.424684e-02 1.105057e-04 
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 2.000458e-01 0.000000e+00 
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 5.383423e-02 
 ];
 Uf = [ 
-2.630485e-01 3.834227e-01 -2.286753e-01 -2.293101e-01 -1.112542e-01 7.791924e-01 1.490663e-01 1.928160e-01 
-2.342793e-01 1.694648e-01 -5.202805e-01 3.851442e-01 -1.293763e-01 -1.561018e-01 3.835930e-01 -5.559833e-01 
-3.009446e-01 -2.776689e-02 -2.385838e-01 -6.926807e-01 5.160826e-01 -2.203839e-01 -5.281333e-02 -2.328798e-01 
-4.670319e-01 4.505680e-01 3.805026e-01 3.632432e-01 4.575581e-01 -1.711215e-01 9.187044e-02 2.346479e-01 
-4.397730e-01 -2.166848e-01 -4.699144e-01 2.052927e-01 -1.460808e-01 -1.311069e-01 -5.527670e-01 3.906442e-01 
-3.786614e-01 -4.531545e-01 3.560796e-01 1.733843e-01 8.848487e-02 4.309356e-01 -2.352114e-01 -4.955323e-01 
-3.720799e-01 -4.925132e-01 1.246324e-01 -1.451726e-01 -2.230632e-01 -1.237276e-01 6.439406e-01 3.203964e-01 
-3.033219e-01 3.544230e-01 3.349118e-01 -3.066698e-01 -6.451800e-01 -2.691409e-01 -2.094541e-01 -2.088131e-01 
 ];
 Vf = [ 
-3.232382e-01 6.700173e-03 -5.179460e-01 -1.896264e-01 -4.219098e-01 -7.662799e-03 5.096451e-01 3.917157e-01 
-3.177619e-01 -6.368799e-01 -7.140105e-02 3.226619e-01 -3.606892e-01 2.385372e-01 -4.440570e-01 4.387743e-03 
-4.178969e-01 -4.104571e-01 2.950856e-02 -3.003321e-01 4.943919e-01 1.989209e-01 3.662432e-01 -3.843046e-01 
-3.343366e-01 1.490211e-01 5.511752e-01 1.397657e-01 -4.900513e-01 -2.881237e-01 2.779974e-01 -3.771406e-01 
-2.759702e-01 1.323785e-02 4.797679e-01 3.049528e-01 3.029940e-01 1.273923e-01 1.520295e-01 6.850873e-01 
-3.480325e-01 5.912583e-01 -5.326370e-02 -3.706634e-02 -5.341519e-02 6.733877e-01 -2.137079e-01 -1.519893e-01 
-3.403486e-01 2.141266e-01 -4.306222e-01 6.027210e-01 3.279437e-01 -3.884398e-01 -3.311543e-02 -1.734226e-01 
-4.415907e-01 8.969370e-02 5.932280e-02 -5.408831e-01 7.152700e-02 -4.470919e-01 -5.116424e-01 1.846451e-01 
 ];

% Computing residuals...
Res. || A - U * T * V' ||_F / || A ||_F = [ 
    1.15101e-15 
 ]; 
Res. || I - U' * U ||_F / || U ||_F = [ 
    7.85750e-16 
 ]; 
Res. || I - V' * V ||_F / || V ||_F = [ 
    5.97083e-16 
 ]; 
Res. || sin.val.A - sin.val.T ||_F / || sin.val.A ||_F = [ 
    2.37709e-16 
 ]; 
% End of Program
```

# Basic installation steps for DM

The steps to install this subpackage are the following:

1. Uncompress the file.

2. Edit the first part of the "Makefile" file.

   The user must define the following variables (maybe replace the current 
   values assigned to these variables): 
     * the Fortran 90 compiler, 
     * the Fortran 90 compiler flags,
     * the Fortran 90 loader or linker (usually the same as the compiler),
     * the Fortran 90 loader flags, and
     * the libraries to be employed: 
         * the ScaLAPACK library (and the PBLAS and BLACS libraries, if not
           included inside it), 
         * the LAPACK library, and 
         * the BLAS library.

   There is no need to define the variable MKLROOT (and install the MKL 
   library) if public-domain ScaLAPACK, LAPACK and BLAS libraries are 
   employed.
   No more changes are required in this file.

3. Type "make" to compile and generate the executable file.

4. Execute the driver with one of the following commands:
     mpirun -np 4 pdttrandutv.x
     mpiexec -n 4 pdttrandutv.x

# Further documentation

Each folder includes a detailed reference manual generated via Doxygen.

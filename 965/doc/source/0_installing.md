# Building and Installing                                               

## Prerequisites

There are no prerequisites for building the base RIDC software and
examples in `explicit/` and `implicit/`.  To build the example in
`brusselator_gsl/`, the GNU Scientific Library and headers need to be
installed.  To build examples in `implicit_mkl/`, `brusselator_mkl/`
and `brusselator\_radau\_mkl/`, the Intel Math Kernel Library needs to
be installed, and appropriate environment variables initialized.


### Required

* A recent C++ compiler that supports (most of) the C++11 standard.  This code has been successfully
tested with GCC 4.1.x and the Intel Compiler 13.0.x.


### Optional


* Intel MKL or GNU Scientific Library are required for some of the examples.


## Obtain the Source Code

The RIDC software is hosted at http://mathgeek.us/software.html.
Users should download the latest `libridc-x.x.tar.gz`, and uncompress
the file to your desired location.


## Configure the Software

In the top level directory, `./configure --help` gives the possible
configuration options.  To configure using standard build
configuration, type `./configure --prefix=/home/user/opt/libridc`.  If
you wish to compile and check the MKL and GSL examples, add the
configuration flags `--with-intel-mkl` or `--with-gsl` respectively.


## Building with Software


The library is built by typing `make && make check && make install`.
By default, only the explicit and implicit examples are part of make
check, unless the `--with-intel-mkl` or `--with-gsl` flags are added
in the configuration step.


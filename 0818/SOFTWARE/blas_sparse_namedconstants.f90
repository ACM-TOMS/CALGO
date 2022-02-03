module blas_sparse_namedconstants

! *** Diagonal entries
      integer, parameter :: blas_non_unit_diag = 0 !DEFAULT
      integer, parameter :: blas_unit_diag = 1

! *** Indices
      integer, parameter :: blas_no_repeated_indices = 0 !DEFAULT
      integer, parameter :: blas_repeated_indices = 2

! *** Use only one half of the matrix: for sym, herm, triang. matrices
      integer, parameter :: blas_upper = 4
      integer, parameter :: blas_lower = 8

! *** structured/unstructured matrix
      integer, parameter :: blas_irregular = 0 !DEFAULT
      integer, parameter :: blas_regular = 16
      integer, parameter :: blas_block_irregular = 0 !DEFAULT
      integer, parameter :: blas_block_regular = 16
      integer, parameter :: blas_unassembled = 32

! *** Index basis of matrix elements
      integer, parameter :: blas_one_base = 0 !DEFAULT
      integer, parameter :: blas_zero_base = 64

! *** Matrix type
      integer, parameter :: blas_general = 0 !DEFAULT
      integer, parameter :: blas_symmetric = 128
      integer, parameter :: blas_upper_symmetric = 132
      integer, parameter :: blas_lower_symmetric = 136
      integer, parameter :: blas_hermitian = 256
      integer, parameter :: blas_upper_hermitian = 260
      integer, parameter :: blas_lower_hermitian = 264
      integer, parameter :: blas_upper_triangular = 516
      integer, parameter :: blas_lower_triangular = 520

! *** For block matrices: specify block-internal storage
      integer, parameter :: blas_col_major = 0 !DEFAULT
      integer, parameter :: blas_row_major = 1024

! *** Other constants
      integer, parameter :: blas_valid_handle = -1
      integer, parameter :: blas_invalid_handle = -10
      integer, parameter :: blas_new_handle = -11
      integer, parameter :: blas_open_handle = -12
      integer, parameter :: blas_real = -2
      integer, parameter :: blas_complex = -3
      integer, parameter :: blas_integer = -4
      integer, parameter :: blas_single_precision = -5
      integer, parameter :: blas_double_precision = -6
      integer, parameter :: blas_num_rows = -7
      integer, parameter :: blas_num_cols = -8
      integer, parameter :: blas_num_nonzeros = -9

! *** Error codes
      integer, parameter :: blas_error_memalloc = -20
      integer, parameter :: blas_error_memdeloc = -21
      integer, parameter :: blas_error_singtria = -22
      integer, parameter :: blas_error_param    = -23

! *** Definition of numerical precisions
      integer, parameter :: sp = SELECTED_REAL_KIND(6,37)
      integer, parameter :: dp = SELECTED_REAL_KIND(15,307)

end module blas_sparse_namedconstants

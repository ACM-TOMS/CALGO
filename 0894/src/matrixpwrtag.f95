! Basic routines for matrix algebra.
!
! This module is intended for internal use only.

module matrixpwrtag

  use floattypes
  use blasinterface
  implicit none


  !    Name : .dgtimes. ( a generic name of an operator )
  !   Usage : (a .dgtimes. b)
  ! Purpose : This binary operator computes a diagonal matrix
  !         : times a diagonal matrix.
  !   Input : A and B are diagonal matrices stored in one dimensional arrays.
  !  Output : The result is the product of A and B. (This is not the Hadamard
  !         : product but the matrix-matrix multiplication.)
  !    Note : The precedence of the operator is lower than + and -;
  !         : hence, it is recommended to enclose the operation by
  !         : a round bracket.

  interface operator(.dgtimes.)
    module procedure rsdg_times               ! Real    in Single precision
    module procedure csdg_times               ! Complex in Single precision
    module procedure rwdg_times               ! Real    in Double precision
    module procedure cwdg_times               ! Complex in Double precision
#ifdef __USE_TPREC
    module procedure rtdg_times               ! Real    in Triple precision
    module procedure ctdg_times               ! Complex in Triple precision
#endif
#ifdef __USE_QPREC
    module procedure rqdg_times               ! Real    in Quadruple precision
    module procedure cqdg_times               ! Complex in Quadruple precision
#endif
  end interface


  !    Name : .sqtimes. ( a generic name of an operator )
  !   Usage : (a .sqtimes. b)
  ! Purpose : This binary operator computes a rectangular matrix
  !         : times a rectangular matrix.
  !   Input : A and B are rectangular matrices stored in two dimensional arrays.
  !  Output : The result is the product of A and B. (This is not the Hadamard
  !         : product but the matrix-matrix multiplication.)
  !    Note : The precedence of the operator is lower than + and -;
  !         : hence, it is recommended to enclose the operation by
  !         : a round bracket.

  interface operator(.sqtimes.)
    module procedure rssq_times               ! Real    in Single    precision
    module procedure cssq_times               ! Complex in Single    precision
    module procedure rwsq_times               ! Real    in Double    precision
    module procedure cwsq_times               ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rtsq_times               ! Real    in Triple    precision
    module procedure ctsq_times               ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqsq_times               ! Real    in Quadruple precision
    module procedure cqsq_times               ! Complex in Quadruple precision
#endif
  end interface


  !    Name : .trtimes. ( a generic name of an operator )
  !   Usage : (a .trtimes. b)
  ! Purpose : This binary operator computes an upper quasi-triangular matrix
  !         : times an upper quasi-triangular matrix.
  !   Input : A and B are upper quasi-triangular matrices stored in two
  !         : dimensional arrays. A and B must have the same block structure.
  !  Output : The result is the product of A and B. (This is not the Hadamard
  !         : product but the matrix-matrix multiplication.)
  !    Note : The precedence of the operator is lower than + and -;
  !         : hence, it is recommended to enclose the operation by
  !         : a round bracket.

  interface operator(.trtimes.)
    module procedure rstr_times               ! Real    in Single    precision
    module procedure cstr_times               ! Complex in Single    precision
    module procedure rwtr_times               ! Real    in Double    precision
    module procedure cwtr_times               ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rttr_times               ! Real    in Triple    precision
    module procedure cttr_times               ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqtr_times               ! Real    in Quadruple precision
    module procedure cqtr_times               ! Complex in Quadruple precision
#endif
  end interface


  !    Name : .urtimes. ( a generic name of an operator )
  !   Usage : (a .urtimes. b)
  ! Purpose : This binary operator computes an upper quasi-triangular matrix
  !         : times a rectangular matrix.
  !   Input : A is an upper quasi-triangular matrix, and B is a rectangular
  !         : matrix, both of which are stored in two dimensional arrays.
  !  Output : The result is the product of A and B. (This is not the Hadamard
  !         : product but the matrix-matrix multiplication.)
  !    Note : The precedence of the operator is lower than + and -;
  !         : hence, it is recommended to enclose the operation by
  !         : a round bracket.

  interface operator(.urtimes.)
    module procedure rsur_times               ! Real    in Single    precision
    module procedure csur_times               ! Complex in Single    precision
    module procedure rwur_times               ! Real    in Double    precision
    module procedure cwur_times               ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rtur_times               ! Real    in Triple    precision
    module procedure ctur_times               ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqur_times               ! Real    in Quadruple precision
    module procedure cqur_times               ! Complex in Quadruple precision
#endif
  end interface


  !    Name : .rutimes. ( a generic name of an operator )
  !   Usage : (a .rutimes. b)
  ! Purpose : This binary operator computes a rectangular matrix
  !         : times an upper quasi-triangular matrix.
  !   Input : A is a rectangular matrix, and B is an upper quasi-triangular
  !         : matrix, both of which are stored in two dimensional arrays.
  !  Output : The result is the product of A and B. (This is not the Hadamard
  !         : product but the matrix-matrix multiplication.)
  !    Note : The precedence of the operator is lower than + and -;
  !         : hence, it is recommended to enclose the operation by
  !         : a round bracket.

  interface operator(.rutimes.)
    module procedure rsru_times               ! Real in Single precision
    module procedure csru_times               ! Complex in Single    precision
    module procedure rwru_times               ! Real    in Double    precision
    module procedure cwru_times               ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rtru_times               ! Real    in Triple    precision
    module procedure ctru_times               ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqru_times               ! Real    in Quadruple precision
    module procedure cqru_times               ! Complex in Quadruple precision
#endif
  end interface

  integer, parameter :: mpt_blksize = 50      ! is the block size used for all
                                              ! block algorithms in the library.
  public

contains

!    Name : rsdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rsdg_norm1(arg)
  real   (kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: rsdg_norm1

  rsdg_norm1 = maxval(abs(arg))
end function rsdg_norm1


!    Name : rsdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rsdg_normi(arg)
  real   (kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: rsdg_normi

  rsdg_normi = maxval(abs(arg))
end function rsdg_normi


!    Name : rsdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rsdg_normf(arg)
  real   (kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: rsdg_normf
  real   (kind=sp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_sp)
  rsdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rsdg_normf


!    Name : rssq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rssq_norm1(arg)
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rssq_norm1

  integer       :: j
  real(kind=sp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  rssq_norm1 = maxval(colsum)
end function rssq_norm1


!    Name : rssq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rssq_normi(arg)
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rssq_normi

  integer       :: i
  real(kind=sp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  rssq_normi = maxval(rowsum)
end function rssq_normi


!    Name : rssq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rssq_normf(arg)
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rssq_normf
  real   (kind=sp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_sp)
  rssq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rssq_normf


!    Name : rstr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rstr_norm1(arg)
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rstr_norm1

  rstr_norm1 = rssq_norm1(arg)
end function rstr_norm1


!    Name : rstr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rstr_normi(arg)
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rstr_normi

  rstr_normi = rssq_normi(arg)
end function rstr_normi


!    Name : rstr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rstr_normf(arg)
  real   (kind=sp), intent(in) :: arg(:,:)
  real   (kind=sp)             :: rstr_normf

  rstr_normf = rssq_normf(arg)
end function rstr_normf


!    Name : rsdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function rsdg_eye(n)
  integer, intent(in) :: n
  real   (kind=sp)    :: rsdg_eye(1:n)

  rsdg_eye = 1.0_sp
end function rsdg_eye


!    Name : rssq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function rssq_eye(n)
  integer, intent(in) :: n
  real   (kind=sp)    :: rssq_eye(1:n,1:n)

  integer :: i

  rssq_eye = 0.0_sp
  do i=1,n
    rssq_eye(i,i) = 1.0_sp
  end do
end function rssq_eye


!    Name : rstr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function rstr_eye(n)
  integer, intent(in) :: n
  real   (kind=sp)    :: rstr_eye(1:n,1:n)

  integer :: i

  rstr_eye = 0.0_sp
  do i=1,n
    rstr_eye(i,i) = 1.0_sp
  end do
end function rstr_eye


!    Name : rsdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rsdg_times(a,b)
  real   (kind=sp), intent(in) :: a(1:), b(1:)
  real   (kind=sp)             :: rsdg_times(size(a))

  rsdg_times = a * b
end function rsdg_times


!    Name : rssq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rssq_times(a,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: rssq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,rssq_times,1,0)
#else
  rssq_times = matmul(a,b)
#endif
end function rssq_times


!    Name : rstr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rstr_times(a,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: rstr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  rstr_times = 0.0_sp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rstr_times,1,0)
#else
    rstr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            rstr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            rstr_times(h(i):t(i), h(j):t(j))                           &
          = rstr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      rstr_times(j+1,j) = 0.0_sp
    end if
  end do

end function rstr_times


!    Name : rsur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rsur_times(a,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: rsur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  rsur_times = 0.0_sp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rsur_times,1,0)
#else
    rsur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),rsur_times(hi:ti,:),1,1)
#else
          rsur_times(hi:ti,:)                                                 &
        = rsur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function rsur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rsru_times(a,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: rsru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  rsru_times = 0.0_sp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rsru_times,1,0)
#else
    rsru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),rsru_times(:,hj:tj),1,1)
#else
          rsru_times(:,hj:tj)                                                 &
        = rsru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function rsru_times


!    Name : rsdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rsdg_ludcmp(a,p)
  real   (kind=sp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rsdg_ludcmp


!    Name : rstr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rstr_ludcmp(a,p)
  real   (kind=sp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rstr_ludcmp


!    Name : rssq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``rssq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rssq_ludcmp(a,p)
  real   (kind=sp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK

  integer  :: i, n, info
  external :: sgetrf

  n = size(a,1); p = (/(i,i=1,n)/)
  call sgetrf(n, n, a, n, p, info)

#else

  real   (kind=sp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine rssq_ludcmp


!    Name : rsdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rsdg_ainvbf(a,p,b)
  real   (kind=sp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  real   (kind=sp)            :: rsdg_ainvbf(1:size(a,1))

  call rsdg_ainvb(a,p,b,rsdg_ainvbf)
end function rsdg_ainvbf


!    Name : rsdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rsdg_ainvb(a,p,b,v)
  real   (kind=sp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  real   (kind=sp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine rsdg_ainvb


!    Name : rstr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rstr_ainvbf(a,p,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=sp)             :: rstr_ainvbf(1:size(a,1),1:size(a,1))

  call rstr_ainvb(a,p,b,rstr_ainvbf)
end function rstr_ainvbf


!    Name : rstr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rstr_ainvb(a,p,b,c)
  real   (kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=sp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  real   (kind=sp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call rstr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = 0.0_sp
    d = 0.0_sp
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call rstr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = 0.0_sp
    end if
  end do
end subroutine rstr_ainvb


!    Name : rstr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rstr_ainvbblock(a,p,b,c)
  real   (kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=sp), intent(out) :: c(1:,1:)

  real   (kind=sp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = 0.0_sp; m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine rstr_ainvbblock


!    Name : rssq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rssq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rssq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function rssq_ainvbf(a,p,b)
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=sp)             :: rssq_ainvbf(1:size(a,1),1:size(b,2))

  call rssq_ainvb(a,p,b,rssq_ainvbf)
end function rssq_ainvbf


!    Name : rssq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rssq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rssq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rssq_ainvb(a,p,b,c)
  real   (kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=sp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK

  integer          :: m, n, info, ptmp(size(p))
  real   (kind=sp) :: atmp(size(a,1), size(a,2))
  external         :: sgetrs

  n = size(a,1); m = size(b,2)
  atmp = a; c = b; ptmp = p
  call sgetrs('N', n, m, atmp, n, ptmp, c, n, info)

#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  real   (kind=sp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = 0.0_sp
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = 0.0_sp
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    real   (kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=sp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    real   (kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=sp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine rssq_ainvb


!    Name : csdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function csdg_norm1(arg)
  complex(kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: csdg_norm1

  csdg_norm1 = maxval(abs(arg))
end function csdg_norm1


!    Name : csdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function csdg_normi(arg)
  complex(kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: csdg_normi

  csdg_normi = maxval(abs(arg))
end function csdg_normi


!    Name : csdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function csdg_normf(arg)
  complex(kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: csdg_normf
  real   (kind=sp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_sp)
  csdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function csdg_normf


!    Name : cssq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cssq_norm1(arg)
  complex(kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: cssq_norm1

  integer       :: j
  real(kind=sp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  cssq_norm1 = maxval(colsum)
end function cssq_norm1


!    Name : cssq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cssq_normi(arg)
  complex(kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: cssq_normi

  integer       :: i
  real(kind=sp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  cssq_normi = maxval(rowsum)
end function cssq_normi


!    Name : cssq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cssq_normf(arg)
  complex(kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: cssq_normf
  real   (kind=sp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_sp)
  cssq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function cssq_normf


!    Name : cstr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cstr_norm1(arg)
  complex(kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: cstr_norm1

  cstr_norm1 = cssq_norm1(arg)
end function cstr_norm1


!    Name : cstr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cstr_normi(arg)
  complex(kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: cstr_normi

  cstr_normi = cssq_normi(arg)
end function cstr_normi


!    Name : cstr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cstr_normf(arg)
  complex(kind=sp), intent(in) :: arg(:,:)
  real   (kind=sp)             :: cstr_normf

  cstr_normf = cssq_normf(arg)
end function cstr_normf


!    Name : csdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function csdg_eye(n)
  integer, intent(in) :: n
  complex(kind=sp)    :: csdg_eye(1:n)

  csdg_eye = cmplx(1.0_sp,0.0_sp,sp)
end function csdg_eye


!    Name : cssq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function cssq_eye(n)
  integer, intent(in) :: n
  complex(kind=sp)    :: cssq_eye(1:n,1:n)

  integer :: i

  cssq_eye = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,n
    cssq_eye(i,i) = cmplx(1.0_sp,0.0_sp,sp)
  end do
end function cssq_eye


!    Name : cstr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function cstr_eye(n)
  integer, intent(in) :: n
  complex(kind=sp)    :: cstr_eye(1:n,1:n)

  integer :: i

  cstr_eye = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,n
    cstr_eye(i,i) = cmplx(1.0_sp,0.0_sp,sp)
  end do
end function cstr_eye


!    Name : csdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function csdg_times(a,b)
  complex(kind=sp), intent(in) :: a(1:), b(1:)
  complex(kind=sp)             :: csdg_times(size(a))

  csdg_times = a * b
end function csdg_times


!    Name : cssq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cssq_times(a,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=sp)             :: cssq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,cssq_times,1,0)
#else
  cssq_times = matmul(a,b)
#endif
end function cssq_times


!    Name : cstr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cstr_times(a,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=sp)             :: cstr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  cstr_times = cmplx(0.0_sp,0.0_sp,sp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cstr_times,1,0)
#else
    cstr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            cstr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            cstr_times(h(i):t(i), h(j):t(j))                           &
          = cstr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      cstr_times(j+1,j) = cmplx(0.0_sp,0.0_sp,sp)
    end if
  end do

end function cstr_times


!    Name : csur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function csur_times(a,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=sp)             :: csur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  csur_times = cmplx(0.0_sp,0.0_sp,sp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,csur_times,1,0)
#else
    csur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),csur_times(hi:ti,:),1,1)
#else
          csur_times(hi:ti,:)                                                 &
        = csur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function csur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function csru_times(a,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=sp)             :: csru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  csru_times = cmplx(0.0_sp,0.0_sp,sp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,csru_times,1,0)
#else
    csru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),csru_times(:,hj:tj),1,1)
#else
          csru_times(:,hj:tj)                                                 &
        = csru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function csru_times


!    Name : csdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine csdg_ludcmp(a,p)
  complex(kind=sp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine csdg_ludcmp


!    Name : cstr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cstr_ludcmp(a,p)
  complex(kind=sp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cstr_ludcmp


!    Name : cssq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``cssq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cssq_ludcmp(a,p)
  complex(kind=sp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK

  integer  :: i, n, info
  external :: cgetrf

  n = size(a,1); p = (/(i,i=1,n)/)
  call cgetrf(n, n, a, n, p, info)

#else

  complex(kind=sp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine cssq_ludcmp


!    Name : csdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function csdg_ainvbf(a,p,b)
  complex(kind=sp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  complex(kind=sp)            :: csdg_ainvbf(1:size(a,1))

  call csdg_ainvb(a,p,b,csdg_ainvbf)
end function csdg_ainvbf


!    Name : csdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine csdg_ainvb(a,p,b,v)
  complex(kind=sp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  complex(kind=sp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine csdg_ainvb


!    Name : cstr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cstr_ainvbf(a,p,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=sp)             :: cstr_ainvbf(1:size(a,1),1:size(a,1))

  call cstr_ainvb(a,p,b,cstr_ainvbf)
end function cstr_ainvbf


!    Name : cstr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cstr_ainvb(a,p,b,c)
  complex(kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=sp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  complex(kind=sp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call cstr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = cmplx(0.0_sp,0.0_sp,sp)
    d = cmplx(0.0_sp,0.0_sp,sp)
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call cstr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = cmplx(0.0_sp,0.0_sp,sp)
    end if
  end do
end subroutine cstr_ainvb


!    Name : cstr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cstr_ainvbblock(a,p,b,c)
  complex(kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=sp), intent(out) :: c(1:,1:)

  complex(kind=sp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = cmplx(0.0_sp,0.0_sp,sp); m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine cstr_ainvbblock


!    Name : cssq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cssq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cssq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function cssq_ainvbf(a,p,b)
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=sp)             :: cssq_ainvbf(1:size(a,1),1:size(b,2))

  call cssq_ainvb(a,p,b,cssq_ainvbf)
end function cssq_ainvbf


!    Name : cssq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cssq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cssq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cssq_ainvb(a,p,b,c)
  complex(kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=sp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK

  integer          :: m, n, info, ptmp(size(p))
  complex(kind=sp) :: atmp(size(a,1), size(a,2))
  external         :: cgetrs

  n = size(a,1); m = size(b,2)
  atmp = a; c = b; ptmp = p
  call cgetrs('N', n, m, atmp, n, ptmp, c, n, info)

#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  complex(kind=sp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = cmplx(0.0_sp,0.0_sp,sp)
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = cmplx(0.0_sp,0.0_sp,sp)
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    complex(kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=sp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    complex(kind=sp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=sp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine cssq_ainvb


!    Name : rwdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rwdg_norm1(arg)
  real   (kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: rwdg_norm1

  rwdg_norm1 = maxval(abs(arg))
end function rwdg_norm1


!    Name : rwdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rwdg_normi(arg)
  real   (kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: rwdg_normi

  rwdg_normi = maxval(abs(arg))
end function rwdg_normi


!    Name : rwdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rwdg_normf(arg)
  real   (kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: rwdg_normf
  real   (kind=wp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_wp)
  rwdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rwdg_normf


!    Name : rwsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rwsq_norm1(arg)
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwsq_norm1

  integer       :: j
  real(kind=wp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  rwsq_norm1 = maxval(colsum)
end function rwsq_norm1


!    Name : rwsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rwsq_normi(arg)
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwsq_normi

  integer       :: i
  real(kind=wp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  rwsq_normi = maxval(rowsum)
end function rwsq_normi


!    Name : rwsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rwsq_normf(arg)
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwsq_normf
  real   (kind=wp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_wp)
  rwsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rwsq_normf


!    Name : rwtr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rwtr_norm1(arg)
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwtr_norm1

  rwtr_norm1 = rwsq_norm1(arg)
end function rwtr_norm1


!    Name : rwtr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rwtr_normi(arg)
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwtr_normi

  rwtr_normi = rwsq_normi(arg)
end function rwtr_normi


!    Name : rwtr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rwtr_normf(arg)
  real   (kind=wp), intent(in) :: arg(:,:)
  real   (kind=wp)             :: rwtr_normf

  rwtr_normf = rwsq_normf(arg)
end function rwtr_normf


!    Name : rwdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function rwdg_eye(n)
  integer, intent(in) :: n
  real   (kind=wp)    :: rwdg_eye(1:n)

  rwdg_eye = 1.0_wp
end function rwdg_eye


!    Name : rwsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function rwsq_eye(n)
  integer, intent(in) :: n
  real   (kind=wp)    :: rwsq_eye(1:n,1:n)

  integer :: i

  rwsq_eye = 0.0_wp
  do i=1,n
    rwsq_eye(i,i) = 1.0_wp
  end do
end function rwsq_eye


!    Name : rwtr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function rwtr_eye(n)
  integer, intent(in) :: n
  real   (kind=wp)    :: rwtr_eye(1:n,1:n)

  integer :: i

  rwtr_eye = 0.0_wp
  do i=1,n
    rwtr_eye(i,i) = 1.0_wp
  end do
end function rwtr_eye


!    Name : rwdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rwdg_times(a,b)
  real   (kind=wp), intent(in) :: a(1:), b(1:)
  real   (kind=wp)             :: rwdg_times(size(a))

  rwdg_times = a * b
end function rwdg_times


!    Name : rwsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rwsq_times(a,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rwsq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,rwsq_times,1,0)
#else
  rwsq_times = matmul(a,b)
#endif
end function rwsq_times


!    Name : rwtr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rwtr_times(a,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rwtr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  rwtr_times = 0.0_wp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rwtr_times,1,0)
#else
    rwtr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            rwtr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            rwtr_times(h(i):t(i), h(j):t(j))                           &
          = rwtr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      rwtr_times(j+1,j) = 0.0_wp
    end if
  end do

end function rwtr_times


!    Name : rwur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rwur_times(a,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rwur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  rwur_times = 0.0_wp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rwur_times,1,0)
#else
    rwur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),rwur_times(hi:ti,:),1,1)
#else
          rwur_times(hi:ti,:)                                                 &
        = rwur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function rwur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rwru_times(a,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rwru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  rwru_times = 0.0_wp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rwru_times,1,0)
#else
    rwru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),rwru_times(:,hj:tj),1,1)
#else
          rwru_times(:,hj:tj)                                                 &
        = rwru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function rwru_times


!    Name : rwdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rwdg_ludcmp(a,p)
  real   (kind=wp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rwdg_ludcmp


!    Name : rwtr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rwtr_ludcmp(a,p)
  real   (kind=wp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rwtr_ludcmp


!    Name : rwsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``rwsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rwsq_ludcmp(a,p)
  real   (kind=wp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK

  integer  :: i, n, info
  external :: dgetrf

  n = size(a,1); p = (/(i,i=1,n)/)
  call dgetrf(n, n, a, n, p, info)

#else

  real   (kind=wp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine rwsq_ludcmp


!    Name : rwdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rwdg_ainvbf(a,p,b)
  real   (kind=wp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  real   (kind=wp)            :: rwdg_ainvbf(1:size(a,1))

  call rwdg_ainvb(a,p,b,rwdg_ainvbf)
end function rwdg_ainvbf


!    Name : rwdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rwdg_ainvb(a,p,b,v)
  real   (kind=wp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  real   (kind=wp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine rwdg_ainvb


!    Name : rwtr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rwtr_ainvbf(a,p,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=wp)             :: rwtr_ainvbf(1:size(a,1),1:size(a,1))

  call rwtr_ainvb(a,p,b,rwtr_ainvbf)
end function rwtr_ainvbf


!    Name : rwtr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rwtr_ainvb(a,p,b,c)
  real   (kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=wp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  real   (kind=wp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call rwtr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = 0.0_wp
    d = 0.0_wp
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call rwtr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = 0.0_wp
    end if
  end do
end subroutine rwtr_ainvb


!    Name : rwtr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rwtr_ainvbblock(a,p,b,c)
  real   (kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=wp), intent(out) :: c(1:,1:)

  real   (kind=wp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = 0.0_wp; m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine rwtr_ainvbblock


!    Name : rwsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rwsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rwsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function rwsq_ainvbf(a,p,b)
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=wp)             :: rwsq_ainvbf(1:size(a,1),1:size(b,2))

  call rwsq_ainvb(a,p,b,rwsq_ainvbf)
end function rwsq_ainvbf


!    Name : rwsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rwsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rwsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rwsq_ainvb(a,p,b,c)
  real   (kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=wp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK

  integer          :: m, n, info, ptmp(size(p))
  real   (kind=wp) :: atmp(size(a,1), size(a,2))
  external         :: dgetrs

  n = size(a,1); m = size(b,2)
  atmp = a; c = b; ptmp = p
  call dgetrs('N', n, m, atmp, n, ptmp, c, n, info)

#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  real   (kind=wp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = 0.0_wp
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = 0.0_wp
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    real   (kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=wp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    real   (kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=wp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine rwsq_ainvb


!    Name : cwdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cwdg_norm1(arg)
  complex(kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: cwdg_norm1

  cwdg_norm1 = maxval(abs(arg))
end function cwdg_norm1


!    Name : cwdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cwdg_normi(arg)
  complex(kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: cwdg_normi

  cwdg_normi = maxval(abs(arg))
end function cwdg_normi


!    Name : cwdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cwdg_normf(arg)
  complex(kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: cwdg_normf
  real   (kind=wp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_wp)
  cwdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function cwdg_normf


!    Name : cwsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cwsq_norm1(arg)
  complex(kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: cwsq_norm1

  integer       :: j
  real(kind=wp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  cwsq_norm1 = maxval(colsum)
end function cwsq_norm1


!    Name : cwsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cwsq_normi(arg)
  complex(kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: cwsq_normi

  integer       :: i
  real(kind=wp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  cwsq_normi = maxval(rowsum)
end function cwsq_normi


!    Name : cwsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cwsq_normf(arg)
  complex(kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: cwsq_normf
  real   (kind=wp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_wp)
  cwsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function cwsq_normf


!    Name : cwtr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cwtr_norm1(arg)
  complex(kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: cwtr_norm1

  cwtr_norm1 = cwsq_norm1(arg)
end function cwtr_norm1


!    Name : cwtr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cwtr_normi(arg)
  complex(kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: cwtr_normi

  cwtr_normi = cwsq_normi(arg)
end function cwtr_normi


!    Name : cwtr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cwtr_normf(arg)
  complex(kind=wp), intent(in) :: arg(:,:)
  real   (kind=wp)             :: cwtr_normf

  cwtr_normf = cwsq_normf(arg)
end function cwtr_normf


!    Name : cwdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function cwdg_eye(n)
  integer, intent(in) :: n
  complex(kind=wp)    :: cwdg_eye(1:n)

  cwdg_eye = cmplx(1.0_wp,0.0_wp,wp)
end function cwdg_eye


!    Name : cwsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function cwsq_eye(n)
  integer, intent(in) :: n
  complex(kind=wp)    :: cwsq_eye(1:n,1:n)

  integer :: i

  cwsq_eye = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,n
    cwsq_eye(i,i) = cmplx(1.0_wp,0.0_wp,wp)
  end do
end function cwsq_eye


!    Name : cwtr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function cwtr_eye(n)
  integer, intent(in) :: n
  complex(kind=wp)    :: cwtr_eye(1:n,1:n)

  integer :: i

  cwtr_eye = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,n
    cwtr_eye(i,i) = cmplx(1.0_wp,0.0_wp,wp)
  end do
end function cwtr_eye


!    Name : cwdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cwdg_times(a,b)
  complex(kind=wp), intent(in) :: a(1:), b(1:)
  complex(kind=wp)             :: cwdg_times(size(a))

  cwdg_times = a * b
end function cwdg_times


!    Name : cwsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cwsq_times(a,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=wp)             :: cwsq_times(1:size(a,1),1:size(b,2))

#if 0

  ! to avoid the problem in __f___pl__pp_zgemm_tn_

  cwsq_times = matmul(a,b)

#else

  integer :: i,j
  do i=1,size(a,1)
    do j=1,size(b,2)
      cwsq_times(i,j) = sum(a(i,:) * b(:,j))
    end do
  end do

#endif
end function cwsq_times


!    Name : cwtr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cwtr_times(a,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=wp)             :: cwtr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  cwtr_times = cmplx(0.0_wp,0.0_wp,wp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cwtr_times,1,0)
#else
    cwtr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            cwtr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            cwtr_times(h(i):t(i), h(j):t(j))                           &
          = cwtr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      cwtr_times(j+1,j) = cmplx(0.0_wp,0.0_wp,wp)
    end if
  end do

end function cwtr_times


!    Name : cwur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cwur_times(a,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=wp)             :: cwur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  cwur_times = cmplx(0.0_wp,0.0_wp,wp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cwur_times,1,0)
#else
    cwur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),cwur_times(hi:ti,:),1,1)
#else
          cwur_times(hi:ti,:)                                                 &
        = cwur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function cwur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cwru_times(a,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=wp)             :: cwru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  cwru_times = cmplx(0.0_wp,0.0_wp,wp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cwru_times,1,0)
#else
    cwru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),cwru_times(:,hj:tj),1,1)
#else
          cwru_times(:,hj:tj)                                                 &
        = cwru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function cwru_times


!    Name : cwdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cwdg_ludcmp(a,p)
  complex(kind=wp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cwdg_ludcmp


!    Name : cwtr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cwtr_ludcmp(a,p)
  complex(kind=wp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cwtr_ludcmp


!    Name : cwsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``cwsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cwsq_ludcmp(a,p)
  complex(kind=wp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK

  integer  :: i, n, info
  external :: zgetrf

  n = size(a,1); p = (/(i,i=1,n)/)
  call zgetrf(n, n, a, n, p, info)

#else

  complex(kind=wp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine cwsq_ludcmp


!    Name : cwdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cwdg_ainvbf(a,p,b)
  complex(kind=wp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  complex(kind=wp)            :: cwdg_ainvbf(1:size(a,1))

  call cwdg_ainvb(a,p,b,cwdg_ainvbf)
end function cwdg_ainvbf


!    Name : cwdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cwdg_ainvb(a,p,b,v)
  complex(kind=wp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  complex(kind=wp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine cwdg_ainvb


!    Name : cwtr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cwtr_ainvbf(a,p,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=wp)             :: cwtr_ainvbf(1:size(a,1),1:size(a,1))

  call cwtr_ainvb(a,p,b,cwtr_ainvbf)
end function cwtr_ainvbf


!    Name : cwtr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cwtr_ainvb(a,p,b,c)
  complex(kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=wp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  complex(kind=wp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call cwtr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = cmplx(0.0_wp,0.0_wp,wp)
    d = cmplx(0.0_wp,0.0_wp,wp)
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call cwtr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = cmplx(0.0_wp,0.0_wp,wp)
    end if
  end do
end subroutine cwtr_ainvb


!    Name : cwtr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cwtr_ainvbblock(a,p,b,c)
  complex(kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=wp), intent(out) :: c(1:,1:)

  complex(kind=wp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = cmplx(0.0_wp,0.0_wp,wp); m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine cwtr_ainvbblock


!    Name : cwsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cwsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cwsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function cwsq_ainvbf(a,p,b)
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=wp)             :: cwsq_ainvbf(1:size(a,1),1:size(b,2))

  call cwsq_ainvb(a,p,b,cwsq_ainvbf)
end function cwsq_ainvbf


!    Name : cwsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cwsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cwsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cwsq_ainvb(a,p,b,c)
  complex(kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=wp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK

  integer          :: m, n, info, ptmp(size(p))
  complex(kind=wp) :: atmp(size(a,1), size(a,2))
  external         :: zgetrs

  n = size(a,1); m = size(b,2)
  atmp = a; c = b; ptmp = p
  call zgetrs('N', n, m, atmp, n, ptmp, c, n, info)

#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  complex(kind=wp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = cmplx(0.0_wp,0.0_wp,wp)
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = cmplx(0.0_wp,0.0_wp,wp)
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    complex(kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=wp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    complex(kind=wp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=wp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine cwsq_ainvb

#ifdef __USE_TPREC

!    Name : rtdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rtdg_norm1(arg)
  real   (kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: rtdg_norm1

  rtdg_norm1 = maxval(abs(arg))
end function rtdg_norm1


!    Name : rtdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rtdg_normi(arg)
  real   (kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: rtdg_normi

  rtdg_normi = maxval(abs(arg))
end function rtdg_normi


!    Name : rtdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rtdg_normf(arg)
  real   (kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: rtdg_normf
  real   (kind=tp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_tp)
  rtdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rtdg_normf


!    Name : rtsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rtsq_norm1(arg)
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rtsq_norm1

  integer       :: j
  real(kind=tp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  rtsq_norm1 = maxval(colsum)
end function rtsq_norm1


!    Name : rtsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rtsq_normi(arg)
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rtsq_normi

  integer       :: i
  real(kind=tp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  rtsq_normi = maxval(rowsum)
end function rtsq_normi


!    Name : rtsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rtsq_normf(arg)
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rtsq_normf
  real   (kind=tp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_tp)
  rtsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rtsq_normf


!    Name : rttr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rttr_norm1(arg)
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rttr_norm1

  rttr_norm1 = rtsq_norm1(arg)
end function rttr_norm1


!    Name : rttr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rttr_normi(arg)
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rttr_normi

  rttr_normi = rtsq_normi(arg)
end function rttr_normi


!    Name : rttr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rttr_normf(arg)
  real   (kind=tp), intent(in) :: arg(:,:)
  real   (kind=tp)             :: rttr_normf

  rttr_normf = rtsq_normf(arg)
end function rttr_normf


!    Name : rtdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function rtdg_eye(n)
  integer, intent(in) :: n
  real   (kind=tp)    :: rtdg_eye(1:n)

  rtdg_eye = 1.0_tp
end function rtdg_eye


!    Name : rtsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function rtsq_eye(n)
  integer, intent(in) :: n
  real   (kind=tp)    :: rtsq_eye(1:n,1:n)

  integer :: i

  rtsq_eye = 0.0_tp
  do i=1,n
    rtsq_eye(i,i) = 1.0_tp
  end do
end function rtsq_eye


!    Name : rttr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function rttr_eye(n)
  integer, intent(in) :: n
  real   (kind=tp)    :: rttr_eye(1:n,1:n)

  integer :: i

  rttr_eye = 0.0_tp
  do i=1,n
    rttr_eye(i,i) = 1.0_tp
  end do
end function rttr_eye


!    Name : rtdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rtdg_times(a,b)
  real   (kind=tp), intent(in) :: a(1:), b(1:)
  real   (kind=tp)             :: rtdg_times(size(a))

  rtdg_times = a * b
end function rtdg_times


!    Name : rtsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rtsq_times(a,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rtsq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,rtsq_times,1,0)
#else
  rtsq_times = matmul(a,b)
#endif
end function rtsq_times


!    Name : rttr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rttr_times(a,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rttr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  rttr_times = 0.0_tp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rttr_times,1,0)
#else
    rttr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            rttr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            rttr_times(h(i):t(i), h(j):t(j))                           &
          = rttr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      rttr_times(j+1,j) = 0.0_tp
    end if
  end do

end function rttr_times


!    Name : rtur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rtur_times(a,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rtur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  rtur_times = 0.0_tp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rtur_times,1,0)
#else
    rtur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),rtur_times(hi:ti,:),1,1)
#else
          rtur_times(hi:ti,:)                                                 &
        = rtur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function rtur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rtru_times(a,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rtru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  rtru_times = 0.0_tp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rtru_times,1,0)
#else
    rtru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),rtru_times(:,hj:tj),1,1)
#else
          rtru_times(:,hj:tj)                                                 &
        = rtru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function rtru_times


!    Name : rtdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rtdg_ludcmp(a,p)
  real   (kind=tp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rtdg_ludcmp


!    Name : rttr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rttr_ludcmp(a,p)
  real   (kind=tp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rttr_ludcmp


!    Name : rtsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``rtsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rtsq_ludcmp(a,p)
  real   (kind=tp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK


#else

  real   (kind=tp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine rtsq_ludcmp


!    Name : rtdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rtdg_ainvbf(a,p,b)
  real   (kind=tp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  real   (kind=tp)            :: rtdg_ainvbf(1:size(a,1))

  call rtdg_ainvb(a,p,b,rtdg_ainvbf)
end function rtdg_ainvbf


!    Name : rtdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rtdg_ainvb(a,p,b,v)
  real   (kind=tp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  real   (kind=tp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine rtdg_ainvb


!    Name : rttr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rttr_ainvbf(a,p,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=tp)             :: rttr_ainvbf(1:size(a,1),1:size(a,1))

  call rttr_ainvb(a,p,b,rttr_ainvbf)
end function rttr_ainvbf


!    Name : rttr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rttr_ainvb(a,p,b,c)
  real   (kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=tp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  real   (kind=tp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call rttr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = 0.0_tp
    d = 0.0_tp
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call rttr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = 0.0_tp
    end if
  end do
end subroutine rttr_ainvb


!    Name : rttr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rttr_ainvbblock(a,p,b,c)
  real   (kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=tp), intent(out) :: c(1:,1:)

  real   (kind=tp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = 0.0_tp; m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine rttr_ainvbblock


!    Name : rtsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rtsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rtsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function rtsq_ainvbf(a,p,b)
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=tp)             :: rtsq_ainvbf(1:size(a,1),1:size(b,2))

  call rtsq_ainvb(a,p,b,rtsq_ainvbf)
end function rtsq_ainvbf


!    Name : rtsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rtsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rtsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rtsq_ainvb(a,p,b,c)
  real   (kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=tp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK


#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  real   (kind=tp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = 0.0_tp
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = 0.0_tp
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    real   (kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=tp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    real   (kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=tp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine rtsq_ainvb


!    Name : ctdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function ctdg_norm1(arg)
  complex(kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: ctdg_norm1

  ctdg_norm1 = maxval(abs(arg))
end function ctdg_norm1


!    Name : ctdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function ctdg_normi(arg)
  complex(kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: ctdg_normi

  ctdg_normi = maxval(abs(arg))
end function ctdg_normi


!    Name : ctdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function ctdg_normf(arg)
  complex(kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: ctdg_normf
  real   (kind=tp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_tp)
  ctdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function ctdg_normf


!    Name : ctsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function ctsq_norm1(arg)
  complex(kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: ctsq_norm1

  integer       :: j
  real(kind=tp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  ctsq_norm1 = maxval(colsum)
end function ctsq_norm1


!    Name : ctsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function ctsq_normi(arg)
  complex(kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: ctsq_normi

  integer       :: i
  real(kind=tp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  ctsq_normi = maxval(rowsum)
end function ctsq_normi


!    Name : ctsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function ctsq_normf(arg)
  complex(kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: ctsq_normf
  real   (kind=tp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_tp)
  ctsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function ctsq_normf


!    Name : cttr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cttr_norm1(arg)
  complex(kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: cttr_norm1

  cttr_norm1 = ctsq_norm1(arg)
end function cttr_norm1


!    Name : cttr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cttr_normi(arg)
  complex(kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: cttr_normi

  cttr_normi = ctsq_normi(arg)
end function cttr_normi


!    Name : cttr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cttr_normf(arg)
  complex(kind=tp), intent(in) :: arg(:,:)
  real   (kind=tp)             :: cttr_normf

  cttr_normf = ctsq_normf(arg)
end function cttr_normf


!    Name : ctdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function ctdg_eye(n)
  integer, intent(in) :: n
  complex(kind=tp)    :: ctdg_eye(1:n)

  ctdg_eye = cmplx(1.0_tp,0.0_tp,tp)
end function ctdg_eye


!    Name : ctsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function ctsq_eye(n)
  integer, intent(in) :: n
  complex(kind=tp)    :: ctsq_eye(1:n,1:n)

  integer :: i

  ctsq_eye = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,n
    ctsq_eye(i,i) = cmplx(1.0_tp,0.0_tp,tp)
  end do
end function ctsq_eye


!    Name : cttr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function cttr_eye(n)
  integer, intent(in) :: n
  complex(kind=tp)    :: cttr_eye(1:n,1:n)

  integer :: i

  cttr_eye = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,n
    cttr_eye(i,i) = cmplx(1.0_tp,0.0_tp,tp)
  end do
end function cttr_eye


!    Name : ctdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function ctdg_times(a,b)
  complex(kind=tp), intent(in) :: a(1:), b(1:)
  complex(kind=tp)             :: ctdg_times(size(a))

  ctdg_times = a * b
end function ctdg_times


!    Name : ctsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function ctsq_times(a,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=tp)             :: ctsq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,ctsq_times,1,0)
#else
  ctsq_times = matmul(a,b)
#endif
end function ctsq_times


!    Name : cttr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cttr_times(a,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=tp)             :: cttr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  cttr_times = cmplx(0.0_tp,0.0_tp,tp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cttr_times,1,0)
#else
    cttr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            cttr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            cttr_times(h(i):t(i), h(j):t(j))                           &
          = cttr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      cttr_times(j+1,j) = cmplx(0.0_tp,0.0_tp,tp)
    end if
  end do

end function cttr_times


!    Name : ctur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function ctur_times(a,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=tp)             :: ctur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  ctur_times = cmplx(0.0_tp,0.0_tp,tp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,ctur_times,1,0)
#else
    ctur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),ctur_times(hi:ti,:),1,1)
#else
          ctur_times(hi:ti,:)                                                 &
        = ctur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function ctur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function ctru_times(a,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=tp)             :: ctru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  ctru_times = cmplx(0.0_tp,0.0_tp,tp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,ctru_times,1,0)
#else
    ctru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),ctru_times(:,hj:tj),1,1)
#else
          ctru_times(:,hj:tj)                                                 &
        = ctru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function ctru_times


!    Name : ctdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine ctdg_ludcmp(a,p)
  complex(kind=tp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine ctdg_ludcmp


!    Name : cttr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cttr_ludcmp(a,p)
  complex(kind=tp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cttr_ludcmp


!    Name : ctsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``ctsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine ctsq_ludcmp(a,p)
  complex(kind=tp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK


#else

  complex(kind=tp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine ctsq_ludcmp


!    Name : ctdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function ctdg_ainvbf(a,p,b)
  complex(kind=tp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  complex(kind=tp)            :: ctdg_ainvbf(1:size(a,1))

  call ctdg_ainvb(a,p,b,ctdg_ainvbf)
end function ctdg_ainvbf


!    Name : ctdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine ctdg_ainvb(a,p,b,v)
  complex(kind=tp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  complex(kind=tp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine ctdg_ainvb


!    Name : cttr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cttr_ainvbf(a,p,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=tp)             :: cttr_ainvbf(1:size(a,1),1:size(a,1))

  call cttr_ainvb(a,p,b,cttr_ainvbf)
end function cttr_ainvbf


!    Name : cttr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cttr_ainvb(a,p,b,c)
  complex(kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=tp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  complex(kind=tp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call cttr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = cmplx(0.0_tp,0.0_tp,tp)
    d = cmplx(0.0_tp,0.0_tp,tp)
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call cttr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = cmplx(0.0_tp,0.0_tp,tp)
    end if
  end do
end subroutine cttr_ainvb


!    Name : cttr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cttr_ainvbblock(a,p,b,c)
  complex(kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=tp), intent(out) :: c(1:,1:)

  complex(kind=tp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = cmplx(0.0_tp,0.0_tp,tp); m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine cttr_ainvbblock


!    Name : ctsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``ctsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``ctsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function ctsq_ainvbf(a,p,b)
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=tp)             :: ctsq_ainvbf(1:size(a,1),1:size(b,2))

  call ctsq_ainvb(a,p,b,ctsq_ainvbf)
end function ctsq_ainvbf


!    Name : ctsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``ctsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``ctsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine ctsq_ainvb(a,p,b,c)
  complex(kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=tp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK


#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  complex(kind=tp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = cmplx(0.0_tp,0.0_tp,tp)
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = cmplx(0.0_tp,0.0_tp,tp)
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    complex(kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=tp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    complex(kind=tp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=tp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine ctsq_ainvb

#endif
#ifdef __USE_QPREC

!    Name : rqdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rqdg_norm1(arg)
  real   (kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: rqdg_norm1

  rqdg_norm1 = maxval(abs(arg))
end function rqdg_norm1


!    Name : rqdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rqdg_normi(arg)
  real   (kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: rqdg_normi

  rqdg_normi = maxval(abs(arg))
end function rqdg_normi


!    Name : rqdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rqdg_normf(arg)
  real   (kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: rqdg_normf
  real   (kind=qp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_qp)
  rqdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rqdg_normf


!    Name : rqsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rqsq_norm1(arg)
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqsq_norm1

  integer       :: j
  real(kind=qp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  rqsq_norm1 = maxval(colsum)
end function rqsq_norm1


!    Name : rqsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rqsq_normi(arg)
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqsq_normi

  integer       :: i
  real(kind=qp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  rqsq_normi = maxval(rowsum)
end function rqsq_normi


!    Name : rqsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rqsq_normf(arg)
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqsq_normf
  real   (kind=qp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_qp)
  rqsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function rqsq_normf


!    Name : rqtr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function rqtr_norm1(arg)
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqtr_norm1

  rqtr_norm1 = rqsq_norm1(arg)
end function rqtr_norm1


!    Name : rqtr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function rqtr_normi(arg)
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqtr_normi

  rqtr_normi = rqsq_normi(arg)
end function rqtr_normi


!    Name : rqtr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function rqtr_normf(arg)
  real   (kind=qp), intent(in) :: arg(:,:)
  real   (kind=qp)             :: rqtr_normf

  rqtr_normf = rqsq_normf(arg)
end function rqtr_normf


!    Name : rqdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function rqdg_eye(n)
  integer, intent(in) :: n
  real   (kind=qp)    :: rqdg_eye(1:n)

  rqdg_eye = 1.0_qp
end function rqdg_eye


!    Name : rqsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function rqsq_eye(n)
  integer, intent(in) :: n
  real   (kind=qp)    :: rqsq_eye(1:n,1:n)

  integer :: i

  rqsq_eye = 0.0_qp
  do i=1,n
    rqsq_eye(i,i) = 1.0_qp
  end do
end function rqsq_eye


!    Name : rqtr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function rqtr_eye(n)
  integer, intent(in) :: n
  real   (kind=qp)    :: rqtr_eye(1:n,1:n)

  integer :: i

  rqtr_eye = 0.0_qp
  do i=1,n
    rqtr_eye(i,i) = 1.0_qp
  end do
end function rqtr_eye


!    Name : rqdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rqdg_times(a,b)
  real   (kind=qp), intent(in) :: a(1:), b(1:)
  real   (kind=qp)             :: rqdg_times(size(a))

  rqdg_times = a * b
end function rqdg_times


!    Name : rqsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rqsq_times(a,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rqsq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,rqsq_times,1,0)
#else
  rqsq_times = matmul(a,b)
#endif
end function rqsq_times


!    Name : rqtr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rqtr_times(a,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rqtr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  rqtr_times = 0.0_qp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rqtr_times,1,0)
#else
    rqtr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            rqtr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            rqtr_times(h(i):t(i), h(j):t(j))                           &
          = rqtr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      rqtr_times(j+1,j) = 0.0_qp
    end if
  end do

end function rqtr_times


!    Name : rqur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rqur_times(a,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rqur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  rqur_times = 0.0_qp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rqur_times,1,0)
#else
    rqur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),rqur_times(hi:ti,:),1,1)
#else
          rqur_times(hi:ti,:)                                                 &
        = rqur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function rqur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function rqru_times(a,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rqru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  rqru_times = 0.0_qp

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,rqru_times,1,0)
#else
    rqru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),rqru_times(:,hj:tj),1,1)
#else
          rqru_times(:,hj:tj)                                                 &
        = rqru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function rqru_times


!    Name : rqdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rqdg_ludcmp(a,p)
  real   (kind=qp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rqdg_ludcmp


!    Name : rqtr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine rqtr_ludcmp(a,p)
  real   (kind=qp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rqtr_ludcmp


!    Name : rqsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``rqsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rqsq_ludcmp(a,p)
  real   (kind=qp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK


#else

  real   (kind=qp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine rqsq_ludcmp


!    Name : rqdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rqdg_ainvbf(a,p,b)
  real   (kind=qp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  real   (kind=qp)            :: rqdg_ainvbf(1:size(a,1))

  call rqdg_ainvb(a,p,b,rqdg_ainvbf)
end function rqdg_ainvbf


!    Name : rqdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rqdg_ainvb(a,p,b,v)
  real   (kind=qp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  real   (kind=qp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine rqdg_ainvb


!    Name : rqtr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function rqtr_ainvbf(a,p,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=qp)             :: rqtr_ainvbf(1:size(a,1),1:size(a,1))

  call rqtr_ainvb(a,p,b,rqtr_ainvbf)
end function rqtr_ainvbf


!    Name : rqtr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rqtr_ainvb(a,p,b,c)
  real   (kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=qp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  real   (kind=qp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call rqtr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = 0.0_qp
    d = 0.0_qp
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call rqtr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = 0.0_qp
    end if
  end do
end subroutine rqtr_ainvb


!    Name : rqtr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine rqtr_ainvbblock(a,p,b,c)
  real   (kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=qp), intent(out) :: c(1:,1:)

  real   (kind=qp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = 0.0_qp; m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine rqtr_ainvbblock


!    Name : rqsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rqsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rqsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function rqsq_ainvbf(a,p,b)
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  real   (kind=qp)             :: rqsq_ainvbf(1:size(a,1),1:size(b,2))

  call rqsq_ainvb(a,p,b,rqsq_ainvbf)
end function rqsq_ainvbf


!    Name : rqsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``rqsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``rqsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine rqsq_ainvb(a,p,b,c)
  real   (kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  real   (kind=qp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK


#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  real   (kind=qp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = 0.0_qp
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = 0.0_qp
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    real   (kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=qp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    real   (kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
    real   (kind=qp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine rqsq_ainvb


!    Name : cqdg_norm1
! Purpose : This function computes the matrix 1-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cqdg_norm1(arg)
  complex(kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: cqdg_norm1

  cqdg_norm1 = maxval(abs(arg))
end function cqdg_norm1


!    Name : cqdg_normi
! Purpose : This function computes the matrix infinity-norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cqdg_normi(arg)
  complex(kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: cqdg_normi

  cqdg_normi = maxval(abs(arg))
end function cqdg_normi


!    Name : cqdg_normf
! Purpose : This function computes Frobenius norm of a diagonal matrix
!         : according to its definition.
!   Input : ``arg'' is a diagonal matrix stored in a one-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cqdg_normf(arg)
  complex(kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: cqdg_normf
  real   (kind=qp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_qp)
  cqdg_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function cqdg_normf


!    Name : cqsq_norm1
! Purpose : This function computes the matrix 1-norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cqsq_norm1(arg)
  complex(kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: cqsq_norm1

  integer       :: j
  real(kind=qp) :: colsum(1:size(arg,2))

  do j=1,size(arg,2)
    colsum(j) = sum(abs(arg(:,j)))
  end do
  cqsq_norm1 = maxval(colsum)
end function cqsq_norm1


!    Name : cqsq_normi
! Purpose : This function computes the matrix infinity-norm of
!         : a rectangular matrix according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cqsq_normi(arg)
  complex(kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: cqsq_normi

  integer       :: i
  real(kind=qp) :: rowsum(1:size(arg,1))

  do i=1,size(arg,1)
    rowsum(i) = sum(abs(arg(i,:)))
  end do
  cqsq_normi = maxval(rowsum)
end function cqsq_normi


!    Name : cqsq_normf
! Purpose : This function computes Frobenius norm of a rectangular matrix
!         : according to its definition.
!   Input : ``arg'' is a rectangular matrix stored in a two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cqsq_normf(arg)
  complex(kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: cqsq_normf
  real   (kind=qp)             :: vmax

  vmax = max(maxval(abs(arg)),1.0_qp)
  cqsq_normf = vmax * sqrt(sum(abs(arg/vmax)**2))
end function cqsq_normf


!    Name : cqtr_norm1
! Purpose : This function computes the matrix 1-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix 1-norm of ``arg''.

pure function cqtr_norm1(arg)
  complex(kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: cqtr_norm1

  cqtr_norm1 = cqsq_norm1(arg)
end function cqtr_norm1


!    Name : cqtr_normi
! Purpose : This function computes the matrix infinity-norm of an upper
!         : quasi-triangular matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored
!         : in a square array.
!  Output : The return value is the matrix infinity-norm of ``arg''.

pure function cqtr_normi(arg)
  complex(kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: cqtr_normi

  cqtr_normi = cqsq_normi(arg)
end function cqtr_normi


!    Name : cqtr_normf
! Purpose : This function computes Frobenius norm of an upper quasi-triangular
!         : matrix according to its definition.
!   Input : ``arg'' is an upper quasi-triangular matrix stored in a
!         : two-dimensional array.
!  Output : The return value is Frobenius norm of ``arg''.

pure function cqtr_normf(arg)
  complex(kind=qp), intent(in) :: arg(:,:)
  real   (kind=qp)             :: cqtr_normf

  cqtr_normf = cqsq_normf(arg)
end function cqtr_normf


!    Name : cqdg_eye
! Purpose : This function generates the diagonal identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the diagonal identity stored
!         : in a one-dimensional array

pure function cqdg_eye(n)
  integer, intent(in) :: n
  complex(kind=qp)    :: cqdg_eye(1:n)

  cqdg_eye = cmplx(1.0_qp,0.0_qp,qp)
end function cqdg_eye


!    Name : cqsq_eye
! Purpose : This function generates the square identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the square identity stored
!         : in a two-dimensional array

pure function cqsq_eye(n)
  integer, intent(in) :: n
  complex(kind=qp)    :: cqsq_eye(1:n,1:n)

  integer :: i

  cqsq_eye = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,n
    cqsq_eye(i,i) = cmplx(1.0_qp,0.0_qp,qp)
  end do
end function cqsq_eye


!    Name : cqtr_eye
! Purpose : This function generates the upper triangular identity.
!   Input : The integer, ``n'', is the size of the identity.
!  Output : The return value is the upper triangular identity stored
!         : in a two-dimensional array

pure function cqtr_eye(n)
  integer, intent(in) :: n
  complex(kind=qp)    :: cqtr_eye(1:n,1:n)

  integer :: i

  cqtr_eye = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,n
    cqtr_eye(i,i) = cmplx(1.0_qp,0.0_qp,qp)
  end do
end function cqtr_eye


!    Name : cqdg_times
!   Usage : This function is invoked through the generic name, ``.dgtimes.''.
! Purpose : This function computes a diagonal matrix
!         : times a diagonal matrix.
!   Input : A and B are diagonal matrices stored in one-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cqdg_times(a,b)
  complex(kind=qp), intent(in) :: a(1:), b(1:)
  complex(kind=qp)             :: cqdg_times(size(a))

  cqdg_times = a * b
end function cqdg_times


!    Name : cqsq_times
!   Usage : This function is invoked through the generic name, ``.sqtimes.''.
! Purpose : This function computes a rectangular matrix
!         : times a rectangular matrix.
!   Input : A and B are rectangular matrices stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cqsq_times(a,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=qp)             :: cqsq_times(1:size(a,1),1:size(b,2))

#ifdef __USE_BLAS
  call ggemm(a,b,cqsq_times,1,0)
#else
  cqsq_times = matmul(a,b)
#endif
end function cqsq_times


!    Name : cqtr_times
!   Usage : This function is invoked through the generic name, ``.trtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A and B are upper quasi-triangular matrices stored in two-
!         : dimensional arrays. A and B must have the same block structure.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cqtr_times(a,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=qp)             :: cqtr_times(1:size(a,1),1:size(a,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))     ! heads and tails

  n = size(a,1)
  cqtr_times = cmplx(0.0_qp,0.0_qp,qp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cqtr_times,1,0)
#else
    cqtr_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product AB is
    ! computed by utilizing the level-3 BLAS.

    do i=1,nblk
      do j=i,nblk
        do k=i,j
#ifdef __USE_BLAS
          call ggemm(                                                  &
            a(h(i):t(i),h(k):t(k)),                                    &
            b(h(k):t(k),h(j):t(j)),                                    &
            cqtr_times(h(i):t(i), h(j):t(j)),1,1                       &
          )
#else
            cqtr_times(h(i):t(i), h(j):t(j))                           &
          = cqtr_times(h(i):t(i), h(j):t(j))                           &
          + matmul( a(h(i):t(i),h(k):t(k)), b(h(k):t(k),h(j):t(j)) )
#endif
        end do
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                               &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                        &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                              &
    ) then
      cqtr_times(j+1,j) = cmplx(0.0_qp,0.0_qp,qp)
    end if
  end do

end function cqtr_times


!    Name : cqur_times
!   Usage : This function is invoked through the generic name, ``.urtimes.''.
! Purpose : This function computes an upper quasi-triangular matrix
!         : times a rectangular matrix.
!   Input : A is an upper quasi-triangular matrix, and B is a rectangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cqur_times(a,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=qp)             :: cqur_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(a,1))
  integer            :: hi, ti, hk, tk, h(1:size(a,1)), t(1:size(a,1))

  n = size(a,1)
  cqur_times = cmplx(0.0_qp,0.0_qp,qp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(a(j+1,j))) <= abs(a(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cqur_times,1,0)
#else
    cqur_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix A. The matrix A
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The rows of the matrix B are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the i-th partition of the rows of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do i=1,nblk
      hi = h(i)
      ti = t(i)
      do k=i,nblk
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(hi:ti,hk:tk),b(hk:tk,:),cqur_times(hi:ti,:),1,1)
#else
          cqur_times(hi:ti,:)                                                 &
        = cqur_times(hi:ti,:) + matmul( a(hi:ti,hk:tk),b(hk:tk,:) )
#endif
      end do
    end do
  end if
end function cqur_times


!    Name : .rutimes. ( a generic name )
!   Usage : This function is invoked through the generic name, ``.rutimes.''.
! Purpose : This function computes a rectangular matrix
!         : times an upper quasi-triangular matrix.
!   Input : A is a rectangular matrix, and B is an upper quasi-triangular
!         : matrix, both of which are stored in two-dimensional arrays.
!  Output : The result is the product of A and B. (This is not the Hadamard
!         : product but the matrix-matrix multiplication.)

pure function cqru_times(a,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  complex(kind=qp)             :: cqru_times(1:size(a,1),1:size(b,2))

  integer, parameter :: blksize = mpt_blksize
  integer            :: i, j, k, m, n, acts, nblk, mask(1:size(b,1))
  integer            :: hj, tj, hk, tk, h(1:size(b,1)), t(1:size(b,1))

  n = size(b,1)
  cqru_times = cmplx(0.0_qp,0.0_qp,qp)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of B is not 0.

  mask = 0
  do j=1,n-1
    if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product is computed without blocking.

#ifdef __USE_BLAS
    call ggemm(a,b,cqru_times,1,0)
#else
    cqru_times = matmul(a,b)
#endif

  else

    ! The product is computed after blocking the matrix B. The matrix B
    ! is block-decomposed so that each 2-by-2 diagonal block is not separated.
    ! The columns of the matrix A are partitioned consistently. The array ``h''
    ! stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the j-th partition of the columns of the product
    ! AB is computed by utilizing the level-3 BLAS.

    do j=1,nblk
      hj = h(j)
      tj = t(j)
      do k=1,j
        hk = h(k)
        tk = t(k)
#ifdef __USE_BLAS
        call ggemm(a(:,hk:tk),b(hk:tk,hj:tj),cqru_times(:,hj:tj),1,1)
#else
          cqru_times(:,hj:tj)                                                 &
        = cqru_times(:,hj:tj) + matmul(a(:,hk:tk),b(hk:tk,hj:tj))
#endif
      end do
    end do
  end if
end function cqru_times


!    Name : cqdg_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a diagonal matrix.
!   Input : - ``A'' is a one-dimensional array that stores a diagonal matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cqdg_ludcmp(a,p)
  complex(kind=qp), intent(inout) :: a(1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cqdg_ludcmp


!    Name : cqtr_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for an upper quasi-triangular matrix.
!   Input : - ``A'' is a two-dimensional array that stores an upper
!         :   quasi-triangular matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting.
!    Note : This routine is formal and does nothing.

pure subroutine cqtr_ludcmp(a,p)
  complex(kind=qp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

  integer :: i
  i = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cqtr_ludcmp


!    Name : cqsq_ludcmp
! Purpose : This subroutine performs an LU factorization
!         : for a square matrix.
!   Input : - ``A'' is a two-dimensional array that stores a square matrix.
!         : - ``p'' is a one-dimensional array of integer.
!  Output : - ``A'' is overwritten by its LU factors.
!         : - ``p'' stores the information on pivoting, which is used
!         :   by ``cqsq_ainvbf'' defined later in this file.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cqsq_ludcmp(a,p)
  complex(kind=qp), intent(inout) :: a(1:,1:)
  integer,          intent(out  ) :: p(1:)

#ifdef __USE_LAPACK


#else

  complex(kind=qp) :: u(1:size(a,1))
  integer          :: i,j,k,m(1),n

  n = size(a,1)
  p = (/(i,i=1,n)/)

  do i=1,n-1                                     ! Gaussian elimination
    m = maxloc(abs(a(i:n,i))); j = i+m(1)-1      ! determine the pivot
    if (i /= j) then                             ! row permutation if necessary
      u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
      k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
    end if
    a(i+1:n,i) = a(i+1:n,i) / a(i,i)             ! the lower factor
    do k=i+1,n                                   ! update the rest
      a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
    end do
  end do

#endif
end subroutine cqsq_ludcmp


!    Name : cqdg_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cqdg_ainvbf(a,p,b)
  complex(kind=qp),intent(in) :: a(:), b(:)
  integer,         intent(in) :: p(:)
  complex(kind=qp)            :: cqdg_ainvbf(1:size(a,1))

  call cqdg_ainvb(a,p,b,cqdg_ainvbf)
end function cqdg_ainvbf


!    Name : cqdg_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are diagonal matrices of the same size
!         :   stored in one-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``V'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cqdg_ainvb(a,p,b,v)
  complex(kind=qp), intent(in ) :: a(:), b(:)
  integer,          intent(in ) :: p(:)
  complex(kind=qp), intent(out) :: v(:)

  v(1) = p(1)                                    ! to avoid compiler's warning
  v = b / a
end subroutine cqdg_ainvb


!    Name : cqtr_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : The return value is A^{-1}B.
!    Note : ``p'' is not used.

pure function cqtr_ainvbf(a,p,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=qp)             :: cqtr_ainvbf(1:size(a,1),1:size(a,1))

  call cqtr_ainvb(a,p,b,cqtr_ainvbf)
end function cqtr_ainvbf


!    Name : cqtr_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by utilizing a block algorithm.
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices of the same
!         :   size and of the same block structure. Both of them are stored in
!         :   two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cqtr_ainvb(a,p,b,c)
  complex(kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=qp), intent(out) :: c(1:,1:)

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,acts,nblk
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: mask(1:size(a,1))
  integer            :: h(1:size(a,1)), t(1:size(a,1))    ! heads and tails
  complex(kind=qp)   :: d(1:size(a,1),1:size(a,1))

  n = size(a,1)

  ! The value of mask(j) is one, when the (j+1,j)-th entry of A or B is not 0.

  mask = 0
  do j=1,n-1
    if ( tiny(abs(a(j+1,j))) <= abs(a(j+1,j)) ) mask(j) = 1
    if ( tiny(abs(b(j+1,j))) <= abs(b(j+1,j)) ) mask(j) = 1
  end do

  if (n <= blksize) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    call cqtr_ainvbblock(a,p,b,c)

  else

    ! The product A^{-1}B is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed so that each 2-by-2 diagonal block is not
    ! separated. The matrices are assumed to have the same structure. The array
    ! ``h'' stores the heads of blocks and ``t'' stores the tails of blocks.

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    m = 0
    do i=1,nblk-1
      h(i) = m + 1
      m = m + acts - mask(m + acts)
      t(i) = m
    end do
    h(nblk) = m + 1; t(nblk) = n

    ! In the following loop, the (i,j)-th block of the product A^{-1}B is
    ! computed by utilizing the level-3 BLAS and a scalar algorithm for
    ! inversion.

    c = cmplx(0.0_qp,0.0_qp,qp)
    d = cmplx(0.0_qp,0.0_qp,qp)
    do i=nblk,1,-1
      hi = h(i); ti = t(i)
      do j=i,nblk
        hj = h(j); tj = t(j)
        do k=i+1,j
          hk = h(k); tk = t(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),d(hi:ti,hj:tj),-1,1)
#else
            d(hi:ti,hj:tj)                                                  &
          = d(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + b(hi:ti,hj:tj)

        call cqtr_ainvbblock(a(hi:ti,hi:ti),p,d(hi:ti,hj:tj),c(hi:ti,hj:tj))
      end do
    end do
  end if

  ! In order to preserve the block structure, 0 is substituted into the
  ! (j+1,j)-th entry of the product, when both A(j+1,j) and B(j+1,j) are 0.

  do j=1,n-1
    if (                                                                    &
      abs(a(j+1,j)) < tiny(abs(a(j+1,j))) .and.                             &
      abs(b(j+1,j)) < tiny(abs(b(j+1,j)))                                   &
    ) then
      c(j+1,j) = cmplx(0.0_qp,0.0_qp,qp)
    end if
  end do
end subroutine cqtr_ainvb


!    Name : cqtr_ainvbblock
! Purpose : This subroutine computes an inverse matrix times a matrix
!         : by a scalar algorithm.
!   Input : - ``A'' is an upper quasi-triangular matrix.
!         : - ``B'' is a rectangular matrix.
!         : Both A and B are stored in two-dimensional arrays.
!         : - ``p'' stores the information on pivoting.
!  Output : ``C'' stores A^{-1}B.
!    Note : ``p'' is not used.

pure subroutine cqtr_ainvbblock(a,p,b,c)
  complex(kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=qp), intent(out) :: c(1:,1:)

  complex(kind=qp) :: blk(2,2), vec(2,1:size(b,2))
  integer          :: i, m, subd

  ! To avoid compiler's warning,

  c = cmplx(0.0_qp,0.0_qp,qp); m = p(1); m = size(a,1)

  ! In the following loop, the solution, ``C'', is computed
  ! in descending order of rows.

  i = m
  do while (1 <= i)

    subd = 0
    if (1 < i) then
      if ( tiny(abs(a(i,i-1))) <= abs(a(i,i-1)) ) subd = 1
    end if

    if ( i == 1 .or. subd == 0 ) then             ! for a 1x1 diagonal block

      if (i == m) then
        c(i,:) = b(i,:) / a(i,i)
      else
#ifdef __USE_BLAS
        c(i:i,:) = b(i:i,:)    
        call ggemm(a(i:i,i+1:m),c(i+1:m,:),c(i:i,:),-1,1)
        c(i:i,:) = c(i:i,:) / a(i,i)
#else
        c(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), c(i+1:m,:) ) ) / a(i,i)
#endif
      end if
      i = i - 1

    else                                          ! for a 2x2 diagonal block

      blk(1,1) =  a(i,i  ); blk(1,2) = -a(i-1,i  )
      blk(2,1) = -a(i,i-1); blk(2,2) =  a(i-1,i-1)
      blk = blk / (a(i-1,i-1)*a(i,i) - a(i-1,i)*a(i,i-1))

      vec(1,:) = b(i-1,:)
      vec(2,:) = b(i  ,:)
      if (i /= m) then
#ifdef __USE_BLAS
        call ggemm(a(i-1:i-1,i+1:m),c(i+1:m,:),vec(1:1,:),-1,1)
        call ggemm(a(i  :i  ,i+1:m),c(i+1:m,:),vec(2:2,:),-1,1)
#else
        vec(1:1,:) = vec(1:1,:) - matmul( a(i-1:i-1,i+1:m), c(i+1:m,:) )
        vec(2:2,:) = vec(2:2,:) - matmul( a(i  :i  ,i+1:m), c(i+1:m,:) )
#endif
      end if
      c(i-1,:) = blk(1,1) * vec(1,:) + blk(1,2) * vec(2,:)
      c(i  ,:) = blk(2,1) * vec(1,:) + blk(2,2) * vec(2,:)
      i = i - 2

    end if
  end do
end subroutine cqtr_ainvbblock


!    Name : cqsq_ainvbf
! Purpose : This function computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cqsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cqsq_ludcmpf'' implemented in this file.
!  Output : The return value is A^{-1}B.

pure function cqsq_ainvbf(a,p,b)
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: p(1:)
  complex(kind=qp)             :: cqsq_ainvbf(1:size(a,1),1:size(b,2))

  call cqsq_ainvb(a,p,b,cqsq_ainvbf)
end function cqsq_ainvbf


!    Name : cqsq_ainvb
! Purpose : This subroutine computes an inverse matrix times a matrix.
!   Input : - ``A'' is a square matrix and stores LU factors computed by
!         :   ``cqsq_ludcmpf'' implemented in this file.
!         : - ``B'' is a rectangular matrix.
!         : - ``p'' stores the information on pivoting provided
!         :   by ``cqsq_ludcmpf'' implemented in this file.
!  Output : ``C'' stores A^{-1}B.
!    Note : If you define ``-Dpure= -D__USE_LAPACK'' on the command line of
!         : your FORTRAN compiler, then, LAPACK routines are used instead of
!         : the implementation here.

pure subroutine cqsq_ainvb(a,p,b,c)
  complex(kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
  integer,          intent(in ) :: p(1:)
  complex(kind=qp), intent(out) :: c(1:,1:)

#ifdef __USE_LAPACK


#else

  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,ell,m,mblk,n,nblk,acts
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))      ! heads and tails
  integer            :: hb(1:size(b,2)),tb(1:size(b,2))      ! heads and tails
  complex(kind=qp)   :: d (1:size(a,1),1:size(b,2))          ! working area
 
  m = size(a,1); n = size(b,2)
  do i=1,m
    d(i,:) = b(p(i),:)                                       ! undo permutation
  end do

  if ( m <= blksize ) then

    ! For small matrices, the product A^{-1}B is computed without blocking.

    c = d; call slower(a,c,d)
    c = d; call supper(a,c,d)
    c = d

  else

    ! The product is computed after blocking the matrices. The matrices,
    ! A and B, are block-decomposed consistently. Since LU factors are
    ! strictly lower/upper triangular, it is not necessary to take into
    ! account the structure of diagonal blocks. The array ``ha'' stores
    ! the heads of blocks of A and ``ta'' stores the tails of blocks of A.
    ! The arrays, ``hb'' and ``tb'', define a column partition of B in the
    ! same way.

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! L^{-1}B, is computed and stored in C by utilizing the level-3 BLAS.

    c = cmplx(0.0_qp,0.0_qp,qp)
    do i=1,mblk
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call slower(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

    ! When A is factorized as A = LU, the (i,j)-th block of the product,
    ! U^{-1}(L^{-1}B), is computed and stored in C by utilizing the
    ! level-3 BLAS.

    d = c; c = cmplx(0.0_qp,0.0_qp,qp)
    do i=mblk,1,-1
      hi = ha(i); ti = ta(i)
      do j=1,nblk
        hj = hb(j); tj = tb(j)
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),c(hk:tk,hj:tj),c(hi:ti,hj:tj),-1,1)
#else
            c(hi:ti,hj:tj)                                                   &
          = c(hi:ti,hj:tj) - matmul( a(hi:ti,hk:tk), c(hk:tk,hj:tj) )
#endif
        end do  
        c(hi:ti,hj:tj) = c(hi:ti,hj:tj) + d(hi:ti,hj:tj)
        call supper(a(hi:ti,hi:ti),c(hi:ti,hj:tj),d(hi:ti,hj:tj))
        c(hi:ti,hj:tj) = d(hi:ti,hj:tj)
      end do
    end do

  end if

contains

  !    Name : slower
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is a lower triangular matrix
  !         :   (the lower factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine slower(a,b,x)
    complex(kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=qp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(1,:) = b(1,:)
    do i=2,m
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,1:i-1),x(1:i-1,:),x(i:i,:),-1,1)
#else
      x(i:i,:) = b(i:i,:) - matmul( a(i:i,1:i-1), x(1:i-1,:) )
#endif
    end do
  end subroutine slower


  !    Name : supper
  ! Purpose : This subroutine computes an inverse matrix times a matrix.
  !   Input : - ``A'' is an upper triangular matrix
  !         :   (the upper factor of LU-factorization).
  !         : - ``B'' is a rectangular matrix.
  !  Output : ``X'' stores A^{-1}B.

  pure subroutine supper(a,b,x)
    complex(kind=qp), intent(in ) :: a(1:,1:), b(1:,1:)
    complex(kind=qp), intent(out) :: x(1:,1:)
    integer :: i, m
    m = size(a,1); x(m,:) = b(m,:) / a(m,m)
    do i=m-1,1,-1
#ifdef __USE_BLAS
      x(i:i,:) = b(i:i,:)
      call ggemm(a(i:i,i+1:m),x(i+1:m,:),x(i:i,:),-1,1)
      x(i:i,:) = x(i:i,:) / a(i,i)
#else
      x(i:i,:) = ( b(i:i,:) - matmul( a(i:i,i+1:m), x(i+1:m,:) ) ) / a(i,i)
#endif
    end do
  end subroutine supper

#endif

end subroutine cqsq_ainvb

#endif

end module matrixpwrtag


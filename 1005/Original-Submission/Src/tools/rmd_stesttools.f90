MODULE RMD_STESTTOOLS
  ! A general module for testing BLAS-RMD functions.
  ! Note: must have lowercase name for rmd_convert_type to work.
  !private
  !public :: fdata, rmd_stestf, rmd_stestrandom, rmd_strishape_m2v, rmd_strishape_v2m,&
  !  & rmd_sbandshape_m2v, rmd_sbandshape_v2m, rmd_srandom, rmd_stolerance, F_interf

  type fdata
    ! A custom datatype that has to be defined for each function to be
    ! tested. An instance of fdata contains the data needed by F and
    ! F_rmd hidden from rmd_stestf, eg. constants and matrix sizes.
    real, allocatable      :: reals(:)
    integer, allocatable   :: integers(:)
    character, allocatable :: characters(:)
  end type fdata

  abstract interface
    subroutine F_dp_interf(testcase, u, v)
      import fdata
      type(fdata), intent(inout) :: testcase
      double precision, intent(in)  :: u(*)
      double precision, intent(out) :: v(*)
    end subroutine F_dp_interf
    subroutine F_interf(testcase, u, v)
      import fdata
      type(fdata), intent(inout) :: testcase
      real, intent(in)           :: u(*)
      real, intent(out)          :: v(*)
    end subroutine F_INTERF
    subroutine F_rmd_interf(testcase, u, v, ua, va)
      import fdata
      type(fdata), intent(inout) :: testcase
      real, intent(in)        :: u(*), v(*), va(*)
      real, intent(inout)     :: ua(*)
    end subroutine F_rmd_interf
  end interface

  abstract interface
    subroutine set_KJ_interf(flag, K, J)
      integer, intent(in) :: flag
      integer, intent(out), allocatable :: K(:), J(:)
    end subroutine set_KJ_interf
  end interface

  interface rmd_sapprox_eq
    module procedure eqs, eqv, eqm, eqsv, eqvs, eqsm, eqms
  end interface rmd_sapprox_eq

  ! Finite difference step for numerical gradients, 1.8e-5 for IEEE real,
  ! 0.0148 for IEEE double:
  real, parameter :: rmd_sdelta = epsilon(0.0)**(1./3)*3
  ! Tolerance for Jacobian comparison, 0.01 for real, 1.5e-8 for double
  real, parameter :: rmd_stolerance = epsilon(0.0)**(2./3)*400
  
CONTAINS

  ! -----APPROXIMATELY EQUAL FUNCTIONS------
  logical function eqs(x, y) result(aeq)
    real :: x, y, maxv
    maxv = max(1.0, abs(x), abs(y))
    aeq = abs(x - y) < maxv*rmd_stolerance
  end function eqs
    
  logical function eqv(x, y) result(aeq)
    real :: x(:), y(:), maxv
    maxv = max(1.0, maxval(abs(x)), maxval(abs(y)))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqv
    
  logical function eqvs(x, y) result(aeq)
    real :: x(:), y, maxv
    maxv = max(1.0, maxval(abs(x)), abs(y))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqvs
    
  logical function eqsv(x, y) result(aeq)
    real :: x, y(:), maxv
    maxv = max(1.0, abs(x), maxval(abs(y)))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqsv
    
  logical function eqm(x, y) result(aeq)
    real :: x(:,:), y(:,:), maxv
    maxv = max(1.0, maxval(abs(x)), maxval(abs(y)))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqm

  logical function eqsm(x, y) result(aeq)
    real :: x, y(:,:), maxv
    maxv = max(1.0, abs(x), maxval(abs(y)))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqsm

  logical function eqms(x, y) result(aeq)
    real :: x(:,:), y, maxv
    maxv = max(1.0, maxval(abs(x)), abs(y))
    aeq = all(abs(x - y) < maxv*rmd_stolerance)
  end function eqms
  ! -----------------------------------------

  INTEGER FUNCTION RMD_DIMBAND(m, n, kl, ku) result(db)
    integer, intent(in) :: m, n, kl, ku
    integer :: j
    db = 0
    do j = 1,n
      db = db + max(0, min(m, j + kl) - max(1, j - ku) + 1)
    enddo
  END FUNCTION RMD_DIMBAND

  SUBROUTINE RMD_STESTF(F, F_rmd, u, dimu, dimv, tol, testcase, F_dp, print, &
    &                   accum, step, testcols, fourpoint, mask)
    implicit none
    procedure(F_interf)        :: F
    procedure(F_rmd_interf)    :: F_rmd
    real, intent(in)           :: tol
    integer, intent(in)        :: dimu, dimv
    real, intent(in)           :: u(dimu)
    type(fdata), intent(inout) :: testcase
    procedure(F_dp_interf), optional :: F_dp
    logical, intent(in), optional :: print        ! = .true. to print Jacobians
    integer, intent(in), optional :: accum(:)     ! indices to check updating for
    real,    intent(in), optional :: step(dimu)   ! discretization step
    integer, intent(in), optional :: testcols(:)  ! Jacobian columns to check
    logical, intent(in), optional :: fourpoint(:) ! use four-point stencil
    logical, intent(in), optional :: mask(dimu)   ! mask to select ua elements

    ! PURPOSE
    !    Subroutine for testing a function and its corresponding rmd operation
    !    against numerical derivatives
    !
    ! ARGUMENTS
    !    F
    !       (input, subroutine)
    !       A wrapper routine for the BLAS routine that corresponds to the BLAS-RMD
    !       routine being tested. F has parameters:
    !          testcase: (input, type(fdata))
    !                    Structure with scalar inputs to the BLAS routine, see below.
    !          u:        (input, real vector)
    !                    All elements of the vector and matrix inputs to the
    !                    BLAS function
    !          v:        (output, real vector)
    !                    All elements of the outputs from the BLAS routine
    !                    (including u-values that it changes)
    !
    !    F_rmd (input, subroutine)
    !       A wrapper routine for the BLAS-RMD routine to be tested, with parameters:
    !          testcase: (input, type(fdata))
    !                    Structure with the same information as on the call to F.
    !          u, v:     (input, real vectors)
    !                    These parameters will will have the same values as they had
    !                    on entry and exit from the corresponding call to F.
    !          va:       (input, real vector with the same dimension as v)
    !                    The adjoint of v
    !          ua:       (input/output, real vector with the same dimension as u)
    !                       On entry: The adjoint already accumulated for u.
    !                                 (always 0 except when accum is used)
    !                       On exit:  The updated adjoint as returned by the
    !                                 wrapped BLAS-RMD function
    !
    !    u
    !       (input, real vector of dimension dimu)
    !       Vector with values for which the test will be carried out. All the
    !       input vector and matrix elements for the BLAS routine concatenated
    !       into one vector. Contains neither the zero elements of triangular and
    !       band matrices nor repeated elements of symmetric matrices
    !   
    !    dimu
    !       (input, integer)
    !       Dimension of u; total number of input elements of the BLAS routine
    !   
    !    dimv
    !       (input, integer)
    !       Total number vector and matrix elements computed by the BLAS routine,
    !       excluding zero elements of triangular and band matrices and repeated
    !       elements of symmetric matrices.
    !   
    !    tol
    !       (input, real)
    !       Tolerance. If the absolute difference or relative difference between
    !       the results of numerical and analytical differentiation is greater
    !       than tol the test will fail.
    !   
    !    testcase
    !       (input, type(fdata))
    !       Structure with two components: integers and characters, containing
    !       scalar inputs to F for which adjoints are not computed (dimensions
    !       and integers, transp, uplo, etc.). This structure is passed on to
    !       both F and F_rmd unchanged.
    !
    !    step
    !       (optional input, real vector of dimension dimu)
    !       If present, use step(i) as the discretization step for approximating
    !       the derivative w.r.t. u(i) instead of the default, (machine
    !       epsilon)**(1/3)
    !
    !    print
    !       (optional input, logical)
    !       If present and true, the analytic and numerical Jacobians are always
    !       printed, instead of only when they don't match.
    !
    !    accum
    !       (optional input, integer vector of assumed shape, i.e. dimension(:))
    !       If present, indices in u for which accumulation is tested. If k is
    !       in accum then rmd_stestf checks that:
    !            ua_out(k) = ua_in(k) + ua_0(k)
    !       where:
    !            ua_out is the value of ua on exit from F_rmd
    !            ua_in is the value of ua on entry to F_rmd
    !            ua_0 is the value of ua computed by F_rmd when ua_in = 0
    !
    !    testcols
    !       (optional input, integer vector of assumed shape)
    !       If present, indices of Jacobian columns to check, i.e. outputs to
    !       check adjoints for, The default is to check all columns.
    !
    !    fourpoint
    !       (optiomal input, logical vector of dimension dimu)
    !       If present, use a four-point finite difference sencil for the ua
    !       corresponding to true entries. The default is to use a two-point
    !       stencil.
    !
    !    mask
    !       (optional input, logical vector of dimension dimu)
    !       If present, the entries in ua corresponding to false elements in
    !       mask are not tested against numerically computed adjoints.
    
    ! Local variables
    !   Assisting variables
    
    integer :: i, j, l, jdim
    integer, allocatable :: jlist(:)
    real, parameter :: delta = rmd_sdelta ! For numerical derivative approximation
    real maxdiff, del
    real, allocatable, dimension(:) :: v, Fm, Fp, Fmm, Fpp, eu, ev
    double precision, allocatable, dimension(:) :: Fm_dp, Fp_dp
    !   Important variables
    real, allocatable :: J_rmd(:,:), J_num(:,:)      ! Jacobians
    real, allocatable :: J_diff(:,:), J_reldiff(:,:) ! Jacobian differences
    real, allocatable :: ua(:), va(:), ua_accum(:)   ! Adjoints
    logical :: err, printon
    !
    if (present(testcols)) then
      allocate(jlist(size(testcols)))
      jlist = testcols(1:size(testcols))
    else
      allocate(jlist(dimv)) ! needed because of bug in gfortran v7 & v8
      jlist = [(l, l=1,dimv)]
    endif
    if (present(print))          printon = print
    if (.not. present(print))    printon = .false.

    jdim = size(jlist)
    
    ! Interface for the subroutine arguments
    allocate(v(dimv), Fm(dimv), Fp(dimv), eu(dimu), ev(dimv))
    allocate(J_rmd(dimu,jdim), J_num(dimu,dimv), ua(dimu), va(dimv))
    allocate(Fmm(dimv), Fpp(dimv))
    if (present(F_dp)) allocate(Fm_dp(dimv), Fp_dp(dimv))
    !
    call F(testcase,u,v) ! v = F(u)
    !
    eu = 0.0
    ev = 0.0
    !
    ! Approximate Jacobian of F using finite differences
    do i = 1,dimu
      del = delta
      if (present(step)) del = step(i)
      eu(i) = 1.0 ! unit vector
      call F(testcase, u - del*eu, Fm)
      call F(testcase, u + del*eu, Fp)
      if (present(F_dp)) then
        call F_dp(testcase, u - dble(del)*eu, Fm_dp)
        call F_dp(testcase, u + dble(del)*eu, Fp_dp)
        J_num(i,:) = sngl(Fp_dp - Fm_dp)/(2*del)
      elseif (present(fourpoint)) then
        if (fourpoint(i)) then
          call F(testcase, u - 2*del*eu, Fmm)
          call F(testcase, u + 2*del*eu, Fpp)
          J_num(i,:) = (Fmm - Fpp + 8*(Fp - Fm))/(12*del)
        else
          J_num(i,:) = (Fp - Fm)/(2*del)
        endif
      else
        J_num(i,:) = (Fp - Fm)/(2*del)
      endif
      eu(i) = 0.0
    end do
    
    ! We now have J_num = Approximation of the Jacobian for F.
    !
    do l = 1,size(jlist)
      j = jlist(l)
      ev(j) = 1.0
      ua = 0.0
      va = ev
      call F_rmd(testcase, u, v, ua, va)
      J_rmd(:,l) = ua
      ev(j) = 0.0
    end do
    ! We now have J_rmd = The Jacobian for F evaluated using RMD
    !
    err = .FALSE.
    !
    ! Check if the Jacobians contain NaN
    if(.not.all(J_num == J_num)) then ! isnan in Fortran 2003
      err = .TRUE.
      print *, 'Error: Numerical Jacobian contains NaN'
    end if
    if(.not.all(J_rmd == J_rmd)) then 
      err = .TRUE.
      print *, 'Error: Analytic Jacobian contains NaN'
    end if
    !
    ! Calculate maximum relative difference between J_rmd and J_num
    J_diff = J_rmd - J_num(:,jlist)
    allocate(J_reldiff(dimu,jdim))
    do i = 1,dimu
      J_reldiff(i,:) = 0
      if (present(mask)) then
        if (.not.mask(i)) cycle
      endif
      J_reldiff(i,:) = abs(J_diff(i,:)) / &
        &              max(1e-2, abs(J_rmd(i,:)), abs(J_num(i,jlist)))
    enddo
    maxdiff = maxval(J_reldiff)
    !print *, 'maximum relative difference = ', maxdiff
    !
    if (maxdiff > tol) then
      err = .TRUE.
      print *, 'Error: Jacobians do not match.'
      print *, 'maximum relative difference = ', maxdiff
    endif
    !
    ! Show more info if any error has been found
    if (err .or. printon) then
      if (err) then
        print *, 'testcase % integers = ', testcase % integers
        print *, 'testcase % characters = ', testcase % characters
        print *, "See test_X_rmd.f90 where X is the name of the function being"
        print *, "tested (eg sgbmv) to see what testcase%integers and"
        print *, "testcase%characters contain."
      endif
      call rmd_sprintj('Numerical Jacobian = ', J_num(:,jlist))
      call rmd_sprintj('Analytic Jacobian = ', J_rmd)
      call rmd_sprintj('Relative difference = ', J_reldiff)
      print *,'Max difference = ', maxdiff
      if (err) stop 1
    end if
    
    if (present(accum)) then ! Check that F_rmd accumulates correctly (currently
      ! used by rot/rotm and sdsdot)
      ua = 0.0
      call random_number(va)
      call F_rmd(testcase, u, v, ua, va)
      allocate(ua_accum(dimu))
      ua_accum = 0.0
      ua_accum(accum) = [(i+1, i=1, size(accum))]
      call F_rmd(testcase, u, v, ua_accum, va)
      ua(accum) = ua(accum) + [(i+1, i=1, size(accum))]
      if (any(abs(ua - ua_accum) > tol)) stop 'accum-check failed'
    endif
    
  END SUBROUTINE RMD_STESTF

  SUBROUTINE RMD_SPRINTJ(s, J)
    implicit none
    real, intent(in) :: J(:,:)
    character(*), intent(in) :: s
    integer :: m, n, i
    m = size(J, 1)
    n = size(J, 2)
    if (n > 99) then ! max 99 pretty-print columns supported
      print *, s, J
    else
      print '(A)', s
      write(*, '(A)', advance = 'no') 'u      '
      do i = 1, n
        write(*, '("v", i0, t12)', advance = 'no') i
      enddo
      write(*, '()', advance='yes')
      do i = 1, m
        write(*, '(i0,t4)', advance = 'no') i
        write(*, '(99(1x,1pg10.3))') J(i,:)
      enddo
    endif
  END SUBROUTINE RMD_SPRINTJ

  SUBROUTINE RMD_STESTRANDOM(F, F_rmd, dimu, dimv, tol, rep, testcase)
    ! Run rmd_stestf rep times with uniform(-2,2) random values in u
    implicit none
    procedure(F_interf)     :: F
    procedure(F_rmd_interf) :: F_rmd
    real, intent(in) :: tol
    integer, intent(in) :: dimu, dimv, rep
    type(fdata), intent(inout) :: testcase ! Defined in rmd_stesttools
    !
    ! Local variables
    real, allocatable :: u(:)
    integer :: i
    !
    allocate(u(dimu))
    do i = 1,rep
      call random_number(u)
      u = 4.0*u - 2.0
      call rmd_stestf(F,F_rmd,u,dimu,dimv,tol,testcase)
    end do
  END SUBROUTINE RMD_STESTRANDOM

  SUBROUTINE RMD_STEST123(F, F_rmd, dimu, dimv, tol, testcase)
    ! Run rmd_stestf with u = [1, 2, 3...]
    implicit none
    procedure(F_interf)     :: F
    procedure(F_rmd_interf) :: F_rmd
    real, intent(in) :: tol
    integer, intent(in) :: dimu, dimv
    type(fdata), intent(inout) :: testcase ! Defined in rmd_stesttools
    !
    ! Local variables
    real, allocatable :: u(:)
    integer :: i
    !
    allocate(u(dimu))
    do i = 1,dimu
      u(i) = 1 + i
    enddo
    call rmd_stestf(F,F_rmd,u,dimu,dimv,tol,testcase)
  END SUBROUTINE RMD_STEST123

  function rmd_s_ivech(uplo, x) result(A)
    ! currently not used
    character, intent(in) :: uplo
    real, intent(in) :: x(:)
    real, allocatable :: A(:,:)
    integer :: nx, n, k, i
    nx = ubound(x, 1)
    n = int(sqrt(nx*2.0))
    allocate(A(n,n))
    k = 1
    do i=1,n
      if (uplo=='l' .or. uplo=='L') then
        A(i:n,i) = x(k:k+n-i)
      else
        A(i,i:n) = x(k:k+n-i)
      endif
      k = k+n-i+1
    enddo
  end function rmd_s_ivech

  function rmd_s_vech(uplo, A) result(x)
    ! currently not used
    character, intent(in) :: uplo
    real, intent(in) :: A(:,:)
    real, allocatable :: x(:)
    integer :: nx, n, k, i
    n = ubound(A, 1)
    nx = n*(n+1)/2
    allocate(x(nx))
    k = 1
    do i=1,n
      if (uplo=='l' .or. uplo=='L') then
        x(k:k+n-i) = A(i:n,i)
      else
        x(k:k+n-i) = A(i,i:n)
      endif
      k = k+n-i+1
    enddo
  end function rmd_s_vech

  SUBROUTINE RMD_STRISHAPE_V2M(uplo, m, u, A, lda)
    ! Reshape a vector into a triangular matrix (vech^(-1))
    ! uplo: (character(*)) Begins with 'u' or 'l' for upper or lower triangular matrix
    ! m:    (integer) Matrix size
    ! u:    (real vector of dimension m*(m+1)/2) Input vector
    ! A:    (real lda by * matrix) Upper left m by m part receives, result
    ! lda:  (integer) Leading dimension of A
    implicit none
    character(*), intent(in) :: uplo
    integer, intent(in)      :: m, lda
    real, intent(in)         :: u(m*(m+1)/2)
    real, intent(out)        :: A(lda,m)
    integer :: i, uax
    uax = 1
    ! Transform a vector into a triangular matrix.
    if(uplo == 'u' .or. uplo == 'U') then
      ! Matrix should be upper triangular.
      do i = 1,m
        call scopy(i, u(uax),1, A(1,i),1)
        uax = uax + i
      end do
    else
      ! Matrix should be lower triangular.
      do i = 1,m
        call scopy(m-i+1, u(uax),1, A(i,i),1)
        uax = uax + m-i+1
      end do
    end if
  END SUBROUTINE RMD_STRISHAPE_V2M
      
  SUBROUTINE RMD_STRISHAPE_M2V(uplo, m, u, A, lda)
    ! Reshape a triangular matrix into a vector (vech)
    ! uplo: (character(*)) Begins with 'u' or 'l' for upper or lower triangular matrix
    ! m:    (integer) Matrix size
    ! u:    (real vector of dimension m*(m+1)/2) Receives result
    ! A:    (real lda by * matrix) Upper left m by m part is input Matrix
    ! lda:  (integer) Leading dimension of A
    implicit none
    character(*), intent(in) :: uplo
    integer, intent(in)      :: m, lda
    real, intent(out)        :: u(m*(m+1)/2)
    real, intent(in)         :: A(lda,m)
    integer :: i, uax
    uax = 1
    if(uplo == 'u' .or. uplo == 'U') then
      ! Matrix is upper triangular. 
      do i = 1,m
        call scopy(i, A(1,i),1, u(uax),1)
        uax = uax + i
      end do
    else
      ! Matrix is lower triangular.
      do i = 1,m
        call scopy(m-i+1, A(i,i),1, u(uax),1)
        uax = uax + m-i+1
      end do
    end if
  END SUBROUTINE RMD_STRISHAPE_M2V

  SUBROUTINE RMD_SCOPY_UNIT_V2M(uplo, n, u, A, lda)
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda
    real, intent(in) :: u(n*(n-1)/2)
    real, intent(out) :: A(lda,n)
    integer k1, i
    k1 = 1
    do i = 1,size(A,2) - 1
      if (uplo == 'U' .or. uplo == 'u') then
        A(1:i, i+1) = u(k1:k1+i-1)
      else
        A(i+1, 1:i) = u(k1:k1+i-1)
      endif
      k1 = k1 + i
    enddo
  END SUBROUTINE RMD_SCOPY_UNIT_V2M

  SUBROUTINE RMD_SCOPY_UNIT_M2V(uplo, n, u, A, lda)
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda
    real, intent(out) :: u(n*(n-1)/2)
    real, intent(in) :: A(lda,n)
    integer k1, i
    k1 = 1
    do i = 1,size(A,2)-1
      if (uplo == 'U' .or. uplo == 'u') then
        u(k1:k1+i-1) = A(1:i, i+1)
      else
        u(k1:k1+i-1) = A(i+1, 1:i)
      endif
      k1 = k1 + i
    enddo
  END SUBROUTINE RMD_SCOPY_UNIT_M2V

  SUBROUTINE RMD_SBANDSHAPE_V2M (m, n, kl, ku, u, A, lda)
    ! Reshape a vector into a banded matrix
    ! n   : Column count of the banded matrix
    ! kl  : Number of subdiagonals in the banded matrix
    ! ku  : Number of superdiagonals in the banded matrix
    ! u   : Vector to be reshaped
    ! A   : Resulting matrix. Should be at least (kl+ku+1) by n elements
    !       A stores an m by n banded matrix
    ! lda : Leading dimension of A
    implicit none
    integer, intent(in) :: m, n, kl, ku, lda
    real, intent(in)    :: u(*)
    real, intent(out)   :: A(lda,*)
    integer :: j, k, p, i1, i2, nrows
    p = 0
    do j = 1,n
      k = ku + 1 - j
      i1 = max(1, j-ku)
      i2 = min(m, j+kl)
      nrows = i2 - i1 + 1
      a(k+i1:k+i2, j) = u(p+1:p+nrows)
      p = p + nrows
    enddo
  END SUBROUTINE RMD_SBANDSHAPE_V2M

  SUBROUTINE RMD_SBANDSHAPE_M2V (m, n, kl, ku, u, A, lda)
    ! Reshape a banded matrix into a vector
    ! n   : Column count of the banded matrix
    ! kl  : Number of subdiagonals in the banded matrix
    ! ku  : Number of superdiagonals in the banded matrix
    ! u   : Resulting vector
    ! A   : Matrix to be reshaped. Should be at least (kl+ku+1) by n elements
    !       A stores an m by n banded matrix
    ! lda : Leading dimension of A
    implicit none
    integer, intent(in) :: m, n, kl, ku, lda
    real, intent(out) :: u(*)
    real, intent(in)    :: A(lda,*)
    integer :: j, k, p, i1, i2, nrows
    p = 0
    do j = 1,n
      k = ku + 1 - j
      i1 = max(1, j-ku)
      i2 = min(m, j+kl)
      nrows = i2 - i1 + 1
      u(p+1:p+nrows) = a(k+i1:k+i2, j)
      p = p + nrows
    enddo
  END SUBROUTINE RMD_SBANDSHAPE_M2V

  FUNCTION RMD_SRANDOM() RESULT(x)
    real :: x
    call random_number(x)
  END FUNCTION RMD_SRANDOM

END MODULE RMD_STESTTOOLS

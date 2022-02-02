! This module computes the matrix values of phi-functions for an overlapped
! block-diagonal matrix that consists of upper quasi-triangular blocks.
! The modified scaling and squaring method is used for computation.
!
! References:
!
!   Nicholas J. Higham,
!   The scaling and squaring method for the matrix exponential revisited,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 26, Number 4, pp. 1179-1193, 2005.
!
!   S. Koikari,
!   An error analysis of the modified scaling and squaring method,
!   Computers and Mathematics with Applications,
!   Volume 53, pp. 1293-1305, 2007.
!
! This module is intended for internal use only.

module ovlpconfgeom

  use floattypes
  use thetamn
  use matrixpwrtag
  use ovlpmatrix
  use ovlpratiofunc
  use mtrcfgphilog

  implicit none

  ! the highest order of rational approximants used in each precision

  integer, parameter :: sp_hord =  7                ! for single    precision
  integer, parameter :: wp_hord = 13                ! for double    precision
  integer, parameter :: tp_hord = 13                ! for triple    precision
  integer, parameter :: qp_hord = 17                ! for quadruple precision

  public

contains


!    Name : rs_ovlp_scsqrf
!         :
!   Usage : values = rs_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a real matrix in single precision.

pure function rs_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: scale
  real   (kind=sp), intent(in) :: arg(:,:)
  real   (kind=sp) :: rs_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rs_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,rs_ovlp_scsqrf,err)

end function rs_ovlp_scsqrf


!    Name : rs_ovlp_scsqr
!         :
!   Usage : call rs_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real matrix in single precision.

pure subroutine rs_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  real   (kind=sp), intent(in ) :: scale
  real   (kind=sp), intent(in ) :: arg(:,:)
  real   (kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_sp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call rs_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rstr_norm1(a), rstr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rs_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    s = real(2**spower,sp)                             ! scaling
    a = a / s
    err % spower = spower

    call rs_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = rs_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call rs_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine rs_ovlp_scsqr


!    Name : rs_ovlp_squaref
!         :
!   Usage : phiv = rs_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a real matrix in single precision.

pure function rs_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: arg(:,:,0:)
  real   (kind=sp) :: rs_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  rs_ovlp_squaref = arg
  call rs_ovlp_square(diag,subd,head,tail,nonz,upto,rs_ovlp_squaref)
end function rs_ovlp_squaref


!    Name : rs_ovlp_square
!         :
!   Usage : call rs_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a real matrix in single precision.

pure subroutine rs_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(inout) :: arg(:,:,0:)

  real   (kind=sp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = 0.0_sp
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + 1.0_sp
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_sp &
                   + arg(i,h:t,4) /  32.0_sp &
                   + arg(i,h:t,3) /  64.0_sp &
                   + arg(i,h:t,2) / 192.0_sp &
                   + arg(i,h:t,1) / 768.0_sp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_sp &
                   + arg(i,h:t,3) /  16.0_sp &
                   + arg(i,h:t,2) /  32.0_sp &
                   + arg(i,h:t,1) /  96.0_sp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_sp &
                   + arg(i,h:t,2) /   8.0_sp &
                   + arg(i,h:t,1) /  16.0_sp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_sp &
                   + arg(i,h:t,1) /   4.0_sp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_sp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=sp), intent(in) :: a(:,:), b(:,:)
    real   (kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rs_ovlp_square


!    Name : cs_ovlp_scsqrf
!         :
!   Usage : values = cs_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a complex matrix in single precision.

pure function cs_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: scale
  complex(kind=sp), intent(in) :: arg(:,:)
  complex(kind=sp) :: cs_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cs_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,cs_ovlp_scsqrf,err)

end function cs_ovlp_scsqrf


!    Name : cs_ovlp_scsqr
!         :
!   Usage : call cs_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex matrix in single precision.

pure subroutine cs_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  complex(kind=sp), intent(in ) :: scale
  complex(kind=sp), intent(in ) :: arg(:,:)
  complex(kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_sp,0.0_sp,sp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call cs_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cstr_norm1(a), cstr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cs_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    s = real(2**spower,sp)                             ! scaling
    a = a / s
    err % spower = spower

    call cs_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = cs_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call cs_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine cs_ovlp_scsqr


!    Name : cs_ovlp_squaref
!         :
!   Usage : phiv = cs_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a complex matrix in single precision.

pure function cs_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: arg(:,:,0:)
  complex(kind=sp) :: cs_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  cs_ovlp_squaref = arg
  call cs_ovlp_square(diag,subd,head,tail,nonz,upto,cs_ovlp_squaref)
end function cs_ovlp_squaref


!    Name : cs_ovlp_square
!         :
!   Usage : call cs_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a complex matrix in single precision.

pure subroutine cs_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(inout) :: arg(:,:,0:)

  complex(kind=sp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + cmplx(1.0_sp,0.0_sp,sp)
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_sp &
                   + arg(i,h:t,4) /  32.0_sp &
                   + arg(i,h:t,3) /  64.0_sp &
                   + arg(i,h:t,2) / 192.0_sp &
                   + arg(i,h:t,1) / 768.0_sp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_sp &
                   + arg(i,h:t,3) /  16.0_sp &
                   + arg(i,h:t,2) /  32.0_sp &
                   + arg(i,h:t,1) /  96.0_sp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_sp &
                   + arg(i,h:t,2) /   8.0_sp &
                   + arg(i,h:t,1) /  16.0_sp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_sp &
                   + arg(i,h:t,1) /   4.0_sp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_sp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=sp), intent(in) :: a(:,:), b(:,:)
    complex(kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cs_ovlp_square


!    Name : rw_ovlp_scsqrf
!         :
!   Usage : values = rw_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a real matrix in double precision.

pure function rw_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: scale
  real   (kind=wp), intent(in) :: arg(:,:)
  real   (kind=wp) :: rw_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rw_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,rw_ovlp_scsqrf,err)

end function rw_ovlp_scsqrf


!    Name : rw_ovlp_scsqr
!         :
!   Usage : call rw_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real matrix in double precision.

pure subroutine rw_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  real   (kind=wp), intent(in ) :: scale
  real   (kind=wp), intent(in ) :: arg(:,:)
  real   (kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_wp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call rw_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rwtr_norm1(a), rwtr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rw_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    s = real(2**spower,wp)                             ! scaling
    a = a / s
    err % spower = spower

    call rw_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = rw_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call rw_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine rw_ovlp_scsqr


!    Name : rw_ovlp_squaref
!         :
!   Usage : phiv = rw_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a real matrix in double precision.

pure function rw_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: arg(:,:,0:)
  real   (kind=wp) :: rw_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  rw_ovlp_squaref = arg
  call rw_ovlp_square(diag,subd,head,tail,nonz,upto,rw_ovlp_squaref)
end function rw_ovlp_squaref


!    Name : rw_ovlp_square
!         :
!   Usage : call rw_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a real matrix in double precision.

pure subroutine rw_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(inout) :: arg(:,:,0:)

  real   (kind=wp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = 0.0_wp
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + 1.0_wp
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_wp &
                   + arg(i,h:t,4) /  32.0_wp &
                   + arg(i,h:t,3) /  64.0_wp &
                   + arg(i,h:t,2) / 192.0_wp &
                   + arg(i,h:t,1) / 768.0_wp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_wp &
                   + arg(i,h:t,3) /  16.0_wp &
                   + arg(i,h:t,2) /  32.0_wp &
                   + arg(i,h:t,1) /  96.0_wp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_wp &
                   + arg(i,h:t,2) /   8.0_wp &
                   + arg(i,h:t,1) /  16.0_wp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_wp &
                   + arg(i,h:t,1) /   4.0_wp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_wp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=wp), intent(in) :: a(:,:), b(:,:)
    real   (kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rw_ovlp_square


!    Name : cw_ovlp_scsqrf
!         :
!   Usage : values = cw_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a complex matrix in double precision.

pure function cw_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: scale
  complex(kind=wp), intent(in) :: arg(:,:)
  complex(kind=wp) :: cw_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cw_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,cw_ovlp_scsqrf,err)

end function cw_ovlp_scsqrf


!    Name : cw_ovlp_scsqr
!         :
!   Usage : call cw_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex matrix in double precision.

pure subroutine cw_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  complex(kind=wp), intent(in ) :: scale
  complex(kind=wp), intent(in ) :: arg(:,:)
  complex(kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_wp,0.0_wp,wp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call cw_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cwtr_norm1(a), cwtr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cw_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    s = real(2**spower,wp)                             ! scaling
    a = a / s
    err % spower = spower

    call cw_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = cw_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call cw_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine cw_ovlp_scsqr


!    Name : cw_ovlp_squaref
!         :
!   Usage : phiv = cw_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a complex matrix in double precision.

pure function cw_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: arg(:,:,0:)
  complex(kind=wp) :: cw_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  cw_ovlp_squaref = arg
  call cw_ovlp_square(diag,subd,head,tail,nonz,upto,cw_ovlp_squaref)
end function cw_ovlp_squaref


!    Name : cw_ovlp_square
!         :
!   Usage : call cw_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a complex matrix in double precision.

pure subroutine cw_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(inout) :: arg(:,:,0:)

  complex(kind=wp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + cmplx(1.0_wp,0.0_wp,wp)
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_wp &
                   + arg(i,h:t,4) /  32.0_wp &
                   + arg(i,h:t,3) /  64.0_wp &
                   + arg(i,h:t,2) / 192.0_wp &
                   + arg(i,h:t,1) / 768.0_wp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_wp &
                   + arg(i,h:t,3) /  16.0_wp &
                   + arg(i,h:t,2) /  32.0_wp &
                   + arg(i,h:t,1) /  96.0_wp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_wp &
                   + arg(i,h:t,2) /   8.0_wp &
                   + arg(i,h:t,1) /  16.0_wp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_wp &
                   + arg(i,h:t,1) /   4.0_wp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_wp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=wp), intent(in) :: a(:,:), b(:,:)
    complex(kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cw_ovlp_square

#ifdef __USE_TPREC

!    Name : rt_ovlp_scsqrf
!         :
!   Usage : values = rt_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a real matrix in triple precision.

pure function rt_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: scale
  real   (kind=tp), intent(in) :: arg(:,:)
  real   (kind=tp) :: rt_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rt_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,rt_ovlp_scsqrf,err)

end function rt_ovlp_scsqrf


!    Name : rt_ovlp_scsqr
!         :
!   Usage : call rt_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real matrix in triple precision.

pure subroutine rt_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  real   (kind=tp), intent(in ) :: scale
  real   (kind=tp), intent(in ) :: arg(:,:)
  real   (kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_tp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call rt_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rttr_norm1(a), rttr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rt_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    s = real(2**spower,tp)                             ! scaling
    a = a / s
    err % spower = spower

    call rt_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = rt_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call rt_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine rt_ovlp_scsqr


!    Name : rt_ovlp_squaref
!         :
!   Usage : phiv = rt_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a real matrix in triple precision.

pure function rt_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: arg(:,:,0:)
  real   (kind=tp) :: rt_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  rt_ovlp_squaref = arg
  call rt_ovlp_square(diag,subd,head,tail,nonz,upto,rt_ovlp_squaref)
end function rt_ovlp_squaref


!    Name : rt_ovlp_square
!         :
!   Usage : call rt_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a real matrix in triple precision.

pure subroutine rt_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(inout) :: arg(:,:,0:)

  real   (kind=tp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = 0.0_tp
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + 1.0_tp
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_tp &
                   + arg(i,h:t,4) /  32.0_tp &
                   + arg(i,h:t,3) /  64.0_tp &
                   + arg(i,h:t,2) / 192.0_tp &
                   + arg(i,h:t,1) / 768.0_tp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_tp &
                   + arg(i,h:t,3) /  16.0_tp &
                   + arg(i,h:t,2) /  32.0_tp &
                   + arg(i,h:t,1) /  96.0_tp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_tp &
                   + arg(i,h:t,2) /   8.0_tp &
                   + arg(i,h:t,1) /  16.0_tp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_tp &
                   + arg(i,h:t,1) /   4.0_tp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_tp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=tp), intent(in) :: a(:,:), b(:,:)
    real   (kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rt_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rt_ovlp_square


!    Name : ct_ovlp_scsqrf
!         :
!   Usage : values = ct_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a complex matrix in triple precision.

pure function ct_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: scale
  complex(kind=tp), intent(in) :: arg(:,:)
  complex(kind=tp) :: ct_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call ct_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,ct_ovlp_scsqrf,err)

end function ct_ovlp_scsqrf


!    Name : ct_ovlp_scsqr
!         :
!   Usage : call ct_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex matrix in triple precision.

pure subroutine ct_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  complex(kind=tp), intent(in ) :: scale
  complex(kind=tp), intent(in ) :: arg(:,:)
  complex(kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_tp,0.0_tp,tp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call ct_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cttr_norm1(a), cttr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call ct_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    s = real(2**spower,tp)                             ! scaling
    a = a / s
    err % spower = spower

    call ct_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = ct_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call ct_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine ct_ovlp_scsqr


!    Name : ct_ovlp_squaref
!         :
!   Usage : phiv = ct_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a complex matrix in triple precision.

pure function ct_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: arg(:,:,0:)
  complex(kind=tp) :: ct_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  ct_ovlp_squaref = arg
  call ct_ovlp_square(diag,subd,head,tail,nonz,upto,ct_ovlp_squaref)
end function ct_ovlp_squaref


!    Name : ct_ovlp_square
!         :
!   Usage : call ct_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a complex matrix in triple precision.

pure subroutine ct_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(inout) :: arg(:,:,0:)

  complex(kind=tp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + cmplx(1.0_tp,0.0_tp,tp)
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_tp &
                   + arg(i,h:t,4) /  32.0_tp &
                   + arg(i,h:t,3) /  64.0_tp &
                   + arg(i,h:t,2) / 192.0_tp &
                   + arg(i,h:t,1) / 768.0_tp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_tp &
                   + arg(i,h:t,3) /  16.0_tp &
                   + arg(i,h:t,2) /  32.0_tp &
                   + arg(i,h:t,1) /  96.0_tp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_tp &
                   + arg(i,h:t,2) /   8.0_tp &
                   + arg(i,h:t,1) /  16.0_tp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_tp &
                   + arg(i,h:t,1) /   4.0_tp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_tp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=tp), intent(in) :: a(:,:), b(:,:)
    complex(kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = ct_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine ct_ovlp_square

#endif
#ifdef __USE_QPREC

!    Name : rq_ovlp_scsqrf
!         :
!   Usage : values = rq_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a real matrix in quadruple precision.

pure function rq_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: scale
  real   (kind=qp), intent(in) :: arg(:,:)
  real   (kind=qp) :: rq_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rq_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,rq_ovlp_scsqrf,err)

end function rq_ovlp_scsqrf


!    Name : rq_ovlp_scsqr
!         :
!   Usage : call rq_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real matrix in quadruple precision.

pure subroutine rq_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  real   (kind=qp), intent(in ) :: scale
  real   (kind=qp), intent(in ) :: arg(:,:)
  real   (kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_qp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call rq_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rqtr_norm1(a), rqtr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rq_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    s = real(2**spower,qp)                             ! scaling
    a = a / s
    err % spower = spower

    call rq_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = rq_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call rq_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine rq_ovlp_scsqr


!    Name : rq_ovlp_squaref
!         :
!   Usage : phiv = rq_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a real matrix in quadruple precision.

pure function rq_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: arg(:,:,0:)
  real   (kind=qp) :: rq_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  rq_ovlp_squaref = arg
  call rq_ovlp_square(diag,subd,head,tail,nonz,upto,rq_ovlp_squaref)
end function rq_ovlp_squaref


!    Name : rq_ovlp_square
!         :
!   Usage : call rq_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a real matrix in quadruple precision.

pure subroutine rq_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(inout) :: arg(:,:,0:)

  real   (kind=qp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = 0.0_qp
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + 1.0_qp
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_qp &
                   + arg(i,h:t,4) /  32.0_qp &
                   + arg(i,h:t,3) /  64.0_qp &
                   + arg(i,h:t,2) / 192.0_qp &
                   + arg(i,h:t,1) / 768.0_qp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_qp &
                   + arg(i,h:t,3) /  16.0_qp &
                   + arg(i,h:t,2) /  32.0_qp &
                   + arg(i,h:t,1) /  96.0_qp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_qp &
                   + arg(i,h:t,2) /   8.0_qp &
                   + arg(i,h:t,1) /  16.0_qp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_qp &
                   + arg(i,h:t,1) /   4.0_qp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_qp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=qp), intent(in) :: a(:,:), b(:,:)
    real   (kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rq_ovlp_square


!    Name : cq_ovlp_scsqrf
!         :
!   Usage : values = cq_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
!         :
! Purpose : This function computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         : for the overlapped block-diagonal region.
!    Note : This function treats a complex matrix in quadruple precision.

pure function cq_ovlp_scsqrf(diag,subd,head,tail,nonz,upto,scale,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: scale
  complex(kind=qp), intent(in) :: arg(:,:)
  complex(kind=qp) :: cq_ovlp_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cq_ovlp_scsqr(diag,subd,head,tail,nonz,upto,scale,arg,cq_ovlp_scsqrf,err)

end function cq_ovlp_scsqrf


!    Name : cq_ovlp_scsqr
!         :
!   Usage : call cq_ovlp_scsqr(                                    &
!         :   diag,subd,head,tail,nonz,upto,scale,arg,values,err   &
!         : )
!         :
! Purpose : This subroutine computes the matrix values of phi-functions
!         : for an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``arg'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``arg'' is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - ``scale'' is a scalar that should multiply the entries of the
!         :   matrix argument, ``arg''.
!         : - ``arg'' is the matrix argument.
!         :
!  Output : - The subarray, values(*,k), stores the value of phi_k(scale*arg(*))
!         :   for the overlapped block-diagonal region.
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in the
!         :   file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex matrix in quadruple precision.

pure subroutine cq_ovlp_scsqr(                                                &
  diag, subd, head, tail, nonz, upto, scale, arg, values, err                 &
)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in ) :: upto
  complex(kind=qp), intent(in ) :: scale
  complex(kind=qp), intent(in ) :: arg(:,:)
  complex(kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs, s
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_qp,0.0_qp,qp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale
  call cq_ovlp_clear(diag,subd,head,tail,nonz,a)       !! to reduce the norm.

  err % upto = upto

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cqtr_norm1(a), cqtr_normi(a))

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cq_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    s = real(2**spower,qp)                             ! scaling
    a = a / s
    err % spower = spower

    call cq_ovlp_ratio(diag,subd,head,tail,nonz,order,upto,a,values)

    do i=1,spower                                      ! repeated squarings
      values = cq_ovlp_squaref(diag,subd,head,tail,nonz,upto,values)
    end do

  end if

  ! Clear unnecessary off-diagonal entries.

  do i=0,upto
    call cq_ovlp_clear(diag,subd,head,tail,nonz,values(:,:,i))
  end do

end subroutine cq_ovlp_scsqr


!    Name : cq_ovlp_squaref
!         :
!   Usage : phiv = cq_ovlp_squaref(diag,subd,head,tail,nonz,upto,phiu)
!         :
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, phiu(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, phiv(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This function treats a complex matrix in quadruple precision.

pure function cq_ovlp_squaref(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: arg(:,:,0:)
  complex(kind=qp) :: cq_ovlp_squaref(size(arg,1),size(arg,2),0:upto)

  cq_ovlp_squaref = arg
  call cq_ovlp_square(diag,subd,head,tail,nonz,upto,cq_ovlp_squaref)
end function cq_ovlp_squaref


!    Name : cq_ovlp_square
!         :
!   Usage : call cq_ovlp_square(diag,subd,head,tail,nonz,upto,values)
!         :
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - ``upto'' is the maximum index of phi-functions to be computed
!         :   simultaneously (0 <= upto <= 5).
!         : - The subarray, values(*,k), stores the value of phi_k(A) for the
!         :   overlapped block-diagonal region (0 <= k <= upto).
!         :
!  Output : The subarray, values(*,k), stores the value of phi_k(2*A) for the
!         : overlapped block-diagonal region (0 <= k <= upto).
!    Note : This subroutine treats a complex matrix in quadruple precision.

pure subroutine cq_ovlp_square(diag,subd,head,tail,nonz,upto,arg)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(inout) :: arg(:,:,0:)

  complex(kind=qp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: i, h, t, band(size(arg,1))

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  Iexp = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,size(arg,1)
    h = max(1,i-1); t = min(band(i),size(arg,2))
    Iexp(i,h:t) = arg(i,h:t,0)
    Iexp(i,i) = Iexp(i,i) + cmplx(1.0_qp,0.0_qp,qp)
  end do

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = matrixmul(arg(:,:,5), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,5) = arg(i,h:t,5) /  32.0_qp &
                   + arg(i,h:t,4) /  32.0_qp &
                   + arg(i,h:t,3) /  64.0_qp &
                   + arg(i,h:t,2) / 192.0_qp &
                   + arg(i,h:t,1) / 768.0_qp
    end do
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = matrixmul(arg(:,:,4), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,4) = arg(i,h:t,4) /  16.0_qp &
                   + arg(i,h:t,3) /  16.0_qp &
                   + arg(i,h:t,2) /  32.0_qp &
                   + arg(i,h:t,1) /  96.0_qp
    end do
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = matrixmul(arg(:,:,3), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,3) = arg(i,h:t,3) /   8.0_qp &
                   + arg(i,h:t,2) /   8.0_qp &
                   + arg(i,h:t,1) /  16.0_qp
    end do
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = matrixmul(arg(:,:,2), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,2) = arg(i,h:t,2) /   4.0_qp &
                   + arg(i,h:t,1) /   4.0_qp
    end do
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = matrixmul(arg(:,:,1), Iexp)
    do i=1,size(arg,1)
      h = max(1,i-1); t = min(band(i),size(arg,2))
      arg(i,h:t,1) = arg(i,h:t,1) /   2.0_qp
    end do
  end if
                                                   ! square the exponential
  arg(:,:,0) = matrixmul(arg(:,:,0), arg(:,:,0))

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices of the same size
  !         : and of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=qp), intent(in) :: a(:,:), b(:,:)
    complex(kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cq_ovlp_square

#endif

end module ovlpconfgeom


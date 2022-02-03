! This module computes phi-functions required in exponential integrators
! by the modified block Schur--Parlett algorithm.
!
! References:
!
!   Beresford N. Parlett,
!   A recurrence among the elements of functions of triangular matrices,
!   Linear Algebra and its Applications,
!   Volume 14, pp. 117-121, 1976.
!
!   Isak Jonsson and Bo K{\aa}gstr{\"{o}}m,
!   Recursive blocked algorithms for solving triangular systems. I.
!   One-sided and coupled Sylvester-type matrix equations",
!   ACM Transactions on Mathematical Software",
!   Volume 28, Number 4, pp. 392-415, 2002,
!
!   Philip I. Davies and Nicholas J. Higham,
!   A Schur--Parlett algorithm for computing matrix functions,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 25, Number 2, pp. 464-485, 2003.
!
!   _. _______,
!   On a sep-inverse estimate and its application to a modified
!   block Schur--Parlett algorithm,
!   submitted to ACM transactions on mathematical software, 2008.
!

module schurparlett

  use floattypes
  use matrixpwrtag
  use ovlpmatrix
  use scalesquare
  use ovlpconfgeom
  use jspsylvester
  use blockdecomp

  implicit none


  !    Name : bspqtrphif ( a generic name )
  !         :
  !   Usage : values = bspqtrphif(upto,hci,a,info)
  !         :
  ! Purpose : This function computes phi-functions by the modified
  !         : block Schur--Parlett algorithm.
  !         :
  !   Input : - ``upto'' is the maximum index of phi-functions to be computed.
  !         : - ``hci'' is a scalar that should multiply the entries of ``A''.
  !         : - The square array ``A'' stores an upper quasi-triangular matrix.
  !         : - ``info'' is of type(bspblock), and stores the blocking
  !         :   information produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
  !         :   The derived type and the functions are defined in the file
  !         :   ``blockdecomp.f95''.
  !         :
  !  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
  !         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)), and
  !         : values(1:N,1:N,0) is the matrix exponential.
  !
  ! The most simple example of usage is:
  !
  !    use blockdecomp,  only : bspblock, bspqtr2dcmpf
  !    use schurparlett, only : bspqtrphif
  !
  !    integer, parameter :: N = 128, upto = 5
  !    real*8             :: hci, a(N,N), values(N,N,0:upto)
  !    type(bspblock)     :: info
  !
  !    ! After defining ``hci'' and ``a'',
  !
  !    info = bspqtr2dcmpf(a)                  ! blocking independent of ``hci''
  !    values = bspqtrphif(upto,hci,a,info)    ! computing phi-functions
  !
  ! The detailed knowledge of the ``info'' is not necessary to use this module.

  interface bspqtrphif
    module procedure rs_blkparlettf            ! Real    in Single    precision
    module procedure cs_blkparlettf            ! Complex in Single    precision
    module procedure rw_blkparlettf            ! Real    in Double    precision
    module procedure cw_blkparlettf            ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_blkparlettf            ! Real    in Triple    precision
    module procedure ct_blkparlettf            ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_blkparlettf            ! Real    in Quadruple precision
    module procedure cq_blkparlettf            ! Complex in Quadruple precision
#endif
  end interface bspqtrphif


  !    Name : bspqtrphi ( a generic name )
  !         :
  !   Usage : call bspqtrphi(upto,hci,a,info,values,errlog)
  !         :
  ! Purpose : This subroutine computes phi-functions by the modified
  !         : block Schur--Parlett algorithm.
  !         :
  !   Input : - ``upto'' is the maximum index of phi-functions to be computed.
  !         : - ``hci'' is a scalar that should multiply the entries of ``A''.
  !         : - The square array ``A'' stores an upper quasi-triangular matrix.
  !         : - ``info'' is of type(bspblock), and stores the blocking
  !         :   information produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
  !         :   The derived type and two functions are defined in the file
  !         :   ``blockdecomp.f95''.
  !         :
  !  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
  !         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
  !         :   values(1:N,1:N,0) is the matrix exponential.
  !         : - ``errlog'' is of type(bspqtrphilog), and stores error
  !         :   information. The derived type, type(bspqtrphilog), is defined
  !         :   later in this file.

  interface bspqtrphi
    module procedure rs_blkparlett             ! Real    in Single    precision
    module procedure cs_blkparlett             ! Complex in Single    precision
    module procedure rw_blkparlett             ! Real    in Double    precision
    module procedure cw_blkparlett             ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_blkparlett             ! Real    in Triple    precision
    module procedure ct_blkparlett             ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_blkparlett             ! Real    in Quadruple precision
    module procedure cq_blkparlett             ! Complex in Quadruple precision
#endif
  end interface bspqtrphi


  ! A type to return an error code.
  ! 0 : OK, no error,
  ! 1 : The index of phi-function is out of the range.
  ! 2 : The size of the matrix is invalid.
  ! 3 : A given matrix is not square.
  ! 4 : A diagonal block larger than 2x2 exists.
  ! 5 : A given matrix is not upper triangular.

  type bspqtrphilog
    integer :: error
  end type bspqtrphilog


  ! ``blockwisefunc''
  ! If you want to compute phi-functions of an overlapped block-diagonal matrix
  ! without repeating equivalent operations, please set
  !     blockwisefunc = .false.
  ! This choice is not recommended for a matrix with a large spectral radius,
  ! while it is recommended for a non-normal matrix with a smaller spectral
  ! radius.

  logical :: blockwisefunc = .true.

  private
  public  :: bspqtrphif, bspqtrphi
  public  :: bspqtrphilog
  public  :: getdiagonal
  public  :: blockwisefunc

contains


!    Name : getdiagonal
!         :
! Purpose : This subroutine extracts the information on diagonal blocks.
!         :
!   Input : - ``info'' is of type(bspblock), and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a block in that depth.
!         : - ``ptr'' is the index of a diagonal block.
!         :
!  Output : - blocks(:,1:2) stores the information on diagonal blocks.
!         :   blocks(k,1) is the row/column number of the head of k-th diagonal
!         :   block, and blocks(k,2) is the row/column number of the tail of
!         :   k-th diagonal block.
!         : - (ptr-1) is equal to the number of diagonal blocks.
!
! The most simple example of usage is:
!
!     integer :: ptr, blocks(1:maxblock,1:2)
!
!     ptr = 1; blocks = 0
!     call getdiagonal(info,1,1,ptr,blocks)
!     ptr = ptr - 1
!
! Then, ``ptr'' stores the number of diagonal blocks. For all k in [1,ptr],
! blocks(k,1) stores the row/column number of the head of k-th block, and
! blocks(k,2) stores the row/column number of the tail of k-th block.

pure recursive &
subroutine getdiagonal(info,depth,partition,ptr,blocks)
  type(bspblock), intent(in   ) :: info
  integer,        intent(in   ) :: depth, partition
  integer,        intent(inout) :: ptr
  integer,        intent(inout) :: blocks(:,:)

  integer :: k, d, p(2)

  k = 2**(depth + 1) + 4*partition - 7

  if      (1 <= info % heads(k+1)  .and. 1 <= info % heads(k+2) ) then

    d = depth + 1
    p(1) = 2 * partition - 1; p(2) = p(1) + 1

    call getdiagonal(info,d,p(1),ptr,blocks)
    call getdiagonal(info,d,p(2),ptr,blocks)

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    blocks(ptr,1) = info % heads(k  )
    blocks(ptr,2) = info % heads(k+3)
    ptr = ptr + 1

  end if
end subroutine getdiagonal


!    Name : rs_blkparlettf
!         :
!   Usage : values = rs_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a real matrix in single precision.

pure function rs_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  real   (kind=sp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  real   (kind=sp)   :: rs_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call rs_blkparlett(upto,hci,a,info,rs_blkparlettf,errlog)
end function rs_blkparlettf


!    Name : rs_blkparlett
!         :
!   Usage : call rs_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a real matrix in single precision.

pure subroutine rs_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  real   (kind=sp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  real   (kind=sp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  real   (kind=sp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  real   (kind=sp) :: rhs(1:size(a,1),1:size(a,2))
  real   (kind=sp) :: tmp(1:size(a,1),1:size(a,2))
  real   (kind=sp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = 0.0_sp
  arg  = 0.0_sp
  rhs  = 0.0_sp
  tmp  = 0.0_sp
  aph  = 0.0_sp

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_sp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_sp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call rs_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call rs_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call rs_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_sp)) phin(j+1,j,:) = 0.0_sp
  end do

end subroutine rs_blkparlett


!    Name : rs_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a real matrix in single precision.

pure recursive &
subroutine rs_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=sp), intent(inout) :: phin(1:,1:,0:)
  real   (kind=sp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rs_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call rs_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call rs_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rs_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call rs_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call rs_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine rs_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  real   (kind=sp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=sp) :: smp, est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine rs_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine rs_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  real   (kind=sp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=sp) :: est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_sp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine rs_checkgrowth2x2

#endif

end subroutine rs_rec3x3nnm1


!    Name : rs_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``rs_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a real matrix in single precision.

pure recursive &
subroutine rs_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=sp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rs_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call rs_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rs_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call rs_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine rs_rec3x3alln


!    Name : rs_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a real matrix in single precision.

pure subroutine rs_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  real   (kind=sp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  real   (kind=sp), intent(out) :: values(:,:,0:)

  real   (kind=sp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = 1.0_sp
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_sp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call rs_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call rs_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine rs_compovlpdiag


!    Name : cs_blkparlettf
!         :
!   Usage : values = cs_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a complex matrix in single precision.

pure function cs_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  complex(kind=sp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  complex(kind=sp)   :: cs_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call cs_blkparlett(upto,hci,a,info,cs_blkparlettf,errlog)
end function cs_blkparlettf


!    Name : cs_blkparlett
!         :
!   Usage : call cs_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a complex matrix in single precision.

pure subroutine cs_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  complex(kind=sp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  complex(kind=sp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  complex(kind=sp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  complex(kind=sp) :: rhs(1:size(a,1),1:size(a,2))
  complex(kind=sp) :: tmp(1:size(a,1),1:size(a,2))
  complex(kind=sp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = cmplx(0.0_sp,0.0_sp,sp)
  arg  = cmplx(0.0_sp,0.0_sp,sp)
  rhs  = cmplx(0.0_sp,0.0_sp,sp)
  tmp  = cmplx(0.0_sp,0.0_sp,sp)
  aph  = cmplx(0.0_sp,0.0_sp,sp)

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_sp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_sp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call cs_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call cs_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call cs_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_sp)) phin(j+1,j,:) = cmplx(0.0_sp,0.0_sp,sp)
  end do

end subroutine cs_blkparlett


!    Name : cs_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a complex matrix in single precision.

pure recursive &
subroutine cs_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=sp), intent(inout) :: phin(1:,1:,0:)
  complex(kind=sp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cs_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call cs_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call cs_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cs_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call cs_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call cs_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine cs_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  complex(kind=sp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=sp) :: smp, est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine cs_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine cs_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  complex(kind=sp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=sp) :: est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_sp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine cs_checkgrowth2x2

#endif

end subroutine cs_rec3x3nnm1


!    Name : cs_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``cs_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a complex matrix in single precision.

pure recursive &
subroutine cs_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=sp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cs_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call cs_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cs_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call cs_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine cs_rec3x3alln


!    Name : cs_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a complex matrix in single precision.

pure subroutine cs_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  complex(kind=sp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  complex(kind=sp), intent(out) :: values(:,:,0:)

  complex(kind=sp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = cmplx(1.0_sp,0.0_sp,sp)
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_sp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call cs_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call cs_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine cs_compovlpdiag


!    Name : rw_blkparlettf
!         :
!   Usage : values = rw_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a real matrix in double precision.

pure function rw_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  real   (kind=wp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  real   (kind=wp)   :: rw_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call rw_blkparlett(upto,hci,a,info,rw_blkparlettf,errlog)
end function rw_blkparlettf


!    Name : rw_blkparlett
!         :
!   Usage : call rw_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a real matrix in double precision.

pure subroutine rw_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  real   (kind=wp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  real   (kind=wp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  real   (kind=wp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  real   (kind=wp) :: rhs(1:size(a,1),1:size(a,2))
  real   (kind=wp) :: tmp(1:size(a,1),1:size(a,2))
  real   (kind=wp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = 0.0_wp
  arg  = 0.0_wp
  rhs  = 0.0_wp
  tmp  = 0.0_wp
  aph  = 0.0_wp

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_wp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_wp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call rw_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call rw_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call rw_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_wp)) phin(j+1,j,:) = 0.0_wp
  end do

end subroutine rw_blkparlett


!    Name : rw_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a real matrix in double precision.

pure recursive &
subroutine rw_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=wp), intent(inout) :: phin(1:,1:,0:)
  real   (kind=wp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rw_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call rw_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call rw_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rw_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call rw_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call rw_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine rw_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  real   (kind=wp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=wp) :: smp, est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine rw_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine rw_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  real   (kind=wp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=wp) :: est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_wp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine rw_checkgrowth2x2

#endif

end subroutine rw_rec3x3nnm1


!    Name : rw_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``rw_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a real matrix in double precision.

pure recursive &
subroutine rw_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=wp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rw_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call rw_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rw_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call rw_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine rw_rec3x3alln


!    Name : rw_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a real matrix in double precision.

pure subroutine rw_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  real   (kind=wp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  real   (kind=wp), intent(out) :: values(:,:,0:)

  real   (kind=wp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = 1.0_wp
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_wp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call rw_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call rw_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine rw_compovlpdiag


!    Name : cw_blkparlettf
!         :
!   Usage : values = cw_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a complex matrix in double precision.

pure function cw_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  complex(kind=wp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  complex(kind=wp)   :: cw_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call cw_blkparlett(upto,hci,a,info,cw_blkparlettf,errlog)
end function cw_blkparlettf


!    Name : cw_blkparlett
!         :
!   Usage : call cw_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a complex matrix in double precision.

pure subroutine cw_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  complex(kind=wp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  complex(kind=wp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  complex(kind=wp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  complex(kind=wp) :: rhs(1:size(a,1),1:size(a,2))
  complex(kind=wp) :: tmp(1:size(a,1),1:size(a,2))
  complex(kind=wp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = cmplx(0.0_wp,0.0_wp,wp)
  arg  = cmplx(0.0_wp,0.0_wp,wp)
  rhs  = cmplx(0.0_wp,0.0_wp,wp)
  tmp  = cmplx(0.0_wp,0.0_wp,wp)
  aph  = cmplx(0.0_wp,0.0_wp,wp)

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_wp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_wp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call cw_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call cw_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call cw_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_wp)) phin(j+1,j,:) = cmplx(0.0_wp,0.0_wp,wp)
  end do

end subroutine cw_blkparlett


!    Name : cw_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a complex matrix in double precision.

pure recursive &
subroutine cw_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=wp), intent(inout) :: phin(1:,1:,0:)
  complex(kind=wp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cw_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call cw_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call cw_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cw_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call cw_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call cw_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine cw_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  complex(kind=wp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=wp) :: smp, est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine cw_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine cw_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  complex(kind=wp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=wp) :: est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_wp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine cw_checkgrowth2x2

#endif

end subroutine cw_rec3x3nnm1


!    Name : cw_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``cw_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a complex matrix in double precision.

pure recursive &
subroutine cw_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=wp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cw_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call cw_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cw_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call cw_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine cw_rec3x3alln


!    Name : cw_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a complex matrix in double precision.

pure subroutine cw_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  complex(kind=wp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  complex(kind=wp), intent(out) :: values(:,:,0:)

  complex(kind=wp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = cmplx(1.0_wp,0.0_wp,wp)
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_wp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call cw_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call cw_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine cw_compovlpdiag

#ifdef __USE_TPREC

!    Name : rt_blkparlettf
!         :
!   Usage : values = rt_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a real matrix in triple precision.

pure function rt_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  real   (kind=tp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  real   (kind=tp)   :: rt_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call rt_blkparlett(upto,hci,a,info,rt_blkparlettf,errlog)
end function rt_blkparlettf


!    Name : rt_blkparlett
!         :
!   Usage : call rt_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a real matrix in triple precision.

pure subroutine rt_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  real   (kind=tp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  real   (kind=tp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  real   (kind=tp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  real   (kind=tp) :: rhs(1:size(a,1),1:size(a,2))
  real   (kind=tp) :: tmp(1:size(a,1),1:size(a,2))
  real   (kind=tp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = 0.0_tp
  arg  = 0.0_tp
  rhs  = 0.0_tp
  tmp  = 0.0_tp
  aph  = 0.0_tp

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_tp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_tp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call rt_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call rt_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call rt_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_tp)) phin(j+1,j,:) = 0.0_tp
  end do

end subroutine rt_blkparlett


!    Name : rt_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a real matrix in triple precision.

pure recursive &
subroutine rt_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=tp), intent(inout) :: phin(1:,1:,0:)
  real   (kind=tp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rt_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call rt_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call rt_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rt_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call rt_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call rt_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine rt_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  real   (kind=tp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=tp) :: smp, est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine rt_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine rt_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  real   (kind=tp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=tp) :: est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_tp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine rt_checkgrowth2x2

#endif

end subroutine rt_rec3x3nnm1


!    Name : rt_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``rt_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a real matrix in triple precision.

pure recursive &
subroutine rt_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=tp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rt_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call rt_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rt_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call rt_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine rt_rec3x3alln


!    Name : rt_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a real matrix in triple precision.

pure subroutine rt_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  real   (kind=tp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  real   (kind=tp), intent(out) :: values(:,:,0:)

  real   (kind=tp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = 1.0_tp
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_tp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call rt_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call rt_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine rt_compovlpdiag


!    Name : ct_blkparlettf
!         :
!   Usage : values = ct_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a complex matrix in triple precision.

pure function ct_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  complex(kind=tp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  complex(kind=tp)   :: ct_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call ct_blkparlett(upto,hci,a,info,ct_blkparlettf,errlog)
end function ct_blkparlettf


!    Name : ct_blkparlett
!         :
!   Usage : call ct_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a complex matrix in triple precision.

pure subroutine ct_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  complex(kind=tp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  complex(kind=tp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  complex(kind=tp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  complex(kind=tp) :: rhs(1:size(a,1),1:size(a,2))
  complex(kind=tp) :: tmp(1:size(a,1),1:size(a,2))
  complex(kind=tp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = cmplx(0.0_tp,0.0_tp,tp)
  arg  = cmplx(0.0_tp,0.0_tp,tp)
  rhs  = cmplx(0.0_tp,0.0_tp,tp)
  tmp  = cmplx(0.0_tp,0.0_tp,tp)
  aph  = cmplx(0.0_tp,0.0_tp,tp)

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_tp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_tp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call ct_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call ct_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call ct_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_tp)) phin(j+1,j,:) = cmplx(0.0_tp,0.0_tp,tp)
  end do

end subroutine ct_blkparlett


!    Name : ct_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a complex matrix in triple precision.

pure recursive &
subroutine ct_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=tp), intent(inout) :: phin(1:,1:,0:)
  complex(kind=tp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call ct_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call ct_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call ct_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call ct_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call ct_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call ct_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine ct_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  complex(kind=tp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=tp) :: smp, est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine ct_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine ct_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  complex(kind=tp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=tp) :: est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_tp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine ct_checkgrowth2x2

#endif

end subroutine ct_rec3x3nnm1


!    Name : ct_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``ct_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a complex matrix in triple precision.

pure recursive &
subroutine ct_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=tp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call ct_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call ct_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call ct_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call ct_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine ct_rec3x3alln


!    Name : ct_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a complex matrix in triple precision.

pure subroutine ct_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  complex(kind=tp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  complex(kind=tp), intent(out) :: values(:,:,0:)

  complex(kind=tp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = cmplx(1.0_tp,0.0_tp,tp)
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_tp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call ct_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call ct_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine ct_compovlpdiag

#endif

#ifdef __USE_QPREC

!    Name : rq_blkparlettf
!         :
!   Usage : values = rq_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a real matrix in quadruple precision.

pure function rq_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  real   (kind=qp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  real   (kind=qp)   :: rq_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call rq_blkparlett(upto,hci,a,info,rq_blkparlettf,errlog)
end function rq_blkparlettf


!    Name : rq_blkparlett
!         :
!   Usage : call rq_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a real matrix in quadruple precision.

pure subroutine rq_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  real   (kind=qp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  real   (kind=qp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  real   (kind=qp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  real   (kind=qp) :: rhs(1:size(a,1),1:size(a,2))
  real   (kind=qp) :: tmp(1:size(a,1),1:size(a,2))
  real   (kind=qp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = 0.0_qp
  arg  = 0.0_qp
  rhs  = 0.0_qp
  tmp  = 0.0_qp
  aph  = 0.0_qp

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_qp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_qp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call rq_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call rq_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call rq_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_qp)) phin(j+1,j,:) = 0.0_qp
  end do

end subroutine rq_blkparlett


!    Name : rq_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a real matrix in quadruple precision.

pure recursive &
subroutine rq_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=qp), intent(inout) :: phin(1:,1:,0:)
  real   (kind=qp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rq_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call rq_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call rq_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rq_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call rq_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call rq_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine rq_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  real   (kind=qp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=qp) :: smp, est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine rq_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine rq_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  real   (kind=qp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=qp) :: est(3), err(4), imp(1)
  real   (kind=RP) :: argref(size(arg,1),size(arg,2))
  real   (kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  real   (kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_qp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine rq_checkgrowth2x2

#endif

end subroutine rq_rec3x3nnm1


!    Name : rq_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``rq_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a real matrix in quadruple precision.

pure recursive &
subroutine rq_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  real   (kind=qp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rq_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call rq_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call rq_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call rq_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine rq_rec3x3alln


!    Name : rq_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a real matrix in quadruple precision.

pure subroutine rq_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  real   (kind=qp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  real   (kind=qp), intent(out) :: values(:,:,0:)

  real   (kind=qp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = 1.0_qp
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_qp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call rq_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call rq_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine rq_compovlpdiag


!    Name : cq_blkparlettf
!         :
!   Usage : values = cq_blkparlettf(upto,hci,a,info)
!         : This function is invoked through the generic name ``bspqtrphif''.
!         :
! Purpose : This function computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : values(1:N,1:N,0:upto) stores the values of phi-functions.
!         : values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         : values(1:N,1:N,0) is the matrix exponential.
!         :
!    Note : This routine treats a complex matrix in quadruple precision.

pure function cq_blkparlettf(upto,hci,a,info)
  integer,          intent(in) :: upto               ! 0 <= upto <= 5.
  complex(kind=qp), intent(in) :: hci, a(1:,1:)      ! must be quasi-upper tri.
  type(bspblock),   intent(in) :: info

  complex(kind=qp)   :: cq_blkparlettf(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog) :: errlog

  call cq_blkparlett(upto,hci,a,info,cq_blkparlettf,errlog)
end function cq_blkparlettf


!    Name : cq_blkparlett
!         :
!   Usage : call cq_blkparlett(upto,hci,a,info,values,errlog)
!         : This subroutine is invoked through the generic name ``bspqtrphi''.
!         :
! Purpose : This subroutine computes phi-functions by the modified
!         : block Schur--Parlett algorithm.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``hci'' is a scalar that should multiply the entries of ``A''.
!         : - The square array ``A'' is an upper quasi-triangular matrix.
!         : - ``info'' is of type(bspblock) and stores the blocking information
!         :   produced by ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions.
!         :   values(1:N,1:N,k) is phi_k(hci*A(1:N,1:N)).
!         :   values(1:N,1:N,0) is the matrix exponential.
!         : - ``errlog'' is of type(bspqtrphilog), and stores error information.
!         :
!    Note : This routine treats a complex matrix in quadruple precision.

pure subroutine cq_blkparlett(upto,hci,a,info,phin,errlog)
  integer,            intent(in ) :: upto
  complex(kind=qp),   intent(in ) :: hci, a(1:,1:)
  type(bspblock),     intent(in ) :: info
  complex(kind=qp),   intent(out) :: phin(1:size(a,1),1:size(a,2),0:upto)
  type(bspqtrphilog), intent(out) :: errlog

  complex(kind=qp) :: arg(1:size(a,1),1:size(a,2))   ! four working areas
  complex(kind=qp) :: rhs(1:size(a,1),1:size(a,2))
  complex(kind=qp) :: tmp(1:size(a,1),1:size(a,2))
  complex(kind=qp) :: aph(1:size(a,1),1:size(a,2))
  integer          :: i, j, m, n, mask(1:size(a,1))

  phin = cmplx(0.0_qp,0.0_qp,qp)
  arg  = cmplx(0.0_qp,0.0_qp,qp)
  rhs  = cmplx(0.0_qp,0.0_qp,qp)
  tmp  = cmplx(0.0_qp,0.0_qp,qp)
  aph  = cmplx(0.0_qp,0.0_qp,qp)

  m = size(a,1); n = size(a,2)

  ! Checking the argument

  errlog % error = 0
  if (upto < 0 .or. 5 < upto) then; errlog % error = 1; return; end if
  if (m <= 0) then; errlog % error = 2; return; end if
  if (m /= n) then; errlog % error = 3; return; end if

  ! The diagonal block of the matrix "a" should be at most 2-by-2.

  mask = 0
  do i=1,m-1; if (tiny(1.0_qp) <= abs(a(i+1,i))) mask(i) = 1; end do
  do i=2,m-1
    if (mask(i-1) == 1 .and. mask(i) == 1) then
      errlog % error = 4
      return
    end if
  end do

  ! The matrix "a" must be upper quasi-triangular.

  do j=1,n-2
    do i=j+2,m
      if (tiny(1.0_qp) <= abs(a(i,j))) then
        errlog % error = 5
        return
      end if
    end do
  end do

  ! Computation starts from here.
  ! The functions of diagonal blocks are computed.

  arg = hci * a                                        ! scale
  call cq_compovlpdiag(upto,arg,info,phin)             ! diagonal blocks

  ! Parlett's recursion.

  if (blkspfuncstride == 1 .or. upto == 0) then
    call cq_rec3x3alln(upto,arg,info,1,1,phin,rhs,tmp)
  else
    call cq_rec3x3nnm1(upto,arg,info,1,1,phin,rhs,tmp,aph)
  end if

  ! When a(j+1,j) is zero, phin(j+1,j,k) should also be zero
  ! for 0 <= k <= upto, because of the quasi-triangular structure.

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_qp)) phin(j+1,j,:) = cmplx(0.0_qp,0.0_qp,qp)
  end do

end subroutine cq_blkparlett


!    Name : cq_rec3x3nnm1
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. Two phi-functions are computed simultaneously
!         : in one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'', ``tmp'', and ``aph'' are used as working areas
!         :   throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion.
!         :
!         : This routine treats a complex matrix in quadruple precision.

pure recursive &
subroutine cq_rec3x3nnm1(upto,arg,info,depth,partition,phin,rhs,tmp,aph)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=qp), intent(inout) :: phin(1:,1:,0:)
  complex(kind=qp), intent(inout) :: rhs(1:,1:), tmp(1:,1:), aph(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then
    
    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks,
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cq_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)  ! upper block
      call cq_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)  ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1},  are computed simultaneously.

      do j=upto,0,-2

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j /= 0) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
      end do

#ifdef __CHECK_GROWTH
      call cq_checkgrowth3x3(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )   
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)
            
        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), aph (h1:t1,h3:t3)                              &
        )
        if (j >= 1) phin(h1:t1,h3:t3,j-1) = aph(h1:t1,h3:t3) + tmp(h1:t1,h3:t3)
        if (j >= 2) then
          phin(h1:t1,h3:t3,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h3:t3,j-1))               &
            +(arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j-1))               &
            +(arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j-1))
        end if
      end do

      end select

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.
 
      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cq_rec3x3nnm1(upto,arg,info,d,p(1),phin,rhs,tmp,aph)   ! upper block
      call cq_rec3x3nnm1(upto,arg,info,d,p(2),phin,rhs,tmp,aph)   ! lower block

      select case (blkspfuncstride)

      case (2)

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.
      ! Two phi-functions, phi_j and phi_{j-1}, are computed simultaenously.

      do j=upto,0,-2

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j /= 0) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
      end do

#ifdef __CHECK_GROWTH
      call cq_checkgrowth2x2(arg,phin,h,t)
#endif

      case (3)

      ! This mode computes three phi-functions in one recursion.
      ! This mode is not recommended to use, because it may provide
      ! inaccurate result.

      do j=upto,0,-3

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), aph (h1:t1,h2:t2)                              &
        )
        if (j >= 1) phin(h1:t1,h2:t2,j-1) = aph(h1:t1,h2:t2) + tmp(h1:t1,h2:t2)
        if (j >= 2) then
          phin(h1:t1,h2:t2,j-2) =                                             &
             (arg(h1:t1,h1:t1) .urtimes. phin(h1:t1,h2:t2,j-1))               &
            +(arg(h1:t1,h2:t2) .rutimes. phin(h2:t2,h2:t2,j-1))
        end if
      end do

      end select

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

#ifdef __CHECK_GROWTH

contains

! Evaluating the growth rate of diagonal perturbation. (3x3)

subroutine cq_checkgrowth3x3(arg,phin,h,t)
  use norm2estimate

  complex(kind=qp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h2, h3, t1, t2, t3
  real   (kind=qp) :: smp, est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h2 = h(2); h3 = h(3)
  t1 = t(1); t2 = t(2); t3 = t(3)

  smp =       sum(abs(arg(h1:t1,h2:t2))**2)
  smp = smp + sum(abs(arg(h1:t1,h3:t3))**2) * 2
  smp = smp + sum(abs(arg(h2:t2,h3:t3))**2)
  est(1) = sqrt(smp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h3:t3,h3:t3))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h2:t2,h2:t2,:)-phiref(h2:t2,h2:t2,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
                + sum(abs(phinrp(h1:t1,h2:t2,:)-phiref(h1:t1,h2:t2,:))**2)    &
                + sum(abs(phinrp(h2:t2,h3:t3,:)-phiref(h2:t2,h3:t3,:))**2)    &
             )

  err(1) = sqrt( sum(abs(phiref(h1:t1,h3:t3,:))**2) )
  err(2) = sqrt( sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h2:t2,h2:t2,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
               + sum(abs(phiref(h1:t1,h2:t2,:))**2)                           &
               + sum(abs(phiref(h2:t2,h3:t3,:))**2) )

  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,2e12.4,e24.16,e12.4)') h1, h2, h3, t3, &
  est(3), err(3)/err(4), err(1)/err(2), imp(1)
end subroutine cq_checkgrowth3x3

! Evaluating the growth rate of diagonal perturbation. (2x2)

subroutine cq_checkgrowth2x2(arg,phin,h,t)
  use norm2estimate

  complex(kind=qp), intent(in) :: arg(:,:), phin(:,:,:)
  integer,          intent(in) :: h(:), t(:)
  integer          :: n, upto, h1, h3, t1, t3
  real   (kind=qp) :: est(3), err(4), imp(1)
  complex(kind=RP) :: argref(size(arg,1),size(arg,2))
  complex(kind=RP) :: phiref(size(arg,1),size(arg,2),0:size(phin,3)-1)
  complex(kind=RP) :: phinrp(size(arg,1),size(arg,2),0:size(phin,3)-1)

  n = size(arg,1); upto = size(phin,3)
  h1 = h(1); h3 = h(2)
  t1 = t(1); t3 = t(2)

  est(1) = sqrt(sum(abs(arg(h1:t1,h2:t2))**2)) * sqrt(2.0_qp)
  est(2) = sepi2estf(arg(h1:t1,h1:t1), arg(h2:t2,h2:t2))
  est(3) = sqrt(1 + (est(1)*est(2))**2)

  open(unit=1,file='phiref.dat',status='old', &
       form='unformatted', access='sequential')
  read(1) argref, phiref
  close(1)

  phinrp = phin
  err(3) = sqrt(  sum(abs(phinrp(h1:t3,h1:t3,:)-phiref(h1:t3,h1:t3,:))**2) )
  err(4) = sqrt(                                                              &
                  sum(abs(phinrp(h1:t1,h1:t1,:)-phiref(h1:t1,h1:t1,:))**2)    &
                + sum(abs(phinrp(h3:t3,h3:t3,:)-phiref(h3:t3,h3:t3,:))**2)    &
             )

  err(1) = sqrt( err(3) / sum(abs(phiref(h1:t3,h1:t3,:))**2) )
  err(2) = sqrt( err(4)                                                       &
               / (                                                            &
                 sum(abs(phiref(h1:t1,h1:t1,:))**2)                           &
               + sum(abs(phiref(h3:t3,h3:t3,:))**2)                           &
             )                                                                &
           )
  imp(1) = err(3) / sqrt(sum(abs(phiref(:,:,:))**2))

  write(0,'(4i4,3e12.4)') h1, t1, h3, t3, est(3), err(3)/err(4), imp(1)

end subroutine cq_checkgrowth2x2

#endif

end subroutine cq_rec3x3nnm1


!    Name : cq_rec3x3alln
!         :
! Purpose : This subroutine performs the top-down recursion described
!         : in the paper. One phi-function is computed per one recursion.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``arg'' is the matrix argument of phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' specifies a diagonal block in that depth.
!         : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!  Output : - ``phin'' stores partially updated values of phi-functions.
!         : - ``rhs'' and ``tmp'' are used as working areas
!         :    throughout a recursion.
!         :
!    Note : ``phin'' stores completely updated values of phi-functions
!         : at the end of the top-down recursion. The accuracy may be
!         : slightly better than ``cq_rec3x3nnm1'' defined above.
!         :
!         : This routine treats a complex matrix in quadruple precision.

pure recursive &
subroutine cq_rec3x3alln(upto,arg,info,depth,partition,phin,rhs,tmp)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(in   ) :: arg(1:,1:)
  type(bspblock),   intent(in   ) :: info
  integer,          intent(in   ) :: depth, partition
  complex(kind=qp), intent(inout) :: phin(1:,1:,0:), rhs(1:,1:), tmp(1:,1:)

  integer :: d, j, k, p(2), h(3), t(3), h1, h2, h3, t1, t2, t3

  k = 2**(depth + 1) + 4*partition - 7

  if (1 <= info % heads(k+1) .and. 1 <= info % heads(k+2)) then

    ! The given block is further decomposed.

    h(1:3) = info % heads(k:k+2)
    t(1:2) = h(2:3) - 1
    t(3:3) = info % heads(k+3:k+3)

    if (t(1) + 1 /= h(3)) then                    ! overlapped blocking (3-by-3)

      h1 = h(1); t1 = t(1)                        ! the first  block
      h2 = h(2); t2 = t(2)                        ! the second block
      h3 = h(3); t3 = t(3)                        ! the third  block

      ! The functions of the upper and the lower diagonal blocks
      !    A(h(1):t(2),h(1):t(2)) and A(h(2):t(3),h(2):t(3)),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cq_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)  ! the upper block
      call cq_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)  ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h3:t3) = ( arg(h1:t1,h2:t2) .sqtimes. phin(h2:t2,h3:t3,j) ) &
                         + ( arg(h1:t1,h3:t3) .rutimes. phin(h3:t3,h3:t3,j) )
        rhs(h1:t1,h3:t3) = ( phin(h1:t1,h1:t1,j) .urtimes. arg(h1:t1,h3:t3) ) &
                         + ( phin(h1:t1,h2:t2,j) .sqtimes. arg(h2:t2,h3:t3) )
        rhs(h1:t1,h3:t3) = rhs(h1:t1,h3:t3) - tmp(h1:t1,h3:t3)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h3:t3,h3:t3),  rhs (h1:t1,h3:t3),         &
          phin(h1:t1,h3:t3,j), tmp (h1:t1,h3:t3)                              &
        )
      end do

    else                                          ! disjoint blocking (2-by-2)

      h1 = h(1); t1 = t(1)                        ! specifies the first  half.
      h2 = h(3); t2 = t(3)                        ! specifies the second half.

      ! The functions of the upper and the lower diagonal blocks
      !    A(h1:t1,h1:t1) and A(h2:t2,h2:t2),
      ! are computed by this subroutine itself.

      d = depth + 1
      p(1) = 2 * partition - 1; p(2) = p(1) + 1
      call cq_rec3x3alln(upto,arg,info,d,p(1),phin,rhs,tmp)   ! the upper block
      call cq_rec3x3alln(upto,arg,info,d,p(2),phin,rhs,tmp)   ! the lower block

      ! Continuous-time Sylvester equations are solved by the algorithm TSYLV_B.

      do j=0,upto

        tmp(h1:t1,h2:t2) = (arg (h1:t1,h2:t2  ) .rutimes. phin(h2:t2,h2:t2,j))
        rhs(h1:t1,h2:t2) = (phin(h1:t1,h1:t1,j) .urtimes. arg (h1:t1,h2:t2  ))
        rhs(h1:t1,h2:t2) =  rhs(h1:t1,h2:t2) - tmp(h1:t1,h2:t2)

        call axuppersylv(                                                     &
          arg (h1:t1,h1:t1),   arg (h2:t2,h2:t2),  rhs (h1:t1,h2:t2),         &
          phin(h1:t1,h2:t2,j), tmp (h1:t1,h2:t2)                              &
        )
      end do

    end if

  else if (info % heads(k+1) == -1 .and. info % heads(k+2) == -1) then

    ! The given block is not decomposed.

  else if (info % heads(k+1) == -2 .and. info % heads(k+2) == -2) then

    ! We do nothing for an already computed partition.

  end if

end subroutine cq_rec3x3alln


!    Name : cq_compovlpdiag
!         :
! Purpose : This subroutine computes phi-functions of diagonal blocks
!         : by the modified scaling and squaring method.
!         :
!   Input : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - The array ``A'' is the matrix argument for phi-functions.
!         : - ``info'' stores the blocking information produced by
!         :   ``bspqtr1dcmpf'' or ``bspqtr2dcmpf''.
!         :
!  Output : - values(1:N,1:N,0:upto) stores the values of phi-functions
!         :   for the overlapped block-diagonal region.
!         :
!    Note : When the module variable ``blockwisefunc'' is .true.,
!         : each diagonal block is treated independently of the other diagonal
!         : blocks. This mode is less efficient but may be more accurate.
!         :
!         : When the module variable ``blockwisefunc'' is .false.,
!         : all diagonal blocks are treated as an overlapped block-diagonal
!         : matrix. The functions are computed without repeating equivalent
!         : operations. This mode is more efficient but may be less accurate.
!         :
!         : If the spectral radius of the argument ``A'' is large,
!         : the block-wise computation is recommended.
!         :
!         : This routine treats a complex matrix in quadruple precision.

pure subroutine cq_compovlpdiag(upto,a,info,values)
  integer,          intent(in ) :: upto
  complex(kind=qp), intent(in ) :: a(:,:)
  type(bspblock),   intent(in ) :: info
  complex(kind=qp), intent(out) :: values(:,:,0:)

  complex(kind=qp) :: one
  integer :: j, mblk, n, nblk, ptr, sumblock, meanblock, blksize, h, t
  integer :: diag(2**(blksplimitdepth-1),2), subd(size(a,1))
  integer :: head(size(a,1)), tail(size(a,1)), nonz(size(a,1), size(a,2))
  type(mcpsqrlog) :: errlog

  one = cmplx(1.0_qp,0.0_qp,qp)
  n = size(a,1)                              ! is the order of the matrix.

  ptr = 1; diag = 0
  call getdiagonal(info,1,1,ptr,diag)        ! extracts diagonal blocks.
  ptr = ptr - 1                              ! is the # of diagonal blocks.

  call ovlp_uniqsort(ptr,diag,mblk)          ! sorts and uniq the blocks.
                                             ! ``mblk'' is the reduced number.

  ! We introduce another partition into the matrix ``A'' in order to
  ! apply block algorithms to the overlapped block-diagonal region.

  sumblock = 0
  do j=1,mblk; sumblock = sumblock + diag(j,2) - diag(j,1) + 1; end do
  meanblock = sumblock / mblk
  blksize = min(meanblock, mpt_blksize / 2)

  subd = 0
  do j=1,n-1
    if (tiny(1.0_qp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  call ovlp_blockdcmp(                                                       &
    blksize, n, diag(1:mblk,:), subd, head, tail, nonz, nblk                 &
  )

  ! The functions of diagonal blocks are computed.

  if (blockwisefunc) then

    ! The functions are computed for each diagonal block.
    ! Although equivalent operations must be repeated, accuracy is
    ! usually better for a matrix with a large spectral radius.

    do j=1,mblk
      h = diag(j,1)                          ! the head of a block
      t = diag(j,2)                          ! the tail of a block
      values(h:t,h:t,:) = sasqtrphif(upto,one,a(h:t,h:t))
    end do

  else

    ! The functions of the overlapped block-diagonal matrix is computed
    ! without repeating equivalent operations.

    call cq_ovlp_scsqr(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      upto, one, a, values, errlog                                           &
    )
  end if

  ! Clear unnecessary off-diagonal entries.

  do j=0,upto
    call cq_ovlp_clear(                                                      &
      diag(1:mblk,:), subd(1:n),                                             &
      head(1:nblk), tail(1:nblk), nonz(1:nblk,1:nblk),                       &
      values(:,:,j)                                                          &
    )
  end do
end subroutine cq_compovlpdiag

#endif

end module schurparlett


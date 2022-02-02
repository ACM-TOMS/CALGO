! This module determines a block decomposition based on a condition estimate.
! The matrix functions will be computed according to the blocking produced by
! ``bspqtr1dcmpf'' or ``bspqtr2dcmpf'' implemented in this module.
!
! References:
!
!   Isak Jonsson and Bo K{\aa}gstr{\"{o}}m,
!   Recursive blocked algorithms for solving triangular systems. I.
!   One-sided and coupled Sylvester-type matrix equations",
!   ACM Transactions on Mathematical Software",
!   Volume 28, Number 4, pp. 392-415, 2002,
!
!   _. _______,
!   On a sep-inverse estimate and its application to a modified
!   block Schur--Parlett algorithm,
!   submitted to ACM transactions on mathematical software, 2008.
!

module blockdecomp

  use floattypes
  use matrixpwrtag
  use ovlpmatrix
  use scalesquare
  use jspsylvester

  implicit none


  !    Name : bspqtr1dcmpf ( a generic name )
  !         :
  !   Usage : info = bspqtr1dcmpf(A)
  !         :
  ! Purpose : This function determines a block decomposition based on a
  !         : condition estimate for the recurrence. The function uses
  !         : the condition numbers defined in terms of the matrix 1-norm.
  !         : This is the 1-norm mode of the modified block Schur--Parlett
  !         : algorithm.
  !         :
  !   Input : The square array ``A'' is an upper quasi-triangular matrix,
  !         : and is a matrix argument for phi-functions.
  !         :
  !  Output : ``info'' is of type(bspblock), and stores the blocking
  !         : information. The derived type is defined later in this file.
  !         :
  !    Note : ``info'' will be used as an argument for ``bspqtrphif'' to
  !         : compute phi-functions. The function ``bspqtrphif'' is defined
  !         : in the file ``schurparlett.f95''.

  interface bspqtr1dcmpf
    module procedure rs_blk1dcmpf              ! Real    in Single    precision
    module procedure cs_blk1dcmpf              ! Complex in Single    precision
    module procedure rw_blk1dcmpf              ! Real    in Double    precision
    module procedure cw_blk1dcmpf              ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_blk1dcmpf              ! Real    in Triple    precision
    module procedure ct_blk1dcmpf              ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_blk1dcmpf              ! Real    in Quadruple precision
    module procedure cq_blk1dcmpf              ! Complex in Quadruple precision
#endif
  end interface bspqtr1dcmpf


  !    Name : cond1sylvf ( a generic name )
  !   Usage : estimate(1:2) = cond1sylvf(A,B)
  ! Purpose : This function estimates the sep-inverse function,
  !         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
  !         : by the recursive algorithm described in the accompanying paper.
  !   Input : A and B are upper quasi-triangular matrices.
  !  Output : estimate(2) is the estimated value of the sep-inverse function.
  !         : estimate(1) is not used in the current version.

  interface cond1sylvf
    module procedure rs_cond1sylvf             ! Real    in Single    precision
    module procedure cs_cond1sylvf             ! Complex in Single    precision
    module procedure rw_cond1sylvf             ! Real    in Double    precision
    module procedure cw_cond1sylvf             ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_cond1sylvf             ! Real    in Triple    precision
    module procedure ct_cond1sylvf             ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_cond1sylvf             ! Real    in Quadruple precision
    module procedure cq_cond1sylvf             ! Complex in Quadruple precision
#endif
  end interface cond1sylvf


  ! Module parameters, variables and a derived type.
  !
  ! ``blkspminblock''
  !   This parameter limits from below the size of a diagonal block.
  !   In terms of efficiency, excessively small diagonal blocks are not
  !   recommended.
  !
  ! ``blkspmaxdepth''
  !   This parameter limits from above the maximal depth of recursion used in
  !   Parlett's recurrence. In general, efficiency improves if we permit deeper
  !   recursion, while accuracy may deteriorate by repeated applications of
  !   the inverse Sylvester operator.
  !
  ! ``blkspfuncstride''
  !   If the value of this parameter is 2, two phi-functions are computed
  !   simultaneously in one recursion, while if the value is 1, the simultaneous
  !   computation is disabled. The former is more efficient but may be slightly
  !   less accurate than the latter.

  integer, parameter :: blksplimitdepth   = 14
  integer, parameter :: blkspminblock     =  9      ! must be greater than 3.
  integer            :: blkspmaxdepth     =  8
  integer            :: blkspfuncstride   =  2
  logical            :: blkspuseestimate  = .true.  ! only for our experiments.


  ! A derived type to define a block decomposition.

  type bspblock
    integer :: heads(4*(2**blksplimitdepth-1))
  end type bspblock


  ! Module variables for the 1-norm mode.
  !
  ! ``blksp1maxcondition''
  !   This variable limits from above the norm-wise condition number for
  !   the recurrence. The function, ``bspqtr1dcmpf'', determines a block
  !   decomposition so that the estimated condition number is less than
  !   the value of this variable for all partitions in the recurrence.
  !
  ! ``cond1sylvmaxdepth''
  !   This variable limits from above the maximal depth of recursion
  !   used in the recursive estimator, ``cond1sylvf''.
  !
  ! ``cond1sylvminblock''
  !   This variable limits from below the size of a diagonal block
  !   used in the recursive estimator, `cond1sylvf''.

  real(kind=wp) :: blksp1maxcondition = 150.0_wp
  integer       :: cond1sylvmaxdepth  = 4
  integer       :: cond1sylvminblock  = 4           ! must be greater than 3.


  private

  public :: bspblock
  public :: blksplimitdepth, blkspmaxdepth, blkspminblock
  public :: blkspfuncstride, blkspuseestimate

  public :: bspqtr1dcmpf, cond1sylvf
  public :: blksp1maxcondition, cond1sylvmaxdepth, cond1sylvminblock

contains


!    Name : rs_blk1dcmpf
!   Usage : info = rs_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a real matrix in single precision.

function rs_blk1dcmpf(a)
  real   (kind=sp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: rs_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_sp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  rs_blk1dcmpf % heads      = 0
  rs_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call rs_blk1dcmpr(a,subd,1,1,0,rs_blk1dcmpf)
end function rs_blk1dcmpf


!    Name : rs_blk1dcmpr
!         :
!   Usage : call rs_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a real matrix in single precision.

recursive subroutine rs_blk1dcmpr(a,subd,depth,partition,covered,info)
  real   (kind=sp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=sp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, sp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call rs_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call rs_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    real   (kind=sp), intent(in) :: a(:,:)
    real   (kind=sp)             :: getanestimate(2)

    real   (kind=sp) :: trm(2)
    real   (kind=sp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_sp); maxdif = 0.0_sp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_sp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_sp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_sp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,sp)
      getanestimate(2) = 1.0_sp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_sp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,sp)
        getanestimate(2) = 1.0_sp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rstr_norm1(a(h1:t1,h2:t2)),   &
            rstr_norm1(a(h1:t1,h3:t3)),   &
            rstr_normi(a(h1:t1,h3:t3)),   &
            rstr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rstr_norm1(a(h1:t1,h3:t3)), rstr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rstr_norm1(a(h1:t1,h2:t2)),   &
            rstr_norm1(a(h1:t1,h3:t3)),   &
            rstr_normi(a(h1:t1,h3:t3)),   &
            rstr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rstr_norm1(a(h1:t1,h3:t3)), rstr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine rs_blk1dcmpr


!    Name : rs_cond1sylvf 
!   Usage : estimate(1:2) = rs_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a real matrix in single precision.

function rs_cond1sylvf(a,b)                            ! top level
  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: rs_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_sp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_sp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  rs_cond1sylvf = rs_cond1sylvr(a,asd,b,bsd,1)
end function rs_cond1sylvf


!    Name : rs_cond1sylvr
!   Usage : call rs_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a real matrix in single precision.

recursive &
function rs_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  real   (kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=sp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=sp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_sp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = rstr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = rstr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rs_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = rs_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = rs_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = rs_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_sp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_sp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rstr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = rstr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rs_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = rs_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_sp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rstr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = rstr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rs_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = rs_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_sp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rstr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function rs_cond1sylvr


!    Name : cs_blk1dcmpf
!   Usage : info = cs_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a complex matrix in single precision.

function cs_blk1dcmpf(a)
  complex(kind=sp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: cs_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_sp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  cs_blk1dcmpf % heads      = 0
  cs_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call cs_blk1dcmpr(a,subd,1,1,0,cs_blk1dcmpf)
end function cs_blk1dcmpf


!    Name : cs_blk1dcmpr
!         :
!   Usage : call cs_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a complex matrix in single precision.

recursive subroutine cs_blk1dcmpr(a,subd,depth,partition,covered,info)
  complex(kind=sp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=sp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, sp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call cs_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call cs_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    complex(kind=sp), intent(in) :: a(:,:)
    real   (kind=sp)             :: getanestimate(2)

    real   (kind=sp) :: trm(2)
    real   (kind=sp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_sp); maxdif = 0.0_sp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_sp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_sp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_sp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,sp)
      getanestimate(2) = 1.0_sp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_sp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,sp)
        getanestimate(2) = 1.0_sp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cstr_norm1(a(h1:t1,h2:t2)),   &
            cstr_norm1(a(h1:t1,h3:t3)),   &
            cstr_normi(a(h1:t1,h3:t3)),   &
            cstr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cstr_norm1(a(h1:t1,h3:t3)), cstr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cstr_norm1(a(h1:t1,h2:t2)),   &
            cstr_norm1(a(h1:t1,h3:t3)),   &
            cstr_normi(a(h1:t1,h3:t3)),   &
            cstr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cstr_norm1(a(h1:t1,h3:t3)), cstr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine cs_blk1dcmpr


!    Name : cs_cond1sylvf 
!   Usage : estimate(1:2) = cs_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a complex matrix in single precision.

function cs_cond1sylvf(a,b)                            ! top level
  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=sp)             :: cs_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_sp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_sp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  cs_cond1sylvf = cs_cond1sylvr(a,asd,b,bsd,1)
end function cs_cond1sylvf


!    Name : cs_cond1sylvr
!   Usage : call cs_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a complex matrix in single precision.

recursive &
function cs_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  complex(kind=sp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=sp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=sp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_sp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = cstr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = cstr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cs_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = cs_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = cs_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = cs_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_sp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_sp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rstr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = cstr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cs_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = cs_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_sp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rstr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = cstr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cs_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = cs_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_sp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rstr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function cs_cond1sylvr


!    Name : rw_blk1dcmpf
!   Usage : info = rw_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a real matrix in double precision.

function rw_blk1dcmpf(a)
  real   (kind=wp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: rw_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_wp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  rw_blk1dcmpf % heads      = 0
  rw_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call rw_blk1dcmpr(a,subd,1,1,0,rw_blk1dcmpf)
end function rw_blk1dcmpf


!    Name : rw_blk1dcmpr
!         :
!   Usage : call rw_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a real matrix in double precision.

recursive subroutine rw_blk1dcmpr(a,subd,depth,partition,covered,info)
  real   (kind=wp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=wp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, wp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call rw_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call rw_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    real   (kind=wp), intent(in) :: a(:,:)
    real   (kind=wp)             :: getanestimate(2)

    real   (kind=wp) :: trm(2)
    real   (kind=wp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_wp); maxdif = 0.0_wp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_wp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_wp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_wp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,wp)
      getanestimate(2) = 1.0_wp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_wp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,wp)
        getanestimate(2) = 1.0_wp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rwtr_norm1(a(h1:t1,h2:t2)),   &
            rwtr_norm1(a(h1:t1,h3:t3)),   &
            rwtr_normi(a(h1:t1,h3:t3)),   &
            rwtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rwtr_norm1(a(h1:t1,h3:t3)), rwtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rwtr_norm1(a(h1:t1,h2:t2)),   &
            rwtr_norm1(a(h1:t1,h3:t3)),   &
            rwtr_normi(a(h1:t1,h3:t3)),   &
            rwtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rwtr_norm1(a(h1:t1,h3:t3)), rwtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine rw_blk1dcmpr


!    Name : rw_cond1sylvf 
!   Usage : estimate(1:2) = rw_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a real matrix in double precision.

function rw_cond1sylvf(a,b)                            ! top level
  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: rw_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_wp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_wp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  rw_cond1sylvf = rw_cond1sylvr(a,asd,b,bsd,1)
end function rw_cond1sylvf


!    Name : rw_cond1sylvr
!   Usage : call rw_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a real matrix in double precision.

recursive &
function rw_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  real   (kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=wp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=wp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_wp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = rwtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = rwtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rw_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = rw_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = rw_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = rw_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_wp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_wp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rwtr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = rwtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rw_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = rw_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_wp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rwtr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = rwtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rw_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = rw_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_wp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rwtr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function rw_cond1sylvr


!    Name : cw_blk1dcmpf
!   Usage : info = cw_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a complex matrix in double precision.

function cw_blk1dcmpf(a)
  complex(kind=wp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: cw_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_wp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  cw_blk1dcmpf % heads      = 0
  cw_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call cw_blk1dcmpr(a,subd,1,1,0,cw_blk1dcmpf)
end function cw_blk1dcmpf


!    Name : cw_blk1dcmpr
!         :
!   Usage : call cw_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a complex matrix in double precision.

recursive subroutine cw_blk1dcmpr(a,subd,depth,partition,covered,info)
  complex(kind=wp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=wp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, wp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call cw_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call cw_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    complex(kind=wp), intent(in) :: a(:,:)
    real   (kind=wp)             :: getanestimate(2)

    real   (kind=wp) :: trm(2)
    real   (kind=wp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_wp); maxdif = 0.0_wp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_wp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_wp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_wp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,wp)
      getanestimate(2) = 1.0_wp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_wp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,wp)
        getanestimate(2) = 1.0_wp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cwtr_norm1(a(h1:t1,h2:t2)),   &
            cwtr_norm1(a(h1:t1,h3:t3)),   &
            cwtr_normi(a(h1:t1,h3:t3)),   &
            cwtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cwtr_norm1(a(h1:t1,h3:t3)), cwtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cwtr_norm1(a(h1:t1,h2:t2)),   &
            cwtr_norm1(a(h1:t1,h3:t3)),   &
            cwtr_normi(a(h1:t1,h3:t3)),   &
            cwtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cwtr_norm1(a(h1:t1,h3:t3)), cwtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine cw_blk1dcmpr


!    Name : cw_cond1sylvf 
!   Usage : estimate(1:2) = cw_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a complex matrix in double precision.

function cw_cond1sylvf(a,b)                            ! top level
  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=wp)             :: cw_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_wp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_wp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  cw_cond1sylvf = cw_cond1sylvr(a,asd,b,bsd,1)
end function cw_cond1sylvf


!    Name : cw_cond1sylvr
!   Usage : call cw_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a complex matrix in double precision.

recursive &
function cw_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  complex(kind=wp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=wp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=wp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_wp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = cwtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = cwtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cw_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = cw_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = cw_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = cw_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_wp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_wp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rwtr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = cwtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cw_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = cw_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_wp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rwtr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = cwtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cw_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = cw_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_wp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rwtr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function cw_cond1sylvr

#ifdef __USE_TPREC

!    Name : rt_blk1dcmpf
!   Usage : info = rt_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a real matrix in triple precision.

function rt_blk1dcmpf(a)
  real   (kind=tp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: rt_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_tp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  rt_blk1dcmpf % heads      = 0
  rt_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call rt_blk1dcmpr(a,subd,1,1,0,rt_blk1dcmpf)
end function rt_blk1dcmpf


!    Name : rt_blk1dcmpr
!         :
!   Usage : call rt_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a real matrix in triple precision.

recursive subroutine rt_blk1dcmpr(a,subd,depth,partition,covered,info)
  real   (kind=tp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=tp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, tp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call rt_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call rt_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    real   (kind=tp), intent(in) :: a(:,:)
    real   (kind=tp)             :: getanestimate(2)

    real   (kind=tp) :: trm(2)
    real   (kind=tp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_tp); maxdif = 0.0_tp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_tp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_tp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_tp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,tp)
      getanestimate(2) = 1.0_tp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_tp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,tp)
        getanestimate(2) = 1.0_tp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rttr_norm1(a(h1:t1,h2:t2)),   &
            rttr_norm1(a(h1:t1,h3:t3)),   &
            rttr_normi(a(h1:t1,h3:t3)),   &
            rttr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rttr_norm1(a(h1:t1,h3:t3)), rttr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rttr_norm1(a(h1:t1,h2:t2)),   &
            rttr_norm1(a(h1:t1,h3:t3)),   &
            rttr_normi(a(h1:t1,h3:t3)),   &
            rttr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rttr_norm1(a(h1:t1,h3:t3)), rttr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine rt_blk1dcmpr


!    Name : rt_cond1sylvf 
!   Usage : estimate(1:2) = rt_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a real matrix in triple precision.

function rt_cond1sylvf(a,b)                            ! top level
  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: rt_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_tp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_tp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  rt_cond1sylvf = rt_cond1sylvr(a,asd,b,bsd,1)
end function rt_cond1sylvf


!    Name : rt_cond1sylvr
!   Usage : call rt_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a real matrix in triple precision.

recursive &
function rt_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  real   (kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=tp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=tp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_tp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = rttr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = rttr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rt_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = rt_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = rt_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = rt_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_tp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_tp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rttr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = rttr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rt_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = rt_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_tp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rttr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = rttr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rt_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = rt_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_tp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rttr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function rt_cond1sylvr


!    Name : ct_blk1dcmpf
!   Usage : info = ct_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a complex matrix in triple precision.

function ct_blk1dcmpf(a)
  complex(kind=tp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: ct_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_tp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  ct_blk1dcmpf % heads      = 0
  ct_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call ct_blk1dcmpr(a,subd,1,1,0,ct_blk1dcmpf)
end function ct_blk1dcmpf


!    Name : ct_blk1dcmpr
!         :
!   Usage : call ct_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a complex matrix in triple precision.

recursive subroutine ct_blk1dcmpr(a,subd,depth,partition,covered,info)
  complex(kind=tp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=tp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, tp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call ct_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call ct_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    complex(kind=tp), intent(in) :: a(:,:)
    real   (kind=tp)             :: getanestimate(2)

    real   (kind=tp) :: trm(2)
    real   (kind=tp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_tp); maxdif = 0.0_tp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_tp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_tp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_tp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,tp)
      getanestimate(2) = 1.0_tp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_tp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,tp)
        getanestimate(2) = 1.0_tp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cttr_norm1(a(h1:t1,h2:t2)),   &
            cttr_norm1(a(h1:t1,h3:t3)),   &
            cttr_normi(a(h1:t1,h3:t3)),   &
            cttr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cttr_norm1(a(h1:t1,h3:t3)), cttr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cttr_norm1(a(h1:t1,h2:t2)),   &
            cttr_norm1(a(h1:t1,h3:t3)),   &
            cttr_normi(a(h1:t1,h3:t3)),   &
            cttr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cttr_norm1(a(h1:t1,h3:t3)), cttr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine ct_blk1dcmpr


!    Name : ct_cond1sylvf 
!   Usage : estimate(1:2) = ct_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a complex matrix in triple precision.

function ct_cond1sylvf(a,b)                            ! top level
  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=tp)             :: ct_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_tp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_tp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  ct_cond1sylvf = ct_cond1sylvr(a,asd,b,bsd,1)
end function ct_cond1sylvf


!    Name : ct_cond1sylvr
!   Usage : call ct_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a complex matrix in triple precision.

recursive &
function ct_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  complex(kind=tp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=tp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=tp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_tp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = cttr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = cttr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = ct_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = ct_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = ct_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = ct_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_tp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_tp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rttr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = cttr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = ct_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = ct_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_tp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rttr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = cttr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = ct_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = ct_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_tp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rttr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function ct_cond1sylvr

#endif

#ifdef __USE_QPREC

!    Name : rq_blk1dcmpf
!   Usage : info = rq_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a real matrix in quadruple precision.

function rq_blk1dcmpf(a)
  real   (kind=qp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: rq_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_qp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  rq_blk1dcmpf % heads      = 0
  rq_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call rq_blk1dcmpr(a,subd,1,1,0,rq_blk1dcmpf)
end function rq_blk1dcmpf


!    Name : rq_blk1dcmpr
!         :
!   Usage : call rq_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a real matrix in quadruple precision.

recursive subroutine rq_blk1dcmpr(a,subd,depth,partition,covered,info)
  real   (kind=qp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=qp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, qp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call rq_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call rq_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    real   (kind=qp), intent(in) :: a(:,:)
    real   (kind=qp)             :: getanestimate(2)

    real   (kind=qp) :: trm(2)
    real   (kind=qp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_qp); maxdif = 0.0_qp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_qp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_qp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_qp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,qp)
      getanestimate(2) = 1.0_qp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_qp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,qp)
        getanestimate(2) = 1.0_qp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rqtr_norm1(a(h1:t1,h2:t2)),   &
            rqtr_norm1(a(h1:t1,h3:t3)),   &
            rqtr_normi(a(h1:t1,h3:t3)),   &
            rqtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rqtr_norm1(a(h1:t1,h3:t3)), rqtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            rqtr_norm1(a(h1:t1,h2:t2)),   &
            rqtr_norm1(a(h1:t1,h3:t3)),   &
            rqtr_normi(a(h1:t1,h3:t3)),   &
            rqtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(rqtr_norm1(a(h1:t1,h3:t3)), rqtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine rq_blk1dcmpr


!    Name : rq_cond1sylvf 
!   Usage : estimate(1:2) = rq_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a real matrix in quadruple precision.

function rq_cond1sylvf(a,b)                            ! top level
  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: rq_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_qp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_qp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  rq_cond1sylvf = rq_cond1sylvr(a,asd,b,bsd,1)
end function rq_cond1sylvf


!    Name : rq_cond1sylvr
!   Usage : call rq_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a real matrix in quadruple precision.

recursive &
function rq_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  real   (kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=qp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=qp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_qp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = rqtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = rqtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rq_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = rq_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = rq_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = rq_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_qp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_qp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rqtr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = rqtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rq_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = rq_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_qp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rqtr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = rqtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = rq_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = rq_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_qp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rqtr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function rq_cond1sylvr


!    Name : cq_blk1dcmpf
!   Usage : info = cq_blk1dcmpf(a)
! Purpose : This function determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!   Input : The array ``A'' is an upper quasi-triangular matrix
!         : and is a matrix argument for phi-functions.
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This function treats a complex matrix in quadruple precision.

function cq_blk1dcmpf(a)
  complex(kind=qp), intent(in) :: a(1:,1:)     ! must be upper quasi-triangular.
  type(bspblock)               :: cq_blk1dcmpf

  integer :: j, subd(1:size(a,1))

  subd = 0
  do j=1,size(a,1)-1
    if (tiny(1.0_qp) <= abs(a(j+1,j))) subd(j) = 1
  end do

  cq_blk1dcmpf % heads      = 0
  cq_blk1dcmpf % heads(1:4) = (/1, -1, -1, size(a,1)/)
  call cq_blk1dcmpr(a,subd,1,1,0,cq_blk1dcmpf)
end function cq_blk1dcmpf


!    Name : cq_blk1dcmpr
!         :
!   Usage : call cq_blk1dcmpr(a,subd,depth,partition,covered,info)
!         :
! Purpose : This subroutine determines a block decomposition based on
!         : a 1-norm estimate for the condition of the recurrence.
!         :
!   Input : - The array ``A'' is an upper quasi-triangular matrix.
!         : - ``subd'' is an array of size(a,1).
!         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``depth'' is the level of recursion.
!         : - ``partition'' is the index of a diagonal block in that depth.
!         : - ``covered'' is the maximum row/column number already treated
!         :   by the algorithm.
!         :
!  Output : ``info'' is of type(bspblock), and stores the blocking information.
!         : ``info'' may be used as an argument for the function ``bspqtrphif''
!         : implemented in ``schurparlett.f95''.
!    Note : This subroutine treats a complex matrix in quadruple precision.

recursive subroutine cq_blk1dcmpr(a,subd,depth,partition,covered,info)
  complex(kind=qp), intent(in   ) :: a(1:,1:)
  integer,          intent(in   ) :: subd(1:), depth, partition, covered
  type(bspblock),   intent(inout) :: info

  logical       :: yes
  integer       :: i, k, m, n, h(3), t(3), idx(size(a,1),2)
  integer       :: ptr, qtr, mest, nest, acts, bsize, nblk, d, p(2)
  real(kind=qp) :: norms(2), cond, maxcond

  k = address(depth,partition)             ! computes the index in ``info''.
  h(1) = info % heads(k  )                 ! the head of the first block
  t(3) = info % heads(k+3)                 ! the tail of the third block

  ! This routine immediately returns when
  ! - the recursion is too deep, or
  ! - the size of the given block, A(h(1):t(3),h(1):t(3)), is too small, or
  ! - the given block is already covered by another block.

  if (blkspmaxdepth <= depth .or.        &
      (t(3)-h(1)+1) <= blkspminblock) then ! The given block is undecomposable. 
    info % heads(k+1:k+2) = -1
    return
  else if (t(3) <= covered) then           ! The given block is already covered.
    info % heads(k+1:k+2) = -2
    return
  end if

  ! The diagonal block, A(h(1):t(3),h(1):t(3)), is decomposed in the following.

  ! A rough partition is determined so that each 2-by-2 diagonal block is not
  ! separated. A 10-by-10 partition is introduced into the submatrix,
  !        A(h(1):t(3), h(1):t(3)),
  ! when the submatrix is sufficiently large. If it is not sufficiently large,
  ! the number of partition may be smaller. Each partition is regarded as an
  ! atomic block in the succeeding process of this subroutine.

  n = t(3) - h(1) + 1
  bsize = max(4, n/10)
  nblk = n / bsize 
  if (mod(n,bsize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = h(1) - 1
  do i=1,nblk-1
    idx(i,1) = m + 1 
    m = m + acts - subd(m + acts)
    idx(i,2) = m
  end do
  idx(nblk,1) = m + 1; idx(nblk,2) = t(3)

  yes = .false.
  maxcond = real(blksp1maxcondition, qp)

  if (.not. blkspuseestimate) then         ! a decomposition without estimation
                                           ! (only for our experimental purpose)
    ptr = 1; qtr = nblk
    do
      norms = getanestimate(1,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then
        if (qtr - ptr <= 2) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
          exit
        else
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else
        if (1 < ptr .and. qtr < nblk) then
          yes = .true.; ptr = ptr - 1; qtr = qtr + 1;
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1
        else
          yes = .false.
        end if
        exit
      end if
    end do

  else                                     ! an estimation-based decomposition

    ! A 3-by-3 blocking,
    !     A(h(i):t(i), h(j):t(j)), 1 <= i, j <= 3,
    ! is determined based on condition estimation. The size of the overlapped
    ! diagonal block, A(h(2):t(2),h(2):t(2)), is reduced as small as possible.
    ! When nest==3, the     recursive estimator is invoked.
    ! When nest==4, the non-recursive estimator is invoked.

    ptr = 1; qtr = nblk
    nest = 2; mest = 4
    do
      norms = getanestimate(nest,a,subd,h(1),idx(ptr,2),idx(qtr,1),t(3))
      cond = norms(1) * norms(2)

      if (cond < maxcond) then             ! The blocking is well-conditioned.
        if (qtr - ptr <= 3) then
          yes = .true.
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1);
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Reduce the size of A_{22}.
          ptr = ptr + 1; qtr = qtr - 1
        end if
      else                                 ! The blocking is ill-conditioned.
        if (nest < mest) then              ! Choose a more accurate estimator.
          nest = nest + 1
        else if (1 < ptr .and. qtr < nblk) then
          yes = .true.                     ! Choose the previous blocking.
          ptr = ptr - 1; qtr = qtr + 1
          h(2) = idx(ptr,2) + 1; h(3) = idx(qtr,1)
          t(1) = idx(ptr,2)    ; t(2) = idx(qtr,1) - 1;
          exit
        else                               ! Blocking is prohibited.
          yes = .false.; exit
        end if
      end if
    end do
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  ! After checking consistency, the 3-by-3 blocking is permitted.

  yes = .false.
  if ( h(1) <= t(1) .and. t(1) < h(3) .and. h(3) <= t(3) ) then
    if ( subd(t(1)) == 0 .and. subd(h(3)-1) == 0) then
      h(2) = t(1) + 1; t(2) = h(3) - 1; yes = .true.
    end if
  end if

  if (.not. yes) then; info % heads(k+1:k+2) = -1; return; end if

  info % heads(k+1) = h(2)                 ! decomposed as this.
  info % heads(k+2) = h(3)

  ! The upper and the lower diagonal blocks,
  !    A(h(1):t(2), h(1):t(2)), and  A(h(2):t(3), h(2):t(3)),
  ! are recursively decomposed by this subroutine itself.
  ! The depth of recursion is incremented by one (``d = depth + 1'').
  ! The integer p(1) is the index of the upper block in the depth ``d''.
  ! The integer p(2) is the index of the lower block in the depth ``d''.

  d = depth + 1 
  p(1) = 2 * partition - 1; p(2) = p(1) + 1

  ! The upper diagonal block, A(h(1):t(2), h(1):t(2)), is decomposed.

  if (t(2) <= covered) then
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-2,-2,t(2)/)
  else
    info % heads(address(d,p(1)):address(d,p(1))+3) = (/h(1),-1,-1,t(2)/)
    call cq_blk1dcmpr(a,subd,d,p(1),covered,info)
  end if

  ! The lower diagonal block, A(h(2):t(3),h(2):t(3)), is decomposed.

  info % heads(address(d,p(2)):address(d,p(2))+3) = (/h(2),-1,-1,t(3)/)
  call cq_blk1dcmpr(a,subd,d,p(2),max(covered,t(2)),info)
  
contains


  !    Name : address
  ! Purpose : This function computes the index of a specified diagonal block
  !         : in a specified depth. The information on the block is found
  !         : in ``address(depth,partition)''-th entry of the integer array
  !         :``heads'' (= the member of the derived type ``bspblock'').
  !   Input : - ``depth'' is the level of recursion.
  !         : - ``partition'' is the index of a block in that depth.
  !  Output : the index of the integer array ``bspblock % heads''.

  function address(depth, partition)
    integer, intent(in) :: depth, partition
    integer             :: address
    address = 2**(depth + 1) + 4*partition - 7
  end function address


  !    Name : getanestimate
  !         :
  ! Purpose : This function estimates the condition number
  !         : for the 2-by-2 or 3-by-3 partition.
  !         :
  !   Input : - nest == 1 :(an eigenvalue-based estimate, not used)
  !         :   nest == 2 :(a theoretical bound, not implemented)
  !         :   nest == 3 : the     recursive condition estimate,
  !         :   nest == 4 : the non-recursive condition estimate.
  !         : - ``A'' is an upper quasi-triangular matrix.
  !         : - ``subd'' is an array of size(a,1).
  !         :   When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
  !         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
  !         : - ``h1'' is the row/column number of the head of the first block.
  !         : - ``t1'' is the row/column number of the tail of the first block.
  !         : - ``h3'' is the row/column number of the head of the last  block.
  !         : - ``t3'' is the row/column number of the tail of the last  block.
  !         :
  !  Output : When the return value is stored in, e.g., values(1:2),
  !         : the product, values(1)*values(2), is the required condition
  !         : number.

  function getanestimate(nest,a,subd,h1,t1,h3,t3)
    use norm1estimate

    integer,          intent(in) :: nest, subd(:), h1, t1, h3, t3
    complex(kind=qp), intent(in) :: a(:,:)
    real   (kind=qp)             :: getanestimate(2)

    real   (kind=qp) :: trm(2)
    real   (kind=qp) :: mindif, maxdif, tmp, smp, vmax
    integer          :: i, j, h2, t2

    ! We first check the minimum distance of the eigenvalues
    ! between diagonal blocks.

    mindif = huge(1.0_qp); maxdif = 0.0_qp; i = h1; j = h3
    do while(i <= t1)
      do while (j <= t3)

        trm(1) = abs(a(i,i) - a(j,j))
        trm(2) = 0.0_qp
        if      (subd(i)==1 .and. subd(j)==0) then
          trm(2) = sqrt(abs(a(i,i+1))) * sqrt(abs(a(i+1,i)))
        else if (subd(i)==0 .and. subd(j)==1) then
          trm(2) = sqrt(abs(a(j,j+1))) * sqrt(abs(a(j+1,j)))
        else if (subd(i)==1 .and. subd(j)==1) then;
          trm(2) =   sqrt(abs(a(i,i+1)))*sqrt(abs(a(i+1,i))) &
                   - sqrt(abs(a(j,j+1)))*sqrt(abs(a(j+1,j)))
        end if

        vmax = max(maxval(abs(trm(1:2))),1.0_qp)
        tmp = vmax * sqrt(sum(abs(trm(1:2)/vmax)**2))

        if (tmp < mindif) mindif = tmp
        if (maxdif < tmp) maxdif = tmp
    
        if (subd(j) == 1) then; j = j + 2; else; j = j + 1; end if
      end do
      if (subd(i) == 1) then; i = i + 2; else; i = i + 1; end if
    end do

    ! When the minimum distance is small, blocking is prohibited.

    if (mindif < sqrt(sqrt(tiny(1.0_qp)))) then
      getanestimate(1) = 2 * real(blksp1maxcondition,qp)
      getanestimate(2) = 1.0_qp
      return
    end if

    ! Condition estimation.

    select case (nest)
      case (1)
        getanestimate(1) = maxdif
        getanestimate(2) = 1.0_qp / mindif

      case (2)      ! If an accurate theoretical bounds are available,...
        getanestimate(1) = 2 * real(blksp1maxcondition,qp)
        getanestimate(2) = 1.0_qp

      case (3)

        ! Recursive estimation.

        getanestimate    = cond1sylvf(a(h1:t1,h1:t1), a(h3:t3,h3:t3))
        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cqtr_norm1(a(h1:t1,h2:t2)),   &
            cqtr_norm1(a(h1:t1,h3:t3)),   &
            cqtr_normi(a(h1:t1,h3:t3)),   &
            cqtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cqtr_norm1(a(h1:t1,h3:t3)), cqtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if

      case (4)

        ! Non-recursive estimation.

        if (t1+1 < h3) then            ! overlapped blocking (3x3)
          h2 = t1+1; t2 = h3-1
          smp = max(                      &
            cqtr_norm1(a(h1:t1,h2:t2)),   &
            cqtr_norm1(a(h1:t1,h3:t3)),   &
            cqtr_normi(a(h1:t1,h3:t3)),   &
            cqtr_normi(a(h2:t2,h3:t3))    &
          )
          getanestimate(1) = smp

        else                           ! disjoint   blocking (2x2)
          smp = max(cqtr_norm1(a(h1:t1,h3:t3)), cqtr_normi(a(h1:t1,h3:t3)))
          getanestimate(1) = smp
        end if
        getanestimate(2) = sepi1estf(a(h1:t1,h1:t1),a(h3:t3,h3:t3))

    end select
  end function getanestimate

end subroutine cq_blk1dcmpr


!    Name : cq_cond1sylvf 
!   Usage : estimate(1:2) = cq_cond1sylvf(a,b)
! Purpose : This function estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!   Input : ``A'' and ``B'' are upper quasi-triangular matrices.
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This function treats a complex matrix in quadruple precision.

function cq_cond1sylvf(a,b)                            ! top level
  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  real   (kind=qp)             :: cq_cond1sylvf(2)

  integer :: i, j, m, n
  integer :: asd(1:size(a,1)), bsd(1:size(b,1))        ! sub-diagonal flags

  m = size(a,1); n = size(b,1);

  asd(1:m) = 0; bsd(1:n) = 0
  do i=1,m-1; if (tiny(1.0_qp) <= abs(a(i+1,i))) asd(i) = 1; end do
  do j=1,n-1; if (tiny(1.0_qp) <= abs(b(j+1,j))) bsd(j) = 1; end do

  cq_cond1sylvf = cq_cond1sylvr(a,asd,b,bsd,1)
end function cq_cond1sylvf


!    Name : cq_cond1sylvr
!   Usage : call cq_cond1sylvr(a,asd,b,bsd,depth)
!         :
! Purpose : This subroutine estimates the sep-inverse function,
!         :   sep_1^{-1}(A,B) := || ( I (x) A - B^T (x) I )^{-1} ||_1,
!         : by the recursive algorithm described in the paper.
!         :
!   Input : - ``A'' and ``B'' are upper quasi-triangular matrices. 
!         : - ``depth'' is the level of recursion.
!         : - ``asd'' is an array of size(a,1).
!         :   When asd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   asd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
!         : - ``bsd'' is an array of size(b,1).
!         :   When bsd(i)==0, the (i+1,i)-th entry of ``B'' is zero, while
!         :   if   bsd(i)==1, the (i+1,i)-th entry of ``B'' is nonzero.
!         :
!  Output : estimate(2) is the estimated sep-inverse function.
!         : estimate(1) is not used in the current version.
!    Note : This subroutine treats a complex matrix in quadruple precision.

recursive &
function cq_cond1sylvr(a,asd,b,bsd,depth) result(values)
  use norm1estimate

  complex(kind=qp), intent(in) :: a(1:,1:), b(1:,1:)
  integer,          intent(in) :: asd(1:), bsd(1:), depth
  real   (kind=qp)             :: values(2)

  integer       :: maxdepth, minblock
  integer       :: m, n, mh, nh
  real(kind=qp) :: a12, b12, nrmtr(4,4), norms(2,4)

  maxdepth = cond1sylvmaxdepth
  minblock = cond1sylvminblock

  m = size(a,1); n = size(b,1)
  nrmtr = 0.0_qp; values = 0.0_wp

  if (maxdepth <= depth) then                          ! terminates the
                                                       ! recursion.
    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse
    return

  end if

  ! Block decomposition and recursion

  if      (minblock <= m .and. minblock <= n) then     ! Decompose both A and B.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1            ! to avoid separating
    nh = n/2; if (bsd(nh) == 1) nh = nh + 1            ! a two-by-two block.

    a12 = cqtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1
    b12 = cqtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cq_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,2) = cq_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(   1:nh,    1:nh), bsd(   1:nh), depth+1)
    norms(:,3) = cq_cond1sylvr(a(   1:mh,    1:mh), asd(   1:mh),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)
    norms(:,4) = cq_cond1sylvr(a(mh+1:m , mh+1:m ), asd(mh+1:m ),             &
                               b(nh+1:n , nh+1:n ), bsd(nh+1:n ), depth+1)

    nrmtr      = 0.0_qp
    nrmtr(1,1) = norms(2,2)
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2)
    nrmtr(3,1) = b12 * norms(2,2) * norms(2,4)
    nrmtr(4,1) = a12 * b12 * norms(2,2) * norms(2,3) * (norms(2,1) + norms(2,4))
    nrmtr(2,2) = norms(2,1)
    nrmtr(3,2) = 0.0_qp
    nrmtr(4,2) = b12 * norms(2,1) * norms(2,3)
    nrmtr(3,3) = norms(2,4)
    nrmtr(4,3) = a12 * norms(2,3) * norms(2,4)
    nrmtr(4,4) = norms(2,3)

    values(2)  = rqtr_norm1(nrmtr(1:4,1:4))

  else if (minblock <= m .and. minblock >  n) then     ! Decompose A only.

    mh = m/2; if (asd(mh) == 1) mh = mh + 1

    a12 = cqtr_norm1(a(1:mh,mh+1:m))                   ! || A_{12} ||_1

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cq_cond1sylvr(a(   1:mh,    1:mh),asd(   1:mh),b,bsd,depth+1)
    norms(:,2) = cq_cond1sylvr(a(mh+1:m , mh+1:m ),asd(mh+1:m ),b,bsd,depth+1)

    nrmtr(1,1) = norms(2,2)                   ; nrmtr(1,2) = 0.0_qp
    nrmtr(2,1) = a12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,1)

    values(2)  = rqtr_norm1(nrmtr(1:2,1:2))

  else if (minblock >  m .and. minblock <= n) then     ! Decompose B only.

    nh = n/2; if (bsd(nh) == 1) nh = nh + 1

    b12 = cqtr_normi(b(1:nh,nh+1:n))                   ! || B_{12} ||_\infty

    ! Recursion to estimate the sep-inverse function between diagonal blocks.

    norms(:,1) = cq_cond1sylvr(a,asd,b(   1:nh,   1:nh),bsd(   1:nh),depth+1)
    norms(:,2) = cq_cond1sylvr(a,asd,b(nh+1:n ,nh+1:n ),bsd(nh+1:n ),depth+1)

    nrmtr(1,1) = norms(2,1)                   ; nrmtr(1,2) = 0.0_qp
    nrmtr(2,1) = b12 * norms(2,1) * norms(2,2); nrmtr(2,2) = norms(2,2)

    values(2)  = rqtr_norm1(nrmtr(1:2,1:2))

  else                                                 ! No decomposition.

    values(2) = sepi1estf(a,b)                         ! 1-norm of the inverse

  end if
end function cq_cond1sylvr

#endif

end module blockdecomp


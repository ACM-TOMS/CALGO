! A test program to confirm that the programs are succesfully compiled.
! The main program is at the end of this file.

module selftestall

  use floattypes
  use randomnumber
  use matrixpwrtag

  use norm1estimate             ! the norm estimators for the Sylvester operator
  use scalesquare               ! the scaling and squaring method
  use blockdecomp               ! the estimation-based blocking
  use schurparlett              ! the modified block Schur--Parlett algorithm

  implicit none
  public

contains

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine sp_test_realdiag
  integer, parameter :: upto = 5
  real(kind=sp) :: a(16), u(16,0:upto), v(16,0:upto)
  real(kind=sp) :: err(0:upto), eps
  integer       :: j
  character     :: judge*4

  do j=1,16
   a(j) = real((j-11)*2, sp)
  end do

  u = sasmtrphif(upto,0.5_sp,a)
  call sp_realref(v)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (1000*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'real    sp ',err,judge
end subroutine sp_test_realdiag

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine sp_test_compdiag
  integer, parameter :: upto = 5
  complex(kind=sp) :: a(32), u(32,0:upto), v(32,0:upto)
  real   (kind=sp) :: err(0:upto), eps
  integer          :: i, j
  character        :: judge*4

  do i=0,1
    do j=1,16
      a(j+i*16) = cmplx((j-11)*2,(i*2-1)*2,sp)
    end do
  end do

  u = sasmtrphif(upto,cmplx(0.5_sp,0.0_sp,sp),a)
  call sp_compref(v)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'complex sp ',err,judge
end subroutine sp_test_compdiag

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine wp_test_realdiag
  integer, parameter :: upto = 5
  real(kind=wp) :: a(16), u(16,0:upto), v(16,0:upto)
  real(kind=wp) :: err(0:upto), eps
  integer       :: j
  character     :: judge*4

  do j=1,16
   a(j) = real((j-11)*2, wp)
  end do

  u = sasmtrphif(upto,0.5_wp,a)
  call wp_realref(v)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (1000*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'real    wp ',err,judge
end subroutine wp_test_realdiag

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine wp_test_compdiag
  integer, parameter :: upto = 5
  complex(kind=wp) :: a(32), u(32,0:upto), v(32,0:upto)
  real   (kind=wp) :: err(0:upto), eps
  integer          :: i, j
  character        :: judge*4

  do i=0,1
    do j=1,16
      a(j+i*16) = cmplx((j-11)*2,(i*2-1)*2,wp)
    end do
  end do

  u = sasmtrphif(upto,cmplx(0.5_wp,0.0_wp,wp),a)
  call wp_compref(v)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'complex wp ',err,judge
end subroutine wp_test_compdiag


#ifdef __USE_TPREC

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine tp_test_realdiag
  integer, parameter :: upto = 5
  real(kind=tp) :: a(16), u(16,0:upto), v(16,0:upto)
  real(kind=tp) :: err(0:upto), eps
  integer       :: j
  character     :: judge*4

  do j=1,16
   a(j) = real((j-11)*2, tp)
  end do

  u = sasmtrphif(upto,0.5_tp,a)
  call tp_realref(v)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (1000*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'real    tp ',err,judge
end subroutine tp_test_realdiag

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine tp_test_compdiag
  integer, parameter :: upto = 5
  complex(kind=tp) :: a(32), u(32,0:upto), v(32,0:upto)
  real   (kind=tp) :: err(0:upto), eps
  integer          :: i, j
  character        :: judge*4

  do i=0,1
    do j=1,16
      a(j+i*16) = cmplx((j-11)*2,(i*2-1)*2,tp)
    end do
  end do

  u = sasmtrphif(upto,cmplx(0.5_tp,0.0_tp,tp),a)
  call tp_compref(v)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'complex tp ',err,judge
end subroutine tp_test_compdiag


#endif

#ifdef __USE_QPREC

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine qp_test_realdiag
  integer, parameter :: upto = 5
  real(kind=qp) :: a(16), u(16,0:upto), v(16,0:upto)
  real(kind=qp) :: err(0:upto), eps
  integer       :: j
  character     :: judge*4

  do j=1,16
   a(j) = real((j-11)*2, qp)
  end do

  u = sasmtrphif(upto,0.5_qp,a)
  call qp_realref(v)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (1000*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'real    qp ',err,judge
end subroutine qp_test_realdiag

! A test of the S & S for a diagonal matrix. We compare the numerical results
! with those obtained by arbitrary precision arithmetic.

subroutine qp_test_compdiag
  integer, parameter :: upto = 5
  complex(kind=qp) :: a(32), u(32,0:upto), v(32,0:upto)
  real   (kind=qp) :: err(0:upto), eps
  integer          :: i, j
  character        :: judge*4

  do i=0,1
    do j=1,16
      a(j+i*16) = cmplx((j-11)*2,(i*2-1)*2,qp)
    end do
  end do

  u = sasmtrphif(upto,cmplx(0.5_qp,0.0_qp,qp),a)
  call qp_compref(v)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = sum(abs(u(:,j)-v(:,j))) / sum(abs(u(:,j)))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do
  
  write(*,'(a11,6e10.3,a6)') 'complex qp ',err,judge
end subroutine qp_test_compdiag


#endif
! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine rs_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  real   (kind=sp) :: a(n), s(n,n), t(n,n)
  real   (kind=sp) :: b(n,  0:upto), c(n,n)
  real   (kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: pi, argument, rpart, ipart
  real   (kind=sp) :: err(0:upto)
  real   (kind=sp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_sp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_sp, 0.5_sp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = rpart + ipart*0
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_sp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = rpart + ipart*0
    end do
  end do
  s = rssq_eye(n) - s / rssq_norm1(s) / (1.0_sp*2)
  t = s; call rssq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, 1.0_sp/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = rssq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, 1.0_sp/2, rssq_ainvbf(t,perm,c))

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rssq_norm1(u(:,:,j)-v(:,:,j)) / rssq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    sp ',err,judge
end subroutine rs_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine rs_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  real   (kind=sp) :: a(n,n)
  real   (kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: pi, argument
  real   (kind=sp) :: err(0:upto)
  real   (kind=sp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_sp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(i) = 1
  end do

  a = 0.0_sp
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = rrandom(10.0_sp,0.0_sp)
    end do
  end do

  i = 1
  do while (i <= n)    ! diagonal entries
    if      (subd(i) == 0) then
      a(i,i) = rrandom(radius, (center-radius)/2)
      i = i + 1
    else if (subd(i) == 1) then
      argument = rrandom(pi/4,pi)
      radii = radius * rrandom(1.0_sp,0.5_sp)
      rpart = radii * cos(argument) + center
      ipart = radii * sin(argument)
      a(i  ,i) =  rpart; a(i  ,i+1) = ipart
      a(i+1,i) = -ipart; a(i+1,i+1) = rpart
      i = i + 2
    end if
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,1.0_sp/2,a)
  v = sasqtrphif(upto,1.0_sp/2,a)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rssq_norm1(u(:,:,j)-v(:,:,j)) / rssq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    sp ',err,judge
end subroutine rs_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine rs_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  real   (kind=sp) :: a(n,n)
  real   (kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = 0.0_sp
  do i=1,n
    do j=i+1,n
      a(i,j) = rrandom(1.0_sp,0.0_sp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(i) = 1
  end do

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * 1.0_sp / (1.0_sp * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  rrandom(1.0_sp,0.0_sp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, 1.0_sp/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, 1.0_sp/2, a, info)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rssq_norm1(u(:,:,j)-v(:,:,j)) / rssq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    sp ',err,judge
end subroutine rs_test_blksp4triangular

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine cs_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  complex(kind=sp) :: a(n), s(n,n), t(n,n)
  complex(kind=sp) :: b(n,  0:upto), c(n,n)
  complex(kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: pi, argument, rpart, ipart
  real   (kind=sp) :: err(0:upto)
  real   (kind=sp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_sp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_sp, 0.5_sp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = cmplx(rpart,ipart,sp)
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_sp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = cmplx(rpart,ipart,sp)
    end do
  end do
  s = cssq_eye(n) - s / cssq_norm1(s) / (cmplx(1.0_sp,0.0_sp,sp)*2)
  t = s; call cssq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, cmplx(1.0_sp,0.0_sp,sp)/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = cssq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, cmplx(1.0_sp,0.0_sp,sp)/2, cssq_ainvbf(t,perm,c))

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cssq_norm1(u(:,:,j)-v(:,:,j)) / cssq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex sp ',err,judge
end subroutine cs_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine cs_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  complex(kind=sp) :: a(n,n)
  complex(kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: pi, argument
  real   (kind=sp) :: err(0:upto)
  real   (kind=sp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_sp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(i) = 1
  end do

  a = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = crandom(10.0_sp,0.0_sp)
    end do
  end do

  do i=1,n
    argument = rrandom(pi/4,pi)
    radii = radius * rrandom(1.0_sp,0.5_sp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i,i) = cmplx(rpart,ipart,sp)
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,cmplx(1.0_sp,0.0_sp,sp)/2,a)
  v = sasqtrphif(upto,cmplx(1.0_sp,0.0_sp,sp)/2,a)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cssq_norm1(u(:,:,j)-v(:,:,j)) / cssq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex sp ',err,judge
end subroutine cs_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine cs_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=sp), parameter :: center = 2.0_sp, radius = 20.0_sp

  complex(kind=sp) :: a(n,n)
  complex(kind=sp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=sp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,n
    do j=i+1,n
      a(i,j) = crandom(1.0_sp,0.0_sp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_sp < rrandom(1.0_sp,0.5_sp)) subd(i) = 1
  end do
  subd = 0

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * cmplx(1.0_sp,0.0_sp,sp) / (cmplx(1.0_sp,0.0_sp,sp) * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  crandom(1.0_sp,0.0_sp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, cmplx(1.0_sp,0.0_sp,sp)/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, cmplx(1.0_sp,0.0_sp,sp)/2, a, info)

  eps = epsilon(1.0_sp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cssq_norm1(u(:,:,j)-v(:,:,j)) / cssq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex sp ',err,judge
end subroutine cs_test_blksp4triangular

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine rw_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  real   (kind=wp) :: a(n), s(n,n), t(n,n)
  real   (kind=wp) :: b(n,  0:upto), c(n,n)
  real   (kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: pi, argument, rpart, ipart
  real   (kind=wp) :: err(0:upto)
  real   (kind=wp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_wp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_wp, 0.5_wp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = rpart + ipart*0
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_wp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = rpart + ipart*0
    end do
  end do
  s = rwsq_eye(n) - s / rwsq_norm1(s) / (1.0_wp*2)
  t = s; call rwsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, 1.0_wp/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = rwsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, 1.0_wp/2, rwsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rwsq_norm1(u(:,:,j)-v(:,:,j)) / rwsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    wp ',err,judge
end subroutine rw_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine rw_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  real   (kind=wp) :: a(n,n)
  real   (kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: pi, argument
  real   (kind=wp) :: err(0:upto)
  real   (kind=wp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_wp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(i) = 1
  end do

  a = 0.0_wp
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = rrandom(10.0_wp,0.0_wp)
    end do
  end do

  i = 1
  do while (i <= n)    ! diagonal entries
    if      (subd(i) == 0) then
      a(i,i) = rrandom(radius, (center-radius)/2)
      i = i + 1
    else if (subd(i) == 1) then
      argument = rrandom(pi/4,pi)
      radii = radius * rrandom(1.0_wp,0.5_wp)
      rpart = radii * cos(argument) + center
      ipart = radii * sin(argument)
      a(i  ,i) =  rpart; a(i  ,i+1) = ipart
      a(i+1,i) = -ipart; a(i+1,i+1) = rpart
      i = i + 2
    end if
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,1.0_wp/2,a)
  v = sasqtrphif(upto,1.0_wp/2,a)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rwsq_norm1(u(:,:,j)-v(:,:,j)) / rwsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    wp ',err,judge
end subroutine rw_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine rw_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  real   (kind=wp) :: a(n,n)
  real   (kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = 0.0_wp
  do i=1,n
    do j=i+1,n
      a(i,j) = rrandom(1.0_wp,0.0_wp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(i) = 1
  end do

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * 1.0_wp / (1.0_wp * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  rrandom(1.0_wp,0.0_wp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, 1.0_wp/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, 1.0_wp/2, a, info)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rwsq_norm1(u(:,:,j)-v(:,:,j)) / rwsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    wp ',err,judge
end subroutine rw_test_blksp4triangular

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine cw_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  complex(kind=wp) :: a(n), s(n,n), t(n,n)
  complex(kind=wp) :: b(n,  0:upto), c(n,n)
  complex(kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: pi, argument, rpart, ipart
  real   (kind=wp) :: err(0:upto)
  real   (kind=wp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_wp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_wp, 0.5_wp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = cmplx(rpart,ipart,wp)
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_wp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = cmplx(rpart,ipart,wp)
    end do
  end do
  s = cwsq_eye(n) - s / cwsq_norm1(s) / (cmplx(1.0_wp,0.0_wp,wp)*2)
  t = s; call cwsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, cmplx(1.0_wp,0.0_wp,wp)/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = cwsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, cmplx(1.0_wp,0.0_wp,wp)/2, cwsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cwsq_norm1(u(:,:,j)-v(:,:,j)) / cwsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex wp ',err,judge
end subroutine cw_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine cw_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  complex(kind=wp) :: a(n,n)
  complex(kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: pi, argument
  real   (kind=wp) :: err(0:upto)
  real   (kind=wp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_wp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(i) = 1
  end do

  a = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = crandom(10.0_wp,0.0_wp)
    end do
  end do

  do i=1,n
    argument = rrandom(pi/4,pi)
    radii = radius * rrandom(1.0_wp,0.5_wp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i,i) = cmplx(rpart,ipart,wp)
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,cmplx(1.0_wp,0.0_wp,wp)/2,a)
  v = sasqtrphif(upto,cmplx(1.0_wp,0.0_wp,wp)/2,a)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cwsq_norm1(u(:,:,j)-v(:,:,j)) / cwsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex wp ',err,judge
end subroutine cw_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine cw_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=wp), parameter :: center = 2.0_wp, radius = 20.0_wp

  complex(kind=wp) :: a(n,n)
  complex(kind=wp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=wp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,n
    do j=i+1,n
      a(i,j) = crandom(1.0_wp,0.0_wp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_wp < rrandom(1.0_wp,0.5_wp)) subd(i) = 1
  end do
  subd = 0

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * cmplx(1.0_wp,0.0_wp,wp) / (cmplx(1.0_wp,0.0_wp,wp) * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  crandom(1.0_wp,0.0_wp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, cmplx(1.0_wp,0.0_wp,wp)/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, cmplx(1.0_wp,0.0_wp,wp)/2, a, info)

  eps = epsilon(1.0_wp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cwsq_norm1(u(:,:,j)-v(:,:,j)) / cwsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex wp ',err,judge
end subroutine cw_test_blksp4triangular


#ifdef __USE_TPREC

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine rt_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  real   (kind=tp) :: a(n), s(n,n), t(n,n)
  real   (kind=tp) :: b(n,  0:upto), c(n,n)
  real   (kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: pi, argument, rpart, ipart
  real   (kind=tp) :: err(0:upto)
  real   (kind=tp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_tp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_tp, 0.5_tp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = rpart + ipart*0
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_tp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = rpart + ipart*0
    end do
  end do
  s = rtsq_eye(n) - s / rtsq_norm1(s) / (1.0_tp*2)
  t = s; call rtsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, 1.0_tp/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = rtsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, 1.0_tp/2, rtsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rtsq_norm1(u(:,:,j)-v(:,:,j)) / rtsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    tp ',err,judge
end subroutine rt_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine rt_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  real   (kind=tp) :: a(n,n)
  real   (kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: pi, argument
  real   (kind=tp) :: err(0:upto)
  real   (kind=tp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_tp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(i) = 1
  end do

  a = 0.0_tp
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = rrandom(10.0_tp,0.0_tp)
    end do
  end do

  i = 1
  do while (i <= n)    ! diagonal entries
    if      (subd(i) == 0) then
      a(i,i) = rrandom(radius, (center-radius)/2)
      i = i + 1
    else if (subd(i) == 1) then
      argument = rrandom(pi/4,pi)
      radii = radius * rrandom(1.0_tp,0.5_tp)
      rpart = radii * cos(argument) + center
      ipart = radii * sin(argument)
      a(i  ,i) =  rpart; a(i  ,i+1) = ipart
      a(i+1,i) = -ipart; a(i+1,i+1) = rpart
      i = i + 2
    end if
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,1.0_tp/2,a)
  v = sasqtrphif(upto,1.0_tp/2,a)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rtsq_norm1(u(:,:,j)-v(:,:,j)) / rtsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    tp ',err,judge
end subroutine rt_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine rt_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  real   (kind=tp) :: a(n,n)
  real   (kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = 0.0_tp
  do i=1,n
    do j=i+1,n
      a(i,j) = rrandom(1.0_tp,0.0_tp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(i) = 1
  end do

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * 1.0_tp / (1.0_tp * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  rrandom(1.0_tp,0.0_tp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, 1.0_tp/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, 1.0_tp/2, a, info)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rtsq_norm1(u(:,:,j)-v(:,:,j)) / rtsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    tp ',err,judge
end subroutine rt_test_blksp4triangular

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine ct_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  complex(kind=tp) :: a(n), s(n,n), t(n,n)
  complex(kind=tp) :: b(n,  0:upto), c(n,n)
  complex(kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: pi, argument, rpart, ipart
  real   (kind=tp) :: err(0:upto)
  real   (kind=tp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_tp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_tp, 0.5_tp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = cmplx(rpart,ipart,tp)
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_tp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = cmplx(rpart,ipart,tp)
    end do
  end do
  s = ctsq_eye(n) - s / ctsq_norm1(s) / (cmplx(1.0_tp,0.0_tp,tp)*2)
  t = s; call ctsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, cmplx(1.0_tp,0.0_tp,tp)/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = ctsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, cmplx(1.0_tp,0.0_tp,tp)/2, ctsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = ctsq_norm1(u(:,:,j)-v(:,:,j)) / ctsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex tp ',err,judge
end subroutine ct_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine ct_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  complex(kind=tp) :: a(n,n)
  complex(kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: pi, argument
  real   (kind=tp) :: err(0:upto)
  real   (kind=tp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_tp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(i) = 1
  end do

  a = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = crandom(10.0_tp,0.0_tp)
    end do
  end do

  do i=1,n
    argument = rrandom(pi/4,pi)
    radii = radius * rrandom(1.0_tp,0.5_tp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i,i) = cmplx(rpart,ipart,tp)
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,cmplx(1.0_tp,0.0_tp,tp)/2,a)
  v = sasqtrphif(upto,cmplx(1.0_tp,0.0_tp,tp)/2,a)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = ctsq_norm1(u(:,:,j)-v(:,:,j)) / ctsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex tp ',err,judge
end subroutine ct_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine ct_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=tp), parameter :: center = 2.0_tp, radius = 20.0_tp

  complex(kind=tp) :: a(n,n)
  complex(kind=tp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=tp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,n
    do j=i+1,n
      a(i,j) = crandom(1.0_tp,0.0_tp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_tp < rrandom(1.0_tp,0.5_tp)) subd(i) = 1
  end do
  subd = 0

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * cmplx(1.0_tp,0.0_tp,tp) / (cmplx(1.0_tp,0.0_tp,tp) * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  crandom(1.0_tp,0.0_tp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, cmplx(1.0_tp,0.0_tp,tp)/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, cmplx(1.0_tp,0.0_tp,tp)/2, a, info)

  eps = epsilon(1.0_tp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = ctsq_norm1(u(:,:,j)-v(:,:,j)) / ctsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex tp ',err,judge
end subroutine ct_test_blksp4triangular

#endif

#ifdef __USE_QPREC

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine rq_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  real   (kind=qp) :: a(n), s(n,n), t(n,n)
  real   (kind=qp) :: b(n,  0:upto), c(n,n)
  real   (kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: pi, argument, rpart, ipart
  real   (kind=qp) :: err(0:upto)
  real   (kind=qp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_qp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_qp, 0.5_qp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = rpart + ipart*0
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_qp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = rpart + ipart*0
    end do
  end do
  s = rqsq_eye(n) - s / rqsq_norm1(s) / (1.0_qp*2)
  t = s; call rqsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, 1.0_qp/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = rqsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, 1.0_qp/2, rqsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rqsq_norm1(u(:,:,j)-v(:,:,j)) / rqsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    qp ',err,judge
end subroutine rq_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine rq_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  real   (kind=qp) :: a(n,n)
  real   (kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: pi, argument
  real   (kind=qp) :: err(0:upto)
  real   (kind=qp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_qp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(i) = 1
  end do

  a = 0.0_qp
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = rrandom(10.0_qp,0.0_qp)
    end do
  end do

  i = 1
  do while (i <= n)    ! diagonal entries
    if      (subd(i) == 0) then
      a(i,i) = rrandom(radius, (center-radius)/2)
      i = i + 1
    else if (subd(i) == 1) then
      argument = rrandom(pi/4,pi)
      radii = radius * rrandom(1.0_qp,0.5_qp)
      rpart = radii * cos(argument) + center
      ipart = radii * sin(argument)
      a(i  ,i) =  rpart; a(i  ,i+1) = ipart
      a(i+1,i) = -ipart; a(i+1,i+1) = rpart
      i = i + 2
    end if
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,1.0_qp/2,a)
  v = sasqtrphif(upto,1.0_qp/2,a)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rqsq_norm1(u(:,:,j)-v(:,:,j)) / rqsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    qp ',err,judge
end subroutine rq_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine rq_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  real   (kind=qp) :: a(n,n)
  real   (kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = 0.0_qp
  do i=1,n
    do j=i+1,n
      a(i,j) = rrandom(1.0_qp,0.0_qp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(i) = 1
  end do

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * 1.0_qp / (1.0_qp * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  rrandom(1.0_qp,0.0_qp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, 1.0_qp/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, 1.0_qp/2, a, info)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = rqsq_norm1(u(:,:,j)-v(:,:,j)) / rqsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'real    qp ',err,judge
end subroutine rq_test_blksp4triangular

! A test of the S & S for square matrices. We compare the two results,
! S^{-1} phi(D) S and phi(S^{-1}D).

subroutine cq_test_sands4square
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  complex(kind=qp) :: a(n), s(n,n), t(n,n)
  complex(kind=qp) :: b(n,  0:upto), c(n,n)
  complex(kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: pi, argument, rpart, ipart
  real   (kind=qp) :: err(0:upto)
  real   (kind=qp) :: radii, eps
  integer          :: i, j, perm(n)
  character        :: judge*4

  pi = acos(-1.0_qp)

  ! Generating a random diagonal matrix

  do i=1,n                                          ! diagonal entries
    argument = rrandom(pi/4, pi)
    radii = radius * rrandom(1.0_qp, 0.5_qp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i) = cmplx(rpart,ipart,qp)
  end do

  ! Generating a similarity transform

  do i=1,n
    do j=1,n
      argument = rrandom(pi*2, 0.0_qp)
      rpart = cos(argument)
      ipart = sin(argument)
      s(i,j) = cmplx(rpart,ipart,qp)
    end do
  end do
  s = cqsq_eye(n) - s / cqsq_norm1(s) / (cmplx(1.0_qp,0.0_qp,qp)*2)
  t = s; call cqsq_ludcmp(t,perm)

  ! Computing phi-functions upto phi_5
  ! ``.sqtimes.'' means the matrix-matrix multiplication.

  b = sasmtrphif(upto, cmplx(1.0_qp,0.0_qp,qp)/2, a)
  do j=0,upto                                       ! the similatrity transform
    do i=1,n
      u(i,:,j) = b(i,j) * s(i,:)
    end do
    u(:,:,j) = cqsq_ainvbf(t,perm,u(:,:,j))
  end do

  do i=1,n                                          ! the similarity transform
    c(i,:) = a(i) * s(i,:)
  end do
  v = sasmtrphif(upto, cmplx(1.0_qp,0.0_qp,qp)/2, cqsq_ainvbf(t,perm,c))

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cqsq_norm1(u(:,:,j)-v(:,:,j)) / cqsq_norm1(u(:,:,j))
    if (200*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex qp ',err,judge
end subroutine cq_test_sands4square

! A test for the S & S specialized to the triangular matrices. We compare the
! results with those obtained by the non-specialized version.

subroutine cq_test_sands4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  complex(kind=qp) :: a(n,n)
  complex(kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: pi, argument
  real   (kind=qp) :: err(0:upto)
  real   (kind=qp) :: radii, eps, rpart, ipart
  integer          :: i, j, subd(n)
  character        :: judge*4

  pi = acos(-1.0_qp)

  ! Generating a random test matrix

  subd = 0                                          ! sub-diagonal pattern
  if (0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(i) = 1
  end do

  a = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,n                                          ! off-diagonal entries
    do j=i,n                                        ! uniform random number
      a(i,j) = crandom(10.0_qp,0.0_qp)
    end do
  end do

  do i=1,n
    argument = rrandom(pi/4,pi)
    radii = radius * rrandom(1.0_qp,0.5_qp)
    rpart = radii * cos(argument) + center
    ipart = radii * sin(argument)
    a(i,i) = cmplx(rpart,ipart,qp)
  end do

  ! Computing phi-functions upto phi_5

  u = sasmtrphif(upto,cmplx(1.0_qp,0.0_qp,qp)/2,a)
  v = sasqtrphif(upto,cmplx(1.0_qp,0.0_qp,qp)/2,a)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cqsq_norm1(u(:,:,j)-v(:,:,j)) / cqsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex qp ',err,judge
end subroutine cq_test_sands4triangular

! A test for the modified block Schur--Parlett algorithm. We compare the results
! with those obtained by the scaling and squaring method specialized to the
! quasi-upper triangular matrices.

subroutine cq_test_blksp4triangular
  integer,       parameter :: n = 64, upto = 5
  real(kind=qp), parameter :: center = 2.0_qp, radius = 20.0_qp

  complex(kind=qp) :: a(n,n)
  complex(kind=qp) :: u(n,n,0:upto), v(n,n,0:upto)
  real   (kind=qp) :: err(0:upto), eps
  integer          :: i, j, subd(n)
  type(bspblock)   :: info
  character        :: judge*4

  ! Generating a random diagonal matrix

  a = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,n
    do j=i+1,n
      a(i,j) = crandom(1.0_qp,0.0_qp)
    end do
  end do

  ! A sub-diagonal pattern

  subd = 0
  if (0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(1) = 1
  do i=2,n-1
    if (subd(i-1) == 0 .and. 0.5_qp < rrandom(1.0_qp,0.5_qp)) subd(i) = 1
  end do
  subd = 0

  i = 1
  do while (i <= n)
    a(i,i) = - (i-1) * cmplx(1.0_qp,0.0_qp,qp) / (cmplx(1.0_qp,0.0_qp,qp) * n) * radius + center
    if (subd(i) == 1) then
      a(i  ,i+1) =  crandom(1.0_qp,0.0_qp)
      a(i+1,i  ) = -a(i  ,i+1)
      a(i+1,i+1) =  a(i  ,i  )
      i = i + 2
    else
      i = i + 1
    end if
  end do

  ! Computing phi-functions upto phi_5
  u = sasmtrphif(upto, cmplx(1.0_qp,0.0_qp,qp)/2, a)

  info = bspqtr1dcmpf(a)
  v = bspqtrphif(upto, cmplx(1.0_qp,0.0_qp,qp)/2, a, info)

  eps = epsilon(1.0_qp)/2                           ! the unit roundoff
  write(judge,'(a4)') '[OK]'
  do j=0,upto
    err(j) = cqsq_norm1(u(:,:,j)-v(:,:,j)) / cqsq_norm1(u(:,:,j))
    if (100*eps < err(j)) write(judge,'(a4)') '[NG]'
  end do

  write(*,'(a11,6e10.3,a6)') 'complex qp ',err,judge
end subroutine cq_test_blksp4triangular

#endif

subroutine sp_realref(values)
  real(kind=sp), intent(out) :: values(16,0:5)

  values( 1, 0) =  4.539992976248485153559151556055061023791e-5_sp
  values( 2, 0) =  1.234098040866795494976366907300338260721e-4_sp
  values( 3, 0) =  3.354626279025118388213891257808610193109e-4_sp
  values( 4, 0) =  9.118819655545162080031360844092826264737e-4_sp
  values( 5, 0) =  2.478752176666358423045167430816667891506e-3_sp
  values( 6, 0) =  6.737946999085467096636048423148424248849e-3_sp
  values( 7, 0) =  1.831563888873418029371802127324124221191e-2_sp
  values( 8, 0) =  4.978706836786394297934241565006177663169e-2_sp
  values( 9, 0) =  1.353352832366126918939994949724844034076e-1_sp
  values(10, 0) =  3.678794411714423215955237701614608674458e-1_sp
  values(11, 0) =  1.000000000000000000000000000000000000000e+0_sp
  values(12, 0) =  2.718281828459045235360287471352662497757e+0_sp
  values(13, 0) =  7.389056098930650227230427460575007813180e+0_sp
  values(14, 0) =  2.008553692318766774092852965458171789698e+1_sp
  values(15, 0) =  5.459815003314423907811026120286087840279e+1_sp
  values(16, 0) =  1.484131591025766034211155800405522796234e+2_sp
  values( 1, 1) =  9.999546000702375151484644084844394493897e-2_sp
  values( 2, 1) =  1.110973989106570356056113737010299962415e-1_sp
  values( 3, 1) =  1.249580671715121860201473263592773923725e-1_sp
  values( 4, 1) =  1.427268740049207833988566948450843881962e-1_sp
  values( 5, 1) =  1.662535413038889402628258054281972220180e-1_sp
  values( 6, 1) =  1.986524106001829065806727903153703151502e-1_sp
  values( 7, 1) =  2.454210902778164549265704946816896894470e-1_sp
  values( 8, 1) =  3.167376438773786856735525281166460744561e-1_sp
  values( 9, 1) =  4.323323583816936540530002525137577982961e-1_sp
  values(10, 1) =  6.321205588285576784044762298385391325541e-1_sp
  values(11, 1) =  1.000000000000000000000000000000000000000e+0_sp
  values(12, 1) =  1.718281828459045235360287471352662497757e+0_sp
  values(13, 1) =  3.194528049465325113615213730287503906590e+0_sp
  values(14, 1) =  6.361845641062555913642843218193905965662e+0_sp
  values(15, 1) =  1.339953750828605976952756530071521960069e+1_sp
  values(16, 1) =  2.948263182051532068422311600811045592469e+1_sp
  values( 1, 2) =  9.000045399929762484851535591515560550610e-2_sp
  values( 2, 2) =  9.876695567659366271048762514433000041760e-2_sp
  values( 3, 2) =  1.093802416035609767474815842050903259534e-1_sp
  values( 4, 2) =  1.224675894278684595144490435935593731148e-1_sp
  values( 5, 2) =  1.389577431160185099561956990953004629969e-1_sp
  values( 6, 2) =  1.602695178799634186838654419369259369699e-1_sp
  values( 7, 2) =  1.886447274305458862683573763295775776382e-1_sp
  values( 8, 2) =  2.277541187075404381088158239611179751812e-1_sp
  values( 9, 2) =  2.838338208091531729734998737431211008519e-1_sp
  values(10, 2) =  3.678794411714423215955237701614608674458e-1_sp
  values(11, 2) =  5.000000000000000000000000000000000000000e-1_sp
  values(12, 2) =  7.182818284590452353602874713526624977572e-1_sp
  values(13, 2) =  1.097264024732662556807606865143751953295e+0_sp
  values(14, 2) =  1.787281880354185304547614406064635321887e+0_sp
  values(15, 2) =  3.099884377071514942381891325178804900174e+0_sp
  values(16, 2) =  5.696526364103064136844623201622091184939e+0_sp
  values( 1, 3) =  4.099995460007023751514846440848443944938e-2_sp
  values( 2, 3) =  4.458144936926737080994581942840777773137e-2_sp
  values( 3, 3) =  4.882746979955487790656480197436370925582e-2_sp
  values( 4, 3) =  5.393320151030450578365013662949151812645e-2_sp
  values( 5, 3) =  6.017370948066358167396738348411658950050e-2_sp
  values( 6, 3) =  6.794609642400731626322691161261481260600e-2_sp
  values( 7, 3) =  7.783881814236352843291065591760560559043e-2_sp
  values( 8, 3) =  9.074862709748652063039472534629400827290e-2_sp
  values( 9, 3) =  1.080830895954234135132500631284394495740e-1_sp
  values(10, 3) =  1.321205588285576784044762298385391325541e-1_sp
  values(11, 3) =  1.666666666666666666666666666666666666666e-1_sp
  values(12, 3) =  2.182818284590452353602874713526624977572e-1_sp
  values(13, 3) =  2.986320123663312784038034325718759766475e-1_sp
  values(14, 3) =  4.290939601180617681825381353548784406291e-1_sp
  values(15, 3) =  6.499710942678787355954728312947012250436e-1_sp
  values(16, 3) =  1.039305272820612827368924640324418236987e+0_sp
  values( 1, 4) =  1.256667120665964291515182022581822272172e-2_sp
  values( 2, 4) =  1.356502414415547731741342747091765432614e-2_sp
  values( 3, 4) =  1.472989960838897359501273308653786967635e-2_sp
  values( 4, 4) =  1.610478073662316584043093286245359264860e-2_sp
  values( 5, 4) =  1.774882619766718083211654719709167952769e-2_sp
  values( 6, 4) =  1.974411404853187008068795101081037081213e-2_sp
  values( 7, 4) =  2.220696213107578455843900268726526526905e-2_sp
  values( 8, 4) =  2.530601318972671534542398044012421946458e-2_sp
  values( 9, 4) =  2.929178853562162657670830176911360854631e-2_sp
  values(10, 4) =  3.454610783810898826219043682812753411247e-2_sp
  values(11, 4) =  4.166666666666666666666666666666666666666e-2_sp
  values(12, 4) =  5.161516179237856869362080468599583109058e-2_sp
  values(13, 4) =  6.598267284983230586856838295260465499043e-2_sp
  values(14, 4) =  8.747576448379836717195715622940392465417e-2_sp
  values(15, 4) =  1.208261069003030172322015411570086395942e-1_sp
  values(16, 4) =  1.745277212307892321404515947315503140642e-1_sp
  values( 1, 5) =  2.909999546000702375151484644084844394493e-3_sp
  values( 2, 5) =  3.122404724723465483250359910638779148947e-3_sp
  values( 3, 5) =  3.367095882284711633956741697516099623788e-3_sp
  values( 4, 5) =  3.651697990006214403747961972030439145437e-3_sp
  values( 5, 5) =  3.986306744833247639091686578262497856495e-3_sp
  values( 6, 5) =  4.384510523626959317195743131171259170907e-3_sp
  values( 7, 5) =  4.864926133897720527056915994850350349402e-3_sp
  values( 8, 5) =  5.453551158979983773747562075514149067359e-3_sp
  values( 9, 5) =  6.187439065522520044979182448776529060178e-3_sp
  values(10, 5) =  7.120558828557678404476229838539132554188e-3_sp
  values(11, 5) =  8.333333333333333333333333333333333333333e-3_sp
  values(12, 5) =  9.948495125711902026954138019329164423913e-3_sp
  values(13, 5) =  1.215800309158281960095085814296899416188e-2_sp
  values(14, 5) =  1.526969927237723350176349652091241932916e-2_sp
  values(15, 5) =  1.978986005840908764138371862258549323189e-2_sp
  values(16, 5) =  2.657221091282451309475698561297672947951e-2_sp
end subroutine sp_realref

subroutine wp_realref(values)
  real(kind=wp), intent(out) :: values(16,0:5)

  values( 1, 0) =  4.539992976248485153559151556055061023791e-5_wp
  values( 2, 0) =  1.234098040866795494976366907300338260721e-4_wp
  values( 3, 0) =  3.354626279025118388213891257808610193109e-4_wp
  values( 4, 0) =  9.118819655545162080031360844092826264737e-4_wp
  values( 5, 0) =  2.478752176666358423045167430816667891506e-3_wp
  values( 6, 0) =  6.737946999085467096636048423148424248849e-3_wp
  values( 7, 0) =  1.831563888873418029371802127324124221191e-2_wp
  values( 8, 0) =  4.978706836786394297934241565006177663169e-2_wp
  values( 9, 0) =  1.353352832366126918939994949724844034076e-1_wp
  values(10, 0) =  3.678794411714423215955237701614608674458e-1_wp
  values(11, 0) =  1.000000000000000000000000000000000000000e+0_wp
  values(12, 0) =  2.718281828459045235360287471352662497757e+0_wp
  values(13, 0) =  7.389056098930650227230427460575007813180e+0_wp
  values(14, 0) =  2.008553692318766774092852965458171789698e+1_wp
  values(15, 0) =  5.459815003314423907811026120286087840279e+1_wp
  values(16, 0) =  1.484131591025766034211155800405522796234e+2_wp
  values( 1, 1) =  9.999546000702375151484644084844394493897e-2_wp
  values( 2, 1) =  1.110973989106570356056113737010299962415e-1_wp
  values( 3, 1) =  1.249580671715121860201473263592773923725e-1_wp
  values( 4, 1) =  1.427268740049207833988566948450843881962e-1_wp
  values( 5, 1) =  1.662535413038889402628258054281972220180e-1_wp
  values( 6, 1) =  1.986524106001829065806727903153703151502e-1_wp
  values( 7, 1) =  2.454210902778164549265704946816896894470e-1_wp
  values( 8, 1) =  3.167376438773786856735525281166460744561e-1_wp
  values( 9, 1) =  4.323323583816936540530002525137577982961e-1_wp
  values(10, 1) =  6.321205588285576784044762298385391325541e-1_wp
  values(11, 1) =  1.000000000000000000000000000000000000000e+0_wp
  values(12, 1) =  1.718281828459045235360287471352662497757e+0_wp
  values(13, 1) =  3.194528049465325113615213730287503906590e+0_wp
  values(14, 1) =  6.361845641062555913642843218193905965662e+0_wp
  values(15, 1) =  1.339953750828605976952756530071521960069e+1_wp
  values(16, 1) =  2.948263182051532068422311600811045592469e+1_wp
  values( 1, 2) =  9.000045399929762484851535591515560550610e-2_wp
  values( 2, 2) =  9.876695567659366271048762514433000041760e-2_wp
  values( 3, 2) =  1.093802416035609767474815842050903259534e-1_wp
  values( 4, 2) =  1.224675894278684595144490435935593731148e-1_wp
  values( 5, 2) =  1.389577431160185099561956990953004629969e-1_wp
  values( 6, 2) =  1.602695178799634186838654419369259369699e-1_wp
  values( 7, 2) =  1.886447274305458862683573763295775776382e-1_wp
  values( 8, 2) =  2.277541187075404381088158239611179751812e-1_wp
  values( 9, 2) =  2.838338208091531729734998737431211008519e-1_wp
  values(10, 2) =  3.678794411714423215955237701614608674458e-1_wp
  values(11, 2) =  5.000000000000000000000000000000000000000e-1_wp
  values(12, 2) =  7.182818284590452353602874713526624977572e-1_wp
  values(13, 2) =  1.097264024732662556807606865143751953295e+0_wp
  values(14, 2) =  1.787281880354185304547614406064635321887e+0_wp
  values(15, 2) =  3.099884377071514942381891325178804900174e+0_wp
  values(16, 2) =  5.696526364103064136844623201622091184939e+0_wp
  values( 1, 3) =  4.099995460007023751514846440848443944938e-2_wp
  values( 2, 3) =  4.458144936926737080994581942840777773137e-2_wp
  values( 3, 3) =  4.882746979955487790656480197436370925582e-2_wp
  values( 4, 3) =  5.393320151030450578365013662949151812645e-2_wp
  values( 5, 3) =  6.017370948066358167396738348411658950050e-2_wp
  values( 6, 3) =  6.794609642400731626322691161261481260600e-2_wp
  values( 7, 3) =  7.783881814236352843291065591760560559043e-2_wp
  values( 8, 3) =  9.074862709748652063039472534629400827290e-2_wp
  values( 9, 3) =  1.080830895954234135132500631284394495740e-1_wp
  values(10, 3) =  1.321205588285576784044762298385391325541e-1_wp
  values(11, 3) =  1.666666666666666666666666666666666666666e-1_wp
  values(12, 3) =  2.182818284590452353602874713526624977572e-1_wp
  values(13, 3) =  2.986320123663312784038034325718759766475e-1_wp
  values(14, 3) =  4.290939601180617681825381353548784406291e-1_wp
  values(15, 3) =  6.499710942678787355954728312947012250436e-1_wp
  values(16, 3) =  1.039305272820612827368924640324418236987e+0_wp
  values( 1, 4) =  1.256667120665964291515182022581822272172e-2_wp
  values( 2, 4) =  1.356502414415547731741342747091765432614e-2_wp
  values( 3, 4) =  1.472989960838897359501273308653786967635e-2_wp
  values( 4, 4) =  1.610478073662316584043093286245359264860e-2_wp
  values( 5, 4) =  1.774882619766718083211654719709167952769e-2_wp
  values( 6, 4) =  1.974411404853187008068795101081037081213e-2_wp
  values( 7, 4) =  2.220696213107578455843900268726526526905e-2_wp
  values( 8, 4) =  2.530601318972671534542398044012421946458e-2_wp
  values( 9, 4) =  2.929178853562162657670830176911360854631e-2_wp
  values(10, 4) =  3.454610783810898826219043682812753411247e-2_wp
  values(11, 4) =  4.166666666666666666666666666666666666666e-2_wp
  values(12, 4) =  5.161516179237856869362080468599583109058e-2_wp
  values(13, 4) =  6.598267284983230586856838295260465499043e-2_wp
  values(14, 4) =  8.747576448379836717195715622940392465417e-2_wp
  values(15, 4) =  1.208261069003030172322015411570086395942e-1_wp
  values(16, 4) =  1.745277212307892321404515947315503140642e-1_wp
  values( 1, 5) =  2.909999546000702375151484644084844394493e-3_wp
  values( 2, 5) =  3.122404724723465483250359910638779148947e-3_wp
  values( 3, 5) =  3.367095882284711633956741697516099623788e-3_wp
  values( 4, 5) =  3.651697990006214403747961972030439145437e-3_wp
  values( 5, 5) =  3.986306744833247639091686578262497856495e-3_wp
  values( 6, 5) =  4.384510523626959317195743131171259170907e-3_wp
  values( 7, 5) =  4.864926133897720527056915994850350349402e-3_wp
  values( 8, 5) =  5.453551158979983773747562075514149067359e-3_wp
  values( 9, 5) =  6.187439065522520044979182448776529060178e-3_wp
  values(10, 5) =  7.120558828557678404476229838539132554188e-3_wp
  values(11, 5) =  8.333333333333333333333333333333333333333e-3_wp
  values(12, 5) =  9.948495125711902026954138019329164423913e-3_wp
  values(13, 5) =  1.215800309158281960095085814296899416188e-2_wp
  values(14, 5) =  1.526969927237723350176349652091241932916e-2_wp
  values(15, 5) =  1.978986005840908764138371862258549323189e-2_wp
  values(16, 5) =  2.657221091282451309475698561297672947951e-2_wp
end subroutine wp_realref

#ifdef __USE_TPREC

subroutine tp_realref(values)
  real(kind=tp), intent(out) :: values(16,0:5)

  values( 1, 0) =  4.539992976248485153559151556055061023791e-5_tp
  values( 2, 0) =  1.234098040866795494976366907300338260721e-4_tp
  values( 3, 0) =  3.354626279025118388213891257808610193109e-4_tp
  values( 4, 0) =  9.118819655545162080031360844092826264737e-4_tp
  values( 5, 0) =  2.478752176666358423045167430816667891506e-3_tp
  values( 6, 0) =  6.737946999085467096636048423148424248849e-3_tp
  values( 7, 0) =  1.831563888873418029371802127324124221191e-2_tp
  values( 8, 0) =  4.978706836786394297934241565006177663169e-2_tp
  values( 9, 0) =  1.353352832366126918939994949724844034076e-1_tp
  values(10, 0) =  3.678794411714423215955237701614608674458e-1_tp
  values(11, 0) =  1.000000000000000000000000000000000000000e+0_tp
  values(12, 0) =  2.718281828459045235360287471352662497757e+0_tp
  values(13, 0) =  7.389056098930650227230427460575007813180e+0_tp
  values(14, 0) =  2.008553692318766774092852965458171789698e+1_tp
  values(15, 0) =  5.459815003314423907811026120286087840279e+1_tp
  values(16, 0) =  1.484131591025766034211155800405522796234e+2_tp
  values( 1, 1) =  9.999546000702375151484644084844394493897e-2_tp
  values( 2, 1) =  1.110973989106570356056113737010299962415e-1_tp
  values( 3, 1) =  1.249580671715121860201473263592773923725e-1_tp
  values( 4, 1) =  1.427268740049207833988566948450843881962e-1_tp
  values( 5, 1) =  1.662535413038889402628258054281972220180e-1_tp
  values( 6, 1) =  1.986524106001829065806727903153703151502e-1_tp
  values( 7, 1) =  2.454210902778164549265704946816896894470e-1_tp
  values( 8, 1) =  3.167376438773786856735525281166460744561e-1_tp
  values( 9, 1) =  4.323323583816936540530002525137577982961e-1_tp
  values(10, 1) =  6.321205588285576784044762298385391325541e-1_tp
  values(11, 1) =  1.000000000000000000000000000000000000000e+0_tp
  values(12, 1) =  1.718281828459045235360287471352662497757e+0_tp
  values(13, 1) =  3.194528049465325113615213730287503906590e+0_tp
  values(14, 1) =  6.361845641062555913642843218193905965662e+0_tp
  values(15, 1) =  1.339953750828605976952756530071521960069e+1_tp
  values(16, 1) =  2.948263182051532068422311600811045592469e+1_tp
  values( 1, 2) =  9.000045399929762484851535591515560550610e-2_tp
  values( 2, 2) =  9.876695567659366271048762514433000041760e-2_tp
  values( 3, 2) =  1.093802416035609767474815842050903259534e-1_tp
  values( 4, 2) =  1.224675894278684595144490435935593731148e-1_tp
  values( 5, 2) =  1.389577431160185099561956990953004629969e-1_tp
  values( 6, 2) =  1.602695178799634186838654419369259369699e-1_tp
  values( 7, 2) =  1.886447274305458862683573763295775776382e-1_tp
  values( 8, 2) =  2.277541187075404381088158239611179751812e-1_tp
  values( 9, 2) =  2.838338208091531729734998737431211008519e-1_tp
  values(10, 2) =  3.678794411714423215955237701614608674458e-1_tp
  values(11, 2) =  5.000000000000000000000000000000000000000e-1_tp
  values(12, 2) =  7.182818284590452353602874713526624977572e-1_tp
  values(13, 2) =  1.097264024732662556807606865143751953295e+0_tp
  values(14, 2) =  1.787281880354185304547614406064635321887e+0_tp
  values(15, 2) =  3.099884377071514942381891325178804900174e+0_tp
  values(16, 2) =  5.696526364103064136844623201622091184939e+0_tp
  values( 1, 3) =  4.099995460007023751514846440848443944938e-2_tp
  values( 2, 3) =  4.458144936926737080994581942840777773137e-2_tp
  values( 3, 3) =  4.882746979955487790656480197436370925582e-2_tp
  values( 4, 3) =  5.393320151030450578365013662949151812645e-2_tp
  values( 5, 3) =  6.017370948066358167396738348411658950050e-2_tp
  values( 6, 3) =  6.794609642400731626322691161261481260600e-2_tp
  values( 7, 3) =  7.783881814236352843291065591760560559043e-2_tp
  values( 8, 3) =  9.074862709748652063039472534629400827290e-2_tp
  values( 9, 3) =  1.080830895954234135132500631284394495740e-1_tp
  values(10, 3) =  1.321205588285576784044762298385391325541e-1_tp
  values(11, 3) =  1.666666666666666666666666666666666666666e-1_tp
  values(12, 3) =  2.182818284590452353602874713526624977572e-1_tp
  values(13, 3) =  2.986320123663312784038034325718759766475e-1_tp
  values(14, 3) =  4.290939601180617681825381353548784406291e-1_tp
  values(15, 3) =  6.499710942678787355954728312947012250436e-1_tp
  values(16, 3) =  1.039305272820612827368924640324418236987e+0_tp
  values( 1, 4) =  1.256667120665964291515182022581822272172e-2_tp
  values( 2, 4) =  1.356502414415547731741342747091765432614e-2_tp
  values( 3, 4) =  1.472989960838897359501273308653786967635e-2_tp
  values( 4, 4) =  1.610478073662316584043093286245359264860e-2_tp
  values( 5, 4) =  1.774882619766718083211654719709167952769e-2_tp
  values( 6, 4) =  1.974411404853187008068795101081037081213e-2_tp
  values( 7, 4) =  2.220696213107578455843900268726526526905e-2_tp
  values( 8, 4) =  2.530601318972671534542398044012421946458e-2_tp
  values( 9, 4) =  2.929178853562162657670830176911360854631e-2_tp
  values(10, 4) =  3.454610783810898826219043682812753411247e-2_tp
  values(11, 4) =  4.166666666666666666666666666666666666666e-2_tp
  values(12, 4) =  5.161516179237856869362080468599583109058e-2_tp
  values(13, 4) =  6.598267284983230586856838295260465499043e-2_tp
  values(14, 4) =  8.747576448379836717195715622940392465417e-2_tp
  values(15, 4) =  1.208261069003030172322015411570086395942e-1_tp
  values(16, 4) =  1.745277212307892321404515947315503140642e-1_tp
  values( 1, 5) =  2.909999546000702375151484644084844394493e-3_tp
  values( 2, 5) =  3.122404724723465483250359910638779148947e-3_tp
  values( 3, 5) =  3.367095882284711633956741697516099623788e-3_tp
  values( 4, 5) =  3.651697990006214403747961972030439145437e-3_tp
  values( 5, 5) =  3.986306744833247639091686578262497856495e-3_tp
  values( 6, 5) =  4.384510523626959317195743131171259170907e-3_tp
  values( 7, 5) =  4.864926133897720527056915994850350349402e-3_tp
  values( 8, 5) =  5.453551158979983773747562075514149067359e-3_tp
  values( 9, 5) =  6.187439065522520044979182448776529060178e-3_tp
  values(10, 5) =  7.120558828557678404476229838539132554188e-3_tp
  values(11, 5) =  8.333333333333333333333333333333333333333e-3_tp
  values(12, 5) =  9.948495125711902026954138019329164423913e-3_tp
  values(13, 5) =  1.215800309158281960095085814296899416188e-2_tp
  values(14, 5) =  1.526969927237723350176349652091241932916e-2_tp
  values(15, 5) =  1.978986005840908764138371862258549323189e-2_tp
  values(16, 5) =  2.657221091282451309475698561297672947951e-2_tp
end subroutine tp_realref

#endif

#ifdef __USE_QPREC

subroutine qp_realref(values)
  real(kind=qp), intent(out) :: values(16,0:5)

  values( 1, 0) =  4.539992976248485153559151556055061023791e-5_qp
  values( 2, 0) =  1.234098040866795494976366907300338260721e-4_qp
  values( 3, 0) =  3.354626279025118388213891257808610193109e-4_qp
  values( 4, 0) =  9.118819655545162080031360844092826264737e-4_qp
  values( 5, 0) =  2.478752176666358423045167430816667891506e-3_qp
  values( 6, 0) =  6.737946999085467096636048423148424248849e-3_qp
  values( 7, 0) =  1.831563888873418029371802127324124221191e-2_qp
  values( 8, 0) =  4.978706836786394297934241565006177663169e-2_qp
  values( 9, 0) =  1.353352832366126918939994949724844034076e-1_qp
  values(10, 0) =  3.678794411714423215955237701614608674458e-1_qp
  values(11, 0) =  1.000000000000000000000000000000000000000e+0_qp
  values(12, 0) =  2.718281828459045235360287471352662497757e+0_qp
  values(13, 0) =  7.389056098930650227230427460575007813180e+0_qp
  values(14, 0) =  2.008553692318766774092852965458171789698e+1_qp
  values(15, 0) =  5.459815003314423907811026120286087840279e+1_qp
  values(16, 0) =  1.484131591025766034211155800405522796234e+2_qp
  values( 1, 1) =  9.999546000702375151484644084844394493897e-2_qp
  values( 2, 1) =  1.110973989106570356056113737010299962415e-1_qp
  values( 3, 1) =  1.249580671715121860201473263592773923725e-1_qp
  values( 4, 1) =  1.427268740049207833988566948450843881962e-1_qp
  values( 5, 1) =  1.662535413038889402628258054281972220180e-1_qp
  values( 6, 1) =  1.986524106001829065806727903153703151502e-1_qp
  values( 7, 1) =  2.454210902778164549265704946816896894470e-1_qp
  values( 8, 1) =  3.167376438773786856735525281166460744561e-1_qp
  values( 9, 1) =  4.323323583816936540530002525137577982961e-1_qp
  values(10, 1) =  6.321205588285576784044762298385391325541e-1_qp
  values(11, 1) =  1.000000000000000000000000000000000000000e+0_qp
  values(12, 1) =  1.718281828459045235360287471352662497757e+0_qp
  values(13, 1) =  3.194528049465325113615213730287503906590e+0_qp
  values(14, 1) =  6.361845641062555913642843218193905965662e+0_qp
  values(15, 1) =  1.339953750828605976952756530071521960069e+1_qp
  values(16, 1) =  2.948263182051532068422311600811045592469e+1_qp
  values( 1, 2) =  9.000045399929762484851535591515560550610e-2_qp
  values( 2, 2) =  9.876695567659366271048762514433000041760e-2_qp
  values( 3, 2) =  1.093802416035609767474815842050903259534e-1_qp
  values( 4, 2) =  1.224675894278684595144490435935593731148e-1_qp
  values( 5, 2) =  1.389577431160185099561956990953004629969e-1_qp
  values( 6, 2) =  1.602695178799634186838654419369259369699e-1_qp
  values( 7, 2) =  1.886447274305458862683573763295775776382e-1_qp
  values( 8, 2) =  2.277541187075404381088158239611179751812e-1_qp
  values( 9, 2) =  2.838338208091531729734998737431211008519e-1_qp
  values(10, 2) =  3.678794411714423215955237701614608674458e-1_qp
  values(11, 2) =  5.000000000000000000000000000000000000000e-1_qp
  values(12, 2) =  7.182818284590452353602874713526624977572e-1_qp
  values(13, 2) =  1.097264024732662556807606865143751953295e+0_qp
  values(14, 2) =  1.787281880354185304547614406064635321887e+0_qp
  values(15, 2) =  3.099884377071514942381891325178804900174e+0_qp
  values(16, 2) =  5.696526364103064136844623201622091184939e+0_qp
  values( 1, 3) =  4.099995460007023751514846440848443944938e-2_qp
  values( 2, 3) =  4.458144936926737080994581942840777773137e-2_qp
  values( 3, 3) =  4.882746979955487790656480197436370925582e-2_qp
  values( 4, 3) =  5.393320151030450578365013662949151812645e-2_qp
  values( 5, 3) =  6.017370948066358167396738348411658950050e-2_qp
  values( 6, 3) =  6.794609642400731626322691161261481260600e-2_qp
  values( 7, 3) =  7.783881814236352843291065591760560559043e-2_qp
  values( 8, 3) =  9.074862709748652063039472534629400827290e-2_qp
  values( 9, 3) =  1.080830895954234135132500631284394495740e-1_qp
  values(10, 3) =  1.321205588285576784044762298385391325541e-1_qp
  values(11, 3) =  1.666666666666666666666666666666666666666e-1_qp
  values(12, 3) =  2.182818284590452353602874713526624977572e-1_qp
  values(13, 3) =  2.986320123663312784038034325718759766475e-1_qp
  values(14, 3) =  4.290939601180617681825381353548784406291e-1_qp
  values(15, 3) =  6.499710942678787355954728312947012250436e-1_qp
  values(16, 3) =  1.039305272820612827368924640324418236987e+0_qp
  values( 1, 4) =  1.256667120665964291515182022581822272172e-2_qp
  values( 2, 4) =  1.356502414415547731741342747091765432614e-2_qp
  values( 3, 4) =  1.472989960838897359501273308653786967635e-2_qp
  values( 4, 4) =  1.610478073662316584043093286245359264860e-2_qp
  values( 5, 4) =  1.774882619766718083211654719709167952769e-2_qp
  values( 6, 4) =  1.974411404853187008068795101081037081213e-2_qp
  values( 7, 4) =  2.220696213107578455843900268726526526905e-2_qp
  values( 8, 4) =  2.530601318972671534542398044012421946458e-2_qp
  values( 9, 4) =  2.929178853562162657670830176911360854631e-2_qp
  values(10, 4) =  3.454610783810898826219043682812753411247e-2_qp
  values(11, 4) =  4.166666666666666666666666666666666666666e-2_qp
  values(12, 4) =  5.161516179237856869362080468599583109058e-2_qp
  values(13, 4) =  6.598267284983230586856838295260465499043e-2_qp
  values(14, 4) =  8.747576448379836717195715622940392465417e-2_qp
  values(15, 4) =  1.208261069003030172322015411570086395942e-1_qp
  values(16, 4) =  1.745277212307892321404515947315503140642e-1_qp
  values( 1, 5) =  2.909999546000702375151484644084844394493e-3_qp
  values( 2, 5) =  3.122404724723465483250359910638779148947e-3_qp
  values( 3, 5) =  3.367095882284711633956741697516099623788e-3_qp
  values( 4, 5) =  3.651697990006214403747961972030439145437e-3_qp
  values( 5, 5) =  3.986306744833247639091686578262497856495e-3_qp
  values( 6, 5) =  4.384510523626959317195743131171259170907e-3_qp
  values( 7, 5) =  4.864926133897720527056915994850350349402e-3_qp
  values( 8, 5) =  5.453551158979983773747562075514149067359e-3_qp
  values( 9, 5) =  6.187439065522520044979182448776529060178e-3_qp
  values(10, 5) =  7.120558828557678404476229838539132554188e-3_qp
  values(11, 5) =  8.333333333333333333333333333333333333333e-3_qp
  values(12, 5) =  9.948495125711902026954138019329164423913e-3_qp
  values(13, 5) =  1.215800309158281960095085814296899416188e-2_qp
  values(14, 5) =  1.526969927237723350176349652091241932916e-2_qp
  values(15, 5) =  1.978986005840908764138371862258549323189e-2_qp
  values(16, 5) =  2.657221091282451309475698561297672947951e-2_qp
end subroutine qp_realref

#endif

subroutine sp_compref(values)
  complex(kind=sp), intent(out) :: values(32,0:5)

  values( 1, 0) = cmplx(2.452968673692215000628885544182439448649e-5_sp,  &
                       -3.820272360744745896412858462619482830985e-5_sp,sp)
  values( 2, 0) = cmplx(6.667860171476833283424963180029090616739e-5_sp,  &
                       -1.038457693797678114167511646823987469663e-4_sp,sp)
  values( 3, 0) = cmplx(1.812512313883128927908472072600232218423e-4_sp,  &
                       -2.822820678673715782795912746517464184101e-4_sp,sp)
  values( 4, 0) = cmplx(4.926919286686766622234621841806266399811e-4_sp,  &
                       -7.673222155837191136432586925770111584867e-4_sp,sp)
  values( 5, 0) = cmplx(1.339275516728503886085557975188675163702e-3_sp,  &
                       -2.085798035194157686322542562196870032728e-3_sp,sp)
  values( 6, 0) = cmplx(3.640528300423190167962660984149051867360e-3_sp,  &
                       -5.669786896903858940476818879173754216394e-3_sp,sp)
  values( 7, 0) = cmplx(9.895981925031249713864719801586350136119e-3_sp,  &
                       -1.541207869308895788150537629843732635202e-2_sp,sp)
  values( 8, 0) = cmplx(2.690006784157160778122578819852981646620e-2_sp,  &
                       -4.189437345020454468781373408358609461698e-2_sp,sp)
  values( 9, 0) = cmplx(7.312196559805963236599139851015598529567e-2_sp,  &
                       -1.138807140643680892286189302968000854906e-1_sp,sp)
  values(10, 0) = cmplx(1.987661103464129406288031913435846982927e-1_sp,  &
                       -3.095598756531121984439128249151294316712e-1_sp,sp)
  values(11, 0) = cmplx(5.403023058681397174009366074429766037323e-1_sp,  &
                       -8.414709848078965066525023216302989996225e-1_sp,sp)
  values(12, 0) = cmplx(1.468693939915885157138967597326604261326e+0_sp,  &
                       -2.287355287178842391208171906700501808955e+0_sp,sp)
  values(13, 0) = cmplx(3.992324048441271426506695498488872542168e+0_sp,  &
                       -6.217676312367968204252850304087010991267e+0_sp,sp)
  values(14, 0) = cmplx(1.085226191419795717634004705578458510082e+1_sp,  &
                       -1.690139653515009430510735354587017631994e+1_sp,sp)
  values(15, 0) = cmplx(2.949950635904248131176182655246028008012e+1_sp,  &
                       -4.594275907707917015245512895344831034661e+1_sp,sp)
  values(16, 0) = cmplx(8.018797208429722826939072854263268713595e+1_sp,  &
                       -1.248853671484961641963775578906393869008e+2_sp,sp)
  values(17, 0) = cmplx(2.452968673692215000628885544182439448649e-5_sp,  &
                        3.820272360744745896412858462619482830985e-5_sp,sp)
  values(18, 0) = cmplx(6.667860171476833283424963180029090616739e-5_sp,  &
                        1.038457693797678114167511646823987469663e-4_sp,sp)
  values(19, 0) = cmplx(1.812512313883128927908472072600232218423e-4_sp,  &
                        2.822820678673715782795912746517464184101e-4_sp,sp)
  values(20, 0) = cmplx(4.926919286686766622234621841806266399811e-4_sp,  &
                        7.673222155837191136432586925770111584867e-4_sp,sp)
  values(21, 0) = cmplx(1.339275516728503886085557975188675163702e-3_sp,  &
                        2.085798035194157686322542562196870032728e-3_sp,sp)
  values(22, 0) = cmplx(3.640528300423190167962660984149051867360e-3_sp,  &
                        5.669786896903858940476818879173754216394e-3_sp,sp)
  values(23, 0) = cmplx(9.895981925031249713864719801586350136119e-3_sp,  &
                        1.541207869308895788150537629843732635202e-2_sp,sp)
  values(24, 0) = cmplx(2.690006784157160778122578819852981646620e-2_sp,  &
                        4.189437345020454468781373408358609461698e-2_sp,sp)
  values(25, 0) = cmplx(7.312196559805963236599139851015598529567e-2_sp,  &
                        1.138807140643680892286189302968000854906e-1_sp,sp)
  values(26, 0) = cmplx(1.987661103464129406288031913435846982927e-1_sp,  &
                        3.095598756531121984439128249151294316712e-1_sp,sp)
  values(27, 0) = cmplx(5.403023058681397174009366074429766037323e-1_sp,  &
                        8.414709848078965066525023216302989996225e-1_sp,sp)
  values(28, 0) = cmplx(1.468693939915885157138967597326604261326e+0_sp,  &
                        2.287355287178842391208171906700501808955e+0_sp,sp)
  values(29, 0) = cmplx(3.992324048441271426506695498488872542168e+0_sp,  &
                        6.217676312367968204252850304087010991267e+0_sp,sp)
  values(30, 0) = cmplx(1.085226191419795717634004705578458510082e+1_sp,  &
                        1.690139653515009430510735354587017631994e+1_sp,sp)
  values(31, 0) = cmplx(2.949950635904248131176182655246028008012e+1_sp,  &
                        4.594275907707917015245512895344831034661e+1_sp,sp)
  values(32, 0) = cmplx(8.018797208429722826939072854263268713595e+1_sp,  &
                        1.248853671484961641963775578906393869008e+2_sp,sp)
  values( 1, 1) = cmplx(9.900785055303206164315743802010106882062e-2_sp,  &
                       -9.896964782942461418419330943547487399231e-3_sp,sp)
  values( 2, 1) = cmplx(1.097500455896822786928769329814448753730e-1_sp,  &
                       -1.218291109114472343127335353519583073623e-2_sp,sp)
  values( 3, 1) = cmplx(1.230589580341040133605531202094857163175e-1_sp,  &
                       -1.534708449577958022278419111685424623739e-2_sp,sp)
  values( 4, 1) = cmplx(1.399463695742980596495615804680662524935e-1_sp,  &
                       -1.988272105124490579084547453935560590501e-2_sp,sp)
  values( 5, 1) = cmplx(1.620013552685087333613461944516504005148e-1_sp,  &
                       -2.665259287221909594583727531490892174702e-2_sp,sp)
  values( 6, 1) = cmplx(1.918256594382610733884870582291703267261e-1_sp,  &
                       -3.723117450827144288960204786999931450194e-2_sp,sp)
  values( 7, 1) = cmplx(2.338722441760567034721203821818877603416e-1_sp,  &
                       -5.461504137074193639765375147086260849739e-2_sp,sp)
  values( 8, 1) = cmplx(2.961194169925489721344136369487996645218e-1_sp,  &
                       -8.474168118078147581553330095507118996828e-2_sp,sp)
  values( 9, 1) = cmplx(3.935273565736497648993272266552976229798e-1_sp,  &
                       -1.398233212546408378353541481792487687445e-1_sp,sp)
  values(10, 1) = cmplx(5.553968826533496289075548167857723666892e-1_sp,  &
                       -2.458370070002374304636419918706429350179e-1_sp,sp)
  values(11, 1) = cmplx(8.414709848078965066525023216302989996225e-1_sp,  &
                       -4.596976941318602825990633925570233962676e-1_sp,sp)
  values(12, 1) = cmplx(1.378024613547363774173569752013553035141e+0_sp,  &
                       -9.093306736314786170346021546869487738143e-1_sp,sp)
  values(13, 1) = cmplx(2.440464881850102211453248260212951215120e+0_sp,  &
                       -1.888605715258932996399801021937029888073e+0_sp,sp)
  values(14, 1) = cmplx(4.645818227774396583412749471322393162240e+0_sp,  &
                       -4.085192769125232573898201358182594385900e+0_sp,sp)
  values(15, 1) = cmplx(9.408281441955829141147202068428790039242e+0_sp,  &
                       -9.133619408780835252826981721254880076844e+0_sp,sp)
  values(16, 1) = cmplx(2.003173952192239636705120002322318548387e+1_sp,  &
                       -2.097072552531475356586527157348324028340e+1_sp,sp)
  values(17, 1) = cmplx(9.900785055303206164315743802010106882062e-2_sp,  &
                        9.896964782942461418419330943547487399231e-3_sp,sp)
  values(18, 1) = cmplx(1.097500455896822786928769329814448753730e-1_sp,  &
                        1.218291109114472343127335353519583073623e-2_sp,sp)
  values(19, 1) = cmplx(1.230589580341040133605531202094857163175e-1_sp,  &
                        1.534708449577958022278419111685424623739e-2_sp,sp)
  values(20, 1) = cmplx(1.399463695742980596495615804680662524935e-1_sp,  &
                        1.988272105124490579084547453935560590501e-2_sp,sp)
  values(21, 1) = cmplx(1.620013552685087333613461944516504005148e-1_sp,  &
                        2.665259287221909594583727531490892174702e-2_sp,sp)
  values(22, 1) = cmplx(1.918256594382610733884870582291703267261e-1_sp,  &
                        3.723117450827144288960204786999931450194e-2_sp,sp)
  values(23, 1) = cmplx(2.338722441760567034721203821818877603416e-1_sp,  &
                        5.461504137074193639765375147086260849739e-2_sp,sp)
  values(24, 1) = cmplx(2.961194169925489721344136369487996645218e-1_sp,  &
                        8.474168118078147581553330095507118996828e-2_sp,sp)
  values(25, 1) = cmplx(3.935273565736497648993272266552976229798e-1_sp,  &
                        1.398233212546408378353541481792487687445e-1_sp,sp)
  values(26, 1) = cmplx(5.553968826533496289075548167857723666892e-1_sp,  &
                        2.458370070002374304636419918706429350179e-1_sp,sp)
  values(27, 1) = cmplx(8.414709848078965066525023216302989996225e-1_sp,  &
                        4.596976941318602825990633925570233962676e-1_sp,sp)
  values(28, 1) = cmplx(1.378024613547363774173569752013553035141e+0_sp,  &
                        9.093306736314786170346021546869487738143e-1_sp,sp)
  values(29, 1) = cmplx(2.440464881850102211453248260212951215120e+0_sp,  &
                        1.888605715258932996399801021937029888073e+0_sp,sp)
  values(30, 1) = cmplx(4.645818227774396583412749471322393162240e+0_sp,  &
                        4.085192769125232573898201358182594385900e+0_sp,sp)
  values(31, 1) = cmplx(9.408281441955829141147202068428790039242e+0_sp,  &
                        9.133619408780835252826981721254880076844e+0_sp,sp)
  values(32, 1) = cmplx(2.003173952192239636705120002322318548387e+1_sp,  &
                        2.097072552531475356586527157348324028340e+1_sp,sp)
  values( 1, 2) = cmplx(8.930513325992694896026579159151026533854e-2_sp,  &
                       -7.940816847698448754184646064796277793931e-3_sp,sp)
  values( 2, 2) = cmplx(9.785893293639029530726074337441697502900e-2_sp,  &
                       -9.519557982805063541776376648802349365863e-3_sp,sp)
  values( 3, 2) = cmplx(1.081673141572761149744362958375533617799e-1_sp,  &
                       -1.160252870768706684395651309008738944281e-2_sp,sp)
  values( 4, 2) = cmplx(1.208051626806231697648782882252578367690e-1_sp,  &
                       -1.441749166133975199629040195512889012342e-2_sp,sp)
  values( 5, 2) = cmplx(1.366120124665180188048043272595947707745e-1_sp,  &
                       -1.832656993238315380982784199078097483791e-2_sp,sp)
  values( 6, 2) = cmplx(1.568501106660371567671987214124672184950e-1_sp,  &
                       -2.392378723155314277551933470849358079861e-2_sp,sp)
  values( 7, 2) = cmplx(1.834780038039126542652454248672536215959e-1_sp,  &
                       -3.221574060829267946689791834909775327463e-2_sp,sp)
  values( 8, 2) = cmplx(2.196383430203134559412292390108672196402e-1_sp,  &
                       -4.496555394651066004189864601859867655733e-2_sp,sp)
  values( 9, 2) = cmplx(2.705537216214682616073399389737307045569e-1_sp,  &
                       -6.536520018341371188599289539724096790618e-2_sp,sp)
  values(10, 2) = cmplx(3.452200621734439007780435875424352841643e-1_sp,  &
                       -9.938305517320647031440159567179234914639e-2_sp,sp)
  values(11, 2) = cmplx(4.596976941318602825990633925570233962676e-1_sp,  &
                       -1.585290151921034933474976783697010003774e-1_sp,sp)
  values(12, 2) = cmplx(6.436776435894211956040859533502509044777e-1_sp,  &
                       -2.656530300420574214305162013366978693365e-1_sp,sp)
  values(13, 2) = cmplx(9.539070957918274838612595084725864636630e-1_sp,  &
                       -4.673493097335527562692707567322217122051e-1_sp,sp)
  values(14, 2) = cmplx(1.502264745244842232413644977214977387262e+0_sp,  &
                       -8.609760079601301138281854603225389995461e-1_sp,sp)
  values(15, 2) = cmplx(2.515690892741420695142105293821767072577e+0_sp,  &
                       -1.654482129009853639421219106858278251066e+0_sp,sp)
  values(16, 2) = cmplx(4.466516274420259053889279680369198757798e+0_sp,  &
                       -3.300841850178898902395198378622808305121e+0_sp,sp)
  values(17, 2) = cmplx(8.930513325992694896026579159151026533854e-2_sp,  &
                        7.940816847698448754184646064796277793931e-3_sp,sp)
  values(18, 2) = cmplx(9.785893293639029530726074337441697502900e-2_sp,  &
                        9.519557982805063541776376648802349365863e-3_sp,sp)
  values(19, 2) = cmplx(1.081673141572761149744362958375533617799e-1_sp,  &
                        1.160252870768706684395651309008738944281e-2_sp,sp)
  values(20, 2) = cmplx(1.208051626806231697648782882252578367690e-1_sp,  &
                        1.441749166133975199629040195512889012342e-2_sp,sp)
  values(21, 2) = cmplx(1.366120124665180188048043272595947707745e-1_sp,  &
                        1.832656993238315380982784199078097483791e-2_sp,sp)
  values(22, 2) = cmplx(1.568501106660371567671987214124672184950e-1_sp,  &
                        2.392378723155314277551933470849358079861e-2_sp,sp)
  values(23, 2) = cmplx(1.834780038039126542652454248672536215959e-1_sp,  &
                        3.221574060829267946689791834909775327463e-2_sp,sp)
  values(24, 2) = cmplx(2.196383430203134559412292390108672196402e-1_sp,  &
                        4.496555394651066004189864601859867655733e-2_sp,sp)
  values(25, 2) = cmplx(2.705537216214682616073399389737307045569e-1_sp,  &
                        6.536520018341371188599289539724096790618e-2_sp,sp)
  values(26, 2) = cmplx(3.452200621734439007780435875424352841643e-1_sp,  &
                        9.938305517320647031440159567179234914639e-2_sp,sp)
  values(27, 2) = cmplx(4.596976941318602825990633925570233962676e-1_sp,  &
                        1.585290151921034933474976783697010003774e-1_sp,sp)
  values(28, 2) = cmplx(6.436776435894211956040859533502509044777e-1_sp,  &
                        2.656530300420574214305162013366978693365e-1_sp,sp)
  values(29, 2) = cmplx(9.539070957918274838612595084725864636630e-1_sp,  &
                        4.673493097335527562692707567322217122051e-1_sp,sp)
  values(30, 2) = cmplx(1.502264745244842232413644977214977387262e+0_sp,  &
                        8.609760079601301138281854603225389995461e-1_sp,sp)
  values(31, 2) = cmplx(2.515690892741420695142105293821767072577e+0_sp,  &
                        1.654482129009853639421219106858278251066e+0_sp,sp)
  values(32, 2) = cmplx(4.466516274420259053889279680369198757798e+0_sp,  &
                        3.300841850178898902395198378622808305121e+0_sp,sp)
  values( 1, 3) = cmplx(4.074148004206365306090620524900686756840e-2_sp,  &
                       -3.280066319436520430672155918421058977446e-3_sp,sp)
  values( 2, 3) = cmplx(4.425352636043039519239548397901279968420e-2_sp,  &
                       -3.859329819736147961179900814467827813148e-3_sp,sp)
  values( 3, 3) = cmplx(4.840406177614581764689947917522554608004e-2_sp,  &
                       -4.600191633557343850367870760642269579653e-3_sp,sp)
  values( 4, 3) = cmplx(5.337562705793955127284284768756648065480e-2_sp,  &
                       -5.565447913799971325221777961776798647340e-3_sp,sp)
  values( 5, 3) = cmplx(5.942309446306148759408113184954627973488e-2_sp,  &
                       -6.849420755113055630708881643127550816161e-3_sp,sp)
  values( 6, 3) = cmplx(6.691050899620643688228945106331374955089e-2_sp,  &
                       -8.597344352930658821354023270964033750455e-3_sp,sp)
  values( 7, 3) = cmplx(7.637080737603776837681860111059313334652e-2_sp,  &
                       -1.103876669193627222748017069037384501797e-2_sp,sp)
  values( 8, 3) = cmplx(8.860505248855702922182109289859970176364e-2_sp,  &
                       -1.454649951401545639330748229333367506877e-2_sp,sp)
  values( 9, 3) = cmplx(1.048515513880954377342626034899559117584e-1_sp,  &
                       -1.974317560234086292413485404635747192613e-2_sp,sp)
  values(10, 3) = cmplx(1.270814964998812847681790040646785324910e-1_sp,  &
                       -2.769844132667481445377740839288618334462e-2_sp,sp)
  values(11, 3) = cmplx(1.585290151921034933474976783697010003774e-1_sp,  &
                       -4.030230586813971740093660744297660373231e-2_sp,sp)
  values(12, 3) = cmplx(2.046653368157393085173010773434743869071e-1_sp,  &
                       -6.098769322631811291321512399322348242936e-2_sp,sp)
  values(13, 3) = cmplx(2.750327002634415447983579547354789279062e-1_sp,  &
                       -9.615830473505560573545640099837139214946e-2_sp,sp)
  values(14, 3) = cmplx(3.867770243694656811069120391967471161333e-1_sp,  &
                       -1.580663278635548109070911403752639611376e-1_sp,sp)
  values(15, 3) = cmplx(5.716026882338550835288023695379615612573e-1_sp,  &
                       -2.707198601939996389731041843300791724523e-1_sp,sp)
  values(16, 3) = cmplx(8.897470470107766989169844915564923882352e-1_sp,  &
                       -4.822189606336244406956427774132631833772e-1_sp,sp)
  values(17, 3) = cmplx(4.074148004206365306090620524900686756840e-2_sp,  &
                        3.280066319436520430672155918421058977446e-3_sp,sp)
  values(18, 3) = cmplx(4.425352636043039519239548397901279968420e-2_sp,  &
                        3.859329819736147961179900814467827813148e-3_sp,sp)
  values(19, 3) = cmplx(4.840406177614581764689947917522554608004e-2_sp,  &
                        4.600191633557343850367870760642269579653e-3_sp,sp)
  values(20, 3) = cmplx(5.337562705793955127284284768756648065480e-2_sp,  &
                        5.565447913799971325221777961776798647340e-3_sp,sp)
  values(21, 3) = cmplx(5.942309446306148759408113184954627973488e-2_sp,  &
                        6.849420755113055630708881643127550816161e-3_sp,sp)
  values(22, 3) = cmplx(6.691050899620643688228945106331374955089e-2_sp,  &
                        8.597344352930658821354023270964033750455e-3_sp,sp)
  values(23, 3) = cmplx(7.637080737603776837681860111059313334652e-2_sp,  &
                        1.103876669193627222748017069037384501797e-2_sp,sp)
  values(24, 3) = cmplx(8.860505248855702922182109289859970176364e-2_sp,  &
                        1.454649951401545639330748229333367506877e-2_sp,sp)
  values(25, 3) = cmplx(1.048515513880954377342626034899559117584e-1_sp,  &
                        1.974317560234086292413485404635747192613e-2_sp,sp)
  values(26, 3) = cmplx(1.270814964998812847681790040646785324910e-1_sp,  &
                        2.769844132667481445377740839288618334462e-2_sp,sp)
  values(27, 3) = cmplx(1.585290151921034933474976783697010003774e-1_sp,  &
                        4.030230586813971740093660744297660373231e-2_sp,sp)
  values(28, 3) = cmplx(2.046653368157393085173010773434743869071e-1_sp,  &
                        6.098769322631811291321512399322348242936e-2_sp,sp)
  values(29, 3) = cmplx(2.750327002634415447983579547354789279062e-1_sp,  &
                        9.615830473505560573545640099837139214946e-2_sp,sp)
  values(30, 3) = cmplx(3.867770243694656811069120391967471161333e-1_sp,  &
                        1.580663278635548109070911403752639611376e-1_sp,sp)
  values(31, 3) = cmplx(5.716026882338550835288023695379615612573e-1_sp,  &
                        2.707198601939996389731041843300791724523e-1_sp,sp)
  values(32, 3) = cmplx(8.897470470107766989169844915564923882352e-1_sp,  &
                        4.822189606336244406956427774132631833772e-1_sp,sp)
  values( 1, 4) = cmplx(1.250031616401452135136907693163385197980e-2_sp,  &
                       -9.220249844578000920696921013212793002356e-4_sp,sp)
  values( 2, 4) = cmplx(1.348265356799832428328805542687015403238e-2_sp,  &
                       -1.069258194251352924678683845822480691025e-3_sp,sp)
  values( 3, 4) = cmplx(1.462616970396498670782315954911032668111e-2_sp,  &
                       -1.253247258800955357181911098558507137682e-3_sp,sp)
  values( 4, 4) = cmplx(1.597205450349779558163977021630956201460e-2_sp,  &
                       -1.486658084242546322345427464933251909609e-3_sp,sp)
  values( 5, 4) = cmplx(1.757596902639849000178978623096891547045e-2_sp,  &
                       -1.787758045214239061846817431306894109049e-3_sp,sp)
  values( 6, 4) = cmplx(1.951454356558583875935538851106648535881e-2_sp,  &
                       -2.183439842531035987600273048020490321672e-3_sp,sp)
  values( 7, 4) = cmplx(2.189542375614422737569837840674517519403e-2_sp,  &
                       -2.714164266051988787054551929092832544015e-3_sp,sp)
  values( 8, 4) = cmplx(2.487313420483443687278442035975345697778e-2_sp,  &
                       -3.442211563606326826492312688806593969669e-3_sp,sp)
  values( 9, 4) = cmplx(2.867468123189666415778859607995579634851e-2_sp,  &
                       -4.465752814777900616826871016799162211191e-3_sp,sp)
  values(10, 4) = cmplx(3.364180574673009817613253549743715876013e-2_sp,  &
                       -5.943364420055283722355127104550975415511e-3_sp,sp)
  values(11, 4) = cmplx(4.030230586813971740093660744297660373231e-2_sp,  &
                       -8.137651474563173319168988296965666289229e-3_sp,sp)
  values(12, 4) = cmplx(4.949318168769537738192476733501560133492e-2_sp,  &
                       -1.149451153862273553129035665820788109443e-2_sp,sp)
  values(13, 4) = cmplx(6.257807438572107239976779542719918292572e-2_sp,  &
                       -1.679011517466726666784430278558610461187e-2_sp,sp)
  values(14, 4) = cmplx(8.183974009719518542278272579655053095375e-2_sp,  &
                       -2.540886258878654182810280485957114339461e-2_sp,sp)
  values(15, 4) = cmplx(1.112037615566325474365674703420740441656e-1_sp,  &
                       -3.987902465934177288413417849700128207168e-2_sp,sp)
  values(16, 4) = cmplx(1.576008023982374846902781500716304535084e-1_sp,  &
                       -6.492363164707739120107292546832654597374e-2_sp,sp)
  values(17, 4) = cmplx(1.250031616401452135136907693163385197980e-2_sp,  &
                        9.220249844578000920696921013212793002356e-4_sp,sp)
  values(18, 4) = cmplx(1.348265356799832428328805542687015403238e-2_sp,  &
                        1.069258194251352924678683845822480691025e-3_sp,sp)
  values(19, 4) = cmplx(1.462616970396498670782315954911032668111e-2_sp,  &
                        1.253247258800955357181911098558507137682e-3_sp,sp)
  values(20, 4) = cmplx(1.597205450349779558163977021630956201460e-2_sp,  &
                        1.486658084242546322345427464933251909609e-3_sp,sp)
  values(21, 4) = cmplx(1.757596902639849000178978623096891547045e-2_sp,  &
                        1.787758045214239061846817431306894109049e-3_sp,sp)
  values(22, 4) = cmplx(1.951454356558583875935538851106648535881e-2_sp,  &
                        2.183439842531035987600273048020490321672e-3_sp,sp)
  values(23, 4) = cmplx(2.189542375614422737569837840674517519403e-2_sp,  &
                        2.714164266051988787054551929092832544015e-3_sp,sp)
  values(24, 4) = cmplx(2.487313420483443687278442035975345697778e-2_sp,  &
                        3.442211563606326826492312688806593969669e-3_sp,sp)
  values(25, 4) = cmplx(2.867468123189666415778859607995579634851e-2_sp,  &
                        4.465752814777900616826871016799162211191e-3_sp,sp)
  values(26, 4) = cmplx(3.364180574673009817613253549743715876013e-2_sp,  &
                        5.943364420055283722355127104550975415511e-3_sp,sp)
  values(27, 4) = cmplx(4.030230586813971740093660744297660373231e-2_sp,  &
                        8.137651474563173319168988296965666289229e-3_sp,sp)
  values(28, 4) = cmplx(4.949318168769537738192476733501560133492e-2_sp,  &
                        1.149451153862273553129035665820788109443e-2_sp,sp)
  values(29, 4) = cmplx(6.257807438572107239976779542719918292572e-2_sp,  &
                        1.679011517466726666784430278558610461187e-2_sp,sp)
  values(30, 4) = cmplx(8.183974009719518542278272579655053095375e-2_sp,  &
                        2.540886258878654182810280485957114339461e-2_sp,sp)
  values(31, 4) = cmplx(1.112037615566325474365674703420740441656e-1_sp,  &
                        3.987902465934177288413417849700128207168e-2_sp,sp)
  values(32, 4) = cmplx(1.576008023982374846902781500716304535084e-1_sp,  &
                        6.492363164707739120107292546832654597374e-2_sp,sp)
  values( 1, 5) = cmplx(2.896886435752269834109362271798509169988e-3_sp,  &
                       -1.974861451294469742039670170477229869753e-4_sp,sp)
  values( 2, 5) = cmplx(3.106407025393493102135197378097452370726e-3_sp,  &
                       -2.263498701269044641618348369194412977445e-4_sp,sp)
  values( 3, 5) = cmplx(3.347341891698683000429691815984757338801e-3_sp,  &
                       -2.617618291122159554059725896782812751398e-4_sp,sp)
  values( 4, 5) = cmplx(3.626978864528492878350674052348659689480e-3_sp,  &
                       -3.057601114694209365721780839164868256958e-4_sp,sp)
  values( 5, 5) = cmplx(3.954917402346575650029948649878200034765e-3_sp,  &
                       -3.611932261887227646971885364285509876193e-4_sp,sp)
  values( 6, 5) = cmplx(4.344002128766737520159871685616207571573e-3_sp,  &
                       -4.321124572471403065119197275191434499801e-4_sp,sp)
  values( 7, 5) = cmplx(4.811713876949514467701629704045811672620e-3_sp,  &
                       -5.243874027243814201617694437382447821514e-4_sp,sp)
  values( 8, 5) = cmplx(5.382280894910301620813905160954622303632e-3_sp,  &
                       -6.466897771013249314405308240493427779874e-4_sp,sp)
  values( 9, 5) = cmplx(6.089944736863581126916602438044180569499e-3_sp,  &
                       -8.120959610428402550448657106225091791540e-4_sp,sp)
  values(10, 5) = cmplx(6.984112669995926106444629136890241661021e-3_sp,  &
                       -1.040748249940642384089502032339266245509e-3_sp,sp)
  values(11, 5) = cmplx(8.137651474563173319168988296965666289229e-3_sp,  &
                       -1.364360798526949265730059223690062934356e-3_sp,sp)
  values(12, 5) = cmplx(9.660513279825723123274228663278407881348e-3_sp,  &
                       -1.833998258797012408016127994929473213087e-3_sp,sp)
  values(13, 5) = cmplx(1.172258612255521562680931206133022742599e-2_sp,  &
                       -2.533764526056025520517495362127938592938e-3_sp,sp)
  values(14, 5) = cmplx(1.459280828803720980964509822492227362558e-2_sp,  &
                       -3.605351433583110672819235544882956589676e-3_sp,sp)
  values(15, 5) = cmplx(1.870749436583560564492572901168416423925e-2_sp,  &
                       -5.292882573376541809802112371329279458106e-3_sp,sp)
  values(16, 5) = cmplx(2.479208885788198005073578240358251846856e-2_sp,  &
                       -8.026308557839082230067428612948805501036e-3_sp,sp)
  values(17, 5) = cmplx(2.896886435752269834109362271798509169988e-3_sp,  &
                        1.974861451294469742039670170477229869753e-4_sp,sp)
  values(18, 5) = cmplx(3.106407025393493102135197378097452370726e-3_sp,  &
                        2.263498701269044641618348369194412977445e-4_sp,sp)
  values(19, 5) = cmplx(3.347341891698683000429691815984757338801e-3_sp,  &
                        2.617618291122159554059725896782812751398e-4_sp,sp)
  values(20, 5) = cmplx(3.626978864528492878350674052348659689480e-3_sp,  &
                        3.057601114694209365721780839164868256958e-4_sp,sp)
  values(21, 5) = cmplx(3.954917402346575650029948649878200034765e-3_sp,  &
                        3.611932261887227646971885364285509876193e-4_sp,sp)
  values(22, 5) = cmplx(4.344002128766737520159871685616207571573e-3_sp,  &
                        4.321124572471403065119197275191434499801e-4_sp,sp)
  values(23, 5) = cmplx(4.811713876949514467701629704045811672620e-3_sp,  &
                        5.243874027243814201617694437382447821514e-4_sp,sp)
  values(24, 5) = cmplx(5.382280894910301620813905160954622303632e-3_sp,  &
                        6.466897771013249314405308240493427779874e-4_sp,sp)
  values(25, 5) = cmplx(6.089944736863581126916602438044180569499e-3_sp,  &
                        8.120959610428402550448657106225091791540e-4_sp,sp)
  values(26, 5) = cmplx(6.984112669995926106444629136890241661021e-3_sp,  &
                        1.040748249940642384089502032339266245509e-3_sp,sp)
  values(27, 5) = cmplx(8.137651474563173319168988296965666289229e-3_sp,  &
                        1.364360798526949265730059223690062934356e-3_sp,sp)
  values(28, 5) = cmplx(9.660513279825723123274228663278407881348e-3_sp,  &
                        1.833998258797012408016127994929473213087e-3_sp,sp)
  values(29, 5) = cmplx(1.172258612255521562680931206133022742599e-2_sp,  &
                        2.533764526056025520517495362127938592938e-3_sp,sp)
  values(30, 5) = cmplx(1.459280828803720980964509822492227362558e-2_sp,  &
                        3.605351433583110672819235544882956589676e-3_sp,sp)
  values(31, 5) = cmplx(1.870749436583560564492572901168416423925e-2_sp,  &
                        5.292882573376541809802112371329279458106e-3_sp,sp)
  values(32, 5) = cmplx(2.479208885788198005073578240358251846856e-2_sp,  &
                        8.026308557839082230067428612948805501036e-3_sp,sp)
end subroutine sp_compref

subroutine wp_compref(values)
  complex(kind=wp), intent(out) :: values(32,0:5)

  values( 1, 0) = cmplx(2.452968673692215000628885544182439448649e-5_wp,  &
                       -3.820272360744745896412858462619482830985e-5_wp,wp)
  values( 2, 0) = cmplx(6.667860171476833283424963180029090616739e-5_wp,  &
                       -1.038457693797678114167511646823987469663e-4_wp,wp)
  values( 3, 0) = cmplx(1.812512313883128927908472072600232218423e-4_wp,  &
                       -2.822820678673715782795912746517464184101e-4_wp,wp)
  values( 4, 0) = cmplx(4.926919286686766622234621841806266399811e-4_wp,  &
                       -7.673222155837191136432586925770111584867e-4_wp,wp)
  values( 5, 0) = cmplx(1.339275516728503886085557975188675163702e-3_wp,  &
                       -2.085798035194157686322542562196870032728e-3_wp,wp)
  values( 6, 0) = cmplx(3.640528300423190167962660984149051867360e-3_wp,  &
                       -5.669786896903858940476818879173754216394e-3_wp,wp)
  values( 7, 0) = cmplx(9.895981925031249713864719801586350136119e-3_wp,  &
                       -1.541207869308895788150537629843732635202e-2_wp,wp)
  values( 8, 0) = cmplx(2.690006784157160778122578819852981646620e-2_wp,  &
                       -4.189437345020454468781373408358609461698e-2_wp,wp)
  values( 9, 0) = cmplx(7.312196559805963236599139851015598529567e-2_wp,  &
                       -1.138807140643680892286189302968000854906e-1_wp,wp)
  values(10, 0) = cmplx(1.987661103464129406288031913435846982927e-1_wp,  &
                       -3.095598756531121984439128249151294316712e-1_wp,wp)
  values(11, 0) = cmplx(5.403023058681397174009366074429766037323e-1_wp,  &
                       -8.414709848078965066525023216302989996225e-1_wp,wp)
  values(12, 0) = cmplx(1.468693939915885157138967597326604261326e+0_wp,  &
                       -2.287355287178842391208171906700501808955e+0_wp,wp)
  values(13, 0) = cmplx(3.992324048441271426506695498488872542168e+0_wp,  &
                       -6.217676312367968204252850304087010991267e+0_wp,wp)
  values(14, 0) = cmplx(1.085226191419795717634004705578458510082e+1_wp,  &
                       -1.690139653515009430510735354587017631994e+1_wp,wp)
  values(15, 0) = cmplx(2.949950635904248131176182655246028008012e+1_wp,  &
                       -4.594275907707917015245512895344831034661e+1_wp,wp)
  values(16, 0) = cmplx(8.018797208429722826939072854263268713595e+1_wp,  &
                       -1.248853671484961641963775578906393869008e+2_wp,wp)
  values(17, 0) = cmplx(2.452968673692215000628885544182439448649e-5_wp,  &
                        3.820272360744745896412858462619482830985e-5_wp,wp)
  values(18, 0) = cmplx(6.667860171476833283424963180029090616739e-5_wp,  &
                        1.038457693797678114167511646823987469663e-4_wp,wp)
  values(19, 0) = cmplx(1.812512313883128927908472072600232218423e-4_wp,  &
                        2.822820678673715782795912746517464184101e-4_wp,wp)
  values(20, 0) = cmplx(4.926919286686766622234621841806266399811e-4_wp,  &
                        7.673222155837191136432586925770111584867e-4_wp,wp)
  values(21, 0) = cmplx(1.339275516728503886085557975188675163702e-3_wp,  &
                        2.085798035194157686322542562196870032728e-3_wp,wp)
  values(22, 0) = cmplx(3.640528300423190167962660984149051867360e-3_wp,  &
                        5.669786896903858940476818879173754216394e-3_wp,wp)
  values(23, 0) = cmplx(9.895981925031249713864719801586350136119e-3_wp,  &
                        1.541207869308895788150537629843732635202e-2_wp,wp)
  values(24, 0) = cmplx(2.690006784157160778122578819852981646620e-2_wp,  &
                        4.189437345020454468781373408358609461698e-2_wp,wp)
  values(25, 0) = cmplx(7.312196559805963236599139851015598529567e-2_wp,  &
                        1.138807140643680892286189302968000854906e-1_wp,wp)
  values(26, 0) = cmplx(1.987661103464129406288031913435846982927e-1_wp,  &
                        3.095598756531121984439128249151294316712e-1_wp,wp)
  values(27, 0) = cmplx(5.403023058681397174009366074429766037323e-1_wp,  &
                        8.414709848078965066525023216302989996225e-1_wp,wp)
  values(28, 0) = cmplx(1.468693939915885157138967597326604261326e+0_wp,  &
                        2.287355287178842391208171906700501808955e+0_wp,wp)
  values(29, 0) = cmplx(3.992324048441271426506695498488872542168e+0_wp,  &
                        6.217676312367968204252850304087010991267e+0_wp,wp)
  values(30, 0) = cmplx(1.085226191419795717634004705578458510082e+1_wp,  &
                        1.690139653515009430510735354587017631994e+1_wp,wp)
  values(31, 0) = cmplx(2.949950635904248131176182655246028008012e+1_wp,  &
                        4.594275907707917015245512895344831034661e+1_wp,wp)
  values(32, 0) = cmplx(8.018797208429722826939072854263268713595e+1_wp,  &
                        1.248853671484961641963775578906393869008e+2_wp,wp)
  values( 1, 1) = cmplx(9.900785055303206164315743802010106882062e-2_wp,  &
                       -9.896964782942461418419330943547487399231e-3_wp,wp)
  values( 2, 1) = cmplx(1.097500455896822786928769329814448753730e-1_wp,  &
                       -1.218291109114472343127335353519583073623e-2_wp,wp)
  values( 3, 1) = cmplx(1.230589580341040133605531202094857163175e-1_wp,  &
                       -1.534708449577958022278419111685424623739e-2_wp,wp)
  values( 4, 1) = cmplx(1.399463695742980596495615804680662524935e-1_wp,  &
                       -1.988272105124490579084547453935560590501e-2_wp,wp)
  values( 5, 1) = cmplx(1.620013552685087333613461944516504005148e-1_wp,  &
                       -2.665259287221909594583727531490892174702e-2_wp,wp)
  values( 6, 1) = cmplx(1.918256594382610733884870582291703267261e-1_wp,  &
                       -3.723117450827144288960204786999931450194e-2_wp,wp)
  values( 7, 1) = cmplx(2.338722441760567034721203821818877603416e-1_wp,  &
                       -5.461504137074193639765375147086260849739e-2_wp,wp)
  values( 8, 1) = cmplx(2.961194169925489721344136369487996645218e-1_wp,  &
                       -8.474168118078147581553330095507118996828e-2_wp,wp)
  values( 9, 1) = cmplx(3.935273565736497648993272266552976229798e-1_wp,  &
                       -1.398233212546408378353541481792487687445e-1_wp,wp)
  values(10, 1) = cmplx(5.553968826533496289075548167857723666892e-1_wp,  &
                       -2.458370070002374304636419918706429350179e-1_wp,wp)
  values(11, 1) = cmplx(8.414709848078965066525023216302989996225e-1_wp,  &
                       -4.596976941318602825990633925570233962676e-1_wp,wp)
  values(12, 1) = cmplx(1.378024613547363774173569752013553035141e+0_wp,  &
                       -9.093306736314786170346021546869487738143e-1_wp,wp)
  values(13, 1) = cmplx(2.440464881850102211453248260212951215120e+0_wp,  &
                       -1.888605715258932996399801021937029888073e+0_wp,wp)
  values(14, 1) = cmplx(4.645818227774396583412749471322393162240e+0_wp,  &
                       -4.085192769125232573898201358182594385900e+0_wp,wp)
  values(15, 1) = cmplx(9.408281441955829141147202068428790039242e+0_wp,  &
                       -9.133619408780835252826981721254880076844e+0_wp,wp)
  values(16, 1) = cmplx(2.003173952192239636705120002322318548387e+1_wp,  &
                       -2.097072552531475356586527157348324028340e+1_wp,wp)
  values(17, 1) = cmplx(9.900785055303206164315743802010106882062e-2_wp,  &
                        9.896964782942461418419330943547487399231e-3_wp,wp)
  values(18, 1) = cmplx(1.097500455896822786928769329814448753730e-1_wp,  &
                        1.218291109114472343127335353519583073623e-2_wp,wp)
  values(19, 1) = cmplx(1.230589580341040133605531202094857163175e-1_wp,  &
                        1.534708449577958022278419111685424623739e-2_wp,wp)
  values(20, 1) = cmplx(1.399463695742980596495615804680662524935e-1_wp,  &
                        1.988272105124490579084547453935560590501e-2_wp,wp)
  values(21, 1) = cmplx(1.620013552685087333613461944516504005148e-1_wp,  &
                        2.665259287221909594583727531490892174702e-2_wp,wp)
  values(22, 1) = cmplx(1.918256594382610733884870582291703267261e-1_wp,  &
                        3.723117450827144288960204786999931450194e-2_wp,wp)
  values(23, 1) = cmplx(2.338722441760567034721203821818877603416e-1_wp,  &
                        5.461504137074193639765375147086260849739e-2_wp,wp)
  values(24, 1) = cmplx(2.961194169925489721344136369487996645218e-1_wp,  &
                        8.474168118078147581553330095507118996828e-2_wp,wp)
  values(25, 1) = cmplx(3.935273565736497648993272266552976229798e-1_wp,  &
                        1.398233212546408378353541481792487687445e-1_wp,wp)
  values(26, 1) = cmplx(5.553968826533496289075548167857723666892e-1_wp,  &
                        2.458370070002374304636419918706429350179e-1_wp,wp)
  values(27, 1) = cmplx(8.414709848078965066525023216302989996225e-1_wp,  &
                        4.596976941318602825990633925570233962676e-1_wp,wp)
  values(28, 1) = cmplx(1.378024613547363774173569752013553035141e+0_wp,  &
                        9.093306736314786170346021546869487738143e-1_wp,wp)
  values(29, 1) = cmplx(2.440464881850102211453248260212951215120e+0_wp,  &
                        1.888605715258932996399801021937029888073e+0_wp,wp)
  values(30, 1) = cmplx(4.645818227774396583412749471322393162240e+0_wp,  &
                        4.085192769125232573898201358182594385900e+0_wp,wp)
  values(31, 1) = cmplx(9.408281441955829141147202068428790039242e+0_wp,  &
                        9.133619408780835252826981721254880076844e+0_wp,wp)
  values(32, 1) = cmplx(2.003173952192239636705120002322318548387e+1_wp,  &
                        2.097072552531475356586527157348324028340e+1_wp,wp)
  values( 1, 2) = cmplx(8.930513325992694896026579159151026533854e-2_wp,  &
                       -7.940816847698448754184646064796277793931e-3_wp,wp)
  values( 2, 2) = cmplx(9.785893293639029530726074337441697502900e-2_wp,  &
                       -9.519557982805063541776376648802349365863e-3_wp,wp)
  values( 3, 2) = cmplx(1.081673141572761149744362958375533617799e-1_wp,  &
                       -1.160252870768706684395651309008738944281e-2_wp,wp)
  values( 4, 2) = cmplx(1.208051626806231697648782882252578367690e-1_wp,  &
                       -1.441749166133975199629040195512889012342e-2_wp,wp)
  values( 5, 2) = cmplx(1.366120124665180188048043272595947707745e-1_wp,  &
                       -1.832656993238315380982784199078097483791e-2_wp,wp)
  values( 6, 2) = cmplx(1.568501106660371567671987214124672184950e-1_wp,  &
                       -2.392378723155314277551933470849358079861e-2_wp,wp)
  values( 7, 2) = cmplx(1.834780038039126542652454248672536215959e-1_wp,  &
                       -3.221574060829267946689791834909775327463e-2_wp,wp)
  values( 8, 2) = cmplx(2.196383430203134559412292390108672196402e-1_wp,  &
                       -4.496555394651066004189864601859867655733e-2_wp,wp)
  values( 9, 2) = cmplx(2.705537216214682616073399389737307045569e-1_wp,  &
                       -6.536520018341371188599289539724096790618e-2_wp,wp)
  values(10, 2) = cmplx(3.452200621734439007780435875424352841643e-1_wp,  &
                       -9.938305517320647031440159567179234914639e-2_wp,wp)
  values(11, 2) = cmplx(4.596976941318602825990633925570233962676e-1_wp,  &
                       -1.585290151921034933474976783697010003774e-1_wp,wp)
  values(12, 2) = cmplx(6.436776435894211956040859533502509044777e-1_wp,  &
                       -2.656530300420574214305162013366978693365e-1_wp,wp)
  values(13, 2) = cmplx(9.539070957918274838612595084725864636630e-1_wp,  &
                       -4.673493097335527562692707567322217122051e-1_wp,wp)
  values(14, 2) = cmplx(1.502264745244842232413644977214977387262e+0_wp,  &
                       -8.609760079601301138281854603225389995461e-1_wp,wp)
  values(15, 2) = cmplx(2.515690892741420695142105293821767072577e+0_wp,  &
                       -1.654482129009853639421219106858278251066e+0_wp,wp)
  values(16, 2) = cmplx(4.466516274420259053889279680369198757798e+0_wp,  &
                       -3.300841850178898902395198378622808305121e+0_wp,wp)
  values(17, 2) = cmplx(8.930513325992694896026579159151026533854e-2_wp,  &
                        7.940816847698448754184646064796277793931e-3_wp,wp)
  values(18, 2) = cmplx(9.785893293639029530726074337441697502900e-2_wp,  &
                        9.519557982805063541776376648802349365863e-3_wp,wp)
  values(19, 2) = cmplx(1.081673141572761149744362958375533617799e-1_wp,  &
                        1.160252870768706684395651309008738944281e-2_wp,wp)
  values(20, 2) = cmplx(1.208051626806231697648782882252578367690e-1_wp,  &
                        1.441749166133975199629040195512889012342e-2_wp,wp)
  values(21, 2) = cmplx(1.366120124665180188048043272595947707745e-1_wp,  &
                        1.832656993238315380982784199078097483791e-2_wp,wp)
  values(22, 2) = cmplx(1.568501106660371567671987214124672184950e-1_wp,  &
                        2.392378723155314277551933470849358079861e-2_wp,wp)
  values(23, 2) = cmplx(1.834780038039126542652454248672536215959e-1_wp,  &
                        3.221574060829267946689791834909775327463e-2_wp,wp)
  values(24, 2) = cmplx(2.196383430203134559412292390108672196402e-1_wp,  &
                        4.496555394651066004189864601859867655733e-2_wp,wp)
  values(25, 2) = cmplx(2.705537216214682616073399389737307045569e-1_wp,  &
                        6.536520018341371188599289539724096790618e-2_wp,wp)
  values(26, 2) = cmplx(3.452200621734439007780435875424352841643e-1_wp,  &
                        9.938305517320647031440159567179234914639e-2_wp,wp)
  values(27, 2) = cmplx(4.596976941318602825990633925570233962676e-1_wp,  &
                        1.585290151921034933474976783697010003774e-1_wp,wp)
  values(28, 2) = cmplx(6.436776435894211956040859533502509044777e-1_wp,  &
                        2.656530300420574214305162013366978693365e-1_wp,wp)
  values(29, 2) = cmplx(9.539070957918274838612595084725864636630e-1_wp,  &
                        4.673493097335527562692707567322217122051e-1_wp,wp)
  values(30, 2) = cmplx(1.502264745244842232413644977214977387262e+0_wp,  &
                        8.609760079601301138281854603225389995461e-1_wp,wp)
  values(31, 2) = cmplx(2.515690892741420695142105293821767072577e+0_wp,  &
                        1.654482129009853639421219106858278251066e+0_wp,wp)
  values(32, 2) = cmplx(4.466516274420259053889279680369198757798e+0_wp,  &
                        3.300841850178898902395198378622808305121e+0_wp,wp)
  values( 1, 3) = cmplx(4.074148004206365306090620524900686756840e-2_wp,  &
                       -3.280066319436520430672155918421058977446e-3_wp,wp)
  values( 2, 3) = cmplx(4.425352636043039519239548397901279968420e-2_wp,  &
                       -3.859329819736147961179900814467827813148e-3_wp,wp)
  values( 3, 3) = cmplx(4.840406177614581764689947917522554608004e-2_wp,  &
                       -4.600191633557343850367870760642269579653e-3_wp,wp)
  values( 4, 3) = cmplx(5.337562705793955127284284768756648065480e-2_wp,  &
                       -5.565447913799971325221777961776798647340e-3_wp,wp)
  values( 5, 3) = cmplx(5.942309446306148759408113184954627973488e-2_wp,  &
                       -6.849420755113055630708881643127550816161e-3_wp,wp)
  values( 6, 3) = cmplx(6.691050899620643688228945106331374955089e-2_wp,  &
                       -8.597344352930658821354023270964033750455e-3_wp,wp)
  values( 7, 3) = cmplx(7.637080737603776837681860111059313334652e-2_wp,  &
                       -1.103876669193627222748017069037384501797e-2_wp,wp)
  values( 8, 3) = cmplx(8.860505248855702922182109289859970176364e-2_wp,  &
                       -1.454649951401545639330748229333367506877e-2_wp,wp)
  values( 9, 3) = cmplx(1.048515513880954377342626034899559117584e-1_wp,  &
                       -1.974317560234086292413485404635747192613e-2_wp,wp)
  values(10, 3) = cmplx(1.270814964998812847681790040646785324910e-1_wp,  &
                       -2.769844132667481445377740839288618334462e-2_wp,wp)
  values(11, 3) = cmplx(1.585290151921034933474976783697010003774e-1_wp,  &
                       -4.030230586813971740093660744297660373231e-2_wp,wp)
  values(12, 3) = cmplx(2.046653368157393085173010773434743869071e-1_wp,  &
                       -6.098769322631811291321512399322348242936e-2_wp,wp)
  values(13, 3) = cmplx(2.750327002634415447983579547354789279062e-1_wp,  &
                       -9.615830473505560573545640099837139214946e-2_wp,wp)
  values(14, 3) = cmplx(3.867770243694656811069120391967471161333e-1_wp,  &
                       -1.580663278635548109070911403752639611376e-1_wp,wp)
  values(15, 3) = cmplx(5.716026882338550835288023695379615612573e-1_wp,  &
                       -2.707198601939996389731041843300791724523e-1_wp,wp)
  values(16, 3) = cmplx(8.897470470107766989169844915564923882352e-1_wp,  &
                       -4.822189606336244406956427774132631833772e-1_wp,wp)
  values(17, 3) = cmplx(4.074148004206365306090620524900686756840e-2_wp,  &
                        3.280066319436520430672155918421058977446e-3_wp,wp)
  values(18, 3) = cmplx(4.425352636043039519239548397901279968420e-2_wp,  &
                        3.859329819736147961179900814467827813148e-3_wp,wp)
  values(19, 3) = cmplx(4.840406177614581764689947917522554608004e-2_wp,  &
                        4.600191633557343850367870760642269579653e-3_wp,wp)
  values(20, 3) = cmplx(5.337562705793955127284284768756648065480e-2_wp,  &
                        5.565447913799971325221777961776798647340e-3_wp,wp)
  values(21, 3) = cmplx(5.942309446306148759408113184954627973488e-2_wp,  &
                        6.849420755113055630708881643127550816161e-3_wp,wp)
  values(22, 3) = cmplx(6.691050899620643688228945106331374955089e-2_wp,  &
                        8.597344352930658821354023270964033750455e-3_wp,wp)
  values(23, 3) = cmplx(7.637080737603776837681860111059313334652e-2_wp,  &
                        1.103876669193627222748017069037384501797e-2_wp,wp)
  values(24, 3) = cmplx(8.860505248855702922182109289859970176364e-2_wp,  &
                        1.454649951401545639330748229333367506877e-2_wp,wp)
  values(25, 3) = cmplx(1.048515513880954377342626034899559117584e-1_wp,  &
                        1.974317560234086292413485404635747192613e-2_wp,wp)
  values(26, 3) = cmplx(1.270814964998812847681790040646785324910e-1_wp,  &
                        2.769844132667481445377740839288618334462e-2_wp,wp)
  values(27, 3) = cmplx(1.585290151921034933474976783697010003774e-1_wp,  &
                        4.030230586813971740093660744297660373231e-2_wp,wp)
  values(28, 3) = cmplx(2.046653368157393085173010773434743869071e-1_wp,  &
                        6.098769322631811291321512399322348242936e-2_wp,wp)
  values(29, 3) = cmplx(2.750327002634415447983579547354789279062e-1_wp,  &
                        9.615830473505560573545640099837139214946e-2_wp,wp)
  values(30, 3) = cmplx(3.867770243694656811069120391967471161333e-1_wp,  &
                        1.580663278635548109070911403752639611376e-1_wp,wp)
  values(31, 3) = cmplx(5.716026882338550835288023695379615612573e-1_wp,  &
                        2.707198601939996389731041843300791724523e-1_wp,wp)
  values(32, 3) = cmplx(8.897470470107766989169844915564923882352e-1_wp,  &
                        4.822189606336244406956427774132631833772e-1_wp,wp)
  values( 1, 4) = cmplx(1.250031616401452135136907693163385197980e-2_wp,  &
                       -9.220249844578000920696921013212793002356e-4_wp,wp)
  values( 2, 4) = cmplx(1.348265356799832428328805542687015403238e-2_wp,  &
                       -1.069258194251352924678683845822480691025e-3_wp,wp)
  values( 3, 4) = cmplx(1.462616970396498670782315954911032668111e-2_wp,  &
                       -1.253247258800955357181911098558507137682e-3_wp,wp)
  values( 4, 4) = cmplx(1.597205450349779558163977021630956201460e-2_wp,  &
                       -1.486658084242546322345427464933251909609e-3_wp,wp)
  values( 5, 4) = cmplx(1.757596902639849000178978623096891547045e-2_wp,  &
                       -1.787758045214239061846817431306894109049e-3_wp,wp)
  values( 6, 4) = cmplx(1.951454356558583875935538851106648535881e-2_wp,  &
                       -2.183439842531035987600273048020490321672e-3_wp,wp)
  values( 7, 4) = cmplx(2.189542375614422737569837840674517519403e-2_wp,  &
                       -2.714164266051988787054551929092832544015e-3_wp,wp)
  values( 8, 4) = cmplx(2.487313420483443687278442035975345697778e-2_wp,  &
                       -3.442211563606326826492312688806593969669e-3_wp,wp)
  values( 9, 4) = cmplx(2.867468123189666415778859607995579634851e-2_wp,  &
                       -4.465752814777900616826871016799162211191e-3_wp,wp)
  values(10, 4) = cmplx(3.364180574673009817613253549743715876013e-2_wp,  &
                       -5.943364420055283722355127104550975415511e-3_wp,wp)
  values(11, 4) = cmplx(4.030230586813971740093660744297660373231e-2_wp,  &
                       -8.137651474563173319168988296965666289229e-3_wp,wp)
  values(12, 4) = cmplx(4.949318168769537738192476733501560133492e-2_wp,  &
                       -1.149451153862273553129035665820788109443e-2_wp,wp)
  values(13, 4) = cmplx(6.257807438572107239976779542719918292572e-2_wp,  &
                       -1.679011517466726666784430278558610461187e-2_wp,wp)
  values(14, 4) = cmplx(8.183974009719518542278272579655053095375e-2_wp,  &
                       -2.540886258878654182810280485957114339461e-2_wp,wp)
  values(15, 4) = cmplx(1.112037615566325474365674703420740441656e-1_wp,  &
                       -3.987902465934177288413417849700128207168e-2_wp,wp)
  values(16, 4) = cmplx(1.576008023982374846902781500716304535084e-1_wp,  &
                       -6.492363164707739120107292546832654597374e-2_wp,wp)
  values(17, 4) = cmplx(1.250031616401452135136907693163385197980e-2_wp,  &
                        9.220249844578000920696921013212793002356e-4_wp,wp)
  values(18, 4) = cmplx(1.348265356799832428328805542687015403238e-2_wp,  &
                        1.069258194251352924678683845822480691025e-3_wp,wp)
  values(19, 4) = cmplx(1.462616970396498670782315954911032668111e-2_wp,  &
                        1.253247258800955357181911098558507137682e-3_wp,wp)
  values(20, 4) = cmplx(1.597205450349779558163977021630956201460e-2_wp,  &
                        1.486658084242546322345427464933251909609e-3_wp,wp)
  values(21, 4) = cmplx(1.757596902639849000178978623096891547045e-2_wp,  &
                        1.787758045214239061846817431306894109049e-3_wp,wp)
  values(22, 4) = cmplx(1.951454356558583875935538851106648535881e-2_wp,  &
                        2.183439842531035987600273048020490321672e-3_wp,wp)
  values(23, 4) = cmplx(2.189542375614422737569837840674517519403e-2_wp,  &
                        2.714164266051988787054551929092832544015e-3_wp,wp)
  values(24, 4) = cmplx(2.487313420483443687278442035975345697778e-2_wp,  &
                        3.442211563606326826492312688806593969669e-3_wp,wp)
  values(25, 4) = cmplx(2.867468123189666415778859607995579634851e-2_wp,  &
                        4.465752814777900616826871016799162211191e-3_wp,wp)
  values(26, 4) = cmplx(3.364180574673009817613253549743715876013e-2_wp,  &
                        5.943364420055283722355127104550975415511e-3_wp,wp)
  values(27, 4) = cmplx(4.030230586813971740093660744297660373231e-2_wp,  &
                        8.137651474563173319168988296965666289229e-3_wp,wp)
  values(28, 4) = cmplx(4.949318168769537738192476733501560133492e-2_wp,  &
                        1.149451153862273553129035665820788109443e-2_wp,wp)
  values(29, 4) = cmplx(6.257807438572107239976779542719918292572e-2_wp,  &
                        1.679011517466726666784430278558610461187e-2_wp,wp)
  values(30, 4) = cmplx(8.183974009719518542278272579655053095375e-2_wp,  &
                        2.540886258878654182810280485957114339461e-2_wp,wp)
  values(31, 4) = cmplx(1.112037615566325474365674703420740441656e-1_wp,  &
                        3.987902465934177288413417849700128207168e-2_wp,wp)
  values(32, 4) = cmplx(1.576008023982374846902781500716304535084e-1_wp,  &
                        6.492363164707739120107292546832654597374e-2_wp,wp)
  values( 1, 5) = cmplx(2.896886435752269834109362271798509169988e-3_wp,  &
                       -1.974861451294469742039670170477229869753e-4_wp,wp)
  values( 2, 5) = cmplx(3.106407025393493102135197378097452370726e-3_wp,  &
                       -2.263498701269044641618348369194412977445e-4_wp,wp)
  values( 3, 5) = cmplx(3.347341891698683000429691815984757338801e-3_wp,  &
                       -2.617618291122159554059725896782812751398e-4_wp,wp)
  values( 4, 5) = cmplx(3.626978864528492878350674052348659689480e-3_wp,  &
                       -3.057601114694209365721780839164868256958e-4_wp,wp)
  values( 5, 5) = cmplx(3.954917402346575650029948649878200034765e-3_wp,  &
                       -3.611932261887227646971885364285509876193e-4_wp,wp)
  values( 6, 5) = cmplx(4.344002128766737520159871685616207571573e-3_wp,  &
                       -4.321124572471403065119197275191434499801e-4_wp,wp)
  values( 7, 5) = cmplx(4.811713876949514467701629704045811672620e-3_wp,  &
                       -5.243874027243814201617694437382447821514e-4_wp,wp)
  values( 8, 5) = cmplx(5.382280894910301620813905160954622303632e-3_wp,  &
                       -6.466897771013249314405308240493427779874e-4_wp,wp)
  values( 9, 5) = cmplx(6.089944736863581126916602438044180569499e-3_wp,  &
                       -8.120959610428402550448657106225091791540e-4_wp,wp)
  values(10, 5) = cmplx(6.984112669995926106444629136890241661021e-3_wp,  &
                       -1.040748249940642384089502032339266245509e-3_wp,wp)
  values(11, 5) = cmplx(8.137651474563173319168988296965666289229e-3_wp,  &
                       -1.364360798526949265730059223690062934356e-3_wp,wp)
  values(12, 5) = cmplx(9.660513279825723123274228663278407881348e-3_wp,  &
                       -1.833998258797012408016127994929473213087e-3_wp,wp)
  values(13, 5) = cmplx(1.172258612255521562680931206133022742599e-2_wp,  &
                       -2.533764526056025520517495362127938592938e-3_wp,wp)
  values(14, 5) = cmplx(1.459280828803720980964509822492227362558e-2_wp,  &
                       -3.605351433583110672819235544882956589676e-3_wp,wp)
  values(15, 5) = cmplx(1.870749436583560564492572901168416423925e-2_wp,  &
                       -5.292882573376541809802112371329279458106e-3_wp,wp)
  values(16, 5) = cmplx(2.479208885788198005073578240358251846856e-2_wp,  &
                       -8.026308557839082230067428612948805501036e-3_wp,wp)
  values(17, 5) = cmplx(2.896886435752269834109362271798509169988e-3_wp,  &
                        1.974861451294469742039670170477229869753e-4_wp,wp)
  values(18, 5) = cmplx(3.106407025393493102135197378097452370726e-3_wp,  &
                        2.263498701269044641618348369194412977445e-4_wp,wp)
  values(19, 5) = cmplx(3.347341891698683000429691815984757338801e-3_wp,  &
                        2.617618291122159554059725896782812751398e-4_wp,wp)
  values(20, 5) = cmplx(3.626978864528492878350674052348659689480e-3_wp,  &
                        3.057601114694209365721780839164868256958e-4_wp,wp)
  values(21, 5) = cmplx(3.954917402346575650029948649878200034765e-3_wp,  &
                        3.611932261887227646971885364285509876193e-4_wp,wp)
  values(22, 5) = cmplx(4.344002128766737520159871685616207571573e-3_wp,  &
                        4.321124572471403065119197275191434499801e-4_wp,wp)
  values(23, 5) = cmplx(4.811713876949514467701629704045811672620e-3_wp,  &
                        5.243874027243814201617694437382447821514e-4_wp,wp)
  values(24, 5) = cmplx(5.382280894910301620813905160954622303632e-3_wp,  &
                        6.466897771013249314405308240493427779874e-4_wp,wp)
  values(25, 5) = cmplx(6.089944736863581126916602438044180569499e-3_wp,  &
                        8.120959610428402550448657106225091791540e-4_wp,wp)
  values(26, 5) = cmplx(6.984112669995926106444629136890241661021e-3_wp,  &
                        1.040748249940642384089502032339266245509e-3_wp,wp)
  values(27, 5) = cmplx(8.137651474563173319168988296965666289229e-3_wp,  &
                        1.364360798526949265730059223690062934356e-3_wp,wp)
  values(28, 5) = cmplx(9.660513279825723123274228663278407881348e-3_wp,  &
                        1.833998258797012408016127994929473213087e-3_wp,wp)
  values(29, 5) = cmplx(1.172258612255521562680931206133022742599e-2_wp,  &
                        2.533764526056025520517495362127938592938e-3_wp,wp)
  values(30, 5) = cmplx(1.459280828803720980964509822492227362558e-2_wp,  &
                        3.605351433583110672819235544882956589676e-3_wp,wp)
  values(31, 5) = cmplx(1.870749436583560564492572901168416423925e-2_wp,  &
                        5.292882573376541809802112371329279458106e-3_wp,wp)
  values(32, 5) = cmplx(2.479208885788198005073578240358251846856e-2_wp,  &
                        8.026308557839082230067428612948805501036e-3_wp,wp)
end subroutine wp_compref

#ifdef __USE_TPREC

subroutine tp_compref(values)
  complex(kind=tp), intent(out) :: values(32,0:5)

  values( 1, 0) = cmplx(2.452968673692215000628885544182439448649e-5_tp,  &
                       -3.820272360744745896412858462619482830985e-5_tp,tp)
  values( 2, 0) = cmplx(6.667860171476833283424963180029090616739e-5_tp,  &
                       -1.038457693797678114167511646823987469663e-4_tp,tp)
  values( 3, 0) = cmplx(1.812512313883128927908472072600232218423e-4_tp,  &
                       -2.822820678673715782795912746517464184101e-4_tp,tp)
  values( 4, 0) = cmplx(4.926919286686766622234621841806266399811e-4_tp,  &
                       -7.673222155837191136432586925770111584867e-4_tp,tp)
  values( 5, 0) = cmplx(1.339275516728503886085557975188675163702e-3_tp,  &
                       -2.085798035194157686322542562196870032728e-3_tp,tp)
  values( 6, 0) = cmplx(3.640528300423190167962660984149051867360e-3_tp,  &
                       -5.669786896903858940476818879173754216394e-3_tp,tp)
  values( 7, 0) = cmplx(9.895981925031249713864719801586350136119e-3_tp,  &
                       -1.541207869308895788150537629843732635202e-2_tp,tp)
  values( 8, 0) = cmplx(2.690006784157160778122578819852981646620e-2_tp,  &
                       -4.189437345020454468781373408358609461698e-2_tp,tp)
  values( 9, 0) = cmplx(7.312196559805963236599139851015598529567e-2_tp,  &
                       -1.138807140643680892286189302968000854906e-1_tp,tp)
  values(10, 0) = cmplx(1.987661103464129406288031913435846982927e-1_tp,  &
                       -3.095598756531121984439128249151294316712e-1_tp,tp)
  values(11, 0) = cmplx(5.403023058681397174009366074429766037323e-1_tp,  &
                       -8.414709848078965066525023216302989996225e-1_tp,tp)
  values(12, 0) = cmplx(1.468693939915885157138967597326604261326e+0_tp,  &
                       -2.287355287178842391208171906700501808955e+0_tp,tp)
  values(13, 0) = cmplx(3.992324048441271426506695498488872542168e+0_tp,  &
                       -6.217676312367968204252850304087010991267e+0_tp,tp)
  values(14, 0) = cmplx(1.085226191419795717634004705578458510082e+1_tp,  &
                       -1.690139653515009430510735354587017631994e+1_tp,tp)
  values(15, 0) = cmplx(2.949950635904248131176182655246028008012e+1_tp,  &
                       -4.594275907707917015245512895344831034661e+1_tp,tp)
  values(16, 0) = cmplx(8.018797208429722826939072854263268713595e+1_tp,  &
                       -1.248853671484961641963775578906393869008e+2_tp,tp)
  values(17, 0) = cmplx(2.452968673692215000628885544182439448649e-5_tp,  &
                        3.820272360744745896412858462619482830985e-5_tp,tp)
  values(18, 0) = cmplx(6.667860171476833283424963180029090616739e-5_tp,  &
                        1.038457693797678114167511646823987469663e-4_tp,tp)
  values(19, 0) = cmplx(1.812512313883128927908472072600232218423e-4_tp,  &
                        2.822820678673715782795912746517464184101e-4_tp,tp)
  values(20, 0) = cmplx(4.926919286686766622234621841806266399811e-4_tp,  &
                        7.673222155837191136432586925770111584867e-4_tp,tp)
  values(21, 0) = cmplx(1.339275516728503886085557975188675163702e-3_tp,  &
                        2.085798035194157686322542562196870032728e-3_tp,tp)
  values(22, 0) = cmplx(3.640528300423190167962660984149051867360e-3_tp,  &
                        5.669786896903858940476818879173754216394e-3_tp,tp)
  values(23, 0) = cmplx(9.895981925031249713864719801586350136119e-3_tp,  &
                        1.541207869308895788150537629843732635202e-2_tp,tp)
  values(24, 0) = cmplx(2.690006784157160778122578819852981646620e-2_tp,  &
                        4.189437345020454468781373408358609461698e-2_tp,tp)
  values(25, 0) = cmplx(7.312196559805963236599139851015598529567e-2_tp,  &
                        1.138807140643680892286189302968000854906e-1_tp,tp)
  values(26, 0) = cmplx(1.987661103464129406288031913435846982927e-1_tp,  &
                        3.095598756531121984439128249151294316712e-1_tp,tp)
  values(27, 0) = cmplx(5.403023058681397174009366074429766037323e-1_tp,  &
                        8.414709848078965066525023216302989996225e-1_tp,tp)
  values(28, 0) = cmplx(1.468693939915885157138967597326604261326e+0_tp,  &
                        2.287355287178842391208171906700501808955e+0_tp,tp)
  values(29, 0) = cmplx(3.992324048441271426506695498488872542168e+0_tp,  &
                        6.217676312367968204252850304087010991267e+0_tp,tp)
  values(30, 0) = cmplx(1.085226191419795717634004705578458510082e+1_tp,  &
                        1.690139653515009430510735354587017631994e+1_tp,tp)
  values(31, 0) = cmplx(2.949950635904248131176182655246028008012e+1_tp,  &
                        4.594275907707917015245512895344831034661e+1_tp,tp)
  values(32, 0) = cmplx(8.018797208429722826939072854263268713595e+1_tp,  &
                        1.248853671484961641963775578906393869008e+2_tp,tp)
  values( 1, 1) = cmplx(9.900785055303206164315743802010106882062e-2_tp,  &
                       -9.896964782942461418419330943547487399231e-3_tp,tp)
  values( 2, 1) = cmplx(1.097500455896822786928769329814448753730e-1_tp,  &
                       -1.218291109114472343127335353519583073623e-2_tp,tp)
  values( 3, 1) = cmplx(1.230589580341040133605531202094857163175e-1_tp,  &
                       -1.534708449577958022278419111685424623739e-2_tp,tp)
  values( 4, 1) = cmplx(1.399463695742980596495615804680662524935e-1_tp,  &
                       -1.988272105124490579084547453935560590501e-2_tp,tp)
  values( 5, 1) = cmplx(1.620013552685087333613461944516504005148e-1_tp,  &
                       -2.665259287221909594583727531490892174702e-2_tp,tp)
  values( 6, 1) = cmplx(1.918256594382610733884870582291703267261e-1_tp,  &
                       -3.723117450827144288960204786999931450194e-2_tp,tp)
  values( 7, 1) = cmplx(2.338722441760567034721203821818877603416e-1_tp,  &
                       -5.461504137074193639765375147086260849739e-2_tp,tp)
  values( 8, 1) = cmplx(2.961194169925489721344136369487996645218e-1_tp,  &
                       -8.474168118078147581553330095507118996828e-2_tp,tp)
  values( 9, 1) = cmplx(3.935273565736497648993272266552976229798e-1_tp,  &
                       -1.398233212546408378353541481792487687445e-1_tp,tp)
  values(10, 1) = cmplx(5.553968826533496289075548167857723666892e-1_tp,  &
                       -2.458370070002374304636419918706429350179e-1_tp,tp)
  values(11, 1) = cmplx(8.414709848078965066525023216302989996225e-1_tp,  &
                       -4.596976941318602825990633925570233962676e-1_tp,tp)
  values(12, 1) = cmplx(1.378024613547363774173569752013553035141e+0_tp,  &
                       -9.093306736314786170346021546869487738143e-1_tp,tp)
  values(13, 1) = cmplx(2.440464881850102211453248260212951215120e+0_tp,  &
                       -1.888605715258932996399801021937029888073e+0_tp,tp)
  values(14, 1) = cmplx(4.645818227774396583412749471322393162240e+0_tp,  &
                       -4.085192769125232573898201358182594385900e+0_tp,tp)
  values(15, 1) = cmplx(9.408281441955829141147202068428790039242e+0_tp,  &
                       -9.133619408780835252826981721254880076844e+0_tp,tp)
  values(16, 1) = cmplx(2.003173952192239636705120002322318548387e+1_tp,  &
                       -2.097072552531475356586527157348324028340e+1_tp,tp)
  values(17, 1) = cmplx(9.900785055303206164315743802010106882062e-2_tp,  &
                        9.896964782942461418419330943547487399231e-3_tp,tp)
  values(18, 1) = cmplx(1.097500455896822786928769329814448753730e-1_tp,  &
                        1.218291109114472343127335353519583073623e-2_tp,tp)
  values(19, 1) = cmplx(1.230589580341040133605531202094857163175e-1_tp,  &
                        1.534708449577958022278419111685424623739e-2_tp,tp)
  values(20, 1) = cmplx(1.399463695742980596495615804680662524935e-1_tp,  &
                        1.988272105124490579084547453935560590501e-2_tp,tp)
  values(21, 1) = cmplx(1.620013552685087333613461944516504005148e-1_tp,  &
                        2.665259287221909594583727531490892174702e-2_tp,tp)
  values(22, 1) = cmplx(1.918256594382610733884870582291703267261e-1_tp,  &
                        3.723117450827144288960204786999931450194e-2_tp,tp)
  values(23, 1) = cmplx(2.338722441760567034721203821818877603416e-1_tp,  &
                        5.461504137074193639765375147086260849739e-2_tp,tp)
  values(24, 1) = cmplx(2.961194169925489721344136369487996645218e-1_tp,  &
                        8.474168118078147581553330095507118996828e-2_tp,tp)
  values(25, 1) = cmplx(3.935273565736497648993272266552976229798e-1_tp,  &
                        1.398233212546408378353541481792487687445e-1_tp,tp)
  values(26, 1) = cmplx(5.553968826533496289075548167857723666892e-1_tp,  &
                        2.458370070002374304636419918706429350179e-1_tp,tp)
  values(27, 1) = cmplx(8.414709848078965066525023216302989996225e-1_tp,  &
                        4.596976941318602825990633925570233962676e-1_tp,tp)
  values(28, 1) = cmplx(1.378024613547363774173569752013553035141e+0_tp,  &
                        9.093306736314786170346021546869487738143e-1_tp,tp)
  values(29, 1) = cmplx(2.440464881850102211453248260212951215120e+0_tp,  &
                        1.888605715258932996399801021937029888073e+0_tp,tp)
  values(30, 1) = cmplx(4.645818227774396583412749471322393162240e+0_tp,  &
                        4.085192769125232573898201358182594385900e+0_tp,tp)
  values(31, 1) = cmplx(9.408281441955829141147202068428790039242e+0_tp,  &
                        9.133619408780835252826981721254880076844e+0_tp,tp)
  values(32, 1) = cmplx(2.003173952192239636705120002322318548387e+1_tp,  &
                        2.097072552531475356586527157348324028340e+1_tp,tp)
  values( 1, 2) = cmplx(8.930513325992694896026579159151026533854e-2_tp,  &
                       -7.940816847698448754184646064796277793931e-3_tp,tp)
  values( 2, 2) = cmplx(9.785893293639029530726074337441697502900e-2_tp,  &
                       -9.519557982805063541776376648802349365863e-3_tp,tp)
  values( 3, 2) = cmplx(1.081673141572761149744362958375533617799e-1_tp,  &
                       -1.160252870768706684395651309008738944281e-2_tp,tp)
  values( 4, 2) = cmplx(1.208051626806231697648782882252578367690e-1_tp,  &
                       -1.441749166133975199629040195512889012342e-2_tp,tp)
  values( 5, 2) = cmplx(1.366120124665180188048043272595947707745e-1_tp,  &
                       -1.832656993238315380982784199078097483791e-2_tp,tp)
  values( 6, 2) = cmplx(1.568501106660371567671987214124672184950e-1_tp,  &
                       -2.392378723155314277551933470849358079861e-2_tp,tp)
  values( 7, 2) = cmplx(1.834780038039126542652454248672536215959e-1_tp,  &
                       -3.221574060829267946689791834909775327463e-2_tp,tp)
  values( 8, 2) = cmplx(2.196383430203134559412292390108672196402e-1_tp,  &
                       -4.496555394651066004189864601859867655733e-2_tp,tp)
  values( 9, 2) = cmplx(2.705537216214682616073399389737307045569e-1_tp,  &
                       -6.536520018341371188599289539724096790618e-2_tp,tp)
  values(10, 2) = cmplx(3.452200621734439007780435875424352841643e-1_tp,  &
                       -9.938305517320647031440159567179234914639e-2_tp,tp)
  values(11, 2) = cmplx(4.596976941318602825990633925570233962676e-1_tp,  &
                       -1.585290151921034933474976783697010003774e-1_tp,tp)
  values(12, 2) = cmplx(6.436776435894211956040859533502509044777e-1_tp,  &
                       -2.656530300420574214305162013366978693365e-1_tp,tp)
  values(13, 2) = cmplx(9.539070957918274838612595084725864636630e-1_tp,  &
                       -4.673493097335527562692707567322217122051e-1_tp,tp)
  values(14, 2) = cmplx(1.502264745244842232413644977214977387262e+0_tp,  &
                       -8.609760079601301138281854603225389995461e-1_tp,tp)
  values(15, 2) = cmplx(2.515690892741420695142105293821767072577e+0_tp,  &
                       -1.654482129009853639421219106858278251066e+0_tp,tp)
  values(16, 2) = cmplx(4.466516274420259053889279680369198757798e+0_tp,  &
                       -3.300841850178898902395198378622808305121e+0_tp,tp)
  values(17, 2) = cmplx(8.930513325992694896026579159151026533854e-2_tp,  &
                        7.940816847698448754184646064796277793931e-3_tp,tp)
  values(18, 2) = cmplx(9.785893293639029530726074337441697502900e-2_tp,  &
                        9.519557982805063541776376648802349365863e-3_tp,tp)
  values(19, 2) = cmplx(1.081673141572761149744362958375533617799e-1_tp,  &
                        1.160252870768706684395651309008738944281e-2_tp,tp)
  values(20, 2) = cmplx(1.208051626806231697648782882252578367690e-1_tp,  &
                        1.441749166133975199629040195512889012342e-2_tp,tp)
  values(21, 2) = cmplx(1.366120124665180188048043272595947707745e-1_tp,  &
                        1.832656993238315380982784199078097483791e-2_tp,tp)
  values(22, 2) = cmplx(1.568501106660371567671987214124672184950e-1_tp,  &
                        2.392378723155314277551933470849358079861e-2_tp,tp)
  values(23, 2) = cmplx(1.834780038039126542652454248672536215959e-1_tp,  &
                        3.221574060829267946689791834909775327463e-2_tp,tp)
  values(24, 2) = cmplx(2.196383430203134559412292390108672196402e-1_tp,  &
                        4.496555394651066004189864601859867655733e-2_tp,tp)
  values(25, 2) = cmplx(2.705537216214682616073399389737307045569e-1_tp,  &
                        6.536520018341371188599289539724096790618e-2_tp,tp)
  values(26, 2) = cmplx(3.452200621734439007780435875424352841643e-1_tp,  &
                        9.938305517320647031440159567179234914639e-2_tp,tp)
  values(27, 2) = cmplx(4.596976941318602825990633925570233962676e-1_tp,  &
                        1.585290151921034933474976783697010003774e-1_tp,tp)
  values(28, 2) = cmplx(6.436776435894211956040859533502509044777e-1_tp,  &
                        2.656530300420574214305162013366978693365e-1_tp,tp)
  values(29, 2) = cmplx(9.539070957918274838612595084725864636630e-1_tp,  &
                        4.673493097335527562692707567322217122051e-1_tp,tp)
  values(30, 2) = cmplx(1.502264745244842232413644977214977387262e+0_tp,  &
                        8.609760079601301138281854603225389995461e-1_tp,tp)
  values(31, 2) = cmplx(2.515690892741420695142105293821767072577e+0_tp,  &
                        1.654482129009853639421219106858278251066e+0_tp,tp)
  values(32, 2) = cmplx(4.466516274420259053889279680369198757798e+0_tp,  &
                        3.300841850178898902395198378622808305121e+0_tp,tp)
  values( 1, 3) = cmplx(4.074148004206365306090620524900686756840e-2_tp,  &
                       -3.280066319436520430672155918421058977446e-3_tp,tp)
  values( 2, 3) = cmplx(4.425352636043039519239548397901279968420e-2_tp,  &
                       -3.859329819736147961179900814467827813148e-3_tp,tp)
  values( 3, 3) = cmplx(4.840406177614581764689947917522554608004e-2_tp,  &
                       -4.600191633557343850367870760642269579653e-3_tp,tp)
  values( 4, 3) = cmplx(5.337562705793955127284284768756648065480e-2_tp,  &
                       -5.565447913799971325221777961776798647340e-3_tp,tp)
  values( 5, 3) = cmplx(5.942309446306148759408113184954627973488e-2_tp,  &
                       -6.849420755113055630708881643127550816161e-3_tp,tp)
  values( 6, 3) = cmplx(6.691050899620643688228945106331374955089e-2_tp,  &
                       -8.597344352930658821354023270964033750455e-3_tp,tp)
  values( 7, 3) = cmplx(7.637080737603776837681860111059313334652e-2_tp,  &
                       -1.103876669193627222748017069037384501797e-2_tp,tp)
  values( 8, 3) = cmplx(8.860505248855702922182109289859970176364e-2_tp,  &
                       -1.454649951401545639330748229333367506877e-2_tp,tp)
  values( 9, 3) = cmplx(1.048515513880954377342626034899559117584e-1_tp,  &
                       -1.974317560234086292413485404635747192613e-2_tp,tp)
  values(10, 3) = cmplx(1.270814964998812847681790040646785324910e-1_tp,  &
                       -2.769844132667481445377740839288618334462e-2_tp,tp)
  values(11, 3) = cmplx(1.585290151921034933474976783697010003774e-1_tp,  &
                       -4.030230586813971740093660744297660373231e-2_tp,tp)
  values(12, 3) = cmplx(2.046653368157393085173010773434743869071e-1_tp,  &
                       -6.098769322631811291321512399322348242936e-2_tp,tp)
  values(13, 3) = cmplx(2.750327002634415447983579547354789279062e-1_tp,  &
                       -9.615830473505560573545640099837139214946e-2_tp,tp)
  values(14, 3) = cmplx(3.867770243694656811069120391967471161333e-1_tp,  &
                       -1.580663278635548109070911403752639611376e-1_tp,tp)
  values(15, 3) = cmplx(5.716026882338550835288023695379615612573e-1_tp,  &
                       -2.707198601939996389731041843300791724523e-1_tp,tp)
  values(16, 3) = cmplx(8.897470470107766989169844915564923882352e-1_tp,  &
                       -4.822189606336244406956427774132631833772e-1_tp,tp)
  values(17, 3) = cmplx(4.074148004206365306090620524900686756840e-2_tp,  &
                        3.280066319436520430672155918421058977446e-3_tp,tp)
  values(18, 3) = cmplx(4.425352636043039519239548397901279968420e-2_tp,  &
                        3.859329819736147961179900814467827813148e-3_tp,tp)
  values(19, 3) = cmplx(4.840406177614581764689947917522554608004e-2_tp,  &
                        4.600191633557343850367870760642269579653e-3_tp,tp)
  values(20, 3) = cmplx(5.337562705793955127284284768756648065480e-2_tp,  &
                        5.565447913799971325221777961776798647340e-3_tp,tp)
  values(21, 3) = cmplx(5.942309446306148759408113184954627973488e-2_tp,  &
                        6.849420755113055630708881643127550816161e-3_tp,tp)
  values(22, 3) = cmplx(6.691050899620643688228945106331374955089e-2_tp,  &
                        8.597344352930658821354023270964033750455e-3_tp,tp)
  values(23, 3) = cmplx(7.637080737603776837681860111059313334652e-2_tp,  &
                        1.103876669193627222748017069037384501797e-2_tp,tp)
  values(24, 3) = cmplx(8.860505248855702922182109289859970176364e-2_tp,  &
                        1.454649951401545639330748229333367506877e-2_tp,tp)
  values(25, 3) = cmplx(1.048515513880954377342626034899559117584e-1_tp,  &
                        1.974317560234086292413485404635747192613e-2_tp,tp)
  values(26, 3) = cmplx(1.270814964998812847681790040646785324910e-1_tp,  &
                        2.769844132667481445377740839288618334462e-2_tp,tp)
  values(27, 3) = cmplx(1.585290151921034933474976783697010003774e-1_tp,  &
                        4.030230586813971740093660744297660373231e-2_tp,tp)
  values(28, 3) = cmplx(2.046653368157393085173010773434743869071e-1_tp,  &
                        6.098769322631811291321512399322348242936e-2_tp,tp)
  values(29, 3) = cmplx(2.750327002634415447983579547354789279062e-1_tp,  &
                        9.615830473505560573545640099837139214946e-2_tp,tp)
  values(30, 3) = cmplx(3.867770243694656811069120391967471161333e-1_tp,  &
                        1.580663278635548109070911403752639611376e-1_tp,tp)
  values(31, 3) = cmplx(5.716026882338550835288023695379615612573e-1_tp,  &
                        2.707198601939996389731041843300791724523e-1_tp,tp)
  values(32, 3) = cmplx(8.897470470107766989169844915564923882352e-1_tp,  &
                        4.822189606336244406956427774132631833772e-1_tp,tp)
  values( 1, 4) = cmplx(1.250031616401452135136907693163385197980e-2_tp,  &
                       -9.220249844578000920696921013212793002356e-4_tp,tp)
  values( 2, 4) = cmplx(1.348265356799832428328805542687015403238e-2_tp,  &
                       -1.069258194251352924678683845822480691025e-3_tp,tp)
  values( 3, 4) = cmplx(1.462616970396498670782315954911032668111e-2_tp,  &
                       -1.253247258800955357181911098558507137682e-3_tp,tp)
  values( 4, 4) = cmplx(1.597205450349779558163977021630956201460e-2_tp,  &
                       -1.486658084242546322345427464933251909609e-3_tp,tp)
  values( 5, 4) = cmplx(1.757596902639849000178978623096891547045e-2_tp,  &
                       -1.787758045214239061846817431306894109049e-3_tp,tp)
  values( 6, 4) = cmplx(1.951454356558583875935538851106648535881e-2_tp,  &
                       -2.183439842531035987600273048020490321672e-3_tp,tp)
  values( 7, 4) = cmplx(2.189542375614422737569837840674517519403e-2_tp,  &
                       -2.714164266051988787054551929092832544015e-3_tp,tp)
  values( 8, 4) = cmplx(2.487313420483443687278442035975345697778e-2_tp,  &
                       -3.442211563606326826492312688806593969669e-3_tp,tp)
  values( 9, 4) = cmplx(2.867468123189666415778859607995579634851e-2_tp,  &
                       -4.465752814777900616826871016799162211191e-3_tp,tp)
  values(10, 4) = cmplx(3.364180574673009817613253549743715876013e-2_tp,  &
                       -5.943364420055283722355127104550975415511e-3_tp,tp)
  values(11, 4) = cmplx(4.030230586813971740093660744297660373231e-2_tp,  &
                       -8.137651474563173319168988296965666289229e-3_tp,tp)
  values(12, 4) = cmplx(4.949318168769537738192476733501560133492e-2_tp,  &
                       -1.149451153862273553129035665820788109443e-2_tp,tp)
  values(13, 4) = cmplx(6.257807438572107239976779542719918292572e-2_tp,  &
                       -1.679011517466726666784430278558610461187e-2_tp,tp)
  values(14, 4) = cmplx(8.183974009719518542278272579655053095375e-2_tp,  &
                       -2.540886258878654182810280485957114339461e-2_tp,tp)
  values(15, 4) = cmplx(1.112037615566325474365674703420740441656e-1_tp,  &
                       -3.987902465934177288413417849700128207168e-2_tp,tp)
  values(16, 4) = cmplx(1.576008023982374846902781500716304535084e-1_tp,  &
                       -6.492363164707739120107292546832654597374e-2_tp,tp)
  values(17, 4) = cmplx(1.250031616401452135136907693163385197980e-2_tp,  &
                        9.220249844578000920696921013212793002356e-4_tp,tp)
  values(18, 4) = cmplx(1.348265356799832428328805542687015403238e-2_tp,  &
                        1.069258194251352924678683845822480691025e-3_tp,tp)
  values(19, 4) = cmplx(1.462616970396498670782315954911032668111e-2_tp,  &
                        1.253247258800955357181911098558507137682e-3_tp,tp)
  values(20, 4) = cmplx(1.597205450349779558163977021630956201460e-2_tp,  &
                        1.486658084242546322345427464933251909609e-3_tp,tp)
  values(21, 4) = cmplx(1.757596902639849000178978623096891547045e-2_tp,  &
                        1.787758045214239061846817431306894109049e-3_tp,tp)
  values(22, 4) = cmplx(1.951454356558583875935538851106648535881e-2_tp,  &
                        2.183439842531035987600273048020490321672e-3_tp,tp)
  values(23, 4) = cmplx(2.189542375614422737569837840674517519403e-2_tp,  &
                        2.714164266051988787054551929092832544015e-3_tp,tp)
  values(24, 4) = cmplx(2.487313420483443687278442035975345697778e-2_tp,  &
                        3.442211563606326826492312688806593969669e-3_tp,tp)
  values(25, 4) = cmplx(2.867468123189666415778859607995579634851e-2_tp,  &
                        4.465752814777900616826871016799162211191e-3_tp,tp)
  values(26, 4) = cmplx(3.364180574673009817613253549743715876013e-2_tp,  &
                        5.943364420055283722355127104550975415511e-3_tp,tp)
  values(27, 4) = cmplx(4.030230586813971740093660744297660373231e-2_tp,  &
                        8.137651474563173319168988296965666289229e-3_tp,tp)
  values(28, 4) = cmplx(4.949318168769537738192476733501560133492e-2_tp,  &
                        1.149451153862273553129035665820788109443e-2_tp,tp)
  values(29, 4) = cmplx(6.257807438572107239976779542719918292572e-2_tp,  &
                        1.679011517466726666784430278558610461187e-2_tp,tp)
  values(30, 4) = cmplx(8.183974009719518542278272579655053095375e-2_tp,  &
                        2.540886258878654182810280485957114339461e-2_tp,tp)
  values(31, 4) = cmplx(1.112037615566325474365674703420740441656e-1_tp,  &
                        3.987902465934177288413417849700128207168e-2_tp,tp)
  values(32, 4) = cmplx(1.576008023982374846902781500716304535084e-1_tp,  &
                        6.492363164707739120107292546832654597374e-2_tp,tp)
  values( 1, 5) = cmplx(2.896886435752269834109362271798509169988e-3_tp,  &
                       -1.974861451294469742039670170477229869753e-4_tp,tp)
  values( 2, 5) = cmplx(3.106407025393493102135197378097452370726e-3_tp,  &
                       -2.263498701269044641618348369194412977445e-4_tp,tp)
  values( 3, 5) = cmplx(3.347341891698683000429691815984757338801e-3_tp,  &
                       -2.617618291122159554059725896782812751398e-4_tp,tp)
  values( 4, 5) = cmplx(3.626978864528492878350674052348659689480e-3_tp,  &
                       -3.057601114694209365721780839164868256958e-4_tp,tp)
  values( 5, 5) = cmplx(3.954917402346575650029948649878200034765e-3_tp,  &
                       -3.611932261887227646971885364285509876193e-4_tp,tp)
  values( 6, 5) = cmplx(4.344002128766737520159871685616207571573e-3_tp,  &
                       -4.321124572471403065119197275191434499801e-4_tp,tp)
  values( 7, 5) = cmplx(4.811713876949514467701629704045811672620e-3_tp,  &
                       -5.243874027243814201617694437382447821514e-4_tp,tp)
  values( 8, 5) = cmplx(5.382280894910301620813905160954622303632e-3_tp,  &
                       -6.466897771013249314405308240493427779874e-4_tp,tp)
  values( 9, 5) = cmplx(6.089944736863581126916602438044180569499e-3_tp,  &
                       -8.120959610428402550448657106225091791540e-4_tp,tp)
  values(10, 5) = cmplx(6.984112669995926106444629136890241661021e-3_tp,  &
                       -1.040748249940642384089502032339266245509e-3_tp,tp)
  values(11, 5) = cmplx(8.137651474563173319168988296965666289229e-3_tp,  &
                       -1.364360798526949265730059223690062934356e-3_tp,tp)
  values(12, 5) = cmplx(9.660513279825723123274228663278407881348e-3_tp,  &
                       -1.833998258797012408016127994929473213087e-3_tp,tp)
  values(13, 5) = cmplx(1.172258612255521562680931206133022742599e-2_tp,  &
                       -2.533764526056025520517495362127938592938e-3_tp,tp)
  values(14, 5) = cmplx(1.459280828803720980964509822492227362558e-2_tp,  &
                       -3.605351433583110672819235544882956589676e-3_tp,tp)
  values(15, 5) = cmplx(1.870749436583560564492572901168416423925e-2_tp,  &
                       -5.292882573376541809802112371329279458106e-3_tp,tp)
  values(16, 5) = cmplx(2.479208885788198005073578240358251846856e-2_tp,  &
                       -8.026308557839082230067428612948805501036e-3_tp,tp)
  values(17, 5) = cmplx(2.896886435752269834109362271798509169988e-3_tp,  &
                        1.974861451294469742039670170477229869753e-4_tp,tp)
  values(18, 5) = cmplx(3.106407025393493102135197378097452370726e-3_tp,  &
                        2.263498701269044641618348369194412977445e-4_tp,tp)
  values(19, 5) = cmplx(3.347341891698683000429691815984757338801e-3_tp,  &
                        2.617618291122159554059725896782812751398e-4_tp,tp)
  values(20, 5) = cmplx(3.626978864528492878350674052348659689480e-3_tp,  &
                        3.057601114694209365721780839164868256958e-4_tp,tp)
  values(21, 5) = cmplx(3.954917402346575650029948649878200034765e-3_tp,  &
                        3.611932261887227646971885364285509876193e-4_tp,tp)
  values(22, 5) = cmplx(4.344002128766737520159871685616207571573e-3_tp,  &
                        4.321124572471403065119197275191434499801e-4_tp,tp)
  values(23, 5) = cmplx(4.811713876949514467701629704045811672620e-3_tp,  &
                        5.243874027243814201617694437382447821514e-4_tp,tp)
  values(24, 5) = cmplx(5.382280894910301620813905160954622303632e-3_tp,  &
                        6.466897771013249314405308240493427779874e-4_tp,tp)
  values(25, 5) = cmplx(6.089944736863581126916602438044180569499e-3_tp,  &
                        8.120959610428402550448657106225091791540e-4_tp,tp)
  values(26, 5) = cmplx(6.984112669995926106444629136890241661021e-3_tp,  &
                        1.040748249940642384089502032339266245509e-3_tp,tp)
  values(27, 5) = cmplx(8.137651474563173319168988296965666289229e-3_tp,  &
                        1.364360798526949265730059223690062934356e-3_tp,tp)
  values(28, 5) = cmplx(9.660513279825723123274228663278407881348e-3_tp,  &
                        1.833998258797012408016127994929473213087e-3_tp,tp)
  values(29, 5) = cmplx(1.172258612255521562680931206133022742599e-2_tp,  &
                        2.533764526056025520517495362127938592938e-3_tp,tp)
  values(30, 5) = cmplx(1.459280828803720980964509822492227362558e-2_tp,  &
                        3.605351433583110672819235544882956589676e-3_tp,tp)
  values(31, 5) = cmplx(1.870749436583560564492572901168416423925e-2_tp,  &
                        5.292882573376541809802112371329279458106e-3_tp,tp)
  values(32, 5) = cmplx(2.479208885788198005073578240358251846856e-2_tp,  &
                        8.026308557839082230067428612948805501036e-3_tp,tp)
end subroutine tp_compref

#endif

#ifdef __USE_QPREC

subroutine qp_compref(values)
  complex(kind=qp), intent(out) :: values(32,0:5)

  values( 1, 0) = cmplx(2.452968673692215000628885544182439448649e-5_qp,  &
                       -3.820272360744745896412858462619482830985e-5_qp,qp)
  values( 2, 0) = cmplx(6.667860171476833283424963180029090616739e-5_qp,  &
                       -1.038457693797678114167511646823987469663e-4_qp,qp)
  values( 3, 0) = cmplx(1.812512313883128927908472072600232218423e-4_qp,  &
                       -2.822820678673715782795912746517464184101e-4_qp,qp)
  values( 4, 0) = cmplx(4.926919286686766622234621841806266399811e-4_qp,  &
                       -7.673222155837191136432586925770111584867e-4_qp,qp)
  values( 5, 0) = cmplx(1.339275516728503886085557975188675163702e-3_qp,  &
                       -2.085798035194157686322542562196870032728e-3_qp,qp)
  values( 6, 0) = cmplx(3.640528300423190167962660984149051867360e-3_qp,  &
                       -5.669786896903858940476818879173754216394e-3_qp,qp)
  values( 7, 0) = cmplx(9.895981925031249713864719801586350136119e-3_qp,  &
                       -1.541207869308895788150537629843732635202e-2_qp,qp)
  values( 8, 0) = cmplx(2.690006784157160778122578819852981646620e-2_qp,  &
                       -4.189437345020454468781373408358609461698e-2_qp,qp)
  values( 9, 0) = cmplx(7.312196559805963236599139851015598529567e-2_qp,  &
                       -1.138807140643680892286189302968000854906e-1_qp,qp)
  values(10, 0) = cmplx(1.987661103464129406288031913435846982927e-1_qp,  &
                       -3.095598756531121984439128249151294316712e-1_qp,qp)
  values(11, 0) = cmplx(5.403023058681397174009366074429766037323e-1_qp,  &
                       -8.414709848078965066525023216302989996225e-1_qp,qp)
  values(12, 0) = cmplx(1.468693939915885157138967597326604261326e+0_qp,  &
                       -2.287355287178842391208171906700501808955e+0_qp,qp)
  values(13, 0) = cmplx(3.992324048441271426506695498488872542168e+0_qp,  &
                       -6.217676312367968204252850304087010991267e+0_qp,qp)
  values(14, 0) = cmplx(1.085226191419795717634004705578458510082e+1_qp,  &
                       -1.690139653515009430510735354587017631994e+1_qp,qp)
  values(15, 0) = cmplx(2.949950635904248131176182655246028008012e+1_qp,  &
                       -4.594275907707917015245512895344831034661e+1_qp,qp)
  values(16, 0) = cmplx(8.018797208429722826939072854263268713595e+1_qp,  &
                       -1.248853671484961641963775578906393869008e+2_qp,qp)
  values(17, 0) = cmplx(2.452968673692215000628885544182439448649e-5_qp,  &
                        3.820272360744745896412858462619482830985e-5_qp,qp)
  values(18, 0) = cmplx(6.667860171476833283424963180029090616739e-5_qp,  &
                        1.038457693797678114167511646823987469663e-4_qp,qp)
  values(19, 0) = cmplx(1.812512313883128927908472072600232218423e-4_qp,  &
                        2.822820678673715782795912746517464184101e-4_qp,qp)
  values(20, 0) = cmplx(4.926919286686766622234621841806266399811e-4_qp,  &
                        7.673222155837191136432586925770111584867e-4_qp,qp)
  values(21, 0) = cmplx(1.339275516728503886085557975188675163702e-3_qp,  &
                        2.085798035194157686322542562196870032728e-3_qp,qp)
  values(22, 0) = cmplx(3.640528300423190167962660984149051867360e-3_qp,  &
                        5.669786896903858940476818879173754216394e-3_qp,qp)
  values(23, 0) = cmplx(9.895981925031249713864719801586350136119e-3_qp,  &
                        1.541207869308895788150537629843732635202e-2_qp,qp)
  values(24, 0) = cmplx(2.690006784157160778122578819852981646620e-2_qp,  &
                        4.189437345020454468781373408358609461698e-2_qp,qp)
  values(25, 0) = cmplx(7.312196559805963236599139851015598529567e-2_qp,  &
                        1.138807140643680892286189302968000854906e-1_qp,qp)
  values(26, 0) = cmplx(1.987661103464129406288031913435846982927e-1_qp,  &
                        3.095598756531121984439128249151294316712e-1_qp,qp)
  values(27, 0) = cmplx(5.403023058681397174009366074429766037323e-1_qp,  &
                        8.414709848078965066525023216302989996225e-1_qp,qp)
  values(28, 0) = cmplx(1.468693939915885157138967597326604261326e+0_qp,  &
                        2.287355287178842391208171906700501808955e+0_qp,qp)
  values(29, 0) = cmplx(3.992324048441271426506695498488872542168e+0_qp,  &
                        6.217676312367968204252850304087010991267e+0_qp,qp)
  values(30, 0) = cmplx(1.085226191419795717634004705578458510082e+1_qp,  &
                        1.690139653515009430510735354587017631994e+1_qp,qp)
  values(31, 0) = cmplx(2.949950635904248131176182655246028008012e+1_qp,  &
                        4.594275907707917015245512895344831034661e+1_qp,qp)
  values(32, 0) = cmplx(8.018797208429722826939072854263268713595e+1_qp,  &
                        1.248853671484961641963775578906393869008e+2_qp,qp)
  values( 1, 1) = cmplx(9.900785055303206164315743802010106882062e-2_qp,  &
                       -9.896964782942461418419330943547487399231e-3_qp,qp)
  values( 2, 1) = cmplx(1.097500455896822786928769329814448753730e-1_qp,  &
                       -1.218291109114472343127335353519583073623e-2_qp,qp)
  values( 3, 1) = cmplx(1.230589580341040133605531202094857163175e-1_qp,  &
                       -1.534708449577958022278419111685424623739e-2_qp,qp)
  values( 4, 1) = cmplx(1.399463695742980596495615804680662524935e-1_qp,  &
                       -1.988272105124490579084547453935560590501e-2_qp,qp)
  values( 5, 1) = cmplx(1.620013552685087333613461944516504005148e-1_qp,  &
                       -2.665259287221909594583727531490892174702e-2_qp,qp)
  values( 6, 1) = cmplx(1.918256594382610733884870582291703267261e-1_qp,  &
                       -3.723117450827144288960204786999931450194e-2_qp,qp)
  values( 7, 1) = cmplx(2.338722441760567034721203821818877603416e-1_qp,  &
                       -5.461504137074193639765375147086260849739e-2_qp,qp)
  values( 8, 1) = cmplx(2.961194169925489721344136369487996645218e-1_qp,  &
                       -8.474168118078147581553330095507118996828e-2_qp,qp)
  values( 9, 1) = cmplx(3.935273565736497648993272266552976229798e-1_qp,  &
                       -1.398233212546408378353541481792487687445e-1_qp,qp)
  values(10, 1) = cmplx(5.553968826533496289075548167857723666892e-1_qp,  &
                       -2.458370070002374304636419918706429350179e-1_qp,qp)
  values(11, 1) = cmplx(8.414709848078965066525023216302989996225e-1_qp,  &
                       -4.596976941318602825990633925570233962676e-1_qp,qp)
  values(12, 1) = cmplx(1.378024613547363774173569752013553035141e+0_qp,  &
                       -9.093306736314786170346021546869487738143e-1_qp,qp)
  values(13, 1) = cmplx(2.440464881850102211453248260212951215120e+0_qp,  &
                       -1.888605715258932996399801021937029888073e+0_qp,qp)
  values(14, 1) = cmplx(4.645818227774396583412749471322393162240e+0_qp,  &
                       -4.085192769125232573898201358182594385900e+0_qp,qp)
  values(15, 1) = cmplx(9.408281441955829141147202068428790039242e+0_qp,  &
                       -9.133619408780835252826981721254880076844e+0_qp,qp)
  values(16, 1) = cmplx(2.003173952192239636705120002322318548387e+1_qp,  &
                       -2.097072552531475356586527157348324028340e+1_qp,qp)
  values(17, 1) = cmplx(9.900785055303206164315743802010106882062e-2_qp,  &
                        9.896964782942461418419330943547487399231e-3_qp,qp)
  values(18, 1) = cmplx(1.097500455896822786928769329814448753730e-1_qp,  &
                        1.218291109114472343127335353519583073623e-2_qp,qp)
  values(19, 1) = cmplx(1.230589580341040133605531202094857163175e-1_qp,  &
                        1.534708449577958022278419111685424623739e-2_qp,qp)
  values(20, 1) = cmplx(1.399463695742980596495615804680662524935e-1_qp,  &
                        1.988272105124490579084547453935560590501e-2_qp,qp)
  values(21, 1) = cmplx(1.620013552685087333613461944516504005148e-1_qp,  &
                        2.665259287221909594583727531490892174702e-2_qp,qp)
  values(22, 1) = cmplx(1.918256594382610733884870582291703267261e-1_qp,  &
                        3.723117450827144288960204786999931450194e-2_qp,qp)
  values(23, 1) = cmplx(2.338722441760567034721203821818877603416e-1_qp,  &
                        5.461504137074193639765375147086260849739e-2_qp,qp)
  values(24, 1) = cmplx(2.961194169925489721344136369487996645218e-1_qp,  &
                        8.474168118078147581553330095507118996828e-2_qp,qp)
  values(25, 1) = cmplx(3.935273565736497648993272266552976229798e-1_qp,  &
                        1.398233212546408378353541481792487687445e-1_qp,qp)
  values(26, 1) = cmplx(5.553968826533496289075548167857723666892e-1_qp,  &
                        2.458370070002374304636419918706429350179e-1_qp,qp)
  values(27, 1) = cmplx(8.414709848078965066525023216302989996225e-1_qp,  &
                        4.596976941318602825990633925570233962676e-1_qp,qp)
  values(28, 1) = cmplx(1.378024613547363774173569752013553035141e+0_qp,  &
                        9.093306736314786170346021546869487738143e-1_qp,qp)
  values(29, 1) = cmplx(2.440464881850102211453248260212951215120e+0_qp,  &
                        1.888605715258932996399801021937029888073e+0_qp,qp)
  values(30, 1) = cmplx(4.645818227774396583412749471322393162240e+0_qp,  &
                        4.085192769125232573898201358182594385900e+0_qp,qp)
  values(31, 1) = cmplx(9.408281441955829141147202068428790039242e+0_qp,  &
                        9.133619408780835252826981721254880076844e+0_qp,qp)
  values(32, 1) = cmplx(2.003173952192239636705120002322318548387e+1_qp,  &
                        2.097072552531475356586527157348324028340e+1_qp,qp)
  values( 1, 2) = cmplx(8.930513325992694896026579159151026533854e-2_qp,  &
                       -7.940816847698448754184646064796277793931e-3_qp,qp)
  values( 2, 2) = cmplx(9.785893293639029530726074337441697502900e-2_qp,  &
                       -9.519557982805063541776376648802349365863e-3_qp,qp)
  values( 3, 2) = cmplx(1.081673141572761149744362958375533617799e-1_qp,  &
                       -1.160252870768706684395651309008738944281e-2_qp,qp)
  values( 4, 2) = cmplx(1.208051626806231697648782882252578367690e-1_qp,  &
                       -1.441749166133975199629040195512889012342e-2_qp,qp)
  values( 5, 2) = cmplx(1.366120124665180188048043272595947707745e-1_qp,  &
                       -1.832656993238315380982784199078097483791e-2_qp,qp)
  values( 6, 2) = cmplx(1.568501106660371567671987214124672184950e-1_qp,  &
                       -2.392378723155314277551933470849358079861e-2_qp,qp)
  values( 7, 2) = cmplx(1.834780038039126542652454248672536215959e-1_qp,  &
                       -3.221574060829267946689791834909775327463e-2_qp,qp)
  values( 8, 2) = cmplx(2.196383430203134559412292390108672196402e-1_qp,  &
                       -4.496555394651066004189864601859867655733e-2_qp,qp)
  values( 9, 2) = cmplx(2.705537216214682616073399389737307045569e-1_qp,  &
                       -6.536520018341371188599289539724096790618e-2_qp,qp)
  values(10, 2) = cmplx(3.452200621734439007780435875424352841643e-1_qp,  &
                       -9.938305517320647031440159567179234914639e-2_qp,qp)
  values(11, 2) = cmplx(4.596976941318602825990633925570233962676e-1_qp,  &
                       -1.585290151921034933474976783697010003774e-1_qp,qp)
  values(12, 2) = cmplx(6.436776435894211956040859533502509044777e-1_qp,  &
                       -2.656530300420574214305162013366978693365e-1_qp,qp)
  values(13, 2) = cmplx(9.539070957918274838612595084725864636630e-1_qp,  &
                       -4.673493097335527562692707567322217122051e-1_qp,qp)
  values(14, 2) = cmplx(1.502264745244842232413644977214977387262e+0_qp,  &
                       -8.609760079601301138281854603225389995461e-1_qp,qp)
  values(15, 2) = cmplx(2.515690892741420695142105293821767072577e+0_qp,  &
                       -1.654482129009853639421219106858278251066e+0_qp,qp)
  values(16, 2) = cmplx(4.466516274420259053889279680369198757798e+0_qp,  &
                       -3.300841850178898902395198378622808305121e+0_qp,qp)
  values(17, 2) = cmplx(8.930513325992694896026579159151026533854e-2_qp,  &
                        7.940816847698448754184646064796277793931e-3_qp,qp)
  values(18, 2) = cmplx(9.785893293639029530726074337441697502900e-2_qp,  &
                        9.519557982805063541776376648802349365863e-3_qp,qp)
  values(19, 2) = cmplx(1.081673141572761149744362958375533617799e-1_qp,  &
                        1.160252870768706684395651309008738944281e-2_qp,qp)
  values(20, 2) = cmplx(1.208051626806231697648782882252578367690e-1_qp,  &
                        1.441749166133975199629040195512889012342e-2_qp,qp)
  values(21, 2) = cmplx(1.366120124665180188048043272595947707745e-1_qp,  &
                        1.832656993238315380982784199078097483791e-2_qp,qp)
  values(22, 2) = cmplx(1.568501106660371567671987214124672184950e-1_qp,  &
                        2.392378723155314277551933470849358079861e-2_qp,qp)
  values(23, 2) = cmplx(1.834780038039126542652454248672536215959e-1_qp,  &
                        3.221574060829267946689791834909775327463e-2_qp,qp)
  values(24, 2) = cmplx(2.196383430203134559412292390108672196402e-1_qp,  &
                        4.496555394651066004189864601859867655733e-2_qp,qp)
  values(25, 2) = cmplx(2.705537216214682616073399389737307045569e-1_qp,  &
                        6.536520018341371188599289539724096790618e-2_qp,qp)
  values(26, 2) = cmplx(3.452200621734439007780435875424352841643e-1_qp,  &
                        9.938305517320647031440159567179234914639e-2_qp,qp)
  values(27, 2) = cmplx(4.596976941318602825990633925570233962676e-1_qp,  &
                        1.585290151921034933474976783697010003774e-1_qp,qp)
  values(28, 2) = cmplx(6.436776435894211956040859533502509044777e-1_qp,  &
                        2.656530300420574214305162013366978693365e-1_qp,qp)
  values(29, 2) = cmplx(9.539070957918274838612595084725864636630e-1_qp,  &
                        4.673493097335527562692707567322217122051e-1_qp,qp)
  values(30, 2) = cmplx(1.502264745244842232413644977214977387262e+0_qp,  &
                        8.609760079601301138281854603225389995461e-1_qp,qp)
  values(31, 2) = cmplx(2.515690892741420695142105293821767072577e+0_qp,  &
                        1.654482129009853639421219106858278251066e+0_qp,qp)
  values(32, 2) = cmplx(4.466516274420259053889279680369198757798e+0_qp,  &
                        3.300841850178898902395198378622808305121e+0_qp,qp)
  values( 1, 3) = cmplx(4.074148004206365306090620524900686756840e-2_qp,  &
                       -3.280066319436520430672155918421058977446e-3_qp,qp)
  values( 2, 3) = cmplx(4.425352636043039519239548397901279968420e-2_qp,  &
                       -3.859329819736147961179900814467827813148e-3_qp,qp)
  values( 3, 3) = cmplx(4.840406177614581764689947917522554608004e-2_qp,  &
                       -4.600191633557343850367870760642269579653e-3_qp,qp)
  values( 4, 3) = cmplx(5.337562705793955127284284768756648065480e-2_qp,  &
                       -5.565447913799971325221777961776798647340e-3_qp,qp)
  values( 5, 3) = cmplx(5.942309446306148759408113184954627973488e-2_qp,  &
                       -6.849420755113055630708881643127550816161e-3_qp,qp)
  values( 6, 3) = cmplx(6.691050899620643688228945106331374955089e-2_qp,  &
                       -8.597344352930658821354023270964033750455e-3_qp,qp)
  values( 7, 3) = cmplx(7.637080737603776837681860111059313334652e-2_qp,  &
                       -1.103876669193627222748017069037384501797e-2_qp,qp)
  values( 8, 3) = cmplx(8.860505248855702922182109289859970176364e-2_qp,  &
                       -1.454649951401545639330748229333367506877e-2_qp,qp)
  values( 9, 3) = cmplx(1.048515513880954377342626034899559117584e-1_qp,  &
                       -1.974317560234086292413485404635747192613e-2_qp,qp)
  values(10, 3) = cmplx(1.270814964998812847681790040646785324910e-1_qp,  &
                       -2.769844132667481445377740839288618334462e-2_qp,qp)
  values(11, 3) = cmplx(1.585290151921034933474976783697010003774e-1_qp,  &
                       -4.030230586813971740093660744297660373231e-2_qp,qp)
  values(12, 3) = cmplx(2.046653368157393085173010773434743869071e-1_qp,  &
                       -6.098769322631811291321512399322348242936e-2_qp,qp)
  values(13, 3) = cmplx(2.750327002634415447983579547354789279062e-1_qp,  &
                       -9.615830473505560573545640099837139214946e-2_qp,qp)
  values(14, 3) = cmplx(3.867770243694656811069120391967471161333e-1_qp,  &
                       -1.580663278635548109070911403752639611376e-1_qp,qp)
  values(15, 3) = cmplx(5.716026882338550835288023695379615612573e-1_qp,  &
                       -2.707198601939996389731041843300791724523e-1_qp,qp)
  values(16, 3) = cmplx(8.897470470107766989169844915564923882352e-1_qp,  &
                       -4.822189606336244406956427774132631833772e-1_qp,qp)
  values(17, 3) = cmplx(4.074148004206365306090620524900686756840e-2_qp,  &
                        3.280066319436520430672155918421058977446e-3_qp,qp)
  values(18, 3) = cmplx(4.425352636043039519239548397901279968420e-2_qp,  &
                        3.859329819736147961179900814467827813148e-3_qp,qp)
  values(19, 3) = cmplx(4.840406177614581764689947917522554608004e-2_qp,  &
                        4.600191633557343850367870760642269579653e-3_qp,qp)
  values(20, 3) = cmplx(5.337562705793955127284284768756648065480e-2_qp,  &
                        5.565447913799971325221777961776798647340e-3_qp,qp)
  values(21, 3) = cmplx(5.942309446306148759408113184954627973488e-2_qp,  &
                        6.849420755113055630708881643127550816161e-3_qp,qp)
  values(22, 3) = cmplx(6.691050899620643688228945106331374955089e-2_qp,  &
                        8.597344352930658821354023270964033750455e-3_qp,qp)
  values(23, 3) = cmplx(7.637080737603776837681860111059313334652e-2_qp,  &
                        1.103876669193627222748017069037384501797e-2_qp,qp)
  values(24, 3) = cmplx(8.860505248855702922182109289859970176364e-2_qp,  &
                        1.454649951401545639330748229333367506877e-2_qp,qp)
  values(25, 3) = cmplx(1.048515513880954377342626034899559117584e-1_qp,  &
                        1.974317560234086292413485404635747192613e-2_qp,qp)
  values(26, 3) = cmplx(1.270814964998812847681790040646785324910e-1_qp,  &
                        2.769844132667481445377740839288618334462e-2_qp,qp)
  values(27, 3) = cmplx(1.585290151921034933474976783697010003774e-1_qp,  &
                        4.030230586813971740093660744297660373231e-2_qp,qp)
  values(28, 3) = cmplx(2.046653368157393085173010773434743869071e-1_qp,  &
                        6.098769322631811291321512399322348242936e-2_qp,qp)
  values(29, 3) = cmplx(2.750327002634415447983579547354789279062e-1_qp,  &
                        9.615830473505560573545640099837139214946e-2_qp,qp)
  values(30, 3) = cmplx(3.867770243694656811069120391967471161333e-1_qp,  &
                        1.580663278635548109070911403752639611376e-1_qp,qp)
  values(31, 3) = cmplx(5.716026882338550835288023695379615612573e-1_qp,  &
                        2.707198601939996389731041843300791724523e-1_qp,qp)
  values(32, 3) = cmplx(8.897470470107766989169844915564923882352e-1_qp,  &
                        4.822189606336244406956427774132631833772e-1_qp,qp)
  values( 1, 4) = cmplx(1.250031616401452135136907693163385197980e-2_qp,  &
                       -9.220249844578000920696921013212793002356e-4_qp,qp)
  values( 2, 4) = cmplx(1.348265356799832428328805542687015403238e-2_qp,  &
                       -1.069258194251352924678683845822480691025e-3_qp,qp)
  values( 3, 4) = cmplx(1.462616970396498670782315954911032668111e-2_qp,  &
                       -1.253247258800955357181911098558507137682e-3_qp,qp)
  values( 4, 4) = cmplx(1.597205450349779558163977021630956201460e-2_qp,  &
                       -1.486658084242546322345427464933251909609e-3_qp,qp)
  values( 5, 4) = cmplx(1.757596902639849000178978623096891547045e-2_qp,  &
                       -1.787758045214239061846817431306894109049e-3_qp,qp)
  values( 6, 4) = cmplx(1.951454356558583875935538851106648535881e-2_qp,  &
                       -2.183439842531035987600273048020490321672e-3_qp,qp)
  values( 7, 4) = cmplx(2.189542375614422737569837840674517519403e-2_qp,  &
                       -2.714164266051988787054551929092832544015e-3_qp,qp)
  values( 8, 4) = cmplx(2.487313420483443687278442035975345697778e-2_qp,  &
                       -3.442211563606326826492312688806593969669e-3_qp,qp)
  values( 9, 4) = cmplx(2.867468123189666415778859607995579634851e-2_qp,  &
                       -4.465752814777900616826871016799162211191e-3_qp,qp)
  values(10, 4) = cmplx(3.364180574673009817613253549743715876013e-2_qp,  &
                       -5.943364420055283722355127104550975415511e-3_qp,qp)
  values(11, 4) = cmplx(4.030230586813971740093660744297660373231e-2_qp,  &
                       -8.137651474563173319168988296965666289229e-3_qp,qp)
  values(12, 4) = cmplx(4.949318168769537738192476733501560133492e-2_qp,  &
                       -1.149451153862273553129035665820788109443e-2_qp,qp)
  values(13, 4) = cmplx(6.257807438572107239976779542719918292572e-2_qp,  &
                       -1.679011517466726666784430278558610461187e-2_qp,qp)
  values(14, 4) = cmplx(8.183974009719518542278272579655053095375e-2_qp,  &
                       -2.540886258878654182810280485957114339461e-2_qp,qp)
  values(15, 4) = cmplx(1.112037615566325474365674703420740441656e-1_qp,  &
                       -3.987902465934177288413417849700128207168e-2_qp,qp)
  values(16, 4) = cmplx(1.576008023982374846902781500716304535084e-1_qp,  &
                       -6.492363164707739120107292546832654597374e-2_qp,qp)
  values(17, 4) = cmplx(1.250031616401452135136907693163385197980e-2_qp,  &
                        9.220249844578000920696921013212793002356e-4_qp,qp)
  values(18, 4) = cmplx(1.348265356799832428328805542687015403238e-2_qp,  &
                        1.069258194251352924678683845822480691025e-3_qp,qp)
  values(19, 4) = cmplx(1.462616970396498670782315954911032668111e-2_qp,  &
                        1.253247258800955357181911098558507137682e-3_qp,qp)
  values(20, 4) = cmplx(1.597205450349779558163977021630956201460e-2_qp,  &
                        1.486658084242546322345427464933251909609e-3_qp,qp)
  values(21, 4) = cmplx(1.757596902639849000178978623096891547045e-2_qp,  &
                        1.787758045214239061846817431306894109049e-3_qp,qp)
  values(22, 4) = cmplx(1.951454356558583875935538851106648535881e-2_qp,  &
                        2.183439842531035987600273048020490321672e-3_qp,qp)
  values(23, 4) = cmplx(2.189542375614422737569837840674517519403e-2_qp,  &
                        2.714164266051988787054551929092832544015e-3_qp,qp)
  values(24, 4) = cmplx(2.487313420483443687278442035975345697778e-2_qp,  &
                        3.442211563606326826492312688806593969669e-3_qp,qp)
  values(25, 4) = cmplx(2.867468123189666415778859607995579634851e-2_qp,  &
                        4.465752814777900616826871016799162211191e-3_qp,qp)
  values(26, 4) = cmplx(3.364180574673009817613253549743715876013e-2_qp,  &
                        5.943364420055283722355127104550975415511e-3_qp,qp)
  values(27, 4) = cmplx(4.030230586813971740093660744297660373231e-2_qp,  &
                        8.137651474563173319168988296965666289229e-3_qp,qp)
  values(28, 4) = cmplx(4.949318168769537738192476733501560133492e-2_qp,  &
                        1.149451153862273553129035665820788109443e-2_qp,qp)
  values(29, 4) = cmplx(6.257807438572107239976779542719918292572e-2_qp,  &
                        1.679011517466726666784430278558610461187e-2_qp,qp)
  values(30, 4) = cmplx(8.183974009719518542278272579655053095375e-2_qp,  &
                        2.540886258878654182810280485957114339461e-2_qp,qp)
  values(31, 4) = cmplx(1.112037615566325474365674703420740441656e-1_qp,  &
                        3.987902465934177288413417849700128207168e-2_qp,qp)
  values(32, 4) = cmplx(1.576008023982374846902781500716304535084e-1_qp,  &
                        6.492363164707739120107292546832654597374e-2_qp,qp)
  values( 1, 5) = cmplx(2.896886435752269834109362271798509169988e-3_qp,  &
                       -1.974861451294469742039670170477229869753e-4_qp,qp)
  values( 2, 5) = cmplx(3.106407025393493102135197378097452370726e-3_qp,  &
                       -2.263498701269044641618348369194412977445e-4_qp,qp)
  values( 3, 5) = cmplx(3.347341891698683000429691815984757338801e-3_qp,  &
                       -2.617618291122159554059725896782812751398e-4_qp,qp)
  values( 4, 5) = cmplx(3.626978864528492878350674052348659689480e-3_qp,  &
                       -3.057601114694209365721780839164868256958e-4_qp,qp)
  values( 5, 5) = cmplx(3.954917402346575650029948649878200034765e-3_qp,  &
                       -3.611932261887227646971885364285509876193e-4_qp,qp)
  values( 6, 5) = cmplx(4.344002128766737520159871685616207571573e-3_qp,  &
                       -4.321124572471403065119197275191434499801e-4_qp,qp)
  values( 7, 5) = cmplx(4.811713876949514467701629704045811672620e-3_qp,  &
                       -5.243874027243814201617694437382447821514e-4_qp,qp)
  values( 8, 5) = cmplx(5.382280894910301620813905160954622303632e-3_qp,  &
                       -6.466897771013249314405308240493427779874e-4_qp,qp)
  values( 9, 5) = cmplx(6.089944736863581126916602438044180569499e-3_qp,  &
                       -8.120959610428402550448657106225091791540e-4_qp,qp)
  values(10, 5) = cmplx(6.984112669995926106444629136890241661021e-3_qp,  &
                       -1.040748249940642384089502032339266245509e-3_qp,qp)
  values(11, 5) = cmplx(8.137651474563173319168988296965666289229e-3_qp,  &
                       -1.364360798526949265730059223690062934356e-3_qp,qp)
  values(12, 5) = cmplx(9.660513279825723123274228663278407881348e-3_qp,  &
                       -1.833998258797012408016127994929473213087e-3_qp,qp)
  values(13, 5) = cmplx(1.172258612255521562680931206133022742599e-2_qp,  &
                       -2.533764526056025520517495362127938592938e-3_qp,qp)
  values(14, 5) = cmplx(1.459280828803720980964509822492227362558e-2_qp,  &
                       -3.605351433583110672819235544882956589676e-3_qp,qp)
  values(15, 5) = cmplx(1.870749436583560564492572901168416423925e-2_qp,  &
                       -5.292882573376541809802112371329279458106e-3_qp,qp)
  values(16, 5) = cmplx(2.479208885788198005073578240358251846856e-2_qp,  &
                       -8.026308557839082230067428612948805501036e-3_qp,qp)
  values(17, 5) = cmplx(2.896886435752269834109362271798509169988e-3_qp,  &
                        1.974861451294469742039670170477229869753e-4_qp,qp)
  values(18, 5) = cmplx(3.106407025393493102135197378097452370726e-3_qp,  &
                        2.263498701269044641618348369194412977445e-4_qp,qp)
  values(19, 5) = cmplx(3.347341891698683000429691815984757338801e-3_qp,  &
                        2.617618291122159554059725896782812751398e-4_qp,qp)
  values(20, 5) = cmplx(3.626978864528492878350674052348659689480e-3_qp,  &
                        3.057601114694209365721780839164868256958e-4_qp,qp)
  values(21, 5) = cmplx(3.954917402346575650029948649878200034765e-3_qp,  &
                        3.611932261887227646971885364285509876193e-4_qp,qp)
  values(22, 5) = cmplx(4.344002128766737520159871685616207571573e-3_qp,  &
                        4.321124572471403065119197275191434499801e-4_qp,qp)
  values(23, 5) = cmplx(4.811713876949514467701629704045811672620e-3_qp,  &
                        5.243874027243814201617694437382447821514e-4_qp,qp)
  values(24, 5) = cmplx(5.382280894910301620813905160954622303632e-3_qp,  &
                        6.466897771013249314405308240493427779874e-4_qp,qp)
  values(25, 5) = cmplx(6.089944736863581126916602438044180569499e-3_qp,  &
                        8.120959610428402550448657106225091791540e-4_qp,qp)
  values(26, 5) = cmplx(6.984112669995926106444629136890241661021e-3_qp,  &
                        1.040748249940642384089502032339266245509e-3_qp,qp)
  values(27, 5) = cmplx(8.137651474563173319168988296965666289229e-3_qp,  &
                        1.364360798526949265730059223690062934356e-3_qp,qp)
  values(28, 5) = cmplx(9.660513279825723123274228663278407881348e-3_qp,  &
                        1.833998258797012408016127994929473213087e-3_qp,qp)
  values(29, 5) = cmplx(1.172258612255521562680931206133022742599e-2_qp,  &
                        2.533764526056025520517495362127938592938e-3_qp,qp)
  values(30, 5) = cmplx(1.459280828803720980964509822492227362558e-2_qp,  &
                        3.605351433583110672819235544882956589676e-3_qp,qp)
  values(31, 5) = cmplx(1.870749436583560564492572901168416423925e-2_qp,  &
                        5.292882573376541809802112371329279458106e-3_qp,qp)
  values(32, 5) = cmplx(2.479208885788198005073578240358251846856e-2_qp,  &
                        8.026308557839082230067428612948805501036e-3_qp,qp)
end subroutine qp_compref

#endif

end module selftestall

program selftestallmain
  use floattypes
  use selftestall

  print '(a)','# This is a post-installation test of this library. All tests'
  print '(a)','# evaluate relative error of numerical results, which should be'
  print '(a)','# less than 100~1000 times the unit roundoff. The test coverage'
  print '(a)','# of this program is 53.51%. Of course, more extensive tests'
  print '(a)','# whose coverage are as wide as 98% is performed on author''s'
  print '(a)','# computers.'
  print '(a)',''
  print '(a)','# The assumed unit roundoffs are'

  print '(a,e10.3)', 'sp (single    precision) = ', 1/2.0_wp**24
  print '(a,e10.3)', 'wp (double    precision) = ', 1/2.0_wp**53
#ifdef __USE_TPREC
  print '(a,e10.3)', 'tp (triple    precision) = ', 1/2.0_wp**64
#endif
#ifdef __USE_QPREC
  print '(a,e10.3)', 'qp (quadruple precision) = ', 1/2.0_wp**113
#endif

  print '(a)','# These assumed values should coincide respectively the unit'
  print '(a)','# roundoffs available with this compiler shown in the following.'

  print '(a,e10.3)', 'sp (single    precision) = ', epsilon(1.0_sp) / 2
  print '(a,e10.3)', 'wp (double    precision) = ', epsilon(1.0_wp) / 2
#ifdef __USE_TPREC
  print '(a,e10.3)', 'tp (triple    precision) = ', epsilon(1.0_tp) / 2
#endif
#ifdef __USE_QPREC
  print '(a,e10.3)', 'qp (quadruple precision) = ', epsilon(1.0_qp) / 2
#endif

  print *, ''
  print '(a)', '# [Test 1] Relative error of the scaling and squaring method'
  print '(a)', '# for diagonal matrices. Numerical results are compared with'
  print '(a)', '# exact values precomputed by arbitrary precision arithmetic.'
  print '(a)', '#'
  print '(a11,6a10)','# type&prec','phi0','phi1','phi2','phi3','phi4','phi5'

  call sp_test_realdiag
  call sp_test_compdiag
  call wp_test_realdiag
  call wp_test_compdiag
#ifdef __USE_TPREC
  call tp_test_realdiag
  call tp_test_compdiag
#endif
#ifdef __USE_QPREC
  call qp_test_realdiag
  call qp_test_compdiag
#endif

  print *, ''
  print '(a)', '# [Test 2] Relative difference of the following two results,'
  print '(a)', '# - phi-functions obtained as S phi(D) S^{-1},'
  print '(a)', '# - phi-functions obtained as phi(SDS^{-1}),'
  print '(a)', '# where D is a diagonal matrix. The scaling and squaring method'
  print '(a)', '# for diagonal matrices is compared with that for square ones.'
  print '(a)', '#'
  print '(a11,6a10)','# type&prec','phi0','phi1','phi2','phi3','phi4','phi5'

  call rs_test_sands4square
  call cs_test_sands4square
  call rw_test_sands4square
  call cw_test_sands4square
#ifdef __USE_TPREC
  call rt_test_sands4square
  call ct_test_sands4square
#endif
#ifdef __USE_QPREC
  call rq_test_sands4square
  call cq_test_sands4square
#endif

  print *, ''
  print '(a)', '# [Test 3] Relative difference of the following two functions:'
  print '(a)', '# - the scaling and squaring method for square     matrices,'
  print '(a)', '# - the scaling and squaring method for triangular matrices.'
  print '(a)', '# The same triangular matrix is given to these two functions.'
  print '(a)', '#'
  print '(a11,6a10)','# type&prec','phi0','phi1','phi2','phi3','phi4','phi5'

  call rs_test_sands4triangular
  call cs_test_sands4triangular
  call rw_test_sands4triangular
  call cw_test_sands4triangular
#ifdef __USE_TPREC
  call rt_test_sands4triangular
  call ct_test_sands4triangular
#endif
#ifdef __USE_QPREC
  call rq_test_sands4triangular
  call cq_test_sands4triangular
#endif

  print *, ''
  print '(a)', '# [Test 4] Relative difference of the following two functions:'
  print '(a)', '# - the scaling and squaring method for triangular matrices,'
  print '(a)', '# - the modified block Schur--Parlett algorithm.'
  print '(a)', '# The same triangular matrix is given to these two functions.'
  print '(a)', '#'
  print '(a11,6a10)','# type&prec','phi0','phi1','phi2','phi3','phi4','phi5'

  call rs_test_blksp4triangular
  call cs_test_blksp4triangular
  call rw_test_blksp4triangular
  call cw_test_blksp4triangular
#ifdef __USE_TPREC
  call rt_test_blksp4triangular
  call ct_test_blksp4triangular
#endif
#ifdef __USE_QPREC
  call rq_test_blksp4triangular
  call cq_test_blksp4triangular
#endif

  print *, ''
  print '(a)', '# The component-wise relative error may be much larger than the'
  print '(a)', '# the above values. Roughly speaking, the printed errors are'
  print '(a)', '# measured relative to the magnitude of exp(mu(A)), where mu is'
  print '(a)', '# the logarithmic norm.'

  print *, ''
  print '(a)', '# Thank you very much for your interest on this library.'

end program selftestallmain

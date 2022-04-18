! This module solves the continuous-time Sylvester equation, AX - XB = C,
! defined by quasi-triangular matrices, A and B.
!
! The product, ``matmul(A,X)'', is computed simultaneously
! without additional cost.
!
! If you want to use the library ``RECSY'', please add a phrase,
! ``-Dpure= -D__USE_RECSY'', literally to the command line arguments of
! your compiler. The pure attribute is lost while an interface to RECSY
! is enabled.
!
! References:
!
!   R. H. Bartels and G. W. Stewart,
!   Algorithm 432 : Solution of the matrix equation AX+XB=C,
!   Communications of the ACM,
!   Volume 15, Number 9, pp. 820-826, 1972.
!
!   Bo K{\aa}gstr{\"{o}}m and P. Poromaa,
!   Distributed and shared memory block algorithms for the triangular
!   Sylvester equations with sep^{-1} estimators,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 13, Number 1, pp. 90-101, 1992.
!
!   Isak Jonsson and Bo K{\aa}gstr{\"{o}}m,
!   Recursive blocked algorithms for solving triangular systems. I.
!   One-sided and coupled Sylvester-type matrix equations,
!   ACM Transactions on Mathematical Software",
!   Volume 28, Number 4, pp. 392-415, 2002.
!
! The module is intended for internal use only.

module jspsylvester

  use floattypes
  use blasinterface
  use matrixpwrtag

  implicit none


  !    Name : axuppersylv ( a generic name )
  !   Usage : call axuppersylv(A,B,C,X,AX)
  ! Purpose : This routine solves the continuous-time Sylvester equation,
  !         : AX - XB = C. The product AX is computed simultaneously.
  !   Input : A and B are upper quasi-triangular matrices.
  !         : C is a rectangular matrix of consistent size.
  !  Output : X and AX are matrices that have the same shape as C.
  !         : The solution is stored in X, and the product, ``matmul(A,X)'',
  !         : is stored in AX.

  interface axuppersylv
    module procedure rs_axuppersylv            ! Real    in Single    precision
    module procedure cs_axuppersylv            ! Complex in Single    precision
    module procedure rw_axuppersylv            ! Real    in Double    precision
    module procedure cw_axuppersylv            ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_axuppersylv            ! Real    in Triple    precision
    module procedure ct_axuppersylv            ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_axuppersylv            ! Real    in Quadruple precision
    module procedure cq_axuppersylv            ! Complex in Quadruple precision
#endif
  end interface axuppersylv


  !    Name : axlowersylv ( a generic name )
  !   Usage : call axlowersylv(A,B,C,X,AX)
  ! Purpose : This routine solves the continuous-time Sylvester equation,
  !         : AX - XB = C. The product AX is computed simultaneously.
  !   Input : A and B are lower quasi-triangular matrices.
  !         : C is a rectangular matrix of consistent size.
  !  Output : X and AX are matrices that have the same shape as C.
  !         : The solution is stored in X, and the product, ``matmul(A,X)'',
  !         : is stored in AX.

  interface axlowersylv
    module procedure rs_axlowersylv            ! Real    in Single    precision
    module procedure cs_axlowersylv            ! Complex in Single    precision
    module procedure rw_axlowersylv            ! Real    in Double    precision
    module procedure cw_axlowersylv            ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_axlowersylv            ! Real    in Triple    precision
    module procedure ct_axlowersylv            ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_axlowersylv            ! Real    in Quadruple precision
    module procedure cq_axlowersylv            ! Complex in Quadruple precision
#endif
  end interface axlowersylv

  logical :: axsylvuseblock = .true.

  private
  public  :: axuppersylv, axlowersylv, axsylvuseblock

contains


!    Name : rs_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rs_axuppersylv(a,b,c,x,ax)
  real   (kind=sp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=sp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=sp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=sp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=sp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=sp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_sp; ax = 0.0_sp
    call rs_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_sp
    x = 0.0_sp; ax = 0.0_sp

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rs_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rs_axuppersylv


!    Name : rs_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rs_uppersolvlock(a,b,c,x,ax)
  real   (kind=sp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=sp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=sp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=sp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=sp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=sp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=sp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=sp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = 0.0_sp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = 0.0_sp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = 0.0_sp
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  0.0_sp
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  0.0_sp
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  0.0_sp
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  0.0_sp
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = 0.0_sp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=sp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=sp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=sp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=sp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=sp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rs_uppersolvlock


!    Name : rs_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rs_axlowersylv(a,b,c,x,ax)
  real   (kind=sp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=sp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=sp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=sp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=sp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=sp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_sp; ax = 0.0_sp
    call rs_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_sp
    x = 0.0_sp; ax = 0.0_sp

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rs_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rs_axlowersylv


!    Name : rs_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rs_lowersolvlock(a,b,c,x,ax)
  real   (kind=sp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=sp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=sp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=sp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=sp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=sp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=sp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=sp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = 0.0_sp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = 0.0_sp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = 0.0_sp
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  0.0_sp
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  0.0_sp
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  0.0_sp
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  0.0_sp
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = 0.0_sp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=sp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=sp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=sp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=sp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=sp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rs_lowersolvlock


!    Name : cs_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cs_axuppersylv(a,b,c,x,ax)
  complex(kind=sp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=sp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=sp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=sp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=sp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=sp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_sp,0.0_sp,sp); ax = cmplx(0.0_sp,0.0_sp,sp)
    call cs_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_sp,0.0_sp,sp)
    x = cmplx(0.0_sp,0.0_sp,sp); ax = cmplx(0.0_sp,0.0_sp,sp)

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cs_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cs_axuppersylv


!    Name : cs_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cs_uppersolvlock(a,b,c,x,ax)
  complex(kind=sp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=sp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=sp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=sp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=sp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=sp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=sp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=sp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_sp,0.0_sp,sp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=sp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=sp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=sp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=sp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=sp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cs_uppersolvlock


!    Name : cs_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cs_axlowersylv(a,b,c,x,ax)
  complex(kind=sp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=sp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=sp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=sp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=sp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=sp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_sp,0.0_sp,sp); ax = cmplx(0.0_sp,0.0_sp,sp)
    call cs_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_sp,0.0_sp,sp)
    x = cmplx(0.0_sp,0.0_sp,sp); ax = cmplx(0.0_sp,0.0_sp,sp)

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cs_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cs_axlowersylv


!    Name : cs_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cs_lowersolvlock(a,b,c,x,ax)
  complex(kind=sp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=sp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=sp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=sp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=sp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=sp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=sp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=sp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_sp,0.0_sp,sp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  cmplx(0.0_sp,0.0_sp,sp)
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = cmplx(0.0_sp,0.0_sp,sp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=sp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=sp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=sp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=sp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=sp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cs_lowersolvlock


!    Name : rw_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rw_axuppersylv(a,b,c,x,ax)
  real   (kind=wp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=wp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=wp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=wp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=wp), intent(out) :: ax(1:,1:)        ! m by n

#ifdef __USE_RECSY

! An interface to RECSY developed by Bo Kagstrom.

  integer       :: m, n, info
  real(kind=wp) :: machine(16)
  external      :: recsyct

  m = size(a,1); n = size(b,1)
  machine = 0.0_wp

  x = c
  call recsyct(1,1.0_wp,m,n,a,m,b,n,x,m,info,machine)
  ax = a .sqtimes. x

#else

  real   (kind=wp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_wp; ax = 0.0_wp
    call rw_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_wp
    x = 0.0_wp; ax = 0.0_wp

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rw_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif

#endif

end subroutine rw_axuppersylv


!    Name : rw_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rw_uppersolvlock(a,b,c,x,ax)
  real   (kind=wp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=wp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=wp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=wp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=wp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=wp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=wp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=wp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = 0.0_wp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = 0.0_wp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = 0.0_wp
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  0.0_wp
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  0.0_wp
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  0.0_wp
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  0.0_wp
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = 0.0_wp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=wp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=wp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=wp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=wp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=wp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rw_uppersolvlock


!    Name : rw_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rw_axlowersylv(a,b,c,x,ax)
  real   (kind=wp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=wp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=wp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=wp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=wp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=wp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_wp; ax = 0.0_wp
    call rw_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_wp
    x = 0.0_wp; ax = 0.0_wp

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rw_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rw_axlowersylv


!    Name : rw_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rw_lowersolvlock(a,b,c,x,ax)
  real   (kind=wp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=wp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=wp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=wp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=wp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=wp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=wp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=wp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = 0.0_wp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = 0.0_wp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = 0.0_wp
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  0.0_wp
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  0.0_wp
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  0.0_wp
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  0.0_wp
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = 0.0_wp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=wp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=wp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=wp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=wp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=wp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rw_lowersolvlock


!    Name : cw_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cw_axuppersylv(a,b,c,x,ax)
  complex(kind=wp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=wp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=wp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=wp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=wp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=wp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_wp,0.0_wp,wp); ax = cmplx(0.0_wp,0.0_wp,wp)
    call cw_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_wp,0.0_wp,wp)
    x = cmplx(0.0_wp,0.0_wp,wp); ax = cmplx(0.0_wp,0.0_wp,wp)

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cw_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cw_axuppersylv


!    Name : cw_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cw_uppersolvlock(a,b,c,x,ax)
  complex(kind=wp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=wp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=wp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=wp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=wp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=wp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=wp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=wp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_wp,0.0_wp,wp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=wp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=wp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=wp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=wp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=wp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cw_uppersolvlock


!    Name : cw_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cw_axlowersylv(a,b,c,x,ax)
  complex(kind=wp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=wp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=wp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=wp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=wp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=wp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_wp,0.0_wp,wp); ax = cmplx(0.0_wp,0.0_wp,wp)
    call cw_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_wp,0.0_wp,wp)
    x = cmplx(0.0_wp,0.0_wp,wp); ax = cmplx(0.0_wp,0.0_wp,wp)

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cw_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cw_axlowersylv


!    Name : cw_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cw_lowersolvlock(a,b,c,x,ax)
  complex(kind=wp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=wp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=wp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=wp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=wp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=wp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=wp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=wp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_wp,0.0_wp,wp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  cmplx(0.0_wp,0.0_wp,wp)
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = cmplx(0.0_wp,0.0_wp,wp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=wp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=wp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=wp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=wp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=wp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cw_lowersolvlock

#ifdef __USE_TPREC

!    Name : rt_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rt_axuppersylv(a,b,c,x,ax)
  real   (kind=tp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=tp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=tp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=tp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=tp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=tp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_tp; ax = 0.0_tp
    call rt_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_tp
    x = 0.0_tp; ax = 0.0_tp

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rt_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rt_axuppersylv


!    Name : rt_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rt_uppersolvlock(a,b,c,x,ax)
  real   (kind=tp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=tp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=tp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=tp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=tp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=tp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=tp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=tp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = 0.0_tp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = 0.0_tp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = 0.0_tp
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  0.0_tp
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  0.0_tp
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  0.0_tp
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  0.0_tp
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = 0.0_tp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=tp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=tp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=tp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=tp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=tp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rt_uppersolvlock


!    Name : rt_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rt_axlowersylv(a,b,c,x,ax)
  real   (kind=tp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=tp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=tp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=tp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=tp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=tp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_tp; ax = 0.0_tp
    call rt_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_tp
    x = 0.0_tp; ax = 0.0_tp

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rt_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rt_axlowersylv


!    Name : rt_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rt_lowersolvlock(a,b,c,x,ax)
  real   (kind=tp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=tp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=tp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=tp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=tp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=tp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=tp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=tp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = 0.0_tp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = 0.0_tp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = 0.0_tp
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  0.0_tp
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  0.0_tp
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  0.0_tp
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  0.0_tp
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = 0.0_tp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=tp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=tp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=tp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=tp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=tp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rt_lowersolvlock


!    Name : ct_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine ct_axuppersylv(a,b,c,x,ax)
  complex(kind=tp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=tp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=tp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=tp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=tp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=tp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_tp,0.0_tp,tp); ax = cmplx(0.0_tp,0.0_tp,tp)
    call ct_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_tp,0.0_tp,tp)
    x = cmplx(0.0_tp,0.0_tp,tp); ax = cmplx(0.0_tp,0.0_tp,tp)

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call ct_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine ct_axuppersylv


!    Name : ct_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine ct_uppersolvlock(a,b,c,x,ax)
  complex(kind=tp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=tp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=tp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=tp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=tp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=tp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=tp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=tp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_tp,0.0_tp,tp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=tp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=tp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=tp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=tp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=tp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine ct_uppersolvlock


!    Name : ct_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine ct_axlowersylv(a,b,c,x,ax)
  complex(kind=tp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=tp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=tp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=tp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=tp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=tp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_tp,0.0_tp,tp); ax = cmplx(0.0_tp,0.0_tp,tp)
    call ct_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_tp,0.0_tp,tp)
    x = cmplx(0.0_tp,0.0_tp,tp); ax = cmplx(0.0_tp,0.0_tp,tp)

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call ct_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine ct_axlowersylv


!    Name : ct_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine ct_lowersolvlock(a,b,c,x,ax)
  complex(kind=tp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=tp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=tp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=tp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=tp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=tp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=tp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=tp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_tp,0.0_tp,tp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  cmplx(0.0_tp,0.0_tp,tp)
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = cmplx(0.0_tp,0.0_tp,tp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=tp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=tp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=tp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=tp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=tp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine ct_lowersolvlock

#endif
#ifdef __USE_QPREC

!    Name : rq_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rq_axuppersylv(a,b,c,x,ax)
  real   (kind=qp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=qp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=qp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=qp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=qp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=qp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_qp; ax = 0.0_qp
    call rq_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_qp
    x = 0.0_qp; ax = 0.0_qp

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rq_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rq_axuppersylv


!    Name : rq_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rq_uppersolvlock(a,b,c,x,ax)
  real   (kind=qp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=qp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=qp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=qp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=qp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=qp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=qp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=qp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = 0.0_qp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = 0.0_qp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = 0.0_qp
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  0.0_qp
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  0.0_qp
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  0.0_qp
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  0.0_qp
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = 0.0_qp
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=qp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=qp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=qp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=qp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=qp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rq_uppersolvlock


!    Name : rq_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine rq_axlowersylv(a,b,c,x,ax)
  real   (kind=qp), intent(in ) ::  a(1:,1:)        ! m by m
  real   (kind=qp), intent(in ) ::  b(1:,1:)        ! n by n
  real   (kind=qp), intent(in ) ::  c(1:,1:)        ! m by n
  real   (kind=qp), intent(out) ::  x(1:,1:)        ! m by n
  real   (kind=qp), intent(out) :: ax(1:,1:)        ! m by n

  real   (kind=qp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = 0.0_qp; ax = 0.0_qp
    call rq_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = 0.0_qp
    x = 0.0_qp; ax = 0.0_qp

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call rq_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine rq_axlowersylv


!    Name : rq_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine rq_lowersolvlock(a,b,c,x,ax)
  real   (kind=qp), intent(in   ) ::  a(1:,1:)      ! m by m
  real   (kind=qp), intent(in   ) ::  b(1:,1:)      ! n by n
  real   (kind=qp), intent(in   ) ::  c(1:,1:)      ! m by n
  real   (kind=qp), intent(out  ) ::  x(1:,1:)      ! m by n
  real   (kind=qp), intent(inout) :: ax(1:,1:)      ! m by n

  real   (kind=qp) :: tmp(4), m22(2,2), m44(4,4)
  real   (kind=qp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=qp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = 0.0_qp
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = 0.0_qp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = 0.0_qp
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  0.0_qp
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  0.0_qp
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  0.0_qp
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  0.0_qp
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = 0.0_qp
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    real   (kind=qp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    real   (kind=qp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    real   (kind=qp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    real   (kind=qp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    real   (kind=qp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine rq_lowersolvlock


!    Name : cq_axuppersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cq_axuppersylv(a,b,c,x,ax)
  complex(kind=qp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=qp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=qp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=qp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=qp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=qp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_qp,0.0_qp,qp); ax = cmplx(0.0_qp,0.0_qp,qp)
    call cq_uppersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i+1,i))) <= abs(a(i+1,i))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j+1,j))) <= abs(b(j+1,j))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_qp,0.0_qp,qp)
    x = cmplx(0.0_qp,0.0_qp,qp); ax = cmplx(0.0_qp,0.0_qp,qp)

    do j=1,nblk
      hj = hb(j); tj = tb(j)
      do i=mblk,1,-1
        hi = ha(i); ti = ta(i)
        do k=1,j-1
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=i+1,mblk
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cq_uppersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cq_axuppersylv


!    Name : cq_uppersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are upper quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cq_uppersolvlock(a,b,c,x,ax)
  complex(kind=qp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=qp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=qp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=qp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=qp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=qp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=qp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=qp) :: eps
  integer          :: i, j, m, n, p4(4), subd, sube

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_qp,0.0_qp,qp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in ascending order
  ! of columns and in descending order of rows.

  j = 1
  do while (j <= n)

    subd = 0
    if (j < n) then
      if (eps <= abs(b(j+1,j))) subd = 1
    end if

    if (j == n .or. subd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (1 < j) then
            x(i,j) = x(i,j) + sum(x(i,1:j-1) * b(1:j-1,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i  ,i  ) - b(j,j)
          m22(1,2) = -a(i-1,i  )
          m22(2,1) = -a(i  ,i-1)
          m22(2,2) =  a(i-1,i-1) - b(j,j)

          det = (a(i-1,i-1)-b(j,j)) * (a(i,i)-b(j,j)) - a(i-1,i) * a(i,i-1)
          m22 = m22 / det

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j))
            ax(i-1:i,j) = ax(i-1:i,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i-1:i,j) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j))
            v2(2,1) = v2(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j))
          end if

           x(i-1,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i  ,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)
          ax(i-1,j) = a(i-1,i-1)*x(i-1,j) + a(i-1,i)*x(i,j) + ax(i-1,j)
          ax(i  ,j) = a(i  ,i-1)*x(i-1,j) + a(i  ,i)*x(i,j) + ax(i  ,j)

          i = i - 2
        end if
      end do
      j = j + 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

      i = m
      do while (1 <= i)

        sube = 0
        if (1 < i) then
          if (eps <= abs(a(i,i-1))) sube = 1
        end if

        if (i == 1 .or. sube == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m22(1,1) =  a(i,i) - b(j+1,j+1)
          m22(1,2) =  b(j+1,j  )
          m22(2,1) =  b(j  ,j+1)
          m22(2,2) =  a(i,i) - b(j  ,j  )

          det = (a(i,i)-b(j,j)) * (a(i,i)-b(j+1,j+1)) - b(j+1,j) * b(j,j+1)
          m22 = m22 / det

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (i < m) then
            tmp(1) = sum(a(i,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i,i+1:m) * x(i+1:m,j+1))
            ax(i,j:j+1) = ax(i,j:j+1) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j:j+1) - tmp(1:2)
          if (1 < j) then
            v2(1,1) = v2(1,1) + sum(x(i,1:j-1) * b(1:j-1,j  ))
            v2(2,1) = v2(2,1) + sum(x(i,1:j-1) * b(1:j-1,j+1))
          end if

           x(i,j  ) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j+1) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )
          ax(i,j+1) = a(i,i) * x(i,j+1) + ax(i,j+1)

          i = i - 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i-1:i,i-1:i), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j:j+1,j:j+1).

          m44(1,1) =  a(i-1,i-1) - b(j  ,j  )
          m44(1,2) =  a(i-1,i  )
          m44(1,3) = -b(j+1,j  )
          m44(1,4) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(2,1) =  a(i  ,i-1)
          m44(2,2) =  a(i  ,i  ) - b(j  ,j  )
          m44(2,3) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(2,4) = -b(j+1,j  )
          m44(3,1) = -b(j  ,j+1)
          m44(3,2) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(3,3) =  a(i-1,i-1) - b(j+1,j+1)
          m44(3,4) =  a(i-1,i  )
          m44(4,1) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(4,2) = -b(j  ,j+1)
          m44(4,3) =  a(i  ,i-1)
          m44(4,4) =  a(i  ,i  ) - b(j+1,j+1)

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (i < m) then
            tmp(1) = sum(a(i-1,i+1:m) * x(i+1:m,j  ))
            tmp(2) = sum(a(i  ,i+1:m) * x(i+1:m,j  ))
            tmp(3) = sum(a(i-1,i+1:m) * x(i+1:m,j+1))
            tmp(4) = sum(a(i  ,i+1:m) * x(i+1:m,j+1))

            ax(i-1,j  ) = ax(i-1,j  ) + tmp(1)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(2)
            ax(i-1,j+1) = ax(i-1,j+1) + tmp(3)
            ax(i  ,j+1) = ax(i  ,j+1) + tmp(4)
          end if

          v4(1,1) = c(i-1,j  ) - tmp(1)
          v4(2,1) = c(i  ,j  ) - tmp(2)
          v4(3,1) = c(i-1,j+1) - tmp(3)
          v4(4,1) = c(i  ,j+1) - tmp(4)

          if (1 < j) then
            v4(1,1) = v4(1,1) + sum(x(i-1,1:j-1) * b(1:j-1,j  ))
            v4(2,1) = v4(2,1) + sum(x(i  ,1:j-1) * b(1:j-1,j  ))
            v4(3,1) = v4(3,1) + sum(x(i-1,1:j-1) * b(1:j-1,j+1))
            v4(4,1) = v4(4,1) + sum(x(i  ,1:j-1) * b(1:j-1,j+1))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i-1,j  ) = w4(1,1)
          x(i  ,j  ) = w4(2,1)
          x(i-1,j+1) = w4(3,1)
          x(i  ,j+1) = w4(4,1)

          ax(i-1,j) = ax(i-1,j)+a(i-1,i-1)*x(i-1,j)+a(i-1,i)*x(i,j)
          ax(i  ,j) = ax(i  ,j)+a(i  ,i-1)*x(i-1,j)+a(i  ,i)*x(i,j)

          ax(i-1,j+1) = ax(i-1,j+1)+a(i-1,i-1)*x(i-1,j+1)+a(i-1,i)*x(i,j+1)
          ax(i  ,j+1) = ax(i  ,j+1)+a(i  ,i-1)*x(i-1,j+1)+a(i  ,i)*x(i,j+1)

          i = i - 2
        end if
      end do
      j = j + 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=qp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=qp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=qp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=qp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=qp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cq_uppersolvlock


!    Name : cq_axlowersylv
!   Usage : This subroutine is invoked through the generic name,``axuppersylv''.
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a block algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of ``TSYLV_B'' proposed by Kagstrom and
!         : Poromaa.

pure subroutine cq_axlowersylv(a,b,c,x,ax)
  complex(kind=qp), intent(in ) ::  a(1:,1:)        ! m by m
  complex(kind=qp), intent(in ) ::  b(1:,1:)        ! n by n
  complex(kind=qp), intent(in ) ::  c(1:,1:)        ! m by n
  complex(kind=qp), intent(out) ::  x(1:,1:)        ! m by n
  complex(kind=qp), intent(out) :: ax(1:,1:)        ! m by n

  complex(kind=qp)   :: d(1:size(a,1),1:size(b,1))
  integer, parameter :: blksize = mpt_blksize
  integer            :: i,j,k,m,n,mblk,nblk,acts,ell
  integer            :: hi,hj,hk,ti,tj,tk
  integer            :: amask(1:size(a,1)), bmask(1:size(b,1))
  integer            :: ha(1:size(a,1)),ta(1:size(a,1))
  integer            :: hb(1:size(b,1)),tb(1:size(b,1))

  m = size(a,1); n = size(b,1)

  if (m < blksize .or. n < blksize .or. (.not. axsylvuseblock)) then

    ! Bartels and Stewart's scalar algorithm is used for smaller matrices.

    x = cmplx(0.0_qp,0.0_qp,qp); ax = cmplx(0.0_qp,0.0_qp,qp)
    call cq_lowersolvlock(a,b,c,x,ax)

  else

    ! The block algorithm is used for larger matrices.

    ha = 0; ta = 0; hb = 0; tb = 0;

    ! The matrix A is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``ha'' stores the heads of blocks and
    ! ``ta'' stores the tails of blocks.

    amask = 0
    do i=1,m-1; if (tiny(abs(a(i,i+1))) <= abs(a(i,i+1))) amask(i) = 1; end do

    mblk = m / blksize
    if (mod(m,blksize) /= 0) mblk = mblk + 1
    acts = m / mblk

    ell = 0
    do i=1,mblk-1
      ha(i) = ell + 1
      ell   = ell + acts - amask(ell + acts)
      ta(i) = ell
    end do
    ha(mblk) = ell + 1; ta(mblk) = m

    ! The matrix B is block-decomposed so that each 2-by-2 diagonal block
    ! is not separated. The array ``hb'' stores the heads of blocks and
    ! ``tb'' stores the tails of blocks.

    bmask = 0
    do j=1,n-1; if (tiny(abs(b(j,j+1))) <= abs(b(j,j+1))) bmask(j) = 1; end do

    nblk = n / blksize
    if (mod(n,blksize) /= 0) nblk = nblk + 1
    acts = n / nblk

    ell = 0
    do i=1,nblk-1
      hb(i) = ell + 1
      ell   = ell + acts - bmask(ell + acts)
      tb(i) = ell
    end do
    hb(nblk) = ell + 1; tb(nblk) = n

    ! In the following loop, the (i,j)-th block of the solution X and the
    ! product AX are respectively computed by utilizing the level-3 BLAS
    ! and Bartels and Stewart's scalar algorithm.

    d = cmplx(0.0_qp,0.0_qp,qp)
    x = cmplx(0.0_qp,0.0_qp,qp); ax = cmplx(0.0_qp,0.0_qp,qp)

    do j=nblk,1,-1
      hj = hb(j); tj = tb(j)
      do i=1,mblk
        hi = ha(i); ti = ta(i)
        do k=j+1,nblk
          hk = hb(k); tk = tb(k)
#ifdef __USE_BLAS
          call ggemm(x(hi:ti,hk:tk),b(hk:tk,hj:tj),d(hi:ti,hj:tj),1,1)
#else
            d (hi:ti,hj:tj)                                                   &
          = d (hi:ti,hj:tj) + matmul(x(hi:ti,hk:tk),b(hk:tk,hj:tj))
#endif
        end do
        do k=1,i-1
          hk = ha(k); tk = ta(k)
#ifdef __USE_BLAS
          call ggemm(a(hi:ti,hk:tk),x(hk:tk,hj:tj),ax(hi:ti,hj:tj),1,1)
#else
            ax(hi:ti,hj:tj)                                                   &
          = ax(hi:ti,hj:tj) + matmul(a(hi:ti,hk:tk),x(hk:tk,hj:tj))
#endif
        end do
        d(hi:ti,hj:tj) = d(hi:ti,hj:tj) + c(hi:ti,hj:tj) - ax(hi:ti,hj:tj)
        call cq_lowersolvlock(                                                &
          a(hi:ti,hi:ti),b (hj:tj,hj:tj),d(hi:ti,hj:tj),                      &
          x(hi:ti,hj:tj),ax(hi:ti,hj:tj)                                      &
        )
      end do
    end do
  endif
end subroutine cq_axlowersylv


!    Name : cq_lowersolvlock
! Purpose : This routine solves the continuous-time Sylvester equation,
!         : AX - XB = C, by applying a scalar algorithm. The product,
!         : ``matmul(A,X)'', is computed simultaneously without additional cost.
!   Input : A and B are lower quasi-triangular matrices, and
!         : C is a rectangular matrix of consistent size.
!  Output : The matrix X stores the solution, and AX stores the product.
!    Note : This is an implementation of Bartels and Stewart's algorithm.

pure subroutine cq_lowersolvlock(a,b,c,x,ax)
  complex(kind=qp), intent(in   ) ::  a(1:,1:)      ! m by m
  complex(kind=qp), intent(in   ) ::  b(1:,1:)      ! n by n
  complex(kind=qp), intent(in   ) ::  c(1:,1:)      ! m by n
  complex(kind=qp), intent(out  ) ::  x(1:,1:)      ! m by n
  complex(kind=qp), intent(inout) :: ax(1:,1:)      ! m by n

  complex(kind=qp) :: tmp(4), m22(2,2), m44(4,4)
  complex(kind=qp) :: det, v2(2,1), v4(4,1), w4(4,1)
  real   (kind=qp) :: eps
  integer          :: i, j, m, n, p4(4), supd, supe

  m = size(a,1); n = size(b,1)
  x = cmplx(0.0_qp,0.0_qp,qp)
  eps = tiny(abs(a(1,1)))

  ! The solution X and the product AX are computed in descending order
  ! of columns and in ascending order of rows.

  j = n
  do while (1 <= j)

    supd = 0
    if (1 < j) then
      if (eps <= abs(b(j-1,j))) supd = 1
    end if

    if (j == 1 .or. supd == 0) then

      ! b_{jj} is included in a 1-by-1 diagonal block.

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! Both a_{ii} and b_{jj} are included in 1-by-1 diagonal blocks.

          x(i,j) = c(i,j)

          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j))
            ax(i,j) = ax(i,j) + tmp(1)
             x(i,j) =  x(i,j) - tmp(1)
          end if

          if (j < n) then
            x(i,j) = x(i,j) + sum(x(i,j+1:n) * b(j+1:n,j))
          end if

           x(i,j) =  x(i,j) / (a(i,i) - b(j,j))
          ax(i,j) = ax(i,j) + a(i,i) * x(i,j)

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), while
          ! b_{jj} is included in a 1-by-1 diagonal block.

          m22(1,1) =  a(i+1,i+1) - b(j,j)
          m22(1,2) = -a(i  ,i+1)
          m22(2,1) = -a(i+1,i  )
          m22(2,2) =  a(i  ,i  ) - b(j,j)

          det = (a(i,i)-b(j,j)) * (a(i+1,i+1)-b(j,j)) - a(i,i+1) * a(i+1,i)
          m22 = m22 / det

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j))
            ax(i:i+1,j) = ax(i:i+1,j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i:i+1,j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j))
            v2(2,1) = v2(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j))
          end if

           x(i  ,j) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i+1,j) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i  ,j) = a(i  ,i)*x(i,j) + a(i  ,i+1)*x(i+1,j) + ax(i  ,j)
          ax(i+1,j) = a(i+1,i)*x(i,j) + a(i+1,i+1)*x(i+1,j) + ax(i+1,j)

          i = i + 2
        end if
      end do
      j = j - 1

    else

      ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

      i = 1
      do while (i <= m)

        supe = 0
        if (i < m) then
          if (eps <= abs(a(i,i+1))) supe = 1
        end if

        if (i == m .or. supe == 0) then

          ! a_{ii} is included in a 1-by-1 diagonal block, while
          ! b_{jj} is included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m22(1,1) =  a(i,i) - b(j  ,j  )
          m22(1,2) =  b(j  ,j-1)
          m22(2,1) =  b(j-1,j  )
          m22(2,2) =  a(i,i) - b(j-1,j-1)

          det = (a(i,i)-b(j-1,j-1)) * (a(i,i)-b(j,j)) - b(j,j-1) * b(j-1,j)
          m22 = m22 / det

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (1 < i) then
            tmp(1) = sum(a(i,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i,1:i-1) * x(1:i-1,j  ))
            ax(i,j-1:j) = ax(i,j-1:j) + tmp(1:2)
          end if

          v2(1:2,1) = c(i,j-1:j) - tmp(1:2)
          if (j < n) then
            v2(1,1) = v2(1,1) + sum(x(i,j+1:n) * b(j+1:n,j-1))
            v2(2,1) = v2(2,1) + sum(x(i,j+1:n) * b(j+1:n,j  ))
          end if

           x(i,j-1) = m22(1,1)*v2(1,1) + m22(1,2)*v2(2,1)
           x(i,j  ) = m22(2,1)*v2(1,1) + m22(2,2)*v2(2,1)

          ax(i,j-1) = a(i,i) * x(i,j-1) + ax(i,j-1)
          ax(i,j  ) = a(i,i) * x(i,j  ) + ax(i,j  )

          i = i + 1

        else

          ! a_{ii} is included in a 2-by-2 diagonal block, (i:i+1,i:i+1), and
          ! b_{jj} is also included in a 2-by-2 diagonal block, (j-1:j,j-1:j).

          m44(1,1) =  a(i  ,i  ) - b(j-1,j-1)
          m44(1,2) =  a(i  ,i+1)
          m44(1,3) = -b(j  ,j-1)
          m44(1,4) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(2,1) =  a(i+1,i  )
          m44(2,2) =  a(i+1,i+1) - b(j-1,j-1)
          m44(2,3) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(2,4) = -b(j  ,j-1)
          m44(3,1) = -b(j-1,j  )
          m44(3,2) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(3,3) =  a(i  ,i  ) - b(j  ,j  )
          m44(3,4) =  a(i  ,i+1)
          m44(4,1) =  cmplx(0.0_qp,0.0_qp,qp)
          m44(4,2) = -b(j-1,j  )
          m44(4,3) =  a(i+1,i  )
          m44(4,4) =  a(i+1,i+1) - b(j  ,j  )

          tmp = cmplx(0.0_qp,0.0_qp,qp)
          if (1 < i) then
            tmp(1) = sum(a(i  ,1:i-1) * x(1:i-1,j-1))
            tmp(2) = sum(a(i+1,1:i-1) * x(1:i-1,j-1))
            tmp(3) = sum(a(i  ,1:i-1) * x(1:i-1,j  ))
            tmp(4) = sum(a(i+1,1:i-1) * x(1:i-1,j  ))

            ax(i  ,j-1) = ax(i  ,j-1) + tmp(1)
            ax(i+1,j-1) = ax(i+1,j-1) + tmp(2)
            ax(i  ,j  ) = ax(i  ,j  ) + tmp(3)
            ax(i+1,j  ) = ax(i+1,j  ) + tmp(4)
          end if

          v4(1,1) = c(i  ,j-1) - tmp(1)
          v4(2,1) = c(i+1,j-1) - tmp(2)
          v4(3,1) = c(i  ,j  ) - tmp(3)
          v4(4,1) = c(i+1,j  ) - tmp(4)

          if (j < n) then
            v4(1,1) = v4(1,1) + sum(x(i  ,j+1:n) * b(j+1:n,j-1))
            v4(2,1) = v4(2,1) + sum(x(i+1,j+1:n) * b(j+1:n,j-1))
            v4(3,1) = v4(3,1) + sum(x(i  ,j+1:n) * b(j+1:n,j  ))
            v4(4,1) = v4(4,1) + sum(x(i+1,j+1:n) * b(j+1:n,j  ))
          end if

          call ludcmp4x4(m44,p4)
          call ainvb4x4(m44,p4,v4(:,1),w4(:,1))

          x(i  ,j-1) = w4(1,1)
          x(i+1,j-1) = w4(2,1)
          x(i  ,j  ) = w4(3,1)
          x(i+1,j  ) = w4(4,1)

          ax(i  ,j-1) = ax(i  ,j-1)+a(i  ,i)*x(i,j-1)+a(i  ,i+1)*x(i+1,j-1)
          ax(i+1,j-1) = ax(i+1,j-1)+a(i+1,i)*x(i,j-1)+a(i+1,i+1)*x(i+1,j-1)

          ax(i  ,j  ) = ax(i  ,j  )+a(i  ,i)*x(i,j  )+a(i  ,i+1)*x(i+1,j  )
          ax(i+1,j  ) = ax(i+1,j  )+a(i+1,i)*x(i,j  )+a(i+1,i+1)*x(i+1,j  )

          i = i + 2
        end if
      end do
      j = j - 2
    end if
  end do

contains


  !    Name : ludcmp4x4
  ! Purpose : This subroutine performs the LU-factorization for a 4-by-4
  !         : square matrix with partial pivoting.
  !   Input : The array A stores a 4-by-4 matrix.
  !  Output : The factored form is stored in A, and the array p stores
  !         : the pivoting information.
  !    Note : The succeeding substitution is done by ``ainvb4x4'', which is
  !         : defined later in this subroutine.

  pure subroutine ludcmp4x4(a,p)
    complex(kind=qp), intent(inout) :: a(1:4,1:4)
    integer,          intent(out  ) :: p(1:4)

    integer, parameter :: n = 4
    complex(kind=qp)   :: u(n)
    integer            :: i,j,k,m(1)

    p = (/(i,i=1,n)/)

    do i=1,n-1                                   ! Gaussian elimination
      m = maxloc(abs(a(i:n,i))); j = i+m(1)-1    ! determine the pivot
      if (i /= j) then                           ! row permutation if necessary
        u = a(i,:); a(i,:) = a(j,:); a(j,:) = u
        k = p(i  ); p(i  ) = p(j  ); p(j  ) = k
      end if
      a(i+1:n,i) = a(i+1:n,i) / a(i,i)           ! the lower factor
      do k=i+1,n                                 ! update the rest
        a(i+1:n,k) = a(i+1:n,k) - a(i+1:n,i) * a(i,k)
      end do
    end do
  end subroutine ludcmp4x4


  !    Name : ainvb4x4
  ! Purpose : This subroutine performs forward and backward substitutions
  !         : for a 4-by-4 linear system.
  !   Input : The array A stores the factored 4-by-4 matrix,
  !         : and the integer array p stores the pivoting information.
  !         : The vector b is the right-hand side of the linear system.
  !  Output : The solution, A^{-1}b, is stored in c.
  !    Note : The LU-decomposition is assumed to be done by ``ludcmp4x4''.

  pure subroutine ainvb4x4(a,p,b,c)
    complex(kind=qp), intent(in ) :: a(1:4,1:4), b(1:4)
    integer,          intent(in ) :: p(1:4)
    complex(kind=qp), intent(out) :: c(1:4)

    integer, parameter :: m = 4
    integer            :: i
    complex(kind=qp)   :: d(1:4)

    do i=1,m
      c(i) = b(p(i))                                         ! undo permutation
    end do
    d(1) = c(1)
    do i=2,m
      d(i) = c(i) - sum(a(i,1:i-1) * d(1:i-1))
    end do

    c = d
    d(m) = c(m) / a(m,m)
    do i=m-1,1,-1
      d(i) = ( c(i) - sum(a(i,i+1:m)*d(i+1:m)) ) / a(i,i)
    end do

    c = d
  end subroutine ainvb4x4

end subroutine cq_lowersolvlock

#endif

end module jspsylvester


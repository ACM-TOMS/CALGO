! -- A Fully Portable High Performance Minimal --
! -- Storage Hybrid Format Cholesky Algorithm  --
!    F.G. Gustavson - IBM, USA;
!    J.K. Reid - RAL, UK;
!    J. Wasniewski - DTU, Denmark
!    March 20, 2006.

module block_hybrid_Cholesky

  implicit none

! Specific names of the BLAS called
  external dgemm, dtrsm, dsyrk, dcopy, dtpsv, dgemv

! Default block sizes
  integer, parameter :: default_nb = 100

  interface PpHpp
  ! Rearrange a packed triangular matrix to blocked hybrid format
    subroutine dPpHpp(UpLo, n, ap, info, nb )
  ! For the character argument, the case is insignificant. 
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper triangular packed format ('U') or
                ! lower triangular packed format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(inout) :: ap(n*(n+1)/2) ! Holds the matrix.
                ! It is rearranged to the new format if info==0.
                ! Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                !   0 Successful rearrangement.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size to be used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
    end subroutine dPpHpp
  end interface PpHpp

  interface HppPp
  ! Rearrange a blocked hybrid matrix to packed triangular format
    subroutine dHppPp(UpLo, n, ap, info, nb )
  ! For the character argument, the case is insignificant. 
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(inout) :: ap(n*(n+1)/2) ! Holds the matrix.
                ! It is rearranged to the packed triangular format if info==0.
                ! Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                !  0 Successful rearrangement.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
    end subroutine dHppPp
  end interface HppPp

  interface HppTf
  ! Cholesky factorize a matrix in blocked hybrid format.
    subroutine dHppTf(UpLo, n, ap, info, nb )
  ! For the character argument, the case is insignificant. 
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(inout) :: ap(n*(n+1)/2) ! Holds the matrix.
                ! It is replaced by its Cholesky factor if info==0.
                ! It is unaltered if info<0, but altered if info>0.
      integer, intent(out) :: info ! set to one of these values:
                ! 0 Successful factorization.
                ! /= 0 Failure:
                  ! > 0. The leading minor of order info is not positive
                  ! definite, and the factorization could not be completed.
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
    end subroutine dHppTf
  end interface HppTf

  interface HppTs
  ! Solve one or more sets of equations, given the Cholesky factorization
  ! of its matrix in blocked hybrid format.
  ! For the character argument, the case is insignificant. 
    subroutine dHppTs(UpLo, n, nrhs, ap, b, ldb, info, nb, mb )
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      integer , intent(in) :: nrhs ! The number of right-hand sides.
      real(wp), intent(in) :: ap(n*(n+1)/2) ! Holds the matrix.
      integer , intent(in) :: ldb ! The first dimension of the array b.
      real(wp), intent(inout) :: b(ldb,*) ! b(1:n,1:nrhs) holds the 
                ! right-hand sides on entry and is overwritten by the 
                ! solution if info==0. Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                ! 0 Successful solution.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
      integer, intent(in), optional :: mb ! If present, specifies the block
                ! size to be used for the right-hand sides. If absent,
                ! the value of nb is used.
    end subroutine dHppTs
  end interface HppTs

  interface HppTs1
  ! Solve one set of equations, given the Cholesky factorization
  ! of its matrix in blocked hybrid format.
  ! For the character argument, the case is insignificant. 
    subroutine dHppTs1(UpLo, n, ap, b, info, nb )
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(in) :: ap(n*(n+1)/2) ! Holds the matrix.
      real(wp), intent(inout) :: b(n) ! Holds the right-hand side on entry
                ! and is overwritten by the solution if info==0. 
                ! Otherwise, it is unaltered..
      integer, intent(out) :: info ! set to one of these values:
                ! 0 Successful solution.
                ! < 0 The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
    end subroutine dHppTs1
  end interface HppTs1

   interface kcf
     subroutine dkcf(n, a, lda, info)
! Kernel subroutine for Cholesky factorization
       integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
       integer, intent(in) :: n ! Specifies the matrix order.
       integer, intent(in) :: lda  ! Leading extent of array a.
               ! The inequality lda >= n must hold.
               ! For efficiency, the value n for lda is preferable.
       real(wp), intent(inout) :: a(lda,*) ! The upper-triangular part,
               ! a(i,j), i<=j<=n, must be set to hold the upper-triangular part
               ! of the matrix and is overwritten by the Cholesky factor if 
               ! info==0. It is altered if info>0.
               ! For efficient execution, a(:,:) should fit into
               ! level-1 cache.
       integer, intent(out) :: info ! set to one of these values:
                 ! 0 Successful solution.
                 ! > 0. The leading minor of order info is not positive
                 ! definite, and the factorization could not be completed.
     end subroutine dkcf
   end interface kcf
end module block_hybrid_Cholesky


  subroutine dPpHpp(UpLo, n, ap, info, nb )
! Rearrange a packed triangular matrix to blocked hybrid format
  ! .. Use statements ..
  use block_hybrid_Cholesky, only: default_nb, copy => dcopy
  implicit none
  ! For the character argument, the case is insignificant. 
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper triangular packed format ('U') or
                ! lower triangular packed format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(inout) :: ap(n*(n+1)/2) ! Holds the matrix.
                ! It is rearranged to the new format if info==0.
                ! Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                !   0 Successful rearrangement.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size to be used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
  ! .. Locals ..
  integer :: lnb ! contains either min(n,nb) or min(n,default_nb)
  real (wp), allocatable :: buf(:) ! Buffer for block column

  ! .. Executable Statements ..
  info = 0
  if( UpLo /= 'l' .and. UpLo /= 'L' .and. &
    UpLo /= 'u' .and. UpLo /= 'U' ) then
    info = -1
  else if( n < 0 ) then
    info = -2
  else if( present(nb) )then
    if( nb <= 0 )then
      info = -5
    else
      lnb = min(nb,n)
    end if
  else
    lnb = min(default_nb,n)
  end if

  if( info == 0 .and. n /= 0 )then
    allocate(buf(n*lnb),stat=info)
    if( info /= 0 )then
      deallocate(buf,stat=info)
      info = -20
      return
    end if
  else
    return
  end if

    if( UpLo == 'u' .or. UpLo == 'U' )then
      call PphPu( lnb )
    else
      call PphPl( lnb )
    end if
    deallocate(buf,stat=lnb)

contains

    subroutine PphPl(nb)
! Rearrange packed lower-triangular matrix to lower blocked hybrid format
! .. Arguments ..
      integer, intent (in) :: nb ! Block size
! .. Locals ..
      integer :: ib, jb ! Row, column length of block 
      integer :: i1, j1 ! First row, column of block 
      integer :: i, j ! Row, column index within the block
      integer :: ijap, ijbuf ! Positions of ap(i1+i-1,j1+j-1), buf(i,j)
      integer :: ab, bb ! Positions in a, buf of start of block
      integer :: ac, bc ! Positions in a, buf of start of column
! .. Intrisics ..
      INTRINSIC min
! .. Executable Statements ..
      ijap = 1
! Main loop over block columns
      do j1 = 1, n, nb
        ab = ijap
        jb = min(nb,n-j1+1)
        ac = ab
        ab = ab + jb
        bc = 1
        bb = 1 + (jb*(jb+1))/2
! Copy the diagonal block to the buffer
        do j = 1, jb
          ijap = ac
          ac = ac + n - j - j1 + 2
          ijbuf = bc
          bc = bc + j + 1
          do i = j, jb
            buf(ijbuf) = ap(ijap)
            ijap = ijap + 1
            ijbuf = ijbuf + i
          end do
        end do
! Copy the off-diagonal blocks to the buffer
        do i1 = j1 + jb, n, nb
          ib = min(nb,n-i1+1)
          ac = ab
          ab = ab + ib
          bc = bb
          bb = bb + ib*jb
          do j = 1, jb
            ijap = ac
            ac = ac + n - j - j1 + 1
            ijbuf = bc
            bc = bc + 1
            do i = 1, ib
              buf(ijbuf) = ap(ijap)
              ijap = ijap + 1
              ijbuf = ijbuf + jb
            end do
          end do
        end do
! Copy array back from buffer
        call copy(bb-1,buf,1,ap(ijap-bb+1),1)
      end do
    end subroutine PphPl

    subroutine PphPu(nb)
! Rearrange packed upper-triangular matrix to upper blocked hybrid
! .. Arguments ..
      integer, intent (in) :: nb ! Row and column block size
! .. Locals ..
      integer :: i1 ! First row of current block
      integer :: iis ! Current position in buf
      integer :: j ! Column index within block
      integer :: jb ! Number of columns in current block column 
      integer :: jjb ! Position in ap of start of current block column
      integer :: js ! Current position in ap
      integer :: j1 ! First column of current block column
      integer :: ld ! Length of current column
      integer :: ld1 ! Length of first column of the block
! .. Intrisic ..
      intrinsic min

! .. Executable Statements ..
      jjb = (nb*nb+nb)/2 + 1
      ld1 = 1 + nb
!     First triangle is not moved
      do j1 = nb + 1, n, nb
        jb = min(nb,n-j1+1)
        iis = 1
        do i1 = 1, j1 - nb, nb
          ld = ld1
          js = jjb + i1 - 1
          do j = 1, jb
            call copy(nb,ap(js),1,buf(iis),1)
            js = js + ld
            ld = ld + 1
            iis = iis + nb
          end do
        end do
!       Now process current triangle
        ld = ld1
        js = jjb + j1 - 1
        do j = 1, jb
          call copy(j,ap(js),1,buf(iis),1)
          js = js + ld
          ld = ld + 1
          iis = iis + j
        end do
!       Copy block column back to ap
        call copy(iis-1,buf,1,ap(jjb),1)
        jjb = jjb + iis - 1
        ld1 = ld1 + nb
      end do

    end subroutine PphPu

end subroutine dPpHpp

subroutine dHppTf(UpLo, n, ap, info, nb )
  ! Cholesky factorization of a symmetric positive-definite matrix held
  ! in packed blocked hybrid format.
  ! .. Use statements ..
  use block_hybrid_Cholesky, only: default_nb, copy => dcopy, &
            syrk => dsyrk, gemm => dgemm, trsm => dtrsm, &
            kcf => dkcf
  implicit none
  integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
  ! .. Arguments ..
  character, intent (in) :: UpLo ! Specifies whether the matrix
            ! is held in upper triangular packed format ('U') or
            ! lower triangular packed format ('L').
  integer, intent (in) :: n ! Specifies the matrix order
  real (wp), intent (inout) :: ap(0:(n*(n+1))/2-1) ! Holds the matrix
            ! A in packed blocked hybrid format. 
            ! It is replaced by its Cholesky factor if info==0.
            ! It is unaltered if info<0, but altered if info>0.
  integer, intent(out) :: info ! set to one of these values:
            ! 0 Successful solution.
            ! /= 0 Failure:
              ! > 0. The leading minor of order info is not positive
              ! definite, and the factorization could not be completed.
              ! -20 Allocation of temporary array failed.
              ! Other: The value of argument -info is not valid.
  integer, intent(in), optional :: nb ! If present, specifies the block
            ! size to be used for the blocked hybrid format. If absent,
            ! the value of the module constant default_nb is used.

  ! .. Locals ..
  integer :: lnb ! contains either min(n,nb) or min(n,default_nb)
    real (wp), allocatable :: buf(:) ! Holds a full-format
                                     ! trianglar submatrix
  intrinsic min

  ! .. Executable Statements ..
  info = 0
  if( UpLo /= 'l' .and. UpLo /= 'L' .and. &
    UpLo /= 'u' .and. UpLo /= 'U' ) then
    info = -1
  else if( n < 0 ) then
    info = -2
  else if( present(nb) )then
    if( nb <= 0 )then
      info = -5
    else
      lnb = min(nb,n)
    end if
  else
    lnb = min(default_nb,n)
  end if

  if( info == 0 .and. n /= 0 )then
    allocate(buf(n*lnb),stat=info)
    if( info /= 0 )then
      deallocate(buf,stat=info)
      info = -20
      return
    end if
  else
    return
  end if

  if( UpLo == 'u' .or. UpLo == 'U' )then
    call HppFu( lnb )
  else
    call HppFl( lnb )
  end if
  deallocate(buf,stat=lnb)

contains

    subroutine HppFl( nb )
! Cholesky factorization of a symmetric positive-definite matrix held
! in lower blocked hybrid format.
! .. Arguments ..
      integer, intent (in) :: nb ! Block size

! .. parameters ..
      real (wp), parameter :: one = 1.0_wp
! .. Locals ..
      integer i ! Temporary variable
      integer ii ! Position in ap of diagonal entry
      integer j ! First column of current block column
      integer jb ! Number of columns in current block column 
      integer jd ! size of current trapezoid / increment for jj
      integer jj ! Position in ap of start of current diagonal block
      integer jk ! Position in ap of start of current block
      integer jks ! Starting value for jk
      integer kd ! jk pointer delta
      integer kds ! Starting value for kd
      integer k ! First column of current block
      integer kk ! Position in buf of diagonal entry
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb

! .. Executable Statements ..
      info = 0
      jj = 0
      nb2 = nb*nb
      nbt = (nb2+nb)/2
      jks = nb - nbt
      jd = n*nb + jks
      kds = jd - nb2
      do j = 1, n, nb ! Loop over block columns
        jb = min(n-j+1,nb)
! Move diagonal block to buf in full upper format
        ii = jj
        kk = 1
        do i = 1, jb
! Transfer i-th row of ap to i-th col of buf, if n > nb
          call copy(i,ap(ii),1,buf(kk),1)
          ii = ii + i
          kk = kk + jb
        end do
        kd = kds
        jk = jks
        do k = 1, j - nb, nb
! Update the diagonal block, held in buf
          call syrk('U','T',jb,nb,-one,ap(jk:),nb,one,buf,jb)
! Update the off-diagonal part of the block column
          if (n>j-1+nb) call gemm('T','N',nb,n-j+1-nb,nb,-one,ap(jk),nb, &
            ap(jk+nb2),nb,one,ap(jj+nbt),nb)
          jk = jk + kd
          kd = kd - nb2
        end do
! Cholesky factorization of diagonal block
        call kcf( jb, buf, jb, i)
        if (i>0) then
          info = i + j - 1
          return
        end if
! Calculate the sub-diagonal blocks of L in the current block column
        if (n>j-1+nb) call trsm('L','U','T','N',nb,n-j+1-nb,one,buf,jb, &
          ap(jj+nbt),nb)
! Move diagonal block back from buf to ap
        ii = jj
        kk = 1
        do i = 1, jb
          call copy(i,buf(kk),1,ap(ii),1)
          ii = ii + i
          kk = kk + jb
        end do
        jks = jks + nb2
        jj = jj + jd
        jd = jd - nb2
      end do
    end subroutine HppFl

    subroutine HppFu( nb )
! Cholesky factorization of a symmetric positive-definite
! matrix held in upper packed blocked hybrid format.
! .. Arguments ..
      integer, intent (in) :: nb ! Block size

! .. Locals ..
      integer i ! Temporary variable
      integer ib ! Number of columns in block column starting at i
      integer id ! Size of trapezoid starting at i
      integer ii ! Position in ap of entry (i,i)
      integer ij ! Position in ap of entry (i,j)
      integer ik ! Position in ap of entry (i,k)
      integer i0 ! Position in ap of block column starting at i
      integer j ! First column of current block column
      integer jb ! Number of columns in block column starting at j 
      integer jd ! Size of trapezoid starting at j
      integer jj ! Position in ap of start of current diagonal block
      integer jk ! Position in ap of start of current block
      integer j0 ! Position in ap of block column starting at j
      integer k ! First column of current block
      integer kk ! Position in buf of diagonal entry
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp

! .. Executable Statement ..
      info = 0
      j0 = 0
      nb2 = nb*nb
      nbt = (nb2+nb)/2
      jd = nbt
      do j = 1, n, nb
        jb = min(nb,n-j+1)
        jj = j0 + (j-1)*jb ! -> ap(j,j)
        ii = jj
        kk = 1
        do i = 1, jb
          call copy(i,ap(ii),1,buf(kk),1)
          ii = ii + i
          kk = kk + jb
        end do
        jk = j0
        do k = 1, j - nb, nb
          call syrk('U','T',jb,nb,-one,ap(jk),nb,one,buf,jb)
          i0 = j0
          id = jd
          do i = j + nb, n, nb
            i0 = i0 + id
            ib = min(n-i+1,nb)
            ik = i0 + (k-1)*ib
            ij = i0 + (j-1)*ib
            call gemm('T','N',nb,ib,nb,-one,ap(jk),nb,ap(ik),nb,one,ap(ij), &
              nb)
            id = id + nb2
          end do
          jk = jk + nb*jb
        end do
        call kcf( jb, buf, jb, i)
        if (i>0) then
          info = j + i - 1
          return
        end if
        i0 = j0
        id = jd
        do i = j + nb, n, nb
          i0 = i0 + id
          ib = min(n-i+1,nb)
          ij = i0 + (j-1)*ib
          call trsm('L','U','T','N',nb,ib,one,buf,nb,ap(ij),nb)
          id = id + nb2
        end do
        ii = jj
        kk = 1
        do i = 1, jb
          call copy(i,buf(kk),1,ap(ii),1)
          ii = ii + i
          kk = kk + jb
        end do
        j0 = j0 + jd
        jd = jd + nb2
      end do

    end subroutine HppFu

end subroutine dHppTf


subroutine dHppPp(UpLo, n, ap, info, nb )
! Rearrange blocked hybrid packed matrix to packed triangular format
  ! .. Use statements ..
  use block_hybrid_Cholesky, only: default_nb, copy => dcopy
  implicit none
  ! For the character argument, the case is insignificant. 
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(inout) :: ap(n*(n+1)/2) ! Holds the matrix.
                ! It is rearranged to the packed triangular format if info==0.
                ! Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                !  0 Successful rearrangement.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
  ! .. Locals ..
  integer :: lnb ! contains either min(n,nb) or min(n,default_nb)
  real (wp), allocatable :: buf(:) ! Buffer for block column

  ! .. Executable Statements ..
  info = 0
  if( UpLo /= 'l' .and. UpLo /= 'L' .and. &
    UpLo /= 'u' .and. UpLo /= 'U' ) then
    info = -1
  else if( n < 0 ) then
    info = -2
  else if( present(nb) )then
    if( nb <= 0 )then
      info = -5
    else
      lnb = min(nb,n)
    end if
  else
    lnb = min(default_nb,n)
  end if

  if( info == 0 .and. n /= 0 )then
    allocate(buf(n*lnb),stat=info)
    if( info /= 0 )then
      deallocate(buf,stat=info)
      info = -20
      return
    end if
  else
    return
  end if

    if( UpLo == 'u' .or. UpLo == 'U' )then
      call HppPu( lnb )
    else
      call HppPl( lnb )
    end if
    deallocate(buf,stat=lnb)

contains

    subroutine HppPl(nb)
! Rearrange lower blocked hybrid matrix to packed lower-triangular format 
! .. Arguments ..
      integer, intent (in) :: nb ! Block size
! .. Locals ..
      integer :: ib, jb ! Row, column length of block 
      integer :: jb2 ! (jb*(jb+1))/2
      integer :: i1, j1 ! First row, column of block 
      integer :: i, j ! Row, column index within the block
      integer :: ijap, ijbuf ! Positions of ap(i1+i-1,j1+j-1), buf(i,j)
      integer :: ab, bb ! Positions in a, buf of start of block
      integer :: ac, bc ! Positions in a, buf of start of column
      intrinsic min
! .. Executable Statements ..
     ijap = 1
! Main loop over block columns
      do j1 = 1, n, nb
        jb = min(n-j1+1,nb) ! number of cols in current block
        jb2 = (jb*(jb+1))/2 ! size of current triangle
        ab = ijap
! Copy entire block column ( trapezoid ) to buffer buf
        call copy((n-j1+1-jb)*jb+jb2,ap(ab),1,buf,1)
        ac = ab
        ab = ab + jb
        bc = 1
        bb = 1 + jb2
! Rearrange the diagonal block of ap
        do j = 1, jb
          ijap = ac
          ac = ac + n - j - j1 + 2
          ijbuf = bc
          bc = bc + j + 1
          do i = j, jb
            ap(ijap) = buf(ijbuf)
            ijap = ijap + 1
            ijbuf = ijbuf + i
          end do
        end do
! Rearrange the off-diagonal blocks of ap
        do i1 = j1 + jb, n, nb
          ib = min(nb,n-i1+1)
          ac = ab
          ab = ab + ib
          bc = bb
          bb = bb + ib*jb
          do j = 1, jb
            ijap = ac
            ac = ac + n - j - j1 + 1
            ijbuf = bc
            bc = bc + 1
            do i = 1, ib
              ap(ijap) = buf(ijbuf)
              ijap = ijap + 1
              ijbuf = ijbuf + jb
            end do
          end do
        end do
      end do
    end subroutine HppPl

    subroutine HppPu(nb)
! Rearrange upper blocked hybrid matrix to packed upper-triangular 
! .. Arguments ..
      integer, intent (in) :: nb ! Row and column block size
! .. Locals ..
      integer :: i1 ! First row of current block
      integer :: iis ! Current position in buf
      integer :: j ! Column index within block
      integer :: jb ! Number of columns in current block column 
      integer :: jjb ! Position in ap of start of current block column
      integer :: js ! Current position in ap
      integer :: j1 ! First column of current block column
      integer :: ld ! Length of current column
      integer :: ld1 ! Length of first column of the block
      intrinsic min

! .. Executable Statements ..
      jjb = (nb*nb+nb)/2 + 1
      ld1 = 1 + nb
!     First triangle is not moved
      do j1 = nb + 1, n, nb
!     Copy block column to buf
        jb = min(nb,n-j1+1)
        call copy(((j1+j1+jb-1)*jb)/2,ap(jjb),1,buf,1)
        iis = 1
        do i1 = 1, j1 - nb, nb
          ld = ld1
          js = jjb + i1 - 1
          do j = 1, jb
            call copy(nb,buf(iis),1,ap(js),1)
            js = js + ld
            ld = ld + 1
            iis = iis + nb
          end do
        end do
!       Now process current triangle
        ld = ld1
        js = jjb + j1 - 1
        do j = 1, jb
          call copy(j,buf(iis),1,ap(js),1)
          js = js + ld
          ld = ld + 1
          iis = iis + j
        end do
        jjb = jjb + iis - 1
        ld1 = ld1 + nb
      end do
    end subroutine HppPu

end subroutine dHppPp

subroutine dHppTs1( UpLo, n, ap, b, info, nb )
  ! Solve one set of equations, given the Cholesky factorization
  ! of its matrix in blocked hybrid format.
  ! For the character argument, the case is insignificant.
  ! .. Use statements ..
  use block_hybrid_Cholesky, only: default_nb, tpsv => dtpsv, gemv => dgemv
  implicit none
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
  ! .. Arguments ..
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      real(wp), intent(in) :: ap(0:n*(n+1)/2-1) ! Holds the matrix.
      real(wp), intent(inout) :: b(n) ! Holds the right-hand side on entry
                ! and is overwritten by the solution if info==0. 
                ! Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                ! 0 Successful solution.
                ! < 0 The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
  ! .. Locals ..
  integer :: lnb ! block size used for the blocked hybrid format:
                 ! either min(n,nb) or min(n,default_nb)
  intrinsic min

  ! .. Executable Statements ..
  info = 0
  if( UpLo /= 'l' .and. UpLo /= 'L' .and. &
    UpLo /= 'u' .and. UpLo /= 'U' ) then
    info = -1
  else if( n < 0 ) then
    info = -2
  else if( present(nb) ) then
    if( nb <= 0 )then
       info = -6
    else
       lnb = min(n,nb)
    end if
  else
    lnb = min(default_nb,n)
  end if

  if( info /= 0 .or. n == 0 ) then
    return
  end if

  if( UpLo == 'u' .or. UpLo == 'U' )then
    call HppSu1( lnb )
  else
    call HppSl1( lnb )
  end if

contains

    subroutine HppSl1( nb )
! Solve a set of equations with one right-hand side, given the
! Cholesky factorization of its matrix in lower blocked hybrid format
! .. Arguments ..
      integer, intent (in) :: nb ! Block size

! .. Locals ..
      integer i ! Temporary variable
      integer j ! First column of current block column
      integer jb ! Number of columns in current block column 
      integer jd ! size of current trapezoid / increment for jj
      integer jj ! Position in ap of start of current diagonal block
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp

      jj = 0
      nb2 = nb*nb
      nbt = (nb2+nb)/2
      jd = n*nb - nbt + nb

!     solve U'*Y=B
      do j = 1, n - nb, nb
        call tpsv('U','T','N',nb,ap(jj),b(j),1)
        call gemv('T',nb,n+1-j-nb,-one,ap(jj+nbt),nb,b(j),1,one,b(j+nb),1)
        jj = jj + jd
        jd = jd - nb2
      end do
      jb = n - j + 1
      call tpsv('U','T','N',jb,ap(jj),b(j),1)

!     solve U*X=Y
      call tpsv('U','N','N',jb,ap(jj),b(j),1)
      i = j - nb
      do j = i, 1, -nb
        jd = jd + nb2
        jj = jj - jd
        call gemv('N',nb,n+1-j-nb,-one,ap(jj+nbt),nb,b(j+nb),1,one,b(j),1)
        call tpsv('U','N','N',nb,ap(jj),b(j),1)
      end do

    end subroutine HppSl1

    subroutine HppSu1( nb )
! Solve a set of equations with one right-hand side, given the
! Cholesky factorization of its matrix in upper blocked hybrid format
! .. Arguments ..
      integer, intent (in) :: nb ! Block size

! .. Locals ..
      integer i ! Temporary variable
      integer ij ! Position in ap of start of current block
      integer j ! First column of current block column
      integer k ! Final value of j
      integer jb ! Number of columns in current block column 
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp

      nb2 = nb*nb
      nbt = (nb2+nb)/2

!       Solve U'*Y=B
      ij = -nbt ! -> starting block of AP
      do j = 1, n, nb
        jb = min(n-j+1,nb)
        ij = ij + nbt
        do i = 1, j - nb, nb
          call gemv('T',nb,jb,-one,ap(ij),nb,b(i),1,one,b(j),1)
          ij = ij + nb*jb ! -> next block of AP
        end do ! do i
        call tpsv('U','T','N',jb,ap(ij),b(j),1)
      end do ! do j

!       Solve U*X=Y
      k = j - nb
      do j = k, 1, -nb
        jb = min(n-j+1,nb)
        call tpsv('U','N','N',jb,ap(ij),b(j),1)
        do i = j - nb, 1, -nb
          ij = ij - nb*jb ! -> next block of AP
          call gemv('N',nb,jb,-one,ap(ij),nb,b(j),1,one,b(i),1)
        end do ! do i
        ij = ij - nbt
      end do ! do j

    end subroutine HppSu1

end subroutine dHppTs1

subroutine dHppTs(UpLo, n, nrhs, ap, b, ldb, info, nb, mb )
  ! Solve one or more sets of equations, given the Cholesky factorization
  ! of its matrix in blocked hybrid format.
  ! For the character argument, the case is insignificant.
  ! .. Use statements ..
  use block_hybrid_Cholesky, only: default_nb, copy => dcopy, &
                  gemm => dgemm, trsm => dtrsm, tpsv => dtpsv, &
                  gemv => dgemv
  implicit none
      integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
      character(len=1), intent(in) :: UpLo ! Specifies whether the matrix
                ! is held in upper blocked hybrid format ('U') or
                ! lower blocked hybrid format ('L').
      integer, intent(in) :: n ! Specifies the matrix order
      integer , intent(in) :: nrhs ! The number of right-hand sides.
      real(wp), intent(in) :: ap(0:n*(n+1)/2-1) ! Holds the matrix.
      integer , intent(in) :: ldb ! The first dimension of the array b.
      real(wp), intent(inout) :: b(ldb,*) ! b(1:n,1:nrhs) holds the 
                ! right-hand sides on entry and is overwritten by the 
                ! solution if info==0. Otherwise, it is unaltered.
      integer, intent(out) :: info ! set to one of these values:
                ! 0 Successful solution.
                ! < 0 Failure:
                  ! -20 Allocation of temporary array failed.
                  ! Other: The value of argument -info is not valid.
      integer, intent(in), optional :: nb ! If present, specifies the block
                ! size used for the blocked hybrid format. If absent,
                ! the value of the module constant default_nb is used.
      integer, intent(in), optional :: mb ! If present, specifies the block
                ! size to be used for the right-hand sides. If absent,
                ! the value of nb is used.
  ! .. Locals ..
  integer :: lnb ! block size used for the blocked hybrid format:
                 ! either min(n,nb) or min(n,default_nb)
  integer :: lmb ! block size used for the right-hand sides:
                 ! either min(nrhs,mb) or min(nrhs,lnb)
  real(wp), allocatable :: buf(:)  ! Buffer for the diagonal blocks of ap
  real(wp), allocatable :: bufb(:) ! Buffer for a block column of b
  intrinsic min

  ! .. Executable Statements ..

  info = 0
  if( UpLo /= 'l' .and. UpLo /= 'L' .and. &
    UpLo /= 'u' .and. UpLo /= 'U' ) then
    info = -1
  else if( n < 0 ) then
    info = -2
  else if( nrhs < 0 ) then
    info = -3
  else if( ldb < max(1,n) ) then
    info = -6
  else if( present(nb) )then
    if( nb <= 0 )then
      info = -8
    end if
  end if
  if( present(mb) .and. info == 0 )then
    if( mb <= 0 )then
      info = -9
    end if
  end if

  if( info == 0 .and. n /= 0 .and. nrhs /= 0 ) then
    if( present(nb) )then
      lnb = min(nb,n)
    else
      lnb = min(default_nb,n)
    endif
    if( present(mb) )then
      lmb = min(mb,nrhs)
    else
      lmb = min(lnb,nrhs)
    endif
    if( nrhs >= 4 ) then
      allocate( buf(n*lnb),bufb(n*lmb),stat=info )
      if( info /= 0 )then
        deallocate(buf,bufb,stat=info)
        info = -20
        return
      end if
    end if
  else
    return
  end if

  if( UpLo == 'u' .or. UpLo == 'U' )then
    call HppSu( lnb, lmb )
  else
    call HppSl( lnb, lmb )
  end if
  deallocate( buf,bufb,stat=lnb )

contains

    subroutine HppSu( nb, mb )
! Solve a set of equations with one or many right-hand sides, given the
! Cholesky factorization of its matrix in upper blocked hybrid format

! .. Arguments ..
      integer, intent (in) :: nb ! Block size
      integer, intent (in) :: mb ! Column block size

! .. Locals ..
      integer :: i ! Row index
      integer :: ib ! Number of columns in current block column of b
      integer :: ibp ! Current position in bufb
      integer :: ii ! Position in ap of diagonal entry
      integer :: ij ! Position in ap
      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: jl ! Final value of j
      integer :: jbp ! Current position in bufb
      integer :: jbd ! nb*kb, increment for jbp
      integer :: k ! First column of current block
      integer :: kb ! ! Number of columns in current block column of b
      integer :: kk ! Position in buf of diagonal entry
      integer :: kks ! Current position in buf
      integer :: nb2 ! nb*nb
      integer :: nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp
      INTRINSIC min

! .. Executable Statements ..

! If there are few rhs, solve each separately
      if (nrhs<4) then
        do j = 1, nrhs
          call HppSu1( nb, b(1,j) )
        end do
        return
      end if

      nb2 = nb*nb
      nbt = (nb2+nb)/2

!     Load the diagonal blocks (triangles) of ap into buf.
      kk = 1 ! note that kk is not restored to 1 here.
      ii = 0 ! is not restored to zero here.
      do j = 1, n, nb
        jb = min(n-j+1,nb)
        ii = ii + (j-1)*jb ! -> next triangle in ap
        do i = 1, jb
          call copy(i,ap(ii),1,buf(kk),1)
          ii = ii + i
          kk = kk + jb
        end do
      end do ! do j

!     Main loop over the block columns of b
      do k = 1, nrhs, mb
        kb = min(nrhs-k+1,mb)

!       Copy kb RHS from b(1:n,k:k+kb-1) to buffer bufb
        jbp = 1 ! -> to bufb
        do i = 1, n, nb
          ib = min(n-i+1,nb)
          do j = k, k + kb - 1
            call copy(ib,b(i,j),1,bufb(jbp),1)
            jbp = jbp + ib
          end do
        end do

!       Now solve A*X = B for these kb RHS.
        jbd = nb*kb ! delta for jbp

!       Solve U'*Y=B
        kks = 1
        jbp = 1 ! is changed but restored to 1
        ij = -nbt ! -> starting block of AP
        do j = 1, n, nb
          jb = min(n-j+1,nb)
          ij = ij + nbt
          ibp = 1 ! -> starting block of B
          do i = 1, j - nb, nb
            call gemm('T','N',jb,kb,nb,-one,ap(ij),nb,bufb(ibp),nb,one, &
              bufb(jbp),jb)
            ibp = ibp + jbd ! -> next block of B
            ij = ij + nb*jb ! -> next block of AP
          end do ! do i

          call trsm('L','U','T','N',jb,kb,one,buf(kks),jb,bufb(jbp),jb)
          jbp = jbp + jbd ! -> next block of B
          kks = kks + nb2 ! -> next triangle of ap in BUF
        end do ! do j

!       Solve U*X=Y
        jl = j - nb
        do j = jl, 1, -nb
          jb = min(n-j+1,nb)
          jbp = jbp - jbd ! -> current block of B
          kks = kks - nb2 ! -> current triangle of ap in BUF
          call trsm('L','U','N','N',jb,kb,one,buf(kks),jb,bufb(jbp),jb)
          ibp = jbp
          do i = j - nb, 1, -nb
            ij = ij - nb*jb ! -> next block of AP
            ibp = ibp - jbd ! -> next block of B
            call gemm('N','N',nb,kb,jb,-one,ap(ij),nb,bufb(jbp),jb,one, &
              bufb(ibp),nb)
          end do ! do i
          ij = ij - nbt
        end do ! do j

!       Copy solution to b(1:n,k:k+kb-1) 
        jbp = 1 ! -> to buf
        do i = 1, n, nb
          ib = min(n-i+1,nb)
          do j = k, k + kb - 1
            call copy(ib,bufb(jbp),1,b(i,j),1)
            jbp = jbp + ib
          end do
        end do
      end do ! do k 

    end subroutine HppSu

    subroutine HppSl( nb, mb )
! Solve a set of equations with many right-hand sides, given the
! Cholesky factorization of its matrix in lower blocked hybrid format

! .. Arguments ..
      integer, intent (in) :: nb ! Block size
      integer, intent (in) :: mb ! Column block size

! .. Locals ..
      integer :: i ! Row index
      integer :: ib ! Number of columns in current block column of b
      integer :: ibp ! Current position in bufb
      integer :: ii ! Position in ap of diagonal entry
      integer :: ij ! Position in ap
      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: jl ! Final value of j
      integer :: jbp ! Current position in bufb
      integer :: jbd ! nb*kb, increment for jbp
      integer :: jd ! size of current trapezoid / increment for jj
      integer :: jj ! Position in ap of start of current diagonal block
      integer :: k ! First column of current block
      integer :: kb ! ! Number of columns in current block column of b
      integer :: kk ! Position in buf of diagonal entry
      integer :: kks ! Current position in buf
      integer :: nb2 ! nb*nb
      integer :: nbt ! Size of packed triangle of order nb
      integer :: nb1 ! Size of packed triangle of order nb-1
      real (wp), parameter :: one = 1.0_wp

      INTRINSIC min

! .. Executable Statements ..

! If there are few rhs, solve each separately
      if (nrhs<4) then
        do j = 1, nrhs
          call HppSl1( nb, b(1,j))
        end do
        return
      end if

      nb2 = nb*nb
      nbt = (nb2+nb)/2
      nb1 = nbt - nb

!     Load the diagonal blocks (triangles) of ap into buf.
      kks = 1 ! note that kks is not restored to 1 here.
      jd = n*nb - nb1 ! jd is not restored to this value here. 
      jj = 0 ! is not restored to zero here.
      do j = 1, n, nb
        jb = min(n-j+1,nb)
        ii = jj
        kk = kks
        do i = 1, jb
          call copy(i,ap(ii),1,buf(kk),1)
          ii = ii + i
          kk = kk + jb
        end do
        jj = jj + jd ! -> next triangle in ap
        jd = jd - nb2 ! next delta for jj
        kks = kks + nb2 ! -> next triangle of ap in buf
      end do ! do j

!     Main loop over the block columns of b
      kks = 1 ! is changed but restored to 1
      jd = n*nb - nb1 ! is changed but restored to this value
      jj = 0 ! is changed but restored to 0
      do k = 1, nrhs, mb
        kb = min(nrhs-k+1,mb)

!       Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
        jbp = 1 ! -> to bufb
        do i = 1, n, nb
          ib = min(n-i+1,nb)
          do j = k, k + kb - 1
            call copy(ib,b(i,j),1,bufb(jbp),1)
            jbp = jbp + ib
          end do
        end do

!      Solve U'*Y=B     
        jbd = nb*kb ! delta for jbp
        jbp = 1 ! is changed but restored to 1
        do j = 1, n - nb, nb
          call trsm('L','U','T','N',nb,kb,one,buf(kks),nb,bufb(jbp),nb)
          ibp = jbp + jbd ! -> starting block of B
          ij = jj + nbt ! -> starting block of AP
          do i = j + nb, n, nb
            jb = min(n-i+1,nb)
            call gemm('T','N',jb,kb,nb,-one,ap(ij),nb,bufb(jbp),nb,one, &
              bufb(ibp),jb)
            ibp = ibp + jbd ! -> next block of B
            ij = ij + nb2 ! -> next block of AP
          end do ! do i
          jj = jj + jd
          jd = jd - nb2
          jbp = jbp + jbd ! -> next block of B
          kks = kks + nb2 ! -> next triangle of ap
        end do ! do j
        jb = n - j + 1
        call trsm('L','U','T','N',jb,kb,one,buf(kks),jb,bufb(jbp),jb)

!       Solve U*X=Y
        call trsm('L','U','N','N',jb,kb,one,buf(kks),jb,bufb(jbp),jb)
        jl = j - nb
        do j = jl, 1, -nb
          jd = jd + nb2
          jj = jj - jd ! -> next triangle of AP
          ij = jj + nbt ! -> starting block of AP
          ibp = jbp ! -> starting update block of B
          jbp = jbp - jbd ! -> current block of B
          kks = kks - nb2 ! -> current triangle of ap
          do i = j + nb, n, nb
            jb = min(n-i+1,nb)
            call gemm('N','N',nb,kb,jb,-one,ap(ij),nb,bufb(ibp),jb,one, &
              bufb(jbp),nb)
            ibp = ibp + jbd ! -> next block of B
            ij = ij + nb2 ! -> next block of AP
          end do ! do i
          call trsm('L','U','N','N',nb,kb,one,buf(kks),nb,bufb(jbp),nb)
        end do ! do j

!       Copy solution to b(1:n,k:k+kb+1) 
        jbp = 1 ! -> to buf
        do i = 1, n, nb
          ib = min(n-i+1,nb)
!         copy buf(is:is+ib+k-1) to B(i:i+ib-1,0:k-1); K <= NB
          do j = k, k + kb - 1
            call copy(ib,bufb(jbp),1,b(i,j),1)
            jbp = jbp + ib
          end do
        end do
      end do ! do k 

    end subroutine HppSl

    subroutine HppSl1( nb, b )
! Solve a set of equations with one right-hand side, given the
! Cholesky factorization of its matrix in lower blocked hybrid format
! .. Arguments ..
      integer, intent (in) :: nb ! Block size
      real (wp), intent (inout) :: b(n) ! rhs on input, solution on exit

! .. Locals ..
      integer i ! Temporary variable
      integer j ! First column of current block column
      integer jb ! Number of columns in current block column 
      integer jd ! size of current trapezoid / increment for jj
      integer jj ! Position in ap of start of current diagonal block
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp

      jj = 0
      nb2 = nb*nb
      nbt = (nb2+nb)/2
      jd = n*nb - nbt + nb

!     solve U'*Y=B
      do j = 1, n - nb, nb
        call tpsv('U','T','N',nb,ap(jj),b(j),1)
        call gemv('T',nb,n+1-j-nb,-one,ap(jj+nbt),nb,b(j),1,one,b(j+nb),1)
        jj = jj + jd
        jd = jd - nb2
      end do
      jb = n - j + 1
      call tpsv('U','T','N',jb,ap(jj),b(j),1)

!     solve U*X=Y
      call tpsv('U','N','N',jb,ap(jj),b(j),1)
      i = j - nb
      do j = i, 1, -nb
        jd = jd + nb2
        jj = jj - jd
        call gemv('N',nb,n+1-j-nb,-one,ap(jj+nbt),nb,b(j+nb),1,one,b(j),1)
        call tpsv('U','N','N',nb,ap(jj),b(j),1)
      end do

    end subroutine HppSl1

    subroutine HppSu1( nb, b )
! Solve a set of equations with one right-hand side, given the
! Cholesky factorization of its matrix in upper blocked hybrid format
! .. Arguments ..
      integer, intent (in) :: nb ! Block size
      real (wp), intent (inout) :: b(n) ! rhs on input, solution on exit

! .. Locals ..
      integer i ! Temporary variable
      integer ij ! Position in ap of start of current block
      integer j ! First column of current block column
      integer jb ! Number of columns in current block column 
      integer jl ! Final value of j 
      integer nb2 ! nb*nb
      integer nbt ! Size of packed triangle of order nb
      real (wp), parameter :: one = 1.0_wp

      nb2 = nb*nb
      nbt = (nb2+nb)/2

!       Solve U'*Y=B
      ij = -nbt ! -> starting block of AP
      do j = 1, n, nb
        jb = min(n-j+1,nb)
        ij = ij + nbt
        do i = 1, j - nb, nb
          call gemv('T',nb,jb,-one,ap(ij),nb,b(i),1,one,b(j),1)
          ij = ij + nb*jb ! -> next block of AP
        end do ! do i
        call tpsv('U','T','N',jb,ap(ij),b(j),1)
      end do ! do j

!       Solve U*X=Y
      jl = j - nb
      do j = jl, 1, -nb
        jb = min(n-j+1,nb)
        call tpsv('U','N','N',jb,ap(ij),b(j),1)
        do i = j - nb, 1, -nb
          ij = ij - nb*jb ! -> next block of AP
          call gemv('N',nb,jb,-one,ap(ij),nb,b(j),1,one,b(i),1)
        end do ! do i
        ij = ij - nbt
      end do ! do j

    end subroutine HppSu1

end subroutine dHppTs



subroutine dkcf(n,a,lda,info)
! Kernel subroutine for Cholesky factorization
  ! .. Use statements ..
  use block_hybrid_cholesky, only : trsm => dtrsm, syrk => dsyrk
  implicit none
  integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
! .. Arguments ..
  integer, intent (in) :: n ! Specifies the matrix order.
  integer, intent (in) :: lda ! Leading extent of array a.
          ! The inequality lda >= n must hold.
          ! For efficiency, the value n for lda is preferable.
  real (wp), intent (inout) :: a(0:lda-1,0:n-1) ! The upper-triangular part,
               ! a(i,j), i<=j<=n-1, must be set to hold the upper-triangular part
               ! of the matrix and is overwritten by the Cholesky factor if 
               ! info==0. It is altered if info>0.
               ! For efficient execution, a(:,:) should fit into
               ! level-1 cache.
  integer, intent (out) :: info ! set to one of these values:
          ! 0 Successful solution.
          ! > 0. The leading minor of order info is not positive
          ! definite, and the factorization could not be completed.
! .. Locals ..
  real (wp), parameter :: zero = 0.0_wp, one = 1.0_wp
  integer :: i ! Temporary variable
  integer :: j ! Temporary variable
  integer :: k ! Temporary variable
  integer :: nr ! Temporary variable
  real (wp) :: rd0 ! scalar to hold reciprocal of diagonal a(k,k) for 0<=k<n 
  real (wp) :: t00, t01, t02, t03 ! regs to hold a(j,i:i+3)
  real (wp) :: ai0, ai1, ai2, ai3 ! regs to hold a(k,i:i+3)
!     or a(k:k+3,j) or a(k:k+3,j+1)                       

!     this code uses 8 FP regs



!     This is a 1 x 4 blocking kernel
!     It is ONLY optimized for the case MOD(N,4) = 0

  info = 0
  nr = mod(n,4)

!     Main loop: factor A(0:n-nr-1,0:n-nr-1) = U'*U

  do j = 0, n - nr - 4, 4

!       Process row j = a(j,j:n-1)

!       Process a(j  ,j  )

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j)
      ai1 = a(k+1,j)
      ai2 = a(k+2,j)
      ai3 = a(k+3,j)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 2
    ai0 = sqrt(ai0)
    a(j,j) = ai0
    rd0 = one/ai0


!       process 1 by 3 rectangle a(j,j+1:j+3)

    t00 = a(j,j+1)
    t01 = a(j,j+2)
    t02 = a(j,j+3)
    do k = 0, j - 1
      ai0 = a(k,j+1)
      ai1 = a(k,j+2)
      ai2 = a(k,j+3)
      ai0 = ai0*a(k,j)
      ai1 = ai1*a(k,j)
      ai2 = ai2*a(k,j)
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
    end do

!       scale and store a(j,j+1:j+3)

    t00 = t00*rd0
    t01 = t01*rd0
    t02 = t02*rd0
    a(j,j+1) = t00
    a(j,j+2) = t01
    a(j,j+3) = t02

!       Process the remainder of rows j = a(j,j+4:n-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j,i)
      t01 = a(j,i+1)
      t02 = a(j,i+2)
      t03 = a(j,i+3)
      do k = 0, j - 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j)
        ai1 = ai1*a(k,j)
        ai2 = ai2*a(k,j)
        ai3 = ai3*a(k,j)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store a(j,i+0:i+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j,i) = t00
      a(j,i+1) = t01
      a(j,i+2) = t02
      a(j,i+3) = t03
    end do

!       Process row j+1 = a(j+1,j+1:n-1)

!       Process a(j+1,j+1)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+1)
      ai1 = a(k+1,j+1)
      ai2 = a(k+2,j+1)
      ai3 = a(k+3,j+1)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+1)
    ai0 = ai0*ai0
    t00 = t00 - ai0

    ai0 = a(j+1,j+1)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 3
    ai0 = sqrt(ai0)
    a(j+1,j+1) = ai0
    rd0 = one/ai0


!       process 1 by 2 rectangle a(j+1,j+2:j+3)

    t00 = a(j+1,j+2)
    t01 = a(j+1,j+3)
    do k = 0, j
      ai0 = a(k,j+2)
      ai1 = a(k,j+3)
      ai0 = ai0*a(k,j+1)
      ai1 = ai1*a(k,j+1)
      t00 = t00 - ai0
      t01 = t01 - ai1
    end do

!       scale and store a(j+1,j+2:j+3)

    t00 = t00*rd0
    t01 = t01*rd0
    a(j+1,j+2) = t00
    a(j+1,j+3) = t01

!       Process the remainder of row j+1 = a(j+1,j+4:n-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+1,i)
      t01 = a(j+1,i+1)
      t02 = a(j+1,i+2)
      t03 = a(j+1,i+3)
      do k = 0, j
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+1)
        ai1 = ai1*a(k,j+1)
        ai2 = ai2*a(k,j+1)
        ai3 = ai3*a(k,j+1)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store a(j+1,i+0:i+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+1,i) = t00
      a(j+1,i+1) = t01
      a(j+1,i+2) = t02
      a(j+1,i+3) = t03
    end do

!       Process row j+2 = a(j+2,j+2:n-1)

!       Process a(j+2,j+2)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+2)
      ai1 = a(k+1,j+2)
      ai2 = a(k+2,j+2)
      ai3 = a(k+3,j+2)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+2)
    ai1 = a(j+1,j+2)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    t00 = t00 - ai0
    t01 = t01 - ai1

    ai0 = a(j+2,j+2)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 4
    ai0 = sqrt(ai0)
    a(j+2,j+2) = ai0
    rd0 = one/ai0


!       process 1 by 1 rectangle a(j+2,j+3)

    t00 = a(j+2,j+3)
    do k = 0, j + 1
      ai0 = a(k,j+3)
      ai0 = ai0*a(k,j+2)
      t00 = t00 - ai0
    end do

!       scale and store A(J+2,J+3)

    t00 = t00*rd0
    a(j+2,j+3) = t00

!       Process the remainder of row j+2 = a(j+2,j+4:n-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+2,i)
      t01 = a(j+2,i+1)
      t02 = a(j+2,i+2)
      t03 = a(j+2,i+3)
      do k = 0, j + 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+2)
        ai1 = ai1*a(k,j+2)
        ai2 = ai2*a(k,j+2)
        ai3 = ai3*a(k,j+2)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store a(j+2,i+0:i+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+2,i) = t00
      a(j+2,i+1) = t01
      a(j+2,i+2) = t02
      a(j+2,i+3) = t03
    end do

!       Process row j+3 = a(j+3,j+3:n-1)

!       Process a(j+3,j+3)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+3)
      ai1 = a(k+1,j+3)
      ai2 = a(k+2,j+3)
      ai3 = a(k+3,j+3)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+3)
    ai1 = a(j+1,j+3)
    ai2 = a(j+2,j+3)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    ai2 = ai2*ai2
    t00 = t00 - ai0
    t01 = t01 - ai1
    t02 = t02 - ai2

    ai0 = a(j+3,j+3)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 5
    ai0 = sqrt(ai0)
    a(j+3,j+3) = ai0
    rd0 = one/ai0


!       Process the remainder of row j+3 = a(j+3,j+4:n-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+3,i)
      t01 = a(j+3,i+1)
      t02 = a(j+3,i+2)
      t03 = a(j+3,i+3)
      do k = 0, j + 2
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+3)
        ai1 = ai1*a(k,j+3)
        ai2 = ai2*a(k,j+3)
        ai3 = ai3*a(k,j+3)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store a(j+3,i+0:i+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+3,i) = t00
      a(j+3,i+1) = t01
      a(j+3,i+2) = t02
      a(j+3,i+3) = t03
    end do
  end do
  if (nr==0) return

!     here 0 < nr < 4 
!     scale last nr cols of a 

  call trsm('l','u','t','n',n-nr,nr,one,a,lda,a(0,n-nr),lda)

!     update triangle a(n-nr:n-1,n-nr,n-1) with a(0:n-nr-1,n-nr:n-1) 

  call syrk('u','t',nr,n-nr,-one,a(0,n-nr),lda,one,a(n-nr,n-nr),lda)

!     Cholesky factor triangle a(n-nr:n-1,n-nr,n-1)

  j = n - nr
  if (nr==3) then ! factor row j and update a(j+1:j+2,j+1:j+2)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = one/t00
    ai0 = a(j,j+1)
    ai1 = a(j,j+2)
    ai0 = ai0*rd0
    ai1 = ai1*rd0
    a(j,j+1) = ai0
    a(j,j+2) = ai1
    t00 = a(j+1,j+1)
    t01 = a(j+1,j+2)
    t02 = a(j+2,j+2)
    ai2 = ai0*ai0
    ai3 = ai0*ai1
    ai1 = ai1*ai1
    t00 = t00 - ai2
    t01 = t01 - ai3
    t02 = t02 - ai1
    a(j+1,j+1) = t00
    a(j+1,j+2) = t01
    a(j+2,j+2) = t02
    j = j + 1
  end if
  if (nr>=2) then ! factor row J and update a(j+1,j+1)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = one/t00
    ai0 = a(j,j+1)
    ai0 = ai0*rd0
    a(j,j+1) = ai0
    t00 = a(j+1,j+1)
    ai1 = ai0*ai0
    t00 = t00 - ai1
    a(j+1,j+1) = t00
    j = j + 1
  end if
  if (nr>=1) then ! a(j,j) = sqrt( a(j,j) )
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
  end if

  return

2 info = j + 1
  return
3 info = j + 2
  return
4 info = j + 3
  return
5 info = j + 4
  return
end subroutine dkcf

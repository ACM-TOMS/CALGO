! This module computes a polynomial of an overlapped block-diagonal matrix.
! The polynomial consists only of even order terms.
!
! References:
!
!   Nicholas J. Higham,
!   The scaling and squaring method for the matrix exponential revisited,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 26, Number 4, pp. 1179-1193, 2005.
!
!   _. _______,
!   On a sep-inverse estimate and its application to a modified
!   block Schur--Parlett algorithm,
!   submitted to ACM transactions on mathematical software, 2008.
!
! This module is intended for internal use only.

module ovlppolyfunc

  use floattypes
  use matrixpwrtag
  use ovlpmatrix

  implicit none

  public

contains


!    Name : rs_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats  a real matrix in single precision.

pure function rs_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  real   (kind=sp), intent(in) :: z(:,:,0:)
  real   (kind=sp)             :: rs_ovlp_poly(size(z,1),size(z,2))

  real   (kind=sp) :: gx(0:size(g,1)-1)
  real   (kind=sp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  rs_ovlp_poly = 0.0_sp
  temporary = 0.0_sp
  do i=0,size(g,1)-1; gx(i) = g(i); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = rs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = rs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = rs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = rs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call rs_ovlp_clear(diag,subd,head,tail,nonz,rs_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=sp), intent(in) :: a(:,:), b(:,:)
    real   (kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function rs_ovlp_poly


!    Name : cs_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats a complex matrix in single precision.

pure function cs_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=sp), intent(in) :: g(0:)
  complex(kind=sp), intent(in) :: z(:,:,0:)
  complex(kind=sp)             :: cs_ovlp_poly(size(z,1),size(z,2))

  complex(kind=sp) :: gx(0:size(g,1)-1)
  complex(kind=sp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  cs_ovlp_poly = cmplx(0.0_sp,0.0_sp,sp)
  temporary = cmplx(0.0_sp,0.0_sp,sp)
  do i=0,size(g,1)-1; gx(i) = cmplx(g(i),0.0_sp,sp); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = cs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = cs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = cs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = cs_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cs_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call cs_ovlp_clear(diag,subd,head,tail,nonz,cs_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=sp), intent(in) :: a(:,:), b(:,:)
    complex(kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function cs_ovlp_poly


!    Name : rw_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats  a real matrix in double precision.

pure function rw_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  real   (kind=wp), intent(in) :: z(:,:,0:)
  real   (kind=wp)             :: rw_ovlp_poly(size(z,1),size(z,2))

  real   (kind=wp) :: gx(0:size(g,1)-1)
  real   (kind=wp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  rw_ovlp_poly = 0.0_wp
  temporary = 0.0_wp
  do i=0,size(g,1)-1; gx(i) = g(i); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = rw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = rw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = rw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = rw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call rw_ovlp_clear(diag,subd,head,tail,nonz,rw_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=wp), intent(in) :: a(:,:), b(:,:)
    real   (kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function rw_ovlp_poly


!    Name : cw_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats a complex matrix in double precision.

pure function cw_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=wp), intent(in) :: g(0:)
  complex(kind=wp), intent(in) :: z(:,:,0:)
  complex(kind=wp)             :: cw_ovlp_poly(size(z,1),size(z,2))

  complex(kind=wp) :: gx(0:size(g,1)-1)
  complex(kind=wp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  cw_ovlp_poly = cmplx(0.0_wp,0.0_wp,wp)
  temporary = cmplx(0.0_wp,0.0_wp,wp)
  do i=0,size(g,1)-1; gx(i) = cmplx(g(i),0.0_wp,wp); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = cw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = cw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = cw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = cw_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cw_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call cw_ovlp_clear(diag,subd,head,tail,nonz,cw_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=wp), intent(in) :: a(:,:), b(:,:)
    complex(kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function cw_ovlp_poly

#ifdef __USE_TPREC

!    Name : rt_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats  a real matrix in triple precision.

pure function rt_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  real   (kind=tp), intent(in) :: z(:,:,0:)
  real   (kind=tp)             :: rt_ovlp_poly(size(z,1),size(z,2))

  real   (kind=tp) :: gx(0:size(g,1)-1)
  real   (kind=tp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  rt_ovlp_poly = 0.0_tp
  temporary = 0.0_tp
  do i=0,size(g,1)-1; gx(i) = g(i); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = rt_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = rt_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = rt_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = rt_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rt_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call rt_ovlp_clear(diag,subd,head,tail,nonz,rt_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=tp), intent(in) :: a(:,:), b(:,:)
    real   (kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rt_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function rt_ovlp_poly


!    Name : ct_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats a complex matrix in triple precision.

pure function ct_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=tp), intent(in) :: g(0:)
  complex(kind=tp), intent(in) :: z(:,:,0:)
  complex(kind=tp)             :: ct_ovlp_poly(size(z,1),size(z,2))

  complex(kind=tp) :: gx(0:size(g,1)-1)
  complex(kind=tp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  ct_ovlp_poly = cmplx(0.0_tp,0.0_tp,tp)
  temporary = cmplx(0.0_tp,0.0_tp,tp)
  do i=0,size(g,1)-1; gx(i) = cmplx(g(i),0.0_tp,tp); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = ct_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = ct_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = ct_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = ct_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        ct_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call ct_ovlp_clear(diag,subd,head,tail,nonz,ct_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=tp), intent(in) :: a(:,:), b(:,:)
    complex(kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = ct_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function ct_ovlp_poly

#endif
#ifdef __USE_QPREC

!    Name : rq_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats  a real matrix in quadruple precision.

pure function rq_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  real   (kind=qp), intent(in) :: z(:,:,0:)
  real   (kind=qp)             :: rq_ovlp_poly(size(z,1),size(z,2))

  real   (kind=qp) :: gx(0:size(g,1)-1)
  real   (kind=qp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  rq_ovlp_poly = 0.0_qp
  temporary = 0.0_qp
  do i=0,size(g,1)-1; gx(i) = g(i); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = rq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = rq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = rq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = rq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        rq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call rq_ovlp_clear(diag,subd,head,tail,nonz,rq_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=qp), intent(in) :: a(:,:), b(:,:)
    real   (kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function rq_ovlp_poly


!    Name : cq_ovlp_poly
!         :
! Purpose : This function computes a polynomial of an overlapped block-diagonal
!         : matrix. The polynomial consists only of even order terms.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm to the overlapped block-diagonal region.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!         : - 2*n is the order of the polynomial.
!         : - g(i) stores the coefficient of the 2*i-th power term.
!         : - z(*,i) is the (2*i)-th power of a matrix.
!         :
!  Output : - The value of the polynomial.
!    Note : This routine treats a complex matrix in quadruple precision.

pure function cq_ovlp_poly(diag,subd,head,tail,nonz,n,g,z)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in) :: n
  real   (kind=qp), intent(in) :: g(0:)
  complex(kind=qp), intent(in) :: z(:,:,0:)
  complex(kind=qp)             :: cq_ovlp_poly(size(z,1),size(z,2))

  complex(kind=qp) :: gx(0:size(g,1)-1)
  complex(kind=qp) :: temporary(size(z,1),size(z,2))
  integer          :: i, h, t, band(size(z,1))

  ! an even polynomial of 2n-th order

  cq_ovlp_poly = cmplx(0.0_qp,0.0_qp,qp)
  temporary = cmplx(0.0_qp,0.0_qp,qp)
  do i=0,size(g,1)-1; gx(i) = cmplx(g(i),0.0_qp,qp); end do

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  select case (n)                ! z(*,i) is (2*i)-th power of the matrix

    case (8)                     ! a 16-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
        temporary   (i,h:t) = z(i,h:t,0) * gx(4) &
                            + z(i,h:t,1) * gx(5) &
                            + z(i,h:t,2) * gx(6) &
                            + z(i,h:t,3) * gx(7) &
                            + z(i,h:t,4) * gx(8)
      end do

      temporary = matrixmul(temporary, z(:,:,4))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = cq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do
      
    case (6)                     ! a 12-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5) &
                            + z(i,h:t,3) * gx(6)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = cq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (5)                     ! a 10-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
        temporary   (i,h:t) = z(i,h:t,0) * gx(3) &
                            + z(i,h:t,1) * gx(4) &
                            + z(i,h:t,2) * gx(5)
      end do

      temporary = matrixmul(temporary, z(:,:,3))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = cq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (4)                     ! an 8-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
        temporary   (i,h:t) = z(i,h:t,0) * gx(2) &
                            + z(i,h:t,1) * gx(3) &
                            + z(i,h:t,2) * gx(4)
      end do

      temporary = matrixmul(temporary, z(:,:,2))

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = cq_ovlp_poly(i,h:t) + temporary(i,h:t)
      end do

    case (3)                     ! a 6-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2) &
                            + z(i,h:t,3) * gx(3)
      end do

    case (2)                     ! a 4-th order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1) &
                            + z(i,h:t,2) * gx(2)
      end do

    case (1)                     ! a 2-nd order polynomial

      do i=1,size(z,1)
        h = max(1,i-1); t = min(band(i),size(z,2))
        cq_ovlp_poly(i,h:t) = z(i,h:t,0) * gx(0) &
                            + z(i,h:t,1) * gx(1)
      end do

  end select

  ! Clear unnecessary off-diagonal entries.

  call cq_ovlp_clear(diag,subd,head,tail,nonz,cq_ovlp_poly)

contains

  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         : ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=qp), intent(in) :: a(:,:), b(:,:)
    complex(kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end function cq_ovlp_poly

#endif

end module ovlppolyfunc


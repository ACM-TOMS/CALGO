! This module computes rational approximants for phi-functions
! of an overlapped block-diagonal matrix.
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
!   _. _______,
!   On a sep-inverse estimate and its application to a modified
!   block Schur--Parlett algorithm,
!   submitted to ACM transactions on mathematical software, 2008.
!
! This module is intended for internal use only.

module ovlpratiofunc

  use floattypes
  use matrixpwrtag
  use ovlpmatrix
  use ovlppolyfunc
  use mcpcoefficients
  use mtrcfgphilog

  implicit none

  public

contains


!    Name : rs_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a real matrix in single precision.

pure subroutine rs_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  real   (kind=sp), intent(in   ) :: a(:,:)
  real   (kind=sp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=sp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  real   (kind=sp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = 0.0_sp
  do i=1,size(a,1); z2(i,i,0) = 1.0_sp; end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = rs_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = rs_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call rs_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call rs_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = rs_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = rs_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call rs_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call rs_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=sp), intent(in) :: a(:,:), b(:,:)
    real   (kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rs_ovlp_ratio


!    Name : cs_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a complex matrix in single precision.

pure subroutine cs_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  complex(kind=sp), intent(in   ) :: a(:,:)
  complex(kind=sp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=sp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  complex(kind=sp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = cmplx(0.0_sp,0.0_sp,sp)
  do i=1,size(a,1); z2(i,i,0) = cmplx(1.0_sp,0.0_sp,sp); end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = cs_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = cs_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call cs_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call cs_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = cs_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = cs_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call cs_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call cs_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=sp), intent(in) :: a(:,:), b(:,:)
    complex(kind=sp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cs_ovlp_ratio


!    Name : rw_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a real matrix in double precision.

pure subroutine rw_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  real   (kind=wp), intent(in   ) :: a(:,:)
  real   (kind=wp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=wp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  real   (kind=wp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = 0.0_wp
  do i=1,size(a,1); z2(i,i,0) = 1.0_wp; end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = rw_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = rw_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call rw_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call rw_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = rw_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = rw_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call rw_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call rw_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=wp), intent(in) :: a(:,:), b(:,:)
    real   (kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rw_ovlp_ratio


!    Name : cw_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a complex matrix in double precision.

pure subroutine cw_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  complex(kind=wp), intent(in   ) :: a(:,:)
  complex(kind=wp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=wp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  complex(kind=wp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = cmplx(0.0_wp,0.0_wp,wp)
  do i=1,size(a,1); z2(i,i,0) = cmplx(1.0_wp,0.0_wp,wp); end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = cw_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = cw_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call cw_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call cw_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = cw_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = cw_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call cw_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call cw_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=wp), intent(in) :: a(:,:), b(:,:)
    complex(kind=wp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cw_ovlp_ratio

#ifdef __USE_TPREC

!    Name : rt_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a real matrix in triple precision.

pure subroutine rt_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  real   (kind=tp), intent(in   ) :: a(:,:)
  real   (kind=tp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=tp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  real   (kind=tp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = 0.0_tp
  do i=1,size(a,1); z2(i,i,0) = 1.0_tp; end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = rt_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = rt_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call rt_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call rt_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = rt_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = rt_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call rt_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call rt_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=tp), intent(in) :: a(:,:), b(:,:)
    real   (kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rt_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rt_ovlp_ratio


!    Name : ct_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a complex matrix in triple precision.

pure subroutine ct_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  complex(kind=tp), intent(in   ) :: a(:,:)
  complex(kind=tp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=tp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  complex(kind=tp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = cmplx(0.0_tp,0.0_tp,tp)
  do i=1,size(a,1); z2(i,i,0) = cmplx(1.0_tp,0.0_tp,tp); end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = ct_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = ct_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call ct_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call ct_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = ct_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = ct_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call ct_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call ct_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=tp), intent(in) :: a(:,:), b(:,:)
    complex(kind=tp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = ct_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine ct_ovlp_ratio

#endif
#ifdef __USE_QPREC

!    Name : rq_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a real matrix in quadruple precision.

pure subroutine rq_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  real   (kind=qp), intent(in   ) :: a(:,:)
  real   (kind=qp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=qp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  real   (kind=qp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = 0.0_qp
  do i=1,size(a,1); z2(i,i,0) = 1.0_qp; end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = rq_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = rq_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call rq_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call rq_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = rq_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = rq_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call rq_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call rq_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    real   (kind=qp), intent(in) :: a(:,:), b(:,:)
    real   (kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = rq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine rq_ovlp_ratio


!    Name : cq_ovlp_ratio
!         :
! Purpose : This subroutine computes rational approximants
!         : for phi-functions of an overlapped block-diagonal matrix.
!         :
!   Input : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of ``A'' is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of ``A'' is nonzero.
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
!         : - The integer ``m'' specifies the order of approximation.
!         :   An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!         :
!  Output : - v(*,k) stores the value of the rational function that
!         :   approximates the overlapped block-diagonal region of phi_k.
!         : - ``err'' stores error informations.
!    Note : This routine treats a complex matrix in quadruple precision.

pure subroutine cq_ovlp_ratio(diag,subd,head,tail,nonz,m,upto,a,v)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  integer,          intent(in   ) :: m, upto
  complex(kind=qp), intent(in   ) :: a(:,:)
  complex(kind=qp), intent(out  ) :: v(:,:,0:)

  integer          :: i, k, n, h, t, perm(size(a,1)), band(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=qp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))
  complex(kind=qp) :: tm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Band structure of the argument.

  do i=1,size(diag,1)
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8.

  z2(:,:,0) = cmplx(0.0_qp,0.0_qp,qp)
  do i=1,size(a,1); z2(i,i,0) = cmplx(1.0_qp,0.0_qp,qp); end do
  z2(:,:,1) = matrixmul(a,a)

  if ( 5 <= m             ) z2(:,:,2) = matrixmul(z2(:,:,1), z2(:,:,1))
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = matrixmul(z2(:,:,2), z2(:,:,1))
  if (17 <= m             ) z2(:,:,4) = matrixmul(z2(:,:,3), z2(:,:,1))

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  tm = cq_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,1:2*n+1:2),z2)
  qm = matrixmul(tm, a)
  pm = -qm

  v(:,:,0) = cq_ovlp_poly(diag,subd,head,tail,nonz,n,cq(0,0:2*n:2),z2)
  do i=1,size(a,1)
    h = max(1,i-1); t = min(band(i),size(a,2))
    qm(i,h:t) = v(i,h:t,0) + qm(i,h:t)
    pm(i,h:t) = v(i,h:t,0) + pm(i,h:t)
  end do

  call cq_ovlp_ludcmp(diag,subd,head,tail,nonz,qm,perm)
  call cq_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,0))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    tm = cq_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,1:2*n+1:2),z2)
    pm = matrixmul(tm,a)

    tm = cq_ovlp_poly(diag,subd,head,tail,nonz,n,cp(k,0:2*n  :2),z2)
    do i=1,size(a,1)
      h = max(1,i-1); t = min(band(i),size(a,2))
      pm(i,h:t) = tm(i,h:t) + pm(i,h:t)
    end do

    call cq_ovlp_ainvb(diag,subd,head,tail,nonz,qm,perm,pm,v(:,:,k))
  end do

  ! Clear unnecessary off-diagonal entries.

  do k=0,upto
    call cq_ovlp_clear(diag,subd,head,tail,nonz,v(:,:,k))
  end do

contains


  !    Name : matrixmul
  ! Purpose : This function computes the block-wise product of two overlapped
  !         : block-diagonal matrices without repeating equivalent operations.
  !   Input : A and B are overlapped block-diagonal matrices
  !         : of the same block structure.
  !  Output : The return value is equivalent to the block-wise product,
  !         :    ( A_{ii}B_{ii}; 1 <= i <= size(diag,1) ).

  pure function matrixmul(a,b)
    complex(kind=qp), intent(in) :: a(:,:), b(:,:)
    complex(kind=qp)             :: matrixmul(size(a,1),size(b,2))
    matrixmul = cq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  end function matrixmul

end subroutine cq_ovlp_ratio

#endif

end module ovlpratiofunc


! This module computes rational approximants for phi-functions
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

module rationalphi

  use floattypes
  use matrixpwrtag
  use mcpcoefficients
  use polynomial12th
  use mtrcfgphilog

  implicit none

  public

contains


!    Name : rsdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, single precision, diagonal matrix.

pure subroutine rsdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=sp), intent(in   ) :: a(:)
  real   (kind=sp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=sp) :: z2(size(a,1),0:4)
  real   (kind=sp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = rsdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(rsdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(rsdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rsdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(rsdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rsdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(rsdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rsdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = rsdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call rsdg_ludcmp(qm,perm)
  call rsdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(rsdg_normi(qm))
  err % qm_norm(1,0) = real(rsdg_norm1(qm))
  err % pm_norm(0,0) = real(rsdg_normi(pm))
  err % pm_norm(1,0) = real(rsdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rsdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = rsdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rsdg_normi(pm))
    err % pm_norm(1,k) = real(rsdg_norm1(pm))

    call rsdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine rsdg_rationalphi


!    Name : rssq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, single precision, square matrix.

pure subroutine rssq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=sp), intent(in   ) :: a(:,:)
  real   (kind=sp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=sp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rssq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rssq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rssq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rssq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rssq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rssq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rssq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rssq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rssq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rssq_ludcmp(qm,perm)
  call rssq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rssq_normi(qm))
  err % qm_norm(1,0) = real(rssq_norm1(qm))
  err % pm_norm(0,0) = real(rssq_normi(pm))
  err % pm_norm(1,0) = real(rssq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rssq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = rssq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rssq_normi(pm))
    err % pm_norm(1,k) = real(rssq_norm1(pm))

    call rssq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine rssq_rationalphi


!    Name : rstr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, single precision, triangular matrix.

pure subroutine rstr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=sp), intent(in   ) :: a(:,:)
  real   (kind=sp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=sp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rstr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rstr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rstr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rstr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rstr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rstr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rstr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rstr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rstr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rstr_ludcmp(qm,perm)
  call rstr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rstr_normi(qm))
  err % qm_norm(1,0) = real(rstr_norm1(qm))
  err % pm_norm(0,0) = real(rstr_normi(pm))
  err % pm_norm(1,0) = real(rstr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rstr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = rstr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rstr_normi(pm))
    err % pm_norm(1,k) = real(rstr_norm1(pm))

    call rstr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_sp)) v(j+1,j,:) = 0.0_sp
  end do
end subroutine rstr_rationalphi


!    Name : csdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, single precision, diagonal matrix.

pure subroutine csdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=sp), intent(in   ) :: a(:)
  complex(kind=sp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=sp) :: z2(size(a,1),0:4)
  complex(kind=sp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = csdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(csdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(csdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(csdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(csdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(csdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(csdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = csdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = csdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call csdg_ludcmp(qm,perm)
  call csdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(csdg_normi(qm))
  err % qm_norm(1,0) = real(csdg_norm1(qm))
  err % pm_norm(0,0) = real(csdg_normi(pm))
  err % pm_norm(1,0) = real(csdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = csdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = csdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(csdg_normi(pm))
    err % pm_norm(1,k) = real(csdg_norm1(pm))

    call csdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine csdg_rationalphi


!    Name : cssq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, single precision, square matrix.

pure subroutine cssq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=sp), intent(in   ) :: a(:,:)
  complex(kind=sp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=sp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cssq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cssq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cssq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cssq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cssq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cssq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cssq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cssq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cssq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cssq_ludcmp(qm,perm)
  call cssq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cssq_normi(qm))
  err % qm_norm(1,0) = real(cssq_norm1(qm))
  err % pm_norm(0,0) = real(cssq_normi(pm))
  err % pm_norm(1,0) = real(cssq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cssq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = cssq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cssq_normi(pm))
    err % pm_norm(1,k) = real(cssq_norm1(pm))

    call cssq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine cssq_rationalphi


!    Name : cstr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, single precision, triangular matrix.

pure subroutine cstr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=sp), intent(in   ) :: a(:,:)
  complex(kind=sp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=sp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=sp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=sp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_sp; cq = 0.0_sp
  call sp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cstr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cstr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cstr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cstr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cstr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cstr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cstr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cstr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cstr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cstr_ludcmp(qm,perm)
  call cstr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cstr_normi(qm))
  err % qm_norm(1,0) = real(cstr_norm1(qm))
  err % pm_norm(0,0) = real(cstr_normi(pm))
  err % pm_norm(1,0) = real(cstr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cstr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = cstr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cstr_normi(pm))
    err % pm_norm(1,k) = real(cstr_norm1(pm))

    call cstr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_sp)) v(j+1,j,:) = cmplx(0.0_sp,0.0_sp,sp)
  end do
end subroutine cstr_rationalphi


!    Name : rwdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, double precision, diagonal matrix.

pure subroutine rwdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=wp), intent(in   ) :: a(:)
  real   (kind=wp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=wp) :: z2(size(a,1),0:4)
  real   (kind=wp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = rwdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(rwdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(rwdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rwdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(rwdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rwdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(rwdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rwdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = rwdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call rwdg_ludcmp(qm,perm)
  call rwdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(rwdg_normi(qm))
  err % qm_norm(1,0) = real(rwdg_norm1(qm))
  err % pm_norm(0,0) = real(rwdg_normi(pm))
  err % pm_norm(1,0) = real(rwdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rwdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = rwdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rwdg_normi(pm))
    err % pm_norm(1,k) = real(rwdg_norm1(pm))

    call rwdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine rwdg_rationalphi


!    Name : rwsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, double precision, square matrix.

pure subroutine rwsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=wp), intent(in   ) :: a(:,:)
  real   (kind=wp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=wp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rwsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rwsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rwsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rwsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rwsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rwsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rwsq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rwsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rwsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rwsq_ludcmp(qm,perm)
  call rwsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rwsq_normi(qm))
  err % qm_norm(1,0) = real(rwsq_norm1(qm))
  err % pm_norm(0,0) = real(rwsq_normi(pm))
  err % pm_norm(1,0) = real(rwsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rwsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = rwsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rwsq_normi(pm))
    err % pm_norm(1,k) = real(rwsq_norm1(pm))

    call rwsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine rwsq_rationalphi


!    Name : rwtr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, double precision, triangular matrix.

pure subroutine rwtr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=wp), intent(in   ) :: a(:,:)
  real   (kind=wp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=wp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rwtr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rwtr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rwtr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rwtr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rwtr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rwtr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rwtr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rwtr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rwtr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rwtr_ludcmp(qm,perm)
  call rwtr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rwtr_normi(qm))
  err % qm_norm(1,0) = real(rwtr_norm1(qm))
  err % pm_norm(0,0) = real(rwtr_normi(pm))
  err % pm_norm(1,0) = real(rwtr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rwtr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = rwtr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rwtr_normi(pm))
    err % pm_norm(1,k) = real(rwtr_norm1(pm))

    call rwtr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_wp)) v(j+1,j,:) = 0.0_wp
  end do
end subroutine rwtr_rationalphi


!    Name : cwdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, double precision, diagonal matrix.

pure subroutine cwdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=wp), intent(in   ) :: a(:)
  complex(kind=wp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=wp) :: z2(size(a,1),0:4)
  complex(kind=wp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = cwdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(cwdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(cwdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cwdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(cwdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cwdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(cwdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cwdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = cwdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call cwdg_ludcmp(qm,perm)
  call cwdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(cwdg_normi(qm))
  err % qm_norm(1,0) = real(cwdg_norm1(qm))
  err % pm_norm(0,0) = real(cwdg_normi(pm))
  err % pm_norm(1,0) = real(cwdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cwdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = cwdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cwdg_normi(pm))
    err % pm_norm(1,k) = real(cwdg_norm1(pm))

    call cwdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine cwdg_rationalphi


!    Name : cwsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, double precision, square matrix.

pure subroutine cwsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=wp), intent(in   ) :: a(:,:)
  complex(kind=wp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=wp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cwsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cwsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cwsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cwsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cwsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cwsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cwsq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cwsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cwsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cwsq_ludcmp(qm,perm)
  call cwsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cwsq_normi(qm))
  err % qm_norm(1,0) = real(cwsq_norm1(qm))
  err % pm_norm(0,0) = real(cwsq_normi(pm))
  err % pm_norm(1,0) = real(cwsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cwsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = cwsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cwsq_normi(pm))
    err % pm_norm(1,k) = real(cwsq_norm1(pm))

    call cwsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine cwsq_rationalphi


!    Name : cwtr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, double precision, triangular matrix.

pure subroutine cwtr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=wp), intent(in   ) :: a(:,:)
  complex(kind=wp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=wp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=wp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=wp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_wp; cq = 0.0_wp
  call wp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cwtr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cwtr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cwtr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cwtr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cwtr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cwtr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cwtr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cwtr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cwtr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cwtr_ludcmp(qm,perm)
  call cwtr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cwtr_normi(qm))
  err % qm_norm(1,0) = real(cwtr_norm1(qm))
  err % pm_norm(0,0) = real(cwtr_normi(pm))
  err % pm_norm(1,0) = real(cwtr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cwtr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = cwtr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cwtr_normi(pm))
    err % pm_norm(1,k) = real(cwtr_norm1(pm))

    call cwtr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_wp)) v(j+1,j,:) = cmplx(0.0_wp,0.0_wp,wp)
  end do
end subroutine cwtr_rationalphi

#ifdef __USE_TPREC

!    Name : rtdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, triple precision, diagonal matrix.

pure subroutine rtdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=tp), intent(in   ) :: a(:)
  real   (kind=tp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=tp) :: z2(size(a,1),0:4)
  real   (kind=tp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = rtdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(rtdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(rtdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rtdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(rtdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rtdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(rtdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rtdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = rtdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call rtdg_ludcmp(qm,perm)
  call rtdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(rtdg_normi(qm))
  err % qm_norm(1,0) = real(rtdg_norm1(qm))
  err % pm_norm(0,0) = real(rtdg_normi(pm))
  err % pm_norm(1,0) = real(rtdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rtdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = rtdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rtdg_normi(pm))
    err % pm_norm(1,k) = real(rtdg_norm1(pm))

    call rtdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine rtdg_rationalphi


!    Name : rtsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, triple precision, square matrix.

pure subroutine rtsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=tp), intent(in   ) :: a(:,:)
  real   (kind=tp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=tp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rtsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rtsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rtsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rtsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rtsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rtsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rtsq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rtsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rtsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rtsq_ludcmp(qm,perm)
  call rtsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rtsq_normi(qm))
  err % qm_norm(1,0) = real(rtsq_norm1(qm))
  err % pm_norm(0,0) = real(rtsq_normi(pm))
  err % pm_norm(1,0) = real(rtsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rtsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = rtsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rtsq_normi(pm))
    err % pm_norm(1,k) = real(rtsq_norm1(pm))

    call rtsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine rtsq_rationalphi


!    Name : rttr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, triple precision, triangular matrix.

pure subroutine rttr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=tp), intent(in   ) :: a(:,:)
  real   (kind=tp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=tp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rttr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rttr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rttr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rttr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rttr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rttr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rttr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rttr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rttr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rttr_ludcmp(qm,perm)
  call rttr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rttr_normi(qm))
  err % qm_norm(1,0) = real(rttr_norm1(qm))
  err % pm_norm(0,0) = real(rttr_normi(pm))
  err % pm_norm(1,0) = real(rttr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rttr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = rttr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rttr_normi(pm))
    err % pm_norm(1,k) = real(rttr_norm1(pm))

    call rttr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_tp)) v(j+1,j,:) = 0.0_tp
  end do
end subroutine rttr_rationalphi


!    Name : ctdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, triple precision, diagonal matrix.

pure subroutine ctdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=tp), intent(in   ) :: a(:)
  complex(kind=tp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=tp) :: z2(size(a,1),0:4)
  complex(kind=tp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = ctdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(ctdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(ctdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(ctdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(ctdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(ctdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(ctdg_norm1(z2(:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = ctdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = ctdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call ctdg_ludcmp(qm,perm)
  call ctdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(ctdg_normi(qm))
  err % qm_norm(1,0) = real(ctdg_norm1(qm))
  err % pm_norm(0,0) = real(ctdg_normi(pm))
  err % pm_norm(1,0) = real(ctdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = ctdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = ctdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(ctdg_normi(pm))
    err % pm_norm(1,k) = real(ctdg_norm1(pm))

    call ctdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine ctdg_rationalphi


!    Name : ctsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, triple precision, square matrix.

pure subroutine ctsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=tp), intent(in   ) :: a(:,:)
  complex(kind=tp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=tp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = ctsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(ctsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(ctsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(ctsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(ctsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(ctsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(ctsq_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = ctsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = ctsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call ctsq_ludcmp(qm,perm)
  call ctsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(ctsq_normi(qm))
  err % qm_norm(1,0) = real(ctsq_norm1(qm))
  err % pm_norm(0,0) = real(ctsq_normi(pm))
  err % pm_norm(1,0) = real(ctsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = ctsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = ctsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(ctsq_normi(pm))
    err % pm_norm(1,k) = real(ctsq_norm1(pm))

    call ctsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine ctsq_rationalphi


!    Name : cttr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, triple precision, triangular matrix.

pure subroutine cttr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=tp), intent(in   ) :: a(:,:)
  complex(kind=tp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=tp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=tp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=tp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_tp; cq = 0.0_tp
  call tp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cttr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cttr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cttr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cttr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cttr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cttr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cttr_norm1(z2(:,:,3)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cttr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cttr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cttr_ludcmp(qm,perm)
  call cttr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cttr_normi(qm))
  err % qm_norm(1,0) = real(cttr_norm1(qm))
  err % pm_norm(0,0) = real(cttr_normi(pm))
  err % pm_norm(1,0) = real(cttr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cttr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = cttr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cttr_normi(pm))
    err % pm_norm(1,k) = real(cttr_norm1(pm))

    call cttr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_tp)) v(j+1,j,:) = cmplx(0.0_tp,0.0_tp,tp)
  end do
end subroutine cttr_rationalphi

#endif
#ifdef __USE_QPREC

!    Name : rqdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, quadruple precision, diagonal matrix.

pure subroutine rqdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=qp), intent(in   ) :: a(:)
  real   (kind=qp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=qp) :: z2(size(a,1),0:4)
  real   (kind=qp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = rqdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(rqdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(rqdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rqdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(rqdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rqdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(rqdg_norm1(z2(:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(rqdg_normi(z2(:,4)))
    err % zp_norm(1,8) = real(rqdg_norm1(z2(:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rqdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = rqdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call rqdg_ludcmp(qm,perm)
  call rqdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(rqdg_normi(qm))
  err % qm_norm(1,0) = real(rqdg_norm1(qm))
  err % pm_norm(0,0) = real(rqdg_normi(pm))
  err % pm_norm(1,0) = real(rqdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rqdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = rqdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rqdg_normi(pm))
    err % pm_norm(1,k) = real(rqdg_norm1(pm))

    call rqdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine rqdg_rationalphi


!    Name : rqsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, quadruple precision, square matrix.

pure subroutine rqsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=qp), intent(in   ) :: a(:,:)
  real   (kind=qp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=qp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rqsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rqsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rqsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rqsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rqsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rqsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rqsq_norm1(z2(:,:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(rqsq_normi(z2(:,:,4)))
    err % zp_norm(1,8) = real(rqsq_norm1(z2(:,:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rqsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rqsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rqsq_ludcmp(qm,perm)
  call rqsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rqsq_normi(qm))
  err % qm_norm(1,0) = real(rqsq_norm1(qm))
  err % pm_norm(0,0) = real(rqsq_normi(pm))
  err % pm_norm(1,0) = real(rqsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rqsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = rqsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rqsq_normi(pm))
    err % pm_norm(1,k) = real(rqsq_norm1(pm))

    call rqsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine rqsq_rationalphi


!    Name : rqtr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a real, quadruple precision, triangular matrix.

pure subroutine rqtr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  real   (kind=qp), intent(in   ) :: a(:,:)
  real   (kind=qp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  real   (kind=qp) :: z2(size(a,1),size(a,2),0:4)
  real   (kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = rqtr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(rqtr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(rqtr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(rqtr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(rqtr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(rqtr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(rqtr_norm1(z2(:,:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(rqtr_normi(z2(:,:,4)))
    err % zp_norm(1,8) = real(rqtr_norm1(z2(:,:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = rqtr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = rqtr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call rqtr_ludcmp(qm,perm)
  call rqtr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(rqtr_normi(qm))
  err % qm_norm(1,0) = real(rqtr_norm1(qm))
  err % pm_norm(0,0) = real(rqtr_normi(pm))
  err % pm_norm(1,0) = real(rqtr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = rqtr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = rqtr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(rqtr_normi(pm))
    err % pm_norm(1,k) = real(rqtr_norm1(pm))

    call rqtr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_qp)) v(j+1,j,:) = 0.0_qp
  end do
end subroutine rqtr_rationalphi


!    Name : cqdg_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, quadruple precision, diagonal matrix.

pure subroutine cqdg_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=qp), intent(in   ) :: a(:)
  complex(kind=qp), intent(out  ) :: v(size(a,1),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=qp) :: z2(size(a,1),0:4)
  complex(kind=qp) :: qm(size(a,1)), pm(size(a,1))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,0) = cqdg_eye(size(a,1))
  z2(:,1) = a .dgtimes. a

  if ( 5 <= m             ) z2(:,2) = z2(:,1) .dgtimes. z2(:,1)
  if ( 7 <= m .and. m /= 9) z2(:,3) = z2(:,2) .dgtimes. z2(:,1)
  if (17 <= m             ) z2(:,4) = z2(:,3) .dgtimes. z2(:,1)

  err % zp_norm(0,2) = real(cqdg_normi(z2(:,1)))
  err % zp_norm(1,2) = real(cqdg_norm1(z2(:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cqdg_normi(z2(:,2)))
    err % zp_norm(1,4) = real(cqdg_norm1(z2(:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cqdg_normi(z2(:,3)))
    err % zp_norm(1,6) = real(cqdg_norm1(z2(:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(cqdg_normi(z2(:,4)))
    err % zp_norm(1,8) = real(cqdg_norm1(z2(:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cqdg_poly(n,cq(0,1:2*n+1:2),z2) .dgtimes. a     ! the odd part
  pm = -qm

  v(:,0) = cqdg_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,0) + qm
  pm = v(:,0) + pm

  call cqdg_ludcmp(qm,perm)
  call cqdg_ainvb(qm,perm,pm,v(:,0))

  err % qm_norm(0,0) = real(cqdg_normi(qm))
  err % qm_norm(1,0) = real(cqdg_norm1(qm))
  err % pm_norm(0,0) = real(cqdg_normi(pm))
  err % pm_norm(1,0) = real(cqdg_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cqdg_poly(n,cp(k,1:2*n+1:2),z2) .dgtimes. a
    pm = cqdg_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cqdg_normi(pm))
    err % pm_norm(1,k) = real(cqdg_norm1(pm))

    call cqdg_ainvb(qm,perm,pm,v(:,k))
  end do

end subroutine cqdg_rationalphi


!    Name : cqsq_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, quadruple precision, square matrix.

pure subroutine cqsq_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=qp), intent(in   ) :: a(:,:)
  complex(kind=qp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=qp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cqsq_eye(size(a,1))
  z2(:,:,1) = a .sqtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .sqtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .sqtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .sqtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cqsq_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cqsq_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cqsq_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cqsq_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cqsq_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cqsq_norm1(z2(:,:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(cqsq_normi(z2(:,:,4)))
    err % zp_norm(1,8) = real(cqsq_norm1(z2(:,:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cqsq_poly(n,cq(0,1:2*n+1:2),z2) .sqtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cqsq_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cqsq_ludcmp(qm,perm)
  call cqsq_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cqsq_normi(qm))
  err % qm_norm(1,0) = real(cqsq_norm1(qm))
  err % pm_norm(0,0) = real(cqsq_normi(pm))
  err % pm_norm(1,0) = real(cqsq_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cqsq_poly(n,cp(k,1:2*n+1:2),z2) .sqtimes. a
    pm = cqsq_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cqsq_normi(pm))
    err % pm_norm(1,k) = real(cqsq_norm1(pm))

    call cqsq_ainvb(qm,perm,pm,v(:,:,k))
  end do

end subroutine cqsq_rationalphi


!    Name : cqtr_rationalphi
! Purpose : This subroutine computes rational approximants
!         : for phi-functions.
!   Input : - ``m'' must be one of {3,5,7,9,13,17}, and specifies the order
!         :   of approximation. An [m/m] rational approximant is used.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!         : - ``A'' is the matrix argument.
!  Output : - v(*,k) stores the value of the rational function
!         :   that approximates phi_k(A).
!         : - err stores error informations.
!    Note : This routine treats a complex, quadruple precision, triangular matrix.

pure subroutine cqtr_rationalphi(m,upto,a,v,err)
  integer,          intent(in   ) :: m, upto
  complex(kind=qp), intent(in   ) :: a(:,:)
  complex(kind=qp), intent(out  ) :: v(size(a,1),size(a,2),0:upto)
  type(mcpsqrlog),  intent(inout) :: err

  integer          :: j, k, n, perm(size(a,1))
  real   (kind=qp) :: cq(0:0,0:m), cp(0:upto,0:m)
  complex(kind=qp) :: z2(size(a,1),size(a,2),0:4)
  complex(kind=qp) :: qm(size(a,1),size(a,2)), pm(size(a,1),size(a,2))

  ! Coefficients of polynomials that define rational approximations.

  n = m / 2                             ! ``m'' must be an odd integer.
  cp = 0.0_qp; cq = 0.0_qp
  call qp_mcpcoeff(m,n,upto,cq,cp)

  ! Even powers of the matrix, i.e., a^0, a^2, a^4, a^6, and a^8

  z2(:,:,0) = cqtr_eye(size(a,1))
  z2(:,:,1) = a .trtimes. a

  if ( 5 <= m             ) z2(:,:,2) = z2(:,:,1) .trtimes. z2(:,:,1)
  if ( 7 <= m .and. m /= 9) z2(:,:,3) = z2(:,:,2) .trtimes. z2(:,:,1)
  if (17 <= m             ) z2(:,:,4) = z2(:,:,3) .trtimes. z2(:,:,1)

  err % zp_norm(0,2) = real(cqtr_normi(z2(:,:,1)))
  err % zp_norm(1,2) = real(cqtr_norm1(z2(:,:,1)))

  if ( 5 <= m ) then
    err % zp_norm(0,4) = real(cqtr_normi(z2(:,:,2)))
    err % zp_norm(1,4) = real(cqtr_norm1(z2(:,:,2)))
  end if
  if ( 7 <= m .and. m /= 9) then
    err % zp_norm(0,6) = real(cqtr_normi(z2(:,:,3)))
    err % zp_norm(1,6) = real(cqtr_norm1(z2(:,:,3)))
  end if
  if ( 17 <= m ) then
    err % zp_norm(0,8) = real(cqtr_normi(z2(:,:,4)))
    err % zp_norm(1,8) = real(cqtr_norm1(z2(:,:,4)))
  end if

  ! Diagonal Pad{\'{e}} approximant for the matrix exponential.
  ! ``qm'' is the denominator and ``pm'' is the numerator;
  ! hence, qm(A)^{-1} pm(A) is the required rational approximant.
  ! The summations of odd terms and even terms of these polynomials
  ! are computed separately to improve efficiency.

  qm = cqtr_poly(n,cq(0,1:2*n+1:2),z2) .trtimes. a     ! the odd part
  pm = -qm

  v(:,:,0) = cqtr_poly(n,cq(0,0:2*n:2),z2)            ! the even part
  qm = v(:,:,0) + qm
  pm = v(:,:,0) + pm

  call cqtr_ludcmp(qm,perm)
  call cqtr_ainvb(qm,perm,pm,v(:,:,0))

  err % qm_norm(0,0) = real(cqtr_normi(qm))
  err % qm_norm(1,0) = real(cqtr_norm1(qm))
  err % pm_norm(0,0) = real(cqtr_normi(pm))
  err % pm_norm(1,0) = real(cqtr_norm1(pm))

  ! Rational approximations for phi-functions, phi_1, phi_2, ... , and phi_5.
  ! We use the same denominator ``qm'' as that for the exponential.
  ! The summations of odd terms and even terms of each numerator ``pm''
  ! are computed separately to improve efficiency.

  do k=1,upto
    pm = cqtr_poly(n,cp(k,1:2*n+1:2),z2) .trtimes. a
    pm = cqtr_poly(n,cp(k,0:2*n  :2),z2) + pm

    err % pm_norm(0,k) = real(cqtr_normi(pm))
    err % pm_norm(1,k) = real(cqtr_norm1(pm))

    call cqtr_ainvb(qm,perm,pm,v(:,:,k))
  end do

  ! for the block triangular structure to be invariant

  do j=1,size(a,1)-1
    if (abs(a(j+1,j)) < tiny(1.0_qp)) v(j+1,j,:) = cmplx(0.0_qp,0.0_qp,qp)
  end do
end subroutine cqtr_rationalphi

#endif

end module rationalphi


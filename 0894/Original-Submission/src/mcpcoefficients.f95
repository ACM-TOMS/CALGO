! This module computes the coefficients of the polynomials
! that define the numerators and the denominator of the
! rational approximants.
!
! References:
!
!   Nicholas J. Higham,
!   The scaling and squaring method for the matrix
!   exponential revisited,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 26, Number 4, pp. 1179-1193, 2005.
!
!   S. Koikari,
!   An error analysis of the modified scaling and
!   squaring method,
!   Computers and Mathematics with Applications,
!   Volume 53, pp. 1293-1305, 2007.
!
! This module is intended for internal use only.

module mcpcoefficients

  use floattypes

  implicit none

  public

contains


!    Name : sp_mcpcoeff
! Purpose : This subroutine computes the coefficients of the polynomials
!         : that define the numerators and the denominator of the rational
!         : approximants.
!   Input : - ``m'' is the order of the polynomials. [m/m] approximant is used.
!         :   ``m'' must be one of { 3, 5, 7, 9, 13, 17 }
!         : - ``n'' is m/2.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!  Output : - cq(0,j) stores the coefficient of the j-th power term
!         :   in the denominator of the diagonal Pad{\'{e}} approximant
!         :   that approximates the exponential.
!         : - cp(i,j) stores the coefficient of the j-th power term in the
!         :   numerator of the rational function that approximates phi_i.

pure subroutine sp_mcpcoeff(m,n,upto,cq,cp)
  integer,       intent(in ) :: m, n, upto
  real(kind=sp), intent(out) :: cq(0:0,0:m), cp(0:upto,0:m)

  integer       :: k
  real(kind=sp) :: temp, gmk, dnm

  cp = 0.0_sp; cq = 0.0_sp

  do k=0,n                              ! coefficients of the polynomials

    temp = sp_fctr(2*m-2*k-1) / sp_fctr(2*m) * sp_combi(m,2*k+1)

    cp(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! pm(z)
    cp(0,2*k+1) =  temp

    cq(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! qm(z)
    cq(0,2*k+1) = -temp
    if (upto <= 0) cycle

    cp(1,2*k  ) =  temp * 2                                             ! phi_1
    cp(1,2*k+1) =  0.0_sp
    if (upto <= 1) cycle

    cp(2,2*k  ) =  temp                                                 ! phi_2
    cp(2,2*k+1) = -temp  * (m-2*k-1)/(2*(2*k+3)*(m-k-1))
    if (upto <= 2) cycle

    cp(3,2*k  ) =  temp * (2*k*m+2*m-2*k*k-3*k-2)/(2*(2*k+3)*(m-k-1))   ! phi_3
    cp(3,2*k+1) = -temp * (m-2*k-1)/(4*(2*k+3)*(m-k-1))
    if (upto <= 3) cycle

    if (m == 3 .and. k == 1) then                                       ! phi_4
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp *           (4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)        )
    else
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp * (m-2*k-1)*(4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2))
    end if
    if (upto <= 4) cycle

#if 0
    if (m == 3 .and. k == 1) then                                       ! phi_5
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)        

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp *           (2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    else
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    end if
#else
    gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
         -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
    dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

    cp(5,2*k  ) =  temp * gmk / dnm
    cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
#endif
  end do

contains


  !    Name : sp_fctr
  ! Purpose : This function computes the factorial n!.
  !   Input : ``n'' is a non-negative integer.
  !  Output : n!.

  pure function sp_fctr(n)              ! n!, the factorial.
    integer, intent(in) :: n
    real(kind=sp)       :: sp_fctr
    integer             :: k

    sp_fctr = 1.0_sp
    do k=2,n; sp_fctr = sp_fctr * k; end do
  end function sp_fctr


  !    Name : sp_combi
  ! Purpose : This function computes m!/((m-n)!n!).
  !   Input : ``m'' and ``n'' are non-negative integers.
  !  Output : m!/((m-n)!n!).

  pure function sp_combi(m,n)           ! mCn := m!/((m-n)!n!)
    integer, intent(in) :: m,n
    real(kind=sp)       :: sp_combi

    sp_combi = sp_fctr(m) / (sp_fctr(m-n) * sp_fctr(n))
  end function sp_combi

end subroutine sp_mcpcoeff


!    Name : wp_mcpcoeff
! Purpose : This subroutine computes the coefficients of the polynomials
!         : that define the numerators and the denominator of the rational
!         : approximants.
!   Input : - ``m'' is the order of the polynomials. [m/m] approximant is used.
!         :   ``m'' must be one of { 3, 5, 7, 9, 13, 17 }
!         : - ``n'' is m/2.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!  Output : - cq(0,j) stores the coefficient of the j-th power term
!         :   in the denominator of the diagonal Pad{\'{e}} approximant
!         :   that approximates the exponential.
!         : - cp(i,j) stores the coefficient of the j-th power term in the
!         :   numerator of the rational function that approximates phi_i.

pure subroutine wp_mcpcoeff(m,n,upto,cq,cp)
  integer,       intent(in ) :: m, n, upto
  real(kind=wp), intent(out) :: cq(0:0,0:m), cp(0:upto,0:m)

  integer       :: k
  real(kind=wp) :: temp, gmk, dnm

  cp = 0.0_wp; cq = 0.0_wp

  do k=0,n                              ! coefficients of the polynomials

    temp = wp_fctr(2*m-2*k-1) / wp_fctr(2*m) * wp_combi(m,2*k+1)

    cp(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! pm(z)
    cp(0,2*k+1) =  temp

    cq(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! qm(z)
    cq(0,2*k+1) = -temp
    if (upto <= 0) cycle

    cp(1,2*k  ) =  temp * 2                                             ! phi_1
    cp(1,2*k+1) =  0.0_wp
    if (upto <= 1) cycle

    cp(2,2*k  ) =  temp                                                 ! phi_2
    cp(2,2*k+1) = -temp  * (m-2*k-1)/(2*(2*k+3)*(m-k-1))
    if (upto <= 2) cycle

    cp(3,2*k  ) =  temp * (2*k*m+2*m-2*k*k-3*k-2)/(2*(2*k+3)*(m-k-1))   ! phi_3
    cp(3,2*k+1) = -temp * (m-2*k-1)/(4*(2*k+3)*(m-k-1))
    if (upto <= 3) cycle

    if (m == 3 .and. k == 1) then                                       ! phi_4
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp *           (4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)        )
    else
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp * (m-2*k-1)*(4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2))
    end if
    if (upto <= 4) cycle

#if 0
    if (m == 3 .and. k == 1) then                                       ! phi_5
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)        

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp *           (2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    else
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    end if
#else
    gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
         -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
    dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

    cp(5,2*k  ) =  temp * gmk / dnm
    cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
#endif
  end do

contains


  !    Name : wp_fctr
  ! Purpose : This function computes the factorial n!.
  !   Input : ``n'' is a non-negative integer.
  !  Output : n!.

  pure function wp_fctr(n)              ! n!, the factorial.
    integer, intent(in) :: n
    real(kind=wp)       :: wp_fctr
    integer             :: k

    wp_fctr = 1.0_wp
    do k=2,n; wp_fctr = wp_fctr * k; end do
  end function wp_fctr


  !    Name : wp_combi
  ! Purpose : This function computes m!/((m-n)!n!).
  !   Input : ``m'' and ``n'' are non-negative integers.
  !  Output : m!/((m-n)!n!).

  pure function wp_combi(m,n)           ! mCn := m!/((m-n)!n!)
    integer, intent(in) :: m,n
    real(kind=wp)       :: wp_combi

    wp_combi = wp_fctr(m) / (wp_fctr(m-n) * wp_fctr(n))
  end function wp_combi

end subroutine wp_mcpcoeff

#ifdef __USE_TPREC

!    Name : tp_mcpcoeff
! Purpose : This subroutine computes the coefficients of the polynomials
!         : that define the numerators and the denominator of the rational
!         : approximants.
!   Input : - ``m'' is the order of the polynomials. [m/m] approximant is used.
!         :   ``m'' must be one of { 3, 5, 7, 9, 13, 17 }
!         : - ``n'' is m/2.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!  Output : - cq(0,j) stores the coefficient of the j-th power term
!         :   in the denominator of the diagonal Pad{\'{e}} approximant
!         :   that approximates the exponential.
!         : - cp(i,j) stores the coefficient of the j-th power term in the
!         :   numerator of the rational function that approximates phi_i.

pure subroutine tp_mcpcoeff(m,n,upto,cq,cp)
  integer,       intent(in ) :: m, n, upto
  real(kind=tp), intent(out) :: cq(0:0,0:m), cp(0:upto,0:m)

  integer       :: k
  real(kind=tp) :: temp, gmk, dnm

  cp = 0.0_tp; cq = 0.0_tp

  do k=0,n                              ! coefficients of the polynomials

    temp = tp_fctr(2*m-2*k-1) / tp_fctr(2*m) * tp_combi(m,2*k+1)

    cp(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! pm(z)
    cp(0,2*k+1) =  temp

    cq(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! qm(z)
    cq(0,2*k+1) = -temp
    if (upto <= 0) cycle

    cp(1,2*k  ) =  temp * 2                                             ! phi_1
    cp(1,2*k+1) =  0.0_tp
    if (upto <= 1) cycle

    cp(2,2*k  ) =  temp                                                 ! phi_2
    cp(2,2*k+1) = -temp  * (m-2*k-1)/(2*(2*k+3)*(m-k-1))
    if (upto <= 2) cycle

    cp(3,2*k  ) =  temp * (2*k*m+2*m-2*k*k-3*k-2)/(2*(2*k+3)*(m-k-1))   ! phi_3
    cp(3,2*k+1) = -temp * (m-2*k-1)/(4*(2*k+3)*(m-k-1))
    if (upto <= 3) cycle

    if (m == 3 .and. k == 1) then                                       ! phi_4
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp *           (4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)        )
    else
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp * (m-2*k-1)*(4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2))
    end if
    if (upto <= 4) cycle

#if 0
    if (m == 3 .and. k == 1) then                                       ! phi_5
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)        

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp *           (2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    else
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    end if
#else
    gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
         -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
    dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

    cp(5,2*k  ) =  temp * gmk / dnm
    cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
#endif
  end do

contains


  !    Name : tp_fctr
  ! Purpose : This function computes the factorial n!.
  !   Input : ``n'' is a non-negative integer.
  !  Output : n!.

  pure function tp_fctr(n)              ! n!, the factorial.
    integer, intent(in) :: n
    real(kind=tp)       :: tp_fctr
    integer             :: k

    tp_fctr = 1.0_tp
    do k=2,n; tp_fctr = tp_fctr * k; end do
  end function tp_fctr


  !    Name : tp_combi
  ! Purpose : This function computes m!/((m-n)!n!).
  !   Input : ``m'' and ``n'' are non-negative integers.
  !  Output : m!/((m-n)!n!).

  pure function tp_combi(m,n)           ! mCn := m!/((m-n)!n!)
    integer, intent(in) :: m,n
    real(kind=tp)       :: tp_combi

    tp_combi = tp_fctr(m) / (tp_fctr(m-n) * tp_fctr(n))
  end function tp_combi

end subroutine tp_mcpcoeff

#endif

#ifdef __USE_QPREC

!    Name : qp_mcpcoeff
! Purpose : This subroutine computes the coefficients of the polynomials
!         : that define the numerators and the denominator of the rational
!         : approximants.
!   Input : - ``m'' is the order of the polynomials. [m/m] approximant is used.
!         :   ``m'' must be one of { 3, 5, 7, 9, 13, 17 }
!         : - ``n'' is m/2.
!         : - ``upto'' is the maximum index of phi-functions to be computed.
!  Output : - cq(0,j) stores the coefficient of the j-th power term
!         :   in the denominator of the diagonal Pad{\'{e}} approximant
!         :   that approximates the exponential.
!         : - cp(i,j) stores the coefficient of the j-th power term in the
!         :   numerator of the rational function that approximates phi_i.

pure subroutine qp_mcpcoeff(m,n,upto,cq,cp)
  integer,       intent(in ) :: m, n, upto
  real(kind=qp), intent(out) :: cq(0:0,0:m), cp(0:upto,0:m)

  integer       :: k
  real(kind=qp) :: temp, gmk, dnm

  cp = 0.0_qp; cq = 0.0_qp

  do k=0,n                              ! coefficients of the polynomials

    temp = qp_fctr(2*m-2*k-1) / qp_fctr(2*m) * qp_combi(m,2*k+1)

    cp(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! pm(z)
    cp(0,2*k+1) =  temp

    cq(0,2*k  ) =  temp * (2*(m-k)*(2*k+1))/(m-2*k)                     ! qm(z)
    cq(0,2*k+1) = -temp
    if (upto <= 0) cycle

    cp(1,2*k  ) =  temp * 2                                             ! phi_1
    cp(1,2*k+1) =  0.0_qp
    if (upto <= 1) cycle

    cp(2,2*k  ) =  temp                                                 ! phi_2
    cp(2,2*k+1) = -temp  * (m-2*k-1)/(2*(2*k+3)*(m-k-1))
    if (upto <= 2) cycle

    cp(3,2*k  ) =  temp * (2*k*m+2*m-2*k*k-3*k-2)/(2*(2*k+3)*(m-k-1))   ! phi_3
    cp(3,2*k+1) = -temp * (m-2*k-1)/(4*(2*k+3)*(m-k-1))
    if (upto <= 3) cycle

    if (m == 3 .and. k == 1) then                                       ! phi_4
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp *           (4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)        )
    else
      cp(4,2*k  ) =  temp * (4*k*m+3*m-4*k*k-4*k-3)/(12*(2*k+3)*(m-k-1))
      cp(4,2*k+1) = -temp * (m-2*k-1)*(4*k*m+9*m-4*k*k-16*k-18)         &
                          / (24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2))
    end if
    if (upto <= 4) cycle

#if 0
    if (m == 3 .and. k == 1) then                                       ! phi_5
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)        

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp *           (2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    else
      gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
           -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
      dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

      cp(5,2*k  ) =  temp * gmk / dnm
      cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
    end if
#else
    gmk = 4*k*k*m*m + 12*k*m*m + 6*m*m - 8*(k**3)*m - 32*k*k*m        &
         -40*k*m - 18*m + 4*(k**4) + 20*(k**3) + 35*k*k + 25*k + 12
    dnm =24*(2*k+3)*(2*k+5)*(m-k-1)*(m-k-2)

    cp(5,2*k  ) =  temp * gmk / dnm
    cp(5,2*k+1) = -temp * (m-2*k-1)*(2*k*m+4*m-2*k*k-7*k-8) / (2*dnm)
#endif
  end do

contains


  !    Name : qp_fctr
  ! Purpose : This function computes the factorial n!.
  !   Input : ``n'' is a non-negative integer.
  !  Output : n!.

  pure function qp_fctr(n)              ! n!, the factorial.
    integer, intent(in) :: n
    real(kind=qp)       :: qp_fctr
    integer             :: k

    qp_fctr = 1.0_qp
    do k=2,n; qp_fctr = qp_fctr * k; end do
  end function qp_fctr


  !    Name : qp_combi
  ! Purpose : This function computes m!/((m-n)!n!).
  !   Input : ``m'' and ``n'' are non-negative integers.
  !  Output : m!/((m-n)!n!).

  pure function qp_combi(m,n)           ! mCn := m!/((m-n)!n!)
    integer, intent(in) :: m,n
    real(kind=qp)       :: qp_combi

    qp_combi = qp_fctr(m) / (qp_fctr(m-n) * qp_fctr(n))
  end function qp_combi

end subroutine qp_mcpcoeff

#endif

end module mcpcoefficients


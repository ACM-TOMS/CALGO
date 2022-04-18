    MODULE humlicekmod
  USE set_rk, ONLY : rk
   implicit none
    private
    public :: humlicek0

    CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!
!
      ELEMENTAL SUBROUTINE humlicek0(z,prbfct, dVdx, dVdy)
!      ELEMENTAL SUBROUTINE humlicek0(z,prbfct)
!!
!
!!    complex probability function for complex argument Z=X+iY
!
!!    real part = voigt function K(x,y)
!
!!    source:   j. humlicek, JQSRT 27, 437, 1982
!
!!    the stated accuracy is claimed to be 1.0E-04 by the author.
!
!!    r.h.norton has checked the accuracy by comparing values computed
!!    using a program written by b.h.armstrong,
!!    and the accuracy claim seems to be warranted.
!
!!    12/91, converted to f90 may 2009
!!    march 2015, converted to elemental function;
!!
!!                parameterized precision (change r8 to r4 for single
!!                  precision version);
!!
!!                use parameter arrays for rational function
!!                  coefficients
!!
!!                runs against other approximations confirm 4 dp
!!                  accuracy claim
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! .. Use Statements ..
!        USE set_rk, ONLY : rk=>rk
! ..
! .. Parameters ..
        REAL (rk), PARAMETER :: s15 = 15.E0_rk, s55 = 5.5E0_rk
        REAL (rk), PARAMETER :: half = 0.5E0_rk, &
                                recsqrtpi = 5.641895835477563e-001_rk
        REAL (rk), PARAMETER :: a2_1 = 1.410474e0_rk, a2_2 =0.75E0_rk, &
                                a2_3 = 3.0E0_rk
        REAL (rk), PARAMETER :: a3_1 = 16.4955E0_rk,  a3_2 =20.20933E0_rk, &
            a3_3=11.96482E0_rk, a3_4 = 3.778987E0_rk, a3_5 =0.5642236E0_rk, &
            a3_6 =16.4955E0_rk, a3_7 = 38.82363E0_rk, a3_8 =39.27121E0_rk, &
            a3_9 =21.69274E0_rk,a3_10= 6.699398E0_rk
        REAL (rk), PARAMETER :: a4_1 = 36183.31E0_rk, a4_2=3321.99E0_rk, &
        a4_3 = 1540.787E0_rk,    a4_4 =219.031E0_rk, a4_5 =35.7668E0_rk, a4_6 =1.320522E0_rk,&
        a4_7 = 32066.6E0_rk, a4_8 =24322.8E0_rk, a4_9 =9022.23E0_rk, a4_10=2186.18E0_rk,&
        a4_11= 364.219E0_rk, a4_12=61.5704E0_rk, a4_13=1.84144E0_rk


! ..
! .. Scalar Arguments ..
        COMPLEX (rk), INTENT (OUT) :: prbfct
        COMPLEX (rk), INTENT (IN) :: z
        real(rk), intent(out), optional :: dVdx, dVdy
! ..
! .. Local Scalars ..
        COMPLEX (rk) :: t, u
        REAL (rk) :: ax, s, x, y
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, aimag, cmplx, exp, real
! ..
        x = real(z,kind=rk)
        y = aimag(z)
        t = cmplx(y,-x,rk)
        ax = abs(x)
        s = ax + y
        IF (s>=s15) THEN
!          region I
          prbfct = t*recsqrtpi/(half+t*t)
        ELSE IF (s<s15 .AND. s>=s55 ) THEN
!             ELSE IF (s<s15 .AND. s>=s55 .and. y>1.0e-6_rk) THEN   !Mofreh's correction
!          region II
          u = t*t
          prbfct = (t*(a2_1+u*recsqrtpi))/(a2_2+(u*(a2_3+u)))
        ELSE IF (s<s55 .AND. y>=(0.195_rk*ax-0.176_rk)) THEN
!         region III
          prbfct = (a3_1+t*(a3_2+t*(a3_3+ &
            t*(a3_4+a3_5*t))))/(a3_6+t*(a3_7+t*(a3_8+ &
            t*(a3_9+t*(a3_10+t)))))
        ELSE
!         region IV
          u = t*t
          prbfct = exp(u) - (t*(a4_1-u*(a4_2-u*(a4_3-u*(a4_4-u*( &
            a4_5-u*(a4_6-u*recsqrtpi))))))/(a4_7-u*(a4_8- &
            u*(a4_9-u*(a4_10-u*(a4_11-u*(a4_12-u*(a4_13-u))))))))
        END IF

!  --------------- Calculation of the derivatives
            ! Partial derivative of V w.r.t. x       =dLdy
            if ( present(dVdx) ) dVdx=-2.0_rk*real(z*prbfct)
            ! Partial derivative of V w.r.t. y       =-dLdx
            if ( present(dVdy) ) dVdy=2.0_rk*aimag(z*prbfct)-1.12837916709551_rk

      END SUBROUTINE humlicek0
    END MODULE humlicekmod

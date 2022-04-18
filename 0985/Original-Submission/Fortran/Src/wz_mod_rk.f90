MODULE wz_mod_rk
    USE set_rk, ONLY:rk

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: wz_rk

    REAL(rk), PARAMETER :: half   = 0.5_rk
    REAL(rk), PARAMETER :: three_halfs  = 1.5_rk
    REAL(rk), PARAMETER :: two    = 2.0_rk
    REAL(rk), PARAMETER :: four   = 4.0_rk
    REAL(rk), PARAMETER :: one_sqrt_pi = 5.641895835477563e-001_rk   !1/sqrt(pi)

    !:- Hui's p6 coefficients
    REAL(rk), PARAMETER :: aa0 = 122.60793_rk
    REAL(rk), PARAMETER :: aa1 = 214.38239_rk
    REAL(rk), PARAMETER :: aa2 = 181.92853_rk
    REAL(rk), PARAMETER :: aa3 =  93.15558_rk
    REAL(rk), PARAMETER :: aa4 =  30.180142_rk
    REAL(rk), PARAMETER :: aa5 =   5.9126262_rk
    REAL(rk), PARAMETER :: aa6 = one_sqrt_pi

    REAL(rk), PARAMETER :: bb0 = 122.60793_rk
    REAL(rk), PARAMETER :: bb1 = 352.73063_rk
    REAL(rk), PARAMETER :: bb2 = 457.33448_rk
    REAL(rk), PARAMETER :: bb3 = 348.70392_rk
    REAL(rk), PARAMETER :: bb4 = 170.35400_rk
    REAL(rk), PARAMETER :: bb5 =  53.992907_rk
    REAL(rk), PARAMETER :: bb6 =  10.479857_rk

    !:- Humlicek's Region IV coeeficients
    REAL(rk), PARAMETER :: ha0 = 36183.31_rk
    REAL(rk), PARAMETER :: ha1 =  3321.99_rk
    REAL(rk), PARAMETER :: ha2 = 1540.787_rk
    REAL(rk), PARAMETER :: ha3 =  219.031_rk
    REAL(rk), PARAMETER :: ha4 =  35.7668_rk
    REAL(rk), PARAMETER :: ha5 = 1.320522_rk
    REAL(rk), PARAMETER :: ha6=   one_sqrt_pi

    REAL(rk), PARAMETER :: hb0 =  32066.6_rk
    REAL(rk), PARAMETER :: hb1 = 24322.84_rk
    REAL(rk), PARAMETER :: hb2 = 9022.228_rk
    REAL(rk), PARAMETER :: hb3 = 2186.181_rk
    REAL(rk), PARAMETER :: hb4 = 364.2191_rk
    REAL(rk), PARAMETER :: hb5 = 61.57037_rk
    REAL(rk), PARAMETER :: hb6 = 1.841439_rk

    REAL(rk), PARAMETER :: two_sqrt_pi = 1.12837916709551_rk   !2/sqrt(pi)
    COMPLEX(rk), PARAMETER :: j1= (0.0_rk,1.0_rk)

CONTAINS

    elemental subroutine wz_rk( z, w, dVdx, dVdy)
        ! ----------
        !  wz_rk is an elemental Fortran subroutine that receives, as input, a
        !  single dummy scalar complex number (z=x+iy) and returns as output
        !  the Faddeyeva function,w, defined as w(z)=exp(-z^2)*erfc(-i*z) where erfc(z)
        !  is the complex complementary error function. The routine approximates the
        !  the function to accuracy <4.0x10^-5 (for both real and imaginary parts)
        !  In addition, the subroutine returns as optional the derivatives
        !  of the real part of the Faddeyeva function (dVdx and dVdy) with respect to
        !  the real and imaginary parts of z, respectively.
        !------------
        !  This routine depends on partitioning the computational domain into
        !  six regions and the use of Laplace continued fractions, Hui's p6-
        !  approximation [JQSRT, Vol. 19, 509-516 (1978)] together with Humlicek's
        !  approximation for region 4 in Humlicek's w4 algorithm [JQSRT, Vol. 27,
        !  No. 4, 437-444 (1982)].
        !------------
        !  The routine can be run in single precision or in double precision
        !  depending on the choice of the integer parameter rk in the
        !  subsidiary module set_rk
        !-------------
        !  The accompanying driver code “wz_driver_rk.f90” can be run for
        !  computation of the Faddeyeva function w(z)=V+iL and the partial
        !  derivatives of its real part, V(x,y) (optional), on a scalar or an array of
        !  the complex variable z. The partial derivatives of the imaginary part, L(x,y),
        !  are simply given by Cauchy-Riemann relations
        !  Examples of generating arrays of z are included in the driver code.
        !----------
        !  Author: Mofreh R. Zaghloul
        !  United Arab Emirates University, June 05, 2017
        !----------

        COMPLEX(rk), INTENT(IN)  :: z
        COMPLEX(rk), INTENT(OUT) :: w

        real(rk), intent(out), optional :: dVdx, dVdy
        REAL(rk)::x,y,x_sqr,y_sqr,xsqr_plus_ysqr,ysqr_minus_xsqr_3halfs,four_x_sqr_y_sqr

        x=REAL(z, KIND=rk)
        y=AIMAG(z)
        x_sqr=x*x
        y_sqr=y*y
        xsqr_plus_ysqr=x_sqr+y_sqr
!
!-------------Region I ----> Laplace Cont. Fractions, 1 convergent
        IF (xsqr_plus_ysqr>=38000.0_rk) THEN
            !:-  w=((j1* one_sqrt_pi)/z)
            w=(y+j1*x)*(one_sqrt_pi/xsqr_plus_ysqr)
!
!-------------Region II ----> Laplace Cont. Fractions, 2 convergents
! Also identical to the approximation for regions I in Humlicek’s w4 algorithm
        ELSE IF (xsqr_plus_ysqr>=256.0_rk) THEN
            !:-  w=(j1*one_sqrt_pi*z)/(z*z-0.5_rk)
            w=((y*(half+xsqr_plus_ysqr))+j1*(x*(xsqr_plus_ysqr-half)))*&
                (one_sqrt_pi/((xsqr_plus_ysqr*xsqr_plus_ysqr+(y_sqr-x_sqr))+0.25_rk))
!
!-------------Region III ----> Laplace Cont. Fractions, 3 convergents
        ELSE IF (xsqr_plus_ysqr>=62.0_rk) THEN
            !:-  w=(j1*one_sqrt_pi/z)*(z**2-1.0_rk)/(z**2-1.5_rk)
            ysqr_minus_xsqr_3halfs=y_sqr-x_sqr+three_halfs
            four_x_sqr_y_sqr=four*x_sqr*y_sqr
            w=one_sqrt_pi*((y*((ysqr_minus_xsqr_3halfs-half)*(ysqr_minus_xsqr_3halfs)+four_x_sqr_y_sqr+x_sqr)) &
                +j1*(x*((ysqr_minus_xsqr_3halfs-half)*(ysqr_minus_xsqr_3halfs)+four_x_sqr_y_sqr-y_sqr)))/&
                (xsqr_plus_ysqr*((ysqr_minus_xsqr_3halfs)*(ysqr_minus_xsqr_3halfs)+four_x_sqr_y_sqr))
!
!-------------Region IV ----> Laplace Laplace Cont. Fractions, 4 convergents
! Also identical to the approximation for regions II in Humlicek’s w4 algorithm
        ELSE IF (xsqr_plus_ysqr>=30.0_rk .and. y_sqr>=1.0e-13_rk ) THEN
            !:- w=(j1*one_sqrt_pi*z)*(z**2-2.5_rk)/(z**2*(z**2-3.0_rk)+0.75_rk)
            w=(z*z)
            w=(one_sqrt_pi*(-y+j1*x))*(w - 2.5_rk)/(w*(w - 3.0_rk)+0.75_rk)

!
!-------------Region V----> Humlicek's w4 (Region IV)
        ELSE IF ((xsqr_plus_ysqr>2.50_rk .and. y_sqr<5.0e-9_rk)) THEN
            w=-z*z
            w=exp(-x_sqr)+(j1*z*(ha0-w*(ha1-w*(ha2-w*(ha3-w*(ha4-w*(ha5-w*ha6)))))) / &
                (hb0-w*(hb1-w*(hb2-w*(hb3-w*(hb4-w*(hb5-w*(hb6-w))))))))

         ELSE IF (xsqr_plus_ysqr>2.5_rk  .and. y_sqr<0.072_rk ) THEN
            w=-z*z
            w=EXP(w)+(j1*z*(ha0-w*(ha1-w*(ha2-w*(ha3-w*(ha4-w*(ha5-w*ha6)))))) / &
                (hb0-w*(hb1-w*(hb2-w*(hb3-w*(hb4-w*(hb5-w*(hb6-w))))))))
!
!-------------Region VI----> Hui's p-6 Approximation
        ELSE
            w=(y-j1*x)
            w=((((((aa6*w+aa5)*w+aa4)*w+aa3)*w+aa2)*w+aa1)*w+aa0)/ &
                (((((((w+bb6)*w+bb5)*w+bb4)*w+bb3)*w+bb2)*w+bb1)*w+bb0)
        END IF
!
!--------------- Calculation of the derivatives
        ! Partial derivative of V w.r.t. x       =dLdy
        if ( present(dVdx) ) dVdx=-two*real(z*w)
        ! Partial derivative of V w.r.t. y       =-dLdx
        if ( present(dVdy) ) dVdy=two*aimag(z*w)-two_sqrt_pi

    END SUBROUTINE wz_rk
END MODULE wz_mod_rk

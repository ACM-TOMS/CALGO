module Faddeyeva_v2_mod_rk
    use set_rk, only:r4,rk
    use rk_erfcx_Cody, only: erfc_scaled=>rk_erfcx

    implicit none
    private
    public :: Faddeyeva_v2_rk

    real(rk), parameter :: zero   = 0.0_rk
    real(rk), parameter :: Rmin = tiny(1.0_rk)                       !underflow limit
    real(rk), parameter :: Log_Rmin = (1-r4/rk)*708.396418532264e0_rk+(r4/rk)*87.33655_rk        !abs(log(Rmin))
    real(rk), parameter :: sqrt_log_Rmin =(1-r4/rk)*26.615717509251258e0_rk+(r4/rk)*9.3454026_rk !sqrt(Log_Rmin)
    real(rk), parameter :: huntw  = 1.0_rk / 120.0_rk
    real(rk), parameter :: sixth  = 1.0_rk / 6.0_rk
    real(rk), parameter :: quarter= 0.25_rk
    real(rk), parameter :: half   = 0.5_rk
    real(rk), parameter :: one    = 1.0_rk
    real(rk), parameter :: three_halfs  = one+half
    real(rk), parameter :: two    = 2.0_rk
    real(rk), parameter :: four   = 4.0_rk
    real(rk), parameter :: five   = 5.0_rk
    real(rk), parameter :: one_sqrt_pi = 5.641895835477563e-001_rk   !1/sqrt(pi)
    real(rk), parameter :: two_sqrt_pi = 1.128379167095513e+000_rk   !2/sqrt(pi)

    real(rk), parameter :: a      = 0.5_rk+(r4/rk)*.2_rk             !a=0.5 for r8 & 0.7 for r4
    real(rk), parameter :: a_pi   = a*0.318309886183791e0_rk         !a/pi
    real(rk), parameter :: a_sqr  = a**2
    real(rk), parameter :: exp_a_sqr  = (1-r4/rk)*7.788007830714049e-001_rk+(r4/rk)*6.126263941844161e-001_rk
    real(rk), parameter :: inv_a_sqr = 1.0_rk/a_sqr
    real(rk), parameter :: two_a_pi = 2.0_rk * a_pi
    real(rk), parameter :: two_a_pi_half_a = 0.5_rk * two_a_pi * a
    real(rk), parameter :: two_a_sqr= 2.0_rk * a_sqr
    real(rk), parameter :: exp_2a_sqr  =(1-r4/rk)*1.648721270700128e+000_rk+(r4/rk)*2.664456241929417e+000_rk

    integer, parameter :: ncycles_4 = 2*(rk/r4+1)+1
    integer, parameter :: ncycles_5 = 4*(rk/r4)+(r4/rk)
    integer, parameter :: ncycles_6 = ncycles_5+1
    integer, parameter :: ndigits_max =13-7*(r4/rk)
    integer, parameter :: bigz_border(10)=(/107, 110, 180, 380, &
        810, 1750, 3800, 8200, 17500, 38000/)
    integer, parameter :: nCycles(10)=(/ncycles_4, ncycles_5, &
        ncycles_6, 9, 10, 10, 11, 11, 12, 13/)

    complex(rk), parameter :: j1= (zero,one)

    ! exp(-two_a_sqr*(n-1))
    real(rk), parameter :: exp_2a_sqr_n_1(13) = &
        (/ (1.000000000000000E+000_rk*(1+(r4/rk)*0.0_rk)), &
        (6.065306597126334E-001_rk*(1-(r4/rk)*3.812166081938592E-001_rk)),  &
        (3.678794411714423E-001_rk*(1-(r4/rk)*6.171071140248879E-001_rk)),  &
        (2.231301601484298E-001_rk*(1-(r4/rk)*7.630722413178782E-001_rk)),  &
        (1.353352832366127E-001_rk*(1-(r4/rk)*8.533930378696498E-001_rk)),  &
        (8.208499862389880E-002_rk*(1-(r4/rk)*9.092820467105875E-001_rk)),  &
        (4.978706836786394E-002_rk*(1-(r4/rk)*9.438652371658662E-001_rk)),  &
        (3.019738342231850E-002_rk*(1-(r4/rk)*9.652647410552614E-001_rk)),  &
        (1.831563888873418E-002_rk*(1-(r4/rk)*9.785063986549101E-001_rk)),  &
        (1.110899653824231E-002_rk*(1-(r4/rk)*9.867001164575563E-001_rk)),  &
        (6.737946999085467E-003_rk*(1-(r4/rk)*9.917702529509800E-001_rk)),  &
        (4.086771438464067E-003_rk*(1-(r4/rk)*9.949075692073008E-001_rk)),  &
        (2.478752176666359E-003_rk*(1-(r4/rk)*9.968488884015556E-001_rk))    /)

!exp(-a_sqr * nn**2)/a**2
    real(rk), parameter :: inv_asqr_exp_asqr_nsqr(13) = &
       (/inv_a_sqr*7.788007830714049e-001_rk*(1-(r4/rk)*2.133721389334466e-001_rk), &
        inv_a_sqr*3.678794411714423e-001_rk*(1-(r4/rk)*6.171071140248879e-001_rk),  &
        inv_a_sqr*1.053992245618643e-001_rk*(1-(r4/rk)*8.846748789619375e-001_rk),  &
        inv_a_sqr*1.831563888873418e-002_rk*(1-(r4/rk)*9.785063986549101e-001_rk),  &
        inv_a_sqr*1.930454136227709e-003_rk*(1-(r4/rk)*9.975212478233336e-001_rk),  &
        inv_a_sqr*1.234098040866796e-004_rk*(1-(r4/rk)*9.998231130977574e-001_rk),  &
        inv_a_sqr*4.785117392129010e-006_rk*(1-(r4/rk)*9.999921891752664e-001_rk),  &
        inv_a_sqr*1.125351747192591e-007_rk*(1-(r4/rk)*9.999997865791929e-001_rk),  &
        inv_a_sqr*1.605228055185612e-009_rk*(1-(r4/rk)*9.999999963915950e-001_rk),  &
        inv_a_sqr*1.388794386496402e-011_rk*(1-(r4/rk)*9.867001164575563E-001_rk),  &
        inv_a_sqr*7.287724095819692e-014_rk*(1-(r4/rk)*9.999999999997556e-001_rk),  &
        inv_a_sqr*2.319522830243570e-016_rk*(1-(r4/rk)*9.999999999999990e-001_rk),  &
        inv_a_sqr*4.477732441718302e-019_rk*(1-(r4/rk)*1.000000000000000e+000_rk)  /)

contains
    !..
    ! If derivatives are not required, only positive values of the imaginary
    ! part of the input complex number z are considered and when the  optional
    ! "Stat" argument is not needed, efficiency can be improved by
    ! commenting the first line defining the subroutine and uncommenting
    ! the following line (line 89) together with commenting the lines
    ! (153,155,163,165,168, 182-185,249,252)
    !..

    elemental subroutine Faddeyeva_v2_rk ( z, w, sdgts, dVdx, dVdy, Stat )
        !         elemental subroutine Faddeyeva_v2_rk ( z, w, sdgts)

        ! ----------
        !  Faddeyeva_v2_rk is an elemental Fortran subroutine that receives, as input, a
        !  single dummy scalar complex number (z=x+iy) and returns as output
        !  the Faddeyeva function,w, defined as w(z)=exp(-z^2)*erfc(-i*z) where erfc(z)
        !  is the complex complementary error function. An optional integer input (sdgts),
        !  may be used representing the desired number of significant figures in
        !  the calculated Faddeyeva function.
        !  In addition, the subroutine returns as optional the derivatives
        !  of the real part of the Faddeyeva function (dVdx and dVdy) with respect to
        !  the real and imaginary parts of z, respectively.
        !------------
        !  This version uses a reformed form of Humlicek's w4 rational approximation
        !  [Zaghloul, M. R. 2015. A simple reform for treating the loss of accuracy
        !  of Humliceck’s W4 algorithm near the real axis. arXiv:1505.05596v1 [astro-ph.IM]]
        !  for the case of 4 significant figures
        !------------
        !  The routine can be run in default real (single precision) or in double
        !  precision depending on the choice of the integer parameter rk in the
        !  subsidiary module set_rk
        !-------------
        !  The subroutine is set to reproduce the highest possible accuracy obtainable
        !  from algorithm 916 [Zaghloul and Ali,TOMS, Vol. 38, No. 2, article 15:1-22
        !  (2011)] through requesting a number of significant figures (sdgts) equal
        !  to 13 when run in double precision.
        !  The desired number of significant figures can be reduced for marginal
        !  improvements of the efficiency at the expense of accuracy.
        !  The recommended range for "sdgts" is between 4 and 13 for double precision
        !  and between 4 and 6 for single precision.
        !  A number of significant figures smaller than 4 is not recommended for accuracy
        !  concerns, particularly regarding the computations of the derivatives, ant it
        !  will be directly changed to 4.
        !
        !  Values of sdgts >13 for double precision are not recommended for performance
        !  concerns and will be automatically changed to 13.
        !  Similarly, values of sdgts > 6 for single precision are not recommended
        !  for performance concerns and will be automatically changed to 6.
        !
        !  An optional "Stat" argument is used because "WRITE" statements are not permitted
        !  in elemental subroutines. The "Stat" returns integer values 0,1,2 and -1
        !  corresponding to normal, too few significant figures <4, too many significant
        !  figures (>13 for double precision and >6 for default real or single precision)
        !  and overflow, respectively.
        !----------
        !  The accompanying driver code “Faddeyeva_driver_rk.f90” can be run for
        !  computation of the Faddeyeva function w(z)=V+iL and the partial
        !  derivatives of its real part, V(x,y) (optional), on a scalar or an array of
        !  the complex variable z. The partial derivatives of the imaginary part, L(x,y),
        !  are simply given by Eq. (23) in the original article of Algorithm 916,
        !  TOMS. Vol. 38, No. 2, article 15:1-22 (2011) and as commented by the
        !  end of this subroutine.
        !  An example of generating an array of z is included in the driver code.
        !----------
        !  Author: Mofreh R. Zaghloul
        !  United Arab Emirates University, June 05, 2017
        !----------


        ! Private mathematical and repeatedly used constants
        real(rk) :: xx,yy
        complex(rk), intent(in) :: z
        integer, intent(in), optional :: sdgts
        complex(rk),intent(out) :: w
        real(rk), intent(out), optional :: dVdx, dVdy
        real(rk)::xsqr_plus_ysqr,ysqr_minus_xsqr,x_sqr,y_sqr,s
        integer, intent(out), optional :: Stat  !  0 => Normal,
                                                !  1 => too few (<4)
                                                !  2 => too many sdgts (>6) for r4 &
                                                !                      (>13) for r8,
                                                !  -1 => overflow
        integer :: ndgts

        if (present(sdgts))ndgts=sdgts
        if ( present(stat) ) stat = 0
        if (ndgts<4) then
            if ( present(stat) ) stat = 1
            ndgts=4
        elseif (ndgts>ndigits_max)then
            if ( present(stat) ) stat = 2
            ndgts=ndigits_max
        endif
        if (.not. present(sdgts) ) ndgts = 4



        xx=real(z)
        yy=aimag(z)
        x_sqr=xx*xx
        y_sqr=yy*yy
        xsqr_plus_ysqr=(x_sqr+y_sqr)
        ysqr_minus_xsqr=y_sqr-x_sqr

        if ( yy<zero .and. (ysqr_minus_xsqr)>=log_Rmin ) then
            if ( present(stat) ) stat = -1
            return
        endif


        if (ndgts>4)then

            !-------------
            ! Asymptotic expression for large |z| (|z|^2>bigz_border)
            ! exapded form of (j1*one_sqrt_pi)*two*(z*z - one)/(z*(two*z*z - 3.0_rk))
            if (xsqr_plus_ysqr>=bigz_border(ndgts-3)) then
                !        w=(j1*one_sqrt_pi)*two*(z*z - one)/(z*(two*z*z - 3.0_rk))
                w=one_sqrt_pi*((yy*((ysqr_minus_xsqr+one)*(ysqr_minus_xsqr+three_halfs)+four*x_sqr*y_sqr+x_sqr)) &
                    +j1*(xx*((ysqr_minus_xsqr+one)*(ysqr_minus_xsqr+three_halfs)+four*x_sqr*y_sqr-y_sqr)))/&
                    (xsqr_plus_ysqr*((ysqr_minus_xsqr+three_halfs)*(ysqr_minus_xsqr+three_halfs)+four*x_sqr*y_sqr))

                !--------------
                ! For x=0, use the asymptotic expressions for x--->0 from eq. (6) in the original
                ! article of the algorithm 916, TOMS. Vol. 38, No. 2, article 15:1-22 (2011).
            elseif (xx==zero .and. xsqr_plus_ysqr<bigz_border(ndgts-3)) then
                w = cmplx(erfc_scaled((yy)), 0.0_rk,kind=rk)


                ! -------Calculating Faddeyeva fn for values of 0<|x|<sqrt(-log(Rmin))
                ! -------while |z|^2<bigz_border
            elseif ( xx/=zero .and. abs(xx)<sqrt_log_Rmin .and. &
                    xsqr_plus_ysqr<bigz_border(ndgts-3) ) then
                w = Use_Cycles_1 ( xx, yy, nCycles(ndgts-3) )
            else

                !-------Calculating Faddeyeva for values of x>=sqrt(-log(Rmin))
                !-------while |z|^2<bigz_border
                w = Use_Cycles_2 ( xx, yy, nCycles(ndgts-3) )
            endif

        else

            !.. Reformed Humlicek Routine
            s=abs(xx)+yy
            if (s>=15.0_rk) then
                w=((yy*(half+xsqr_plus_ysqr))+j1*(xx*(xsqr_plus_ysqr-half)))*&
                    (one_sqrt_pi/((xsqr_plus_ysqr*xsqr_plus_ysqr+ysqr_minus_xsqr)+quarter))

            elseif (s<15.0_rk .and. s>=5.5_rk .and. y_sqr>1.000e-12_rk ) then
                w=z*z
                w=(one_sqrt_pi*(-yy+j1*xx))*(w - 2.5_rk)/ &
                    (w*(w - 3.0_rk)+0.75_rk)

            elseif (yy >= 0.195_rk*abs(xx) - 0.176_rk ) then
                w=-j1*z
                w = (16.4955_rk + w*(20.20933_rk + w*(11.96482_rk + w*(3.778987_rk + &
                    0.5642236_rk*w)))) / &
                    (16.4955_rk + w*(38.82363_rk + w*(39.27121_rk + w*(21.69274_rk + &
                    w*(6.699398_rk + w)))))
            else
                w=-z*z
                w=exp(w)-(-j1*z*(36183.31_rk-w*(3321.99_rk-w*(1540.787_rk-w*(219.031_rk-&
                    w*(35.7668_rk-w*(1.320522_rk-w*0.56419_rk)))))) / &
                    (32066.6_rk-w*(24322.84_rk-w*(9022.228_rk-w*(2186.181_rk-&
                    w*(364.2191_rk-w*(61.57037_rk-w*(1.841439_rk-w))))))))
            endif

        endif

        !--------------- Calculation of the derivatives
        ! Partial derivative of V w.r.t. x       =dLdy
        if ( present(dVdx) ) dVdx=-two*real(z*w)

        ! Partial derivative of V w.r.t. y       =-dLdx
        if ( present(dVdy) ) dVdy=two*aimag(z*w)-two_sqrt_pi

    end subroutine Faddeyeva_v2_rk


    !==========================================
    elemental function Use_Cycles_1 ( xx, yy, nCycles ) result ( w )

        real(rk), intent(in) :: xx, yy
        integer, intent(in) :: nCycles
        complex(rk) :: w
        real(rk) :: cos_2yx, del2_tmp, del3_tmp, del3_3_tmp, den1,     &
            erfcx_y, exp_x_sqr, exp1, exp2, exp3, exp3_den, exp3_3_den,&
            L_old, n3, n3_3, sigma1, sigma2_3, sigma4_5, two_a_pi_y,   &
            two_a_sqr_n3, two_x, two_yx,two_a_x, V_old, x, x_sqr, y,   &
            y_sqr_a_sqr

        integer :: n

        x=abs(xx)
        y=max(Rmin,abs(yy))
        erfcx_y=erfc_scaled(y)
        x_sqr=x*x
        two_x=two*x
        two_a_x=a*two_x
        two_yx=y*two_x
        cos_2yx=cos(two_yx)
        exp_x_sqr=exp(-x_sqr)
        two_a_pi_y=two_a_pi/y
        n3=real(ceiling(x/a),kind=rk)
        n3_3=n3-one
        two_a_sqr_n3=two_a_sqr*n3
        y_sqr_a_sqr=inv_a_sqr*y*y

        sigma1=zero
        sigma2_3=zero
        sigma4_5=zero
        V_old=exp_x_sqr*(erfcx_y*cos_2yx+two_a_pi_y*sin(two_yx*half)**2)
        L_old=exp_x_sqr*(-erfcx_y+half*two_a_pi_y);

        exp1=exp(-two_a_x)
        exp3=exp(-two_a_sqr_n3+two_a_x+two_a_sqr)
        exp2=exp_2a_sqr/(exp3*exp3)

        del2_tmp=one;
        del3_tmp=exp(-((a_sqr*n3*n3-two_a_x*n3-two_a_sqr_n3)+x_sqr+two_a_x+a_sqr));
        del3_3_tmp=exp_a_sqr*exp3;

        do n = 1, nCycles
            den1=inv_asqr_exp_asqr_nsqr(n)*exp_x_sqr/(n*n+y_sqr_a_sqr);
            del2_tmp=del2_tmp*exp1;
            del3_tmp=del3_tmp*exp3;
            exp3_den=del3_tmp*inv_asqr_exp_asqr_nsqr(n)/((n3_3+n)**2+y_sqr_a_sqr);
            sigma1=sigma1+den1;
            if (n3_3>=n)then
                del3_3_tmp=del3_3_tmp*exp2;
                exp3_3_den=del3_3_tmp*del3_tmp*inv_asqr_exp_asqr_nsqr(n)/((n3-n)**2+y_sqr_a_sqr);
                sigma2_3=sigma2_3+del2_tmp*den1+exp3_3_den+exp3_den;
                sigma4_5=sigma4_5+(n3-n)*exp3_3_den+(n3_3+n)*exp3_den-n*del2_tmp*den1;
            else
                sigma2_3=sigma2_3+del2_tmp*den1+exp3_den;
                if (x>=5.0e-3_rk)then
                    sigma4_5= sigma4_5+(n3_3+n)*exp3_den-n*del2_tmp*den1;
                else
                    sigma4_5=sigma4_5+two*n*n*two_a_x*den1*(one+sixth*n*n*two_a_x*two_a_x+ &
                        huntw*(n*n*two_a_x*two_a_x)**2);
                endif
            endif
        end do

        if ((y < five).and.two_yx>Rmin)then
            w=(V_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+ &
                j1*sign(one,xx)*(sin(two_yx)*(L_old+two_a_pi*y*sigma1)+&
                two_a_pi_half_a*sigma4_5));
        elseif ((y < five) .and. two_yx<=Rmin)then
            w=(V_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+ &
                j1*sign(one,xx)*(y*(two_x*L_old+two_x*two_a_pi*y*sigma1)+&
                two_a_pi_half_a*sigma4_5));
        else
            w=(V_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+ &
                j1*sign(one,xx)*two_a_pi_half_a*sigma4_5);
        endif
        IF (yy<zero) THEN
            w=two*EXP(y*y-x*x+j1*two_yx)-(REAL(w)-j1*AIMAG(w))
        ENDIF

    end function Use_Cycles_1

    !=========================================
    elemental function Use_Cycles_2 ( xx, yy, nCycles ) result ( w )

        real(rk), intent(in) :: xx, yy
        integer, intent(in) :: nCycles
        complex(rk) :: w
        real(rk) :: del3_3_tmp, den1, exp2, exp3_den, exp3_3_den, factor, &
            n3, n3_3, sigma3, sigma5, x, y, y_sqr_a_sqr
        integer :: n

        x=abs(xx)
        y=max(Rmin,abs(yy))
        n3=real(ceiling(x/a),kind=rk) ! ceiling of (x/a)
        n3_3=n3-1
        y_sqr_a_sqr=inv_a_sqr*y*y
        sigma3=zero
        sigma5=zero
        del3_3_tmp=exp((two*a*x-two_a_sqr*n3)+a_sqr)
        exp2=1/(del3_3_tmp*del3_3_tmp)
        factor=del3_3_tmp
        exp3_den=inv_a_sqr*exp(-(a*n3_3-x)**2)

        do n = 1, nCycles
            del3_3_tmp = del3_3_tmp*exp2
            exp3_den = exp3_den*factor*exp_2a_sqr_n_1(n)
            exp3_3_den = exp3_den*del3_3_tmp/((n3-n)**2+y_sqr_a_sqr)
            den1 = exp3_den/((n3_3+n)**2+y_sqr_a_sqr);
            sigma3 = sigma3+exp3_3_den+den1
            sigma5 = sigma5+(n3-n)*exp3_3_den+(n3_3+n)*den1
        end do

        w=y*a_pi*sigma3+j1*sign(one,xx)*two_a_pi_half_a*sigma5

        IF (yy<zero) THEN
            w=two*EXP(y*y-x*x+j1*two*y*x)-(REAL(w)-j1*AIMAG(w))
        ENDIF

    end function Use_Cycles_2
end module Faddeyeva_v2_mod_rk

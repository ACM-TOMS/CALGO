  MODULE Parabolw
  IMPLICIT NONE
  INTEGER, PARAMETER  :: r8 = KIND(0.0d0)
  PRIVATE
  PUBLIC  :: parabw
  CONTAINS           
    SUBROUTINE parabw(a,x,mode,waxx,wamx,ierr)
    ! ---------------------------------------------------------
    ! Calculation of the real parabolic cylinder functions
    ! W(a,x), W(a,-x) and their derivatives, x,a real and x>=0.   
    ! ----------------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   waxx , waxx(1), W(a,x)
    !          waxx(2), W'(a,x)
    !   wamx,  wamx(1), W(a,-x)
    !          wamx(2), W'(a,-x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems
    !          ierr=2, the argument x is out of range 
    ! -----------------------------------------------------------
    !           METHODS OF COMPUTATION
    ! -----------------------------------------------------------
    ! The present code uses different methods of computation
    ! depending on the values of a and x:
    ! For a>0:
    !    a) McLaurin series; b)Numerical solution of ODEs using Taylor 
    !       method; c) Asymptotic expansions in terms of elementary 
    !       functions; d) Airy-type asymptotic expansions.
    ! For a<0:
    !    a) McLaurin series; b)Numerical solution of ODEs using Taylor 
    !       method; c) Asymptotic expansions
    ! ----------------------------------------------------------------
    !                  SCALING
    ! ----------------------------------------------------------------
    ! For a>50, there is the possibility of computing scaled functions.
    ! The scaling factors are:
    ! 
    ! exp(a*arcsin(x/(2*sqrt(a)))+x/2*sqrt(a-x**2/4)), if x**2 <= 2a
    ! exp(a*pi/2),  if  x > sqrt(2a)
    ! exp(-a*pi/2), if  x < -sqrt(2a)
    !
    ! If mode=1 and a<=50, unscaled functions are computed
    ! ----------------------------------------------------------------  
    !                  ACCURACY  
    !-----------------------------------------------------------------
    !  The aimed relative accuracy for scaled functions is better than 
    !  5.0e-14 and better than 5.0e-13 in the computable range of 
    !  the unscaled PCFs W(a,x). 
    ! ----------------------------------------------------------------
    ! Authors:
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! -------------------------------------------------------------
    !  References:
    !  1) Fast and Accurate Computation of the Weber Parabolic Cylinder 
    !      Function W(a,x)
    !     A. Gil, J. Segura, N.M. Temme         
    !     IMA Journal of Numerical Analysis (2010) 
    ! -------------------------------------------------------------
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(OUT) :: waxx(2), wamx(2)   
    REAL(r8) :: xx, wax, waxd, waxm, wamxd, f1,f2,f3,f4,f5,f6
    INTEGER,  INTENT(INOUT) :: ierr
    INTEGER :: modea
    ierr=0
    IF (x<0) THEN
      waxx(1)=0.0_r8
      waxx(2)=0.0_r8
      wamx(1)=0.0_r8
      wamx(2)=0.0_r8
      ierr=2
    ENDIF
    IF (ierr==0) THEN
      xx=x
      modea=mode
      IF (a<50) modea=0 
      IF (a>0) THEN
        IF (xx<2) THEN
          f1=7.0_r8*(xx-2.0_r8)*(xx-2.0_r8)+2.5_r8;
          f2=1.1_r8*xx*xx+30.5_r8;
          IF (a<f1) THEN
          ! McLaurin Series 
            call waxsma(a, xx, wax, waxd, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a<f2) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))<1;
            call expaelem(a,xx,modea,wax,waxd,ierr)
            call expaelem(a,-xx,modea,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (xx<3.8_r8) THEN
          f2=1.1_r8*xx*xx+30.5_r8;     
          IF (a<f2) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))<1;
            call expaelem(a,xx,modea,wax,waxd,ierr)
            call expaelem(a,-xx,modea,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (xx<12.0_r8) THEN 
          f2=1.1_r8*xx*xx+30.5_r8;       
          f3=0.2_r8*(xx-12.0_r8)*(xx-12.0_r8)+33.0_r8;
          IF (a<f3) THEN
          ! Numerical solution of ODEs
            call taylor(a, xx, wax, waxd)
            call taylor(a, -xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a<f2) THEN
          ! Airy-type expansions
            call expair(a,xx,modea,wax,waxd,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))<1;
            call expaelem(a,xx,modea,wax,waxd,ierr)
            call expaelem(a,-xx,modea,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (xx<22.49_r8) THEN
          f2=1.1_r8*xx*xx+30.5_r8; 
          f3=0.2_r8*(xx-12.0_r8)*(xx-12.0_r8)+33.0_r8;
          f4=0.22_r8*xx*xx-40.0_r8;
          IF (a<f4) THEN
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))>1;
            call wammx(a, xx, modea,wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a<f3) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a<f2) THEN
          ! Airy-type expansions
            call expair(a,xx,modea,wax,waxd,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))<1;
            call expaelem(a,xx,mode,wax,waxd,ierr)
            call expaelem(a,-xx,mode,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSE
          f2=1.1_r8*xx*xx+30.5_r8;  
          f4=0.22_r8*xx*xx-40.0_r8;
          IF (a<f4) THEN
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))>1;
            call wammx(a, xx, modea,wax, waxd, waxm, wamxd, ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a<f2) THEN      
          ! Airy-type expansions
            call expair(a, xx, modea,wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a positive
          ! x+ a large, t=x/(2a**(1/2))<1;
            call expaelem(a,xx,modea,wax,waxd,ierr)
            call expaelem(a,-xx,modea,waxm,wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ENDIF
      ELSE
        IF (xx<1.3_r8) THEN
          f6=0.13_r8*xx*xx-25.0_r8;
          IF (a>f6) THEN
          ! McLaurin Series 
            call waxsma(a, xx, wax, waxd, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE  
          ! Asymptotic expansions, a negative
            call waxane(abs(a), xx, wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (xx<3.8_r8) THEN
          f6=0.13_r8*xx*xx-25.0_r8;
          f5=-2.0_r8*(xx-4.6_r8)*(xx-4.6_r8)-3.0_r8;
          IF (a>f5) THEN
          ! Series 
            call waxsma(a, xx, wax, waxd, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a>f6) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a negative
            call waxane(abs(a), xx, wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (x<4.5_r8) THEN
          f6=0.13_r8*xx*xx-25.0_r8;
          IF (a>-4.0_r8) THEN
          ! McLaurin series 
            call waxsma(a, xx, wax, waxd, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSEIF (a>f6) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a negative
            call waxane(abs(a), xx, wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSEIF (x<13.87_r8) THEN
          f6=0.13_r8*xx*xx-25.0_r8;
          IF (a>f6) THEN
          ! Numerical solution of ODEs
            call taylor(a,xx, wax, waxd)
            call taylor(a,-xx, waxm, wamxd)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ELSE
          ! Asymptotic expansions, a negative
            call waxane(abs(a), xx, wax, waxd, waxm, wamxd,ierr)
            waxx(1)=wax
            waxx(2)=waxd
            wamx(1)=waxm
            wamx(2)=wamxd
          ENDIF
        ELSE
        ! Asymptotic expansions, a negative
          call waxane(abs(a), xx, wax, waxd, waxm, wamxd,ierr)
          waxx(1)=wax
          waxx(2)=waxd
          wamx(1)=waxm
          wamx(2)=wamxd
        ENDIF
      ENDIF
    ENDIF
    END SUBROUTINE parabw



 
    FUNCTION xpowy(x,y)
    IMPLICIT NONE
    REAL(r8) :: x,y,xpowy
    xpowy=x**y
    END FUNCTION xpowy

    FUNCTION fractio2(x,n,pk,qk)
    IMPLICIT NONE	
    INTEGER m,n,k
    REAL(r8) x,fractio2,pk(0:8),qk(0:8),p,q
    p= pk(n); q= qk(n); m= n-1
    DO k= m,0,-1 
      p= p*x + pk(k); q= q * x + qk(k)
    ENDDO
    fractio2= p/q
    END FUNCTION fractio2


    FUNCTION polysum(pk, x, n)
    IMPLICIT NONE
    INTEGER n, k
    REAL(r8) x,pk(0:20),polysum,p
    ! {evaluates a polynomial p=p0 + p1x+...pnx^n}  
    p= pk(n); 
    DO k= n-1,0,-1 
      p= p*x + pk(k);
    ENDDO
    polysum= p
    END FUNCTION polysum


    RECURSIVE FUNCTION phase(x,y) RESULT(res)
    USE someconstants
    IMPLICIT NONE
    REAL(r8)x,y,res
    !{computes the phase of the complex number z=x+iy}
    IF ((x==0.0_r8).AND.(y==0.0_r8)) THEN
       res= 0.0_r8
    ELSE
      IF (y < 0) THEN
        res= -phase(x,-y)
      ELSEIF (x >= y) THEN
        res= atan(y/x)
      ELSEIF ((x + y).ge.0) THEN
        res= pi*0.5_r8 - atan(x/y)
      ELSE
        res= pi + atan(y/x)
      ENDIF
    ENDIF
    END FUNCTION phase
    
    FUNCTION rhoa(a)
    USE someconstants
    IMPLICIT NONE
    INTEGER s
    REAL(r8) a, rhoa, xi, beta, phi, u, v, delta
    beta= sqrt(0.25_r8+a*a); xi= atan(2.0_r8*a); 
    CALL gammaplushalf(a,u,v) 
    phi= phase(u,v);
    delta= (xi*(beta*cos(xi)-0.5_r8)+sin(xi)*(beta*log(beta)-beta-1.0_r8/&
           (12.0_r8*beta))-phi)/(2.0_r8*pi);
    s=nint(delta);
    rhoa= phi/2.0_r8+pi*s
    END FUNCTION rhoa     


  
    FUNCTION rhostar(a)
    IMPLICIT NONE
    REAL(r8) a, rhostar, a2, rk(0:20)
    !-----------------------------------------
    !Computes rhostar of (4.11) in Gil (2004);
    !rhoa is as in (4.9), without term pi/8
    !-----------------------------------------
    IF (a==0) THEN
      rhostar= 0.0_r8 
    ELSE
      IF (abs(a)< 1.0e-12_r8) THEN
        rhostar= a*(-0.48175501301071174_r8 -0.5_r8*log(abs(a))) 
      ELSE
        IF (abs(a)> 8) THEN 
          rk(0)= 0.20833333333333333333e-1_r8;
          rk(1)= 0.12152777777777777778e-2_r8;
          rk(2)= 0.38442460317460317460e-3_r8;
          rk(3)= 0.29529389880952380952e-3_r8;
          rk(4)= 0.42005339856902356902e-3_r8;
          rk(5)= 0.95829531254335941836e-3_r8;
          rk(6)= 0.32047369541266025641e-2_r8;
          rk(7)= 0.14774875890195759293e-1_r8;
          rk(8)= 0.89821500895519226285e-1_r8;
          a2= 1.0_r8/(a*a);
          rhostar= polysum(rk, a2,8)/a;
        ELSE
          rhostar= rhoa(a)+a*(1.0_r8-log(abs(a)))*0.5_r8
        ENDIF
      ENDIF
    ENDIF
    END FUNCTION rhostar


    SUBROUTINE gammaplushalf(a, reg, img)
    !-------------------------------------------------------
    !Computes reg and img in Gamma(1/2+i*a)=reg+i img; any a
    !------------------------------------------------------- 
    USE someconstants
    IMPLICIT NONE
    INTEGER n, k
    REAL(r8) a, reg, img
    REAL(r8) a2, g, theta, r, lnr, c, s, sigma, res, ims, rer, imr, u, v, gk(0:10)
    a2= a*a;
    IF (a > 20) THEN
      n= 0 
    ELSE 
      n= 1+nint(sqrt(400.0_r8-a2));
    ENDIF
    gk(0)=  1.00000000000000000000e0_r8;
    gk(1)= -0.41666666666666666667e-1_r8;
    gk(2)=  0.86805555555555555556e-3_r8;
    gk(3)=  0.24184992283950617284e-2_r8;
    gk(4)= -0.10114756140689300412e-3_r8;
    gk(5)= -0.76674039565229705565e-3_r8;
    gk(6)=  0.34959887446996281973e-4_r8;
    gk(7)=  0.58979762398995290743e-3_r8;
    gk(8)= -0.26464724593955018560e-4_r8;
    gk(9)= -0.83951400936616112498e-3_r8;
    gk(10)=  0.36726629864520529433e-4_r8;
    r= sqrt(n*n+a2); lnr= log(r); theta= atan(a/n); g= 1; 
    res= 1.0_r8; ims= 0.0_r8;
    DO k= 1,10  
      g= g/r; 
      res= res+g*gk(k)*cos(k*theta); 
      ims= ims-g*gk(k)*sin(k*theta);
    ENDDO
    sigma= a*lnr+n*theta-a;
    c= cos(sigma); s= sin(sigma);
    g= sqrt(2.0_r8*pi)*exp(n*lnr-theta*a-n);
    img= g*(c*ims+s*res);
    reg= g*(c*res-s*ims);
    rer= 1.0_r8; imr= 0.0_r8;
    DO k= 0,n-1  
      u= rer; v= imr; r= (k+0.5_r8)*(k+0.5_r8)+a2;
      rer= (u*(k+0.5_r8)+v*a)/r; imr= (-u*a+v*(k+0.5_r8))/r
    ENDDO
    u= reg; v= img;
    reg= rer*u-imr*v; img= rer*v+imr*u
    END SUBROUTINE gammaplushalf
  


    SUBROUTINE waxane(a, x, wax, waxp, waxm, waxmp,ierr)
    ! --------------------------------------------------
    ! Calculation of W(a,x), W'(a,x) by using asymptotic
    ! expansions for a negative
    ! --------------------------------------------------
    ! Inputs:
    !   abs(a) , being a, the order of the functions
    !   x ,    argument of the functions
    ! Outputs: 
    !   wax=W(-abs(a),x)
    !   waxp=W'(-abs(a),x) 
    !   waxm= W(-abs(a),-x)
    !   waxmp= W'(-abs(a),-x)
    ! ----------------------------------------------------
    USE someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    REAL(r8), INTENT(OUT) :: waxm
    REAL(r8), INTENT(OUT) :: waxmp
    INTEGER, INTENT(INOUT) :: ierr
    REAL(r8) ka, x2, x24, a2, a2xit, taut, w, &
             sr, si, tr, ti, sqa, p, q, r, s, c, chit
    REAL(r8) phik(0:8), psik(0:8)
    REAL(r8), DIMENSION(0:20) :: pjk1, qjk1, pjk2, qjk2, pjk3, qjk3
    REAL(r8), DIMENSION(0:20) :: pjk4, qjk4, pjk5, qjk5, pjk6, qjk6
    REAL(r8), DIMENSION(0:20) :: pjk7, qjk7, pjk8, qjk8
    INTEGER k
    DATA pjk1/-0.75_r8,-2.5_r8,-1.6666666666666666667_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk1/1.25_r8,3.5_r8,2.3333333333333333333_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk2/3.28125_r8,27.875_r8,67.375_r8,64.166666666666666667_r8,&
              21.388888888888888889_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8/
    DATA qjk2/-4.21875_r8,-33.625_r8,-79.958333333333333333_r8,&
              -75.833333333333333333_r8,-25.277777777777777778_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk3/-27.0703125_r8,-396.328125_r8,-1814.5125_r8,&
              -3861.2291666666666667_r8,-4254.25_r8,&
              -2363.4722222222222222_r8,-525.21604938271604938_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk3/31.9921875_r8,453.109375_r8,2047.1125_r8,&
              4330.5208333333333333_r8,4759.0277777777777778_r8,&
              2641.5277777777777778_r8,587.00617283950617284_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk4/329.91943359375_r8,7077.45703125_r8,49877.23984375_r8,&
              173584.125_r8,344694.42447916666667_r8,&
              411244.16666666666667_r8,292637.25231481481481_r8,&
              114759.70679012345679_r8,19126.617798353909465_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk4/-373.90869140625_r8,-7859.49609375_r8,-54833.22265625_r8,&
              -189793.22083333333333_r8,-375736.00052083333333_r8,&
              -447530.41666666666667_r8,-318189.01311728395062_r8,&
              -124738.81172839506173_r8,-20789.801954732510288_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk5/-5328.1988525390625_r8,-153266.177490234375_r8,&
              -1496065.7216099330357_r8,-7447882.350390625_r8,&
              -22100856.077647569444_r8,-42015571.25390625_r8,&
              -52708852.666377314815_r8,-43565653.690200617284_r8,&
              -22880216.541280864198_r8,-6933398.9519032921811_r8,&
              -924453.19358710562414_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
               0.0_r8,0.0_r8,0.0_r8/
    DATA qjk5/5889.0618896484375_r8,167081.378173828125_r8,&
             1618334.6748744419643_r8,8018883.473046875_r8,&
             23725519.159852430556_r8,45020597.459635416667_r8,&
             56413247.888020833333_r8,46595060.471836419753_r8,&
             24461987.833204732510_r8,7411564.3968621399177_r8,&
             988208.58624828532236_r8,&
             0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
             0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk6/107230.00190734863281_r8,3911743.4463500976562_r8,&
             49592558.348275320871_r8,327327256.24445452009_r8,&
             1320314667.6672144717_r8,3526512665.7394314236_r8,&
             6506075517.7748770255_r8,8449099135.7631655093_r8,&
             7726873009.8828125000_r8,4881459532.0875128601_r8,&
             2031139222.9600694444_r8,501515857.52100480110_r8,&
             55723984.169000533455_r8,&
             0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk6/-116554.34989929199219_r8,-4209790.5294799804688_r8,&
             -53046844.358616420201_r8,-348738575.26079101562_r8,&
             -1402928860.3276453993_r8,-3740348588.5545789931_r8,&
             -6891971112.2487340856_r8,-8942654096.4852912809_r8,&
             -8173643819.4593942901_r8,-5161879657.2836291152_r8,&
             -2147389212.0536479767_r8,-530173906.52220507545_r8,&
             -58908211.835800563938_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk7/-2585008.9745521545410_r8,-115109601.28864288330_r8,&
             -1814958496.4086914062_r8,-15118535470.703432792_r8,&
             -78228459666.499029134_r8,-273482436440.61774631_r8,&
             -677676841014.59605577_r8,-1223138659577.8502062_r8,&
             -1629074225394.3409650_r8,-1602934313802.8913484_r8,&
             -1152270884857.3288805_r8,-588917533590.47791281_r8,&
             -202879881562.49714220_r8,-42266641992.186904626_r8,&
             -4025394475.4463718691_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk7/2776491.1208152770996_r8,122708978.79421234131_r8,&
             1925312041.3459734235_r8,15983788916.836100551_r8,&
             82508725486.563457380_r8,287952809783.26722765_r8,&
             712650296975.96442744_r8,1285107966580.4986798_r8,&
             1710504064563.1047815_r8,1682287248047.2539706_r8,&
             1208933528167.4145608_r8,617751113141.81053884_r8,&
             212787605947.74543705_r8,44328429406.439924364_r8,&
             4221755181.5657070823_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk8/72622595.878824591637_r8,3839296326.6825199127_r8,&
             72972869796.835744858_r8,741010050728.21005118_r8,&
             4728818470068.3977036_r8,20667683341938.663361_r8,&
             65097477007191.150894_r8,152479523376491.21695_r8,&
             270645546827389.17508_r8,367482071900019.28684_r8,&
             382220257297374.97871_r8,302446999586890.60270_r8,&
             179020321106247.51147_r8,76811168492401.261776_r8,&
             22575519201784.001194_r8,4067661117438.5587738_r8,&
             338971759786.54656448_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk8/-77307924.645200371742_r8,-4062968417.3650646210_r8,&
             -76913846048.558694839_r8,-778781085218.40286473_r8,&
             -4959356633780.3453088_r8,-21641132068083.787434_r8,&
             -68082808844916.255468_r8,-159329371436569.03712_r8,&
             -282612499203202.34409_r8,-383535732852395.08975_r8,&
             -398767118978192.21973_r8,-315453595729638.25935_r8,&
             -186682859410003.59807_r8,-80088585766097.861952_r8,&
             -23536984672244.367115_r8,-4240753079882.7527641_r8,&
             -353396089990.22939701_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    IF (2.0_r8*log(x*0.5_r8)>log(giant)) THEN
      ierr=1
      wax=dwarf
      waxp=dwarf
      waxm=giant
      waxmp=giant
    ENDIF
    IF (ierr==0) THEN
      x2= x/2.0_r8; x24= x2*x2; a2= 2.0_r8*a; sqa= sqrt(a);
      p= sqrt(x24+a);
      w= 4.0_r8*p*(x2+p); taut= -a2/w; a2xit= x2*p;
      IF (a>0) THEN
        a2xit= a2xit+a*log((x2+p)/sqa);
      ENDIF
      phik(0)= 1; psik(0)= 1; r= -1.0_r8/w;
      phik(1)= r*polysum(pjk1, taut, 2); psik(1)= r*polysum(qjk1, taut,2); r= -r/w;
      phik(2)= r*polysum(pjk2, taut, 4); psik(2)= r*polysum(qjk2, taut,4); r= -r/w;
      phik(3)= r*polysum(pjk3, taut, 6); psik(3)= r*polysum(qjk3, taut,6); r= -r/w;
      phik(4)= r*polysum(pjk4, taut, 8); psik(4)= r*polysum(qjk4, taut,8); r= -r/w;
      phik(5)= r*polysum(pjk5, taut,10); psik(5)= r*polysum(qjk5,taut,10); r= -r/w;
      phik(6)= r*polysum(pjk6, taut,12); psik(6)= r*polysum(qjk6,taut,12); r= -r/w;
      phik(7)= r*polysum(pjk7, taut,14); psik(7)= r*polysum(qjk7,taut,14); r= -r/w;
      phik(8)= r*polysum(pjk8, taut,16); psik(8)= r*polysum(qjk8, taut,16);
      r= 1.0_r8; sr= 1.0_r8; tr= 1.0_r8; si= 0.0_r8; ti= 0.0_r8;
      DO k= 1,4 
        r= -r;
        sr= sr+r*phik(2*k); si= si+r*phik(2*k-1);  
        tr= tr+r*psik(2*k); ti= ti+r*psik(2*k-1)
      ENDDO
      chit= pi*0.25_r8+rhostar(-a)+a2xit; s= sin(chit); c= cos(chit);
      q= exp(-pi*a); ka= 1.0_r8/(q+sqrt(1.0_r8+q*q)); p= sqrt(p);
      ka=sqrt(ka);
      wax= ka*(c*sr-s*si)/p;
      waxp=-ka*p*(s*tr+c*ti);
      waxm= (s*sr+c*si)/(p*ka);
      waxmp= -p*(c*tr-s*ti)/ka;
    ENDIF
    END SUBROUTINE waxane

    
     SUBROUTINE expaelem(a,x,mode,wax,waxp,ierr)
    ! --------------------------------------------------
    ! Calculation of W(a,x), W'(a,x) by using asymptotic
    ! expansions in terms of elementary functions for 
    ! a positive and t=x/(2a**(1/2))<1
    ! --------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode,  mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   wax,  W(a,x)
    !   waxp, W'(a,x)
    !   ierr, error flag
    ! ---------------------------------------------------------
    USE someconstants
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN)  :: mode
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    INTEGER, INTENT(INOUT) :: ierr
    REAL(r8) mu, mu2, t, t2p1, facu, t2, t3, t4, t5, t6, t7,&
             t8, t9, t10, t11, t12, t13, t14, t15, t16, t17,&
             t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,&
             t28, t29, t30, t31, t32, t33, t34, t35, t36, t37,&
             t38, t39, t40, t41, xargu, delta1,&
             delta2, mu4, tts, ttss, mus, lmus, s1, s2, mus1
    REAL(r8), DIMENSION(0:15) :: lk,u,v
    INTEGER k
    mu=sqrt(2.0_r8*abs(a))
    t=x/(mu*sqrt2)
    t2=t*t
    t2p1=1.0_r8-t2
    ierr=0
    IF (mode==0) THEN
      ! ------------------------------------------------
      ! Check for possible overflow/underflow problems
      ! ------------------------------------------------        	
      IF (a*abs(asin(t)+t*sqrt(t2p1))>log(giant)) THEN
        ierr=1
        IF (x>0) THEN
          wax=dwarf
          waxp=dwarf
        ELSE
          wax=giant
          waxp=giant
        ENDIF
      ENDIF
    ENDIF
    lk(0)=1.0_r8
    lk(1)=-.86805555555555555555555555555555555555555555555556e-3_r8;
    lk(2)=-.40496700585133744855967078189300411522633744855967e-3_r8;
    lk(3)=-.55970973301118646421903237746859146036100768611057e-3_r8;
    lk(4)=-16.943102322943963245388661733269423877173517985977_r8;
    lk(5)=-.94037146647019565951398813260001779563624393715061e-2_r8;
    lk(6)=-.84327434386741291387140061090400783753265461929755e-1_r8;
    lk(7)=-1.1162982841593720168762155424414027479385891257438_r8;
    IF (ierr==0) THEN      
        mu2=mu*mu
        t3=t2*t; t4=t3*t; t5=t4*t; t6=t5*t; t7=t6*t; t8=t7*t; t9=t8*t;
        t10=t9*t; t11=t10*t; t12=t11*t; t13=t12*t; t14=t13*t; t15=t14*t;
        t16=t15*t; t17=t16*t; t18=t17*t; t19=t18*t; t20=t19*t; t21=t20*t;
        t22=t21*t; t23=t22*t; t24=t23*t; t25=t24*t; t26=t25*t; t27=t26*t;
        t28=t27*t; t29=t28*t; t30=t29*t; t31=t30*t; t32=t31*t; t33=t32*t;
        t34=t33*t; t35=t34*t; t36=t35*t; t37=t36*t; t38=t37*t; t39=t38*t;
        t40=t39*t; t41=t40*t  
        u(0)=1.0_r8
        u(1)=.4166666666666666666666667e-1_r8*t*(t2-6.0_r8)
        u(2)=.1258680555555555555555556_r8+.2161458333333333333333333_r8*t2-&
             .7812500000000000000000000e-2_r8*t4
        u(3)=-.6252170138888888888888889_r8*t-.3665002893518518518518519_r8*t3&
             -.6820746527777777777777778e-1_r8*t5+.4385850694444444444444444e-1_r8*t7&
             -.9746334876543209876543210e-2_r8*t9
        u(4)=-.3892736866640946502057613e-2_r8*t6+.3208267274707433127572016_r8&
             -.8071183569637345679012346e-2_r8*t8+.1827437789351851851851852e-2_r8*t10&
             +3.079461293161651234567901_r8*t2+1.279432885440779320987654_r8*t4
        u(5)=.1495080253037428963354889e-9_r8*t*(-34009066266.0_r8-119582875013.0_r8*t2&
             +1994971575.0_r8*t10-3630137104.0_r8*t8+4433574213.0_r8*t6&
             -37370295816.0_r8*t4+82393456.0_r8*t14-617950920.0_r8*t12)
        u(6)=25.84962798338368105790252_r8*t6+.1713039028908536874755046e-1_r8*t14&
            +2.581943865065753536501370_r8-.2309715544595780055849500e-2_r8*t16&
            -.5076538089737495643136799e-1_r8*t12-.3186668406760671571869489_r8*t8&
            +.8601473529456107027269498e-1_r8*t10+62.68920892021549979121503_r8*t2&
            +121.7179460820865799300044*t4;
        u(7)=-.7865531634245733182633045e-15_r8*t*(110065274285017326.0_r8&
            -125776738623632286.0_r8*t10+839442893976461885.0_r8*t2-51589093643118259.0_r8*t14&
            +96739097673864045.0_r8*t12+1014872932457342145.0_r8*t4+18208957214179542.0_r8*t16&
            +80153092162171710.0_r8*t6+114221755795701455.0_r8*t8-3833677186738392.0_r8*t18&
            +365112113022704.0_r8*t20)
        u(8)=5596.879044604074283137768_r8*t6+13.69521272837774834409345_r8*t14&
            +2.638371331228334872355069_r8*t18-7.390056136076261345369326_r8*t16&
            -17.49598518388745077640624_r8*t12+679.5130553150720497530461_r8*t8&
            +.5384626640674465938317402e-1_r8*t22+11.60187790297422517485544_r8*t10&
            -.5608986084035902019080627_r8*t20+1769.339898073369191761939_r8*t2&
            +7072.617883848053965523396_r8*t4+46.00893448266345421442011_r8
        u(9)=.1716369641152583388897846e-22_r8*t*(-.1359234799748243655770368e27_r8&
            +196970983448449506582784.0_r8*t26+.3935158819364260070148776e27_r8*t10&
            -.4824119188126012644493997e27_r8*t12+.4501074357324417329611653e27_r8*t14&
            -.5072207498713186429103949e27_r8*t8-.4586684500765682622356092e28_r8*t4&
            -.2551586343346702497275826e28_r8*t6+.1671990097243978109434037e27_r8*t18&
            -.3176839480721397029728588e27_r8*t16+.1661834830959361591291224e26_r8*t22&
            -.6369835812416149521724992e26_r8*t20-.1825949193391480686660700e28_r8*t2&
            -2659108276554068338867584.0_r8*t24)
        u(10)=914112.2670988274244565704_r8*t6+1457.987476485281604308530_r8*t14&
             +981.4809913805934941962403_r8*t18-52.76053458609560685621252_r8*t24&
             -.6338906553354126033547987_r8*t28+8.504699625750119095010216_r8*t26&
             +1213.983451915325067883233_r8-1375.290270384664889535434_r8*t16&
             -1323.383956656255885190354_r8*t12+415376.7168317029950197364_r8*t8&
             +200.5566577594646000950658_r8*t22+35768.04278735927027177930_r8*t10&
             -521.6682753837636432341761_r8*t20+75167.67289016707645386544_r8*t2&
             +532260.3452982493777355059_r8*t4
        u(11)=-.3359006667469412593393634e-31_r8*t*(.2788947144378624557826141e37_r8+&
             .1212394967494477547839171e34_r8*t32-.2000451696365887953934633e35_r8*t30&
             +.1550285194341124782611328e36_r8*t28-.7492595928004916139721089e36_r8*t26&
             +.2528614856038361343385493e37_r8*t24-.6321380226416724568027542e37_r8*t22&
             +.1211628603883129655678535e38_r8*t20-.1817616151670370488438873e38_r8*t18&
             -.2039447456964012213382125e38_r8*t14+.1526729786684424672364467e38_r8*t12&
             -.4167726542712734865896683e35_r8*t10+.1286775691925214910694996e39_r8*t8&
             +.3375397423010140569245605e39_r8*t6+.2668957183771871703360651e39_r8*t4&
             +.6033537049761229856043034e38_r8*t2+.2158825430879550563836032e38_r8*t16)
        u(12)=153053249.2910557053308385_r8*t6+78424.56411652519477233020_r8*t14+&
             128182.4650008160699499318_r8*t18-38558.08568947489765851297_r8*t24&
             -4637.584888683471835001935_r8*t28+15541.01589815635919400013_r8*t26&
             -119592.7457128667753060960_r8*t16+2857220.354347191785055925_r8*t12&
             +151197224.7859743079793903_r8*t8+965.8070562838161807622995_r8*t30&
             +73303.13746459791327933258_r8*t22+46150931.05589019507584725_r8*t10&
             +7.635830211413084690217302_r8*t34-125.3548793040314736644007_r8*t32&
             -108987.9643972328459888073_r8*t20+4457999.843158161718417669_r8*t2&
             +50231562.58692156905275991_r8*t4+48082.01508844808488819284_r8
        u(13)=.7703374808710407147877427e-40_r8*t*(-.6926180733786124297783135e47_r8&
             -.3490732449110820673555558e50_r8*t6+.8082605068869249980049792e48_r8*t22&
             +.2682009600079171280174385e48_r8*t26-.5172558239822205556645540e48_r8*t24&
             -.1109805333478410099533049e48_r8*t28+.1487945383913154715160031e46_r8*t34&
             +.3580113360465377112428930e47_r8*t30-.8679377325149432331738492e46_r8*t32&
             +.8249366361784484160838564e43_r8*t38-.1608626440547974411363520e45_r8*t36&
             -.9370275546601561256767938e48_r8*t16+.6651854621331235420882608e48_r8*t14&
             -.7836351143702021345383874e48_r8*t12-.7007020293314608957292234e49_r8*t10&
             -.2782407158660428165686684e50_r8*t8-.1553433992040236087729803e50_r8*t4&
             -.2221765150562660149162442e49_r8*t2-.1032888153750967506466895e49_r8*t20&
             +.1084709216891677773429753e49_r8*t18)
        u(14)=29777856511.31303518396996_r8*t6+405725357.0120253215335119_r8*t14&
             +12581472.10822014648877845_r8*t18-11128041.92193217407836428_r8*t24&
             -3746199.514655045078233992_r8*t28+7174530.658451246530535730_r8*t26&
             -9965848.257226286974264551_r8*t16+2313.542955855017656723246_r8*t38&
             +7954048017.175341662787203_r8*t12+52707258002.79912630559216_r8*t8&
             +1560495.910682062257811549_r8*t30+14108384.19513404699749891_r8*t22&
             +35147200819.52248373734473_r8*t10+123538.5886859896169363174_r8*t34&
             -506576.1722262477285863913_r8*t32-21294.73752084727438922035_r8*t36&
             -14687754.72489477423506592_r8*t20+358625041.4599066211836996_r8*t2&
             +6078341745.370762847643443_r8*t4-119.1524269109880338226565_r8*t40&
             +2674542.439978280740741820_r8
        u(15)=-414790432.6994505167392778_r8*t+6945.991168730899603944917_r8+&
             510652.2760188471673646379_r8*t2+3506139.062413368413555602_r8*t18&
             -4424351.480098610189178738_r8*t20+4507540.556879283080170755_r8*t22&
             +2493813.214756674570007723_r8*t26-3725064.570331952034347715_r8*t24&
             +5707621.898105787921162792_r8*t4-19329184319.93965811896088_r8*t3&
             -204822062695.9835996571353_r8*t5-746830449214.8311385338279_r8*t7&
             -1071759390145.197565881790_r8*t9+21556032.59225403242059432_r8*t6&
             +35040737.47008682938237320_r8*t8+24773160.94128538674006004_r8*t10&
             +1461033.332711555204742717_r8*t14+5971318.466177208703908447_r8*t12&
             -2216567.774820235631058887_r8*t16-1343045.511157374057238458_r8*t28&
             -190443.7698884165190432579_r8*t32+574046.2785186452537735790_r8*t30&
             +47314.98337952692225675031_r8*t34+913.5019396509082593070329_r8*t38&
             -8288.707145438840911983684_r8*t36-47.66097076439521352906259_r8*t40&
             +2815957.942354381680464703_r8*t29-4949452.062439341149009949_r8*t27&
             +6883478.737523361808313990_r8*t25-122150797756.7143021925426_r8*t13&
             +16322778.36996041226643312_r8*t17-5600428212.110544205675050_r8*t15&
             -610079501215.4117897742389_r8*t11-4100118.049694015514682209_r8*t19&
             +6373075.827367999613237065_r8*t21-7547338.888154381034617815_r8*t23&
             -110849.0793687075287393032_r8*t35+432724.8534694151550145610_r8*t33&
             -1257756.708472409889487035_r8*t31+119.1524269109880338226565_r8*t41&
             -2244.037373490274636993364_r8*t39+19926.89640129648768991938_r8*t37
        v(0)=u(0)
        v(1)=u(1)+0.5_r8*t*u(0)
        v(2)=-.1241319444444444444444444_r8-.2838541666666666666666667_r8*t2&
            +.1302083333333333333333333e-1_r8*t4
        v(3)=.6252170138888888888888889_r8*t+.5749059606481481481481482_r8*t3&
            -.8773871527777777777777778e-1_r8*t5+.4385850694444444444444444e-1_r8*t7&
            -.9746334876543209876543210e-2_r8*t9
        v(4)=-.3816782005529835390946500e-2_r8*t6-.3043902864181455761316872_r8&
            +.1385806990258487654320987e-1_r8*t8-.3045729648919753086419753e-2_r8*t10&
            -3.334384192949459876543209_r8*t2-1.443856321735146604938272_r8*t4
        v(5)=5.084628339853796939300410_r8*t+19.57347561662252711334020_r8*t3&
            +.3028328551887274489271016_r8*t11-.5607805781707374951009209_r8*t9&
            +.5729826674329610798794825_r8*t7+5.264663972579893261316874_r8*t5&
            +.1231848290451082696453067e-1_r8*t15-.9238862178383120223397999e-1_r8*t13
        v(6)=-28.17555842815603248376414_r8*t6-.2906392060283023236943954e-1_r8*t14&
            -2.502684474788043402799041_r8+.3849525907659633426415834e-2_r8*t16&
            +.9037170913188460136551619e-1_r8*t12+.4309883571133621996252204_r8*t8&
            -.1608534918423846448703581_r8*t10-64.67370051767834022936110_r8*t2&
            -129.7003433719719481799922_r8*t4
        v(7)=86.57218967207372000385268_r8*t+860.9772677404372680912355_r8*t3&
            +1121.756470592261700060500_r8*t5+257.5263468684699442083095_r8*t7&
            -92.32733782717116555772998_r8*t9-76.18059281400228878625290_r8*t13&
            -14.32808718833711404666971_r8*t17+40.61769611078856415222867_r8*t15&
            +99.05203232887895800013400_r8*t11+3.015390918777700925457745_r8*t19&
            -.2871800875026381833769281_r8*t21
        v(8)=-8015.583638506729061155720_r8*t6-24.41008329244471177080736_r8*t14&
            -4.519609259635658228333184_r8*t18+12.88479579364600159304798_r8*t16&
            +32.09421420855616969862385_r8*t12-1043.207983739144800727558_r8*t8&
            -.8974377734457443230529005e-1_r8*t22-23.18565586367109152015149_r8*t10&
            +.9467968509852602608208101_r8*t20-1986.189381518536462868531_r8*t2&
            -9133.569273415523777237354_r8*t4-40.56325518941026578943259_r8
        v(9)=2332.949345485996505888585_r8*t+36760.58162380117473168501_r8*t3&
            +99845.48057991644426620303_r8*t5+65183.71340269782868490893_r8*t7&
            +2622.201227310636710039871_r8*t9+3.380750161788867217892261_r8*t27&
            -45.64012718414970744154552_r8*t25-8306.938556402440897450785_r8*t13&
            -5467.955632719681768619915_r8*t17+7751.066214952751247768557_r8*t15&
            +6712.462538057135629060264_r8*t11+2875.684484702232084974429_r8*t19&
            -1094.638706632625876620833_r8*t21+285.3669009128752146814545_r8*t23
        v(10)=-80675.59761993737478733866_r8*t2-1745.835002504339507610403_r8*t18&
             +911.9476960925337440830110_r8*t20-345.5461064752209510674120_r8*t22&
             -14.31536396632473462576254_r8*t26+89.78156942102429575363385_r8*t24&
             -629192.8778816414355475387_r8*t4-1163750.663354504353033275_r8*t6&
             -560765.7238305186542428340_r8*t8-50992.88361699267893018953_r8*t10&
             -2718.991768401117794132920_r8*t14+2481.071341994842760629068_r8*t12&
             +2498.827701212957733409790_r8*t16+1.056484425559021005591331_r8*t28&
             -1118.965893570671438005353_r8
        v(11)=93680.92053187578373908477_r8*t+2244650.940705589178033889_r8*t3+&
             10624527.77751646968478233_r8*t5+14380781.82932476042057132_r8*t7&
             +5381368.427872546425567057_r8*t9-5209.003030909290742154185_r8*t29&
             +25188.20188878916819094802_r8*t27-85058.16289349076452263097_r8*t25&
             -516900.3318170604402267402_r8*t13-727425.9526984207442109454_r8*t17&
             +687103.9096254180006068751_r8*t15+718764.2117855061019497234_r8*t11&
             +612373.6241206691026703661_r8*t19-408052.1263148301297946975_r8*t21&
             +212774.0620423199828384262_r8*t23-40.72442779420311834782560_r8*t33&
             +671.9530586043514527391227_r8*t31
        v(12)=-4666637.427699311288632484_r8*t2-236576.3442879792221912621_r8*t18&
             +197095.9817518425408469117_r8*t20-130314.1457376034450561984_r8*t22&
             -26890.53644971165967586539_r8*t26+67556.99854178264515084666_r8*t24&
             -56407070.38512263316820071_r8*t4-181816409.6548781888410872_r8*t6&
             -188159903.2775835541339353_r8*t8-59729728.73424203726696958_r8*t10&
             -158804.8694791222203286819_r8*t14-3778065.140617472679475600_r8*t12&
             +226299.4976083004679687682_r8*t16+7937.202199806361407713609_r8*t28&
             +210.6216499981442527051607_r8*t32-1637.030496200573732230987_r8*t30&
             -12.72638368568847448369550_r8*t34-45598.90544342769885089189_r8
       v(13)=5335496.618522339267627189_r8*t+183576348.0194920448725676_r8*t3&
             +1354176851.693565444970589_r8*t5+3188659860.422165550420750_r8*t7&
             +2627237695.599145398984883_r8*t9-8559691.708355332039465864_r8*t29&
             +20693853.48001686619594154_r8*t27-39924064.49900019391213026_r8*t25&
             +10733492.58653441250421101_r8*t13-72327111.51346468735941796_r8*t17&
             +51092798.08238289191185946_r8*t15+717509281.1168328033139334_r8*t11&
             +83750832.87636387843170424_r8*t19-79752878.22653923708993263_r8*t21&
             +62400837.71297169795531112_r8*t23+114641.0994472638480012588_r8*t35&
             -668909.4451481646322592866_r8*t33+2760161.897061770530576782_r8*t31&
             +635.4796101919361803875012_r8*t39-12391.85239874275551755627_r8*t37
       v(14)=-2660954.178544058526885368_r8-360244269.1046558724915919_r8*t2&
             -23905166.63593070422718515_r8*t18+27364022.18954367225325487_r8*t20&
             -25828929.45718207123186239_r8*t22-12760534.41026642324470170_r8*t26&
             +20063135.60588227514388882_r8*t24-6152008132.743840646273932_r8*t4&
             -30330498547.79309275046998_r8*t6-53977141029.83300025224457_r8*t8&
             -36163818479.23891864661839_r8*t10-425356977.9744537092972964_r8*t14&
             -8214360062.046573670722601_r8*t12+18349861.97314849584742475_r8*t16&
             +6582218.353750901954628165_r8*t28+871571.7482998278264261389_r8*t32&
             -2711831.071120029682914462_r8*t30-210623.3475234238287272353_r8*t34&
             -3882.383243516360102054890_r8*t38+36005.76814847969026409722_r8*t36&
             +198.5873781849800563710940_r8*t40
       v(15)=416127703.9194396571096487_r8*t+6945.991168730899603944917_r8&
             +510652.2760188471673646379_r8*t2+3506139.062413368413555602_r8*t18&
             -4424351.480098610189178738_r8*t20+4507540.556879283080170755_r8*t22&
             +2493813.214756674570007723_r8*t26-3725064.570331952034347715_r8*t24&
             +5707621.898105787921162792_r8*t4+19508496840.66961142955273_r8*t3&
             +207861233568.6689810809570_r8*t5+761719377470.4876561258129_r8*t7&
             +1098113019146.597129034586_r8*t9+21556032.59225403242059432_r8*t6&
             +35040737.47008682938237320_r8*t8+24773160.94128538674006004_r8*t10&
             +1461033.332711555204742717_r8*t14+5971318.466177208703908447_r8*t12&
             -2216567.774820235631058887_r8*t16-1343045.511157374057238458_r8*t28&
             -190443.7698884165190432579_r8*t32+574046.2785186452537735790_r8*t30&
             +47314.98337952692225675031_r8*t34+913.5019396509082593070329_r8*t38&
             -8288.707145438840911983684_r8*t36-47.66097076439521352906259_r8*t40&
             -4689057.699681904219581699_r8*t29+8536717.391664964414277814_r8*t27&
             -12447499.69848944884749613_r8*t25+126127821765.3019730239362_r8*t13&
             -21305702.49857355575356540_r8*t17+5803290890.616556866441806_r8*t15&
             +627653101625.1730316429113_r8*t11+10390854.10380408875907143_r8*t19&
             -13716953.18981538673077002_r8*t21+14601530.98572140453336727_r8*t23&
             +172618.3737117023372074619_r8*t35-686012.9395825390193077566_r8*t33&
             +2038004.663813441018392810_r8*t31-178.7286403664820507339848_r8*t41&
             +3400.808851417783465354987_r8*t39-30574.26516172012488452956_r8*t37
        mu4=mu2*mu2
        tts=t2p1*t2p1*t2p1
        ttss=sqrt(tts)
        IF (mode==0) THEN
          xargu=exp(-a*(asin(t)+t*sqrt(t2p1)))
        ELSE
          xargu=1.0_r8          
        ENDIF
        ! ------------------------------
        ! Computation of l(mu) 
        ! -----------------------------
        DO k=0,7
          IF (k==0) THEN
            mus=1.0_r8
            lmus=0.0_r8
          ELSE
            mus=mus*mu4
          ENDIF
          lmus=lmus+lk(k)/mus
        ENDDO
        lmus=sqrt(sqrt2)/sqrt(mu)*lmus
        k=0
        DO WHILE (k<15)
          IF (k==0) THEN
            mus1=1.0_r8
            s1=0.0_r8
            s2=0.0_r8
          ELSE
            mus1=-mus1*mu2*ttss
          ENDIF
          delta1=u(k)/mus1
          delta2=v(k)/mus1
          s1=s1+delta1
          s2=s2+delta2
          k=k+1
        ENDDO
        facu=xpowy(t2p1,0.25_r8)
        wax=lmus*xargu*s1/(sqrt2*facu);
        waxp=-0.5_r8*mu*lmus*xargu*facu*s2
    ENDIF 
    END SUBROUTINE expaelem

    SUBROUTINE expair(a,x,mode,wax,waxp,waxm,waxmp,ierr)
    ! ---------------------------------------------------
    ! Airy-type asymptotic expansions for the W-function
    ! ---------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   wax,   W(a,x)
    !   waxp,  W'(a,x)
    !   waxm,  W(a,-x)
    !   waxmp, W'(a,-x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ----------------------------------------------------
    !  The Airy-type asymptotic expansions are used  
    !  around the turning point x*x/4-a=0, a>0
    ! ----------------------------------------------------
    USE Someconstants
    USE AiryFunction
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    INTEGER, INTENT(IN)  :: mode
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    REAL(r8), INTENT(OUT) :: waxm
    REAL(r8), INTENT(OUT) :: waxmp
    INTEGER, INTENT(INOUT) :: ierr
    REAL(r8) :: eps, mu, mu2, mu4, mu8, mu13, mu23, mu43, mu83, t, t2, &
                t2p1, zeta, y, zmasf, psi, sq, etta, eta, argu, phis,&
                chi, sas, sbs, scs, sds, as, bs, asp, bso, bsp, bspo, etal,&
                twom13, mus, mu2k, f2, air, dair, bir, dbir, dfacu, dfacup,&
                dfacv, dfacvp, ffa, lmus, expia, mulm1, mulm2
    REAL(r8), DIMENSION(0:10,0:40) :: ast, bst
    REAL(r8), DIMENSION(0:30) :: fik,chik
    REAL(r8), DIMENSION(0:18) :: cozmas
    REAL(r8), DIMENSION(0:15) :: lk
    INTEGER, DIMENSION(0:6) :: inda, indb
    INTEGER k, j, l
    DATA inda/40,30,30,30,20,20,20/
    DATA indb/37,30,30,30,20,20,20/
    DATA cozmas/1.885618083164126731735585_r8,.2828427124746190097603378_r8,&
         -.2525381361380526872860159e-1_r8,.4910463758239913363894753e-2_r8,&
         -.1255516301822705121450363e-2_r8,.3718259816936472859679921e-3_r8,&
         -.1208434440504353679395974e-3_r8,.4188900896706268006309575e-4_r8,&
         -.1522610358835666495714500e-4_r8,.5739999368626520519558632e-5_r8,&
         -.2227369320217030245089599e-5_r8,.8848730844862201973674137e-6_r8,&
         -.3584555319099271632854106e-6_r8,.1476133191751092628648807e-6_r8,&
         -.6164726751264643754437699e-7_r8,.2605998126670963041648663e-7_r8,&
         -.1113366163939335549490076e-7_r8,.4801280953394988359287492e-8_r8,&
         -.2087736482939914810074796e-8_r8/
    DATA ast(1,0:39)/-.86458333333333333333e-2_r8,.11118506493506493506e-1_r8,&
       -.11338461538461538462e-1_r8,.10485072951739618406e-1_r8,-.91750414291590762179e-2_r8,&
        .77467458168452212003e-2_r8,-.63785507498774352182e-2_r8,.51549769720023356846e-2_r8,&
       -.41065415400340848415e-2_r8,.32340443813071477202e-2_r8,-.25232263480572751191e-2_r8,&
        .19534081514864215152e-2_r8,-.15023865064299121639e-2_r8,.11490357422627664756e-2_r8,&
       -.87453393997191988595e-3_r8,.66279160411710427603e-3_r8,-.50044124329147178642e-3_r8,&
        .37660561306010746994e-3_r8,-.28257320676903403086e-3_r8,.21145378792045091763e-3_r8,&
       -.15785247292738707697e-3_r8,.11758031706072545009e-3_r8,-.87407457034219086216e-4_r8,&
        .64858335619867878134e-4_r8,-.48045267802574547378e-4_r8,.35535265384237208329e-4_r8,&
       -.26244779856397880734e-4_r8,.19357337318355832102e-4_r8,-.14259616234056385015e-4_r8,&
        .10492181059381593327e-4_r8,-.77117350891421712751e-5_r8,&
        .56623460684682441023388000908606482749399795844647e-5_r8,&
        -.41536087769992473560528795985505977591447090572397e-5_r8,&
        .30441367115808179827127378360809570747317003238776e-5_r8,&
         -.22291254179152811330836670701699386921930200692429e-5_r8,&
        .16310120592387743503983919363881550173527463656721e-5_r8,&
       -.11924797593173726521905994130534586329421615084496e-5_r8,&
         .87123116433670338207099105124965843641863370154514e-6_r8,&
        -.63609261679946347627956411620578214436890625764449e-6_r8,&
        .46411616429358644143012047921855806058625426820438e-6_r8/
    DATA ast(2,0:30)/.56266224308317252266e-2_r8,-.11647291855662719633e-1_r8,&
        .17438935550993679015e-1_r8,-.22248420786275842558e-1_r8,.25681645587606477421e-1_r8,&
       -.27651752813946224558e-1_r8,.28276480334473051867e-1_r8,-.27784472881795235575e-1_r8,&
        .26445276180024289855e-1_r8,-.24523529386345315127e-1_r8,.22252947171162342216e-1_r8,&
       -.19824623260148005537e-1_r8,.17384726703684670525e-1_r8,-.15037729762938088493e-1_r8,&
        .12852408813094941275e-1_r8,-.10868797174027507674e-1_r8,.91049867228358693020e-2_r8,&
       -.75631853088958125607e-2_r8,.62347763146567497230e-2_r8,-.51043369269941506053e-2_r8,&
        .41526903371370890122e-2_r8,-.33591243802274860175e-2_r8,.27029276999493624052e-2_r8,&
       -.21643904392277104271e-2_r8,.17254005528927170693e-2_r8,-.13697460566867932515e-2_r8,&
        .10832120482574887245e-2_r8,-.85354146066446264128e-3_r8,.67031131907963448661e-3_r8,&
       -.52476209811926895365e-3_r8,.40960652212630493071e-3_r8/
    DATA ast(3,0:30)/-.11580906180050927561e-1_r8,.31661887057447902497e-1_r8,&
       -.60537519726707374400e-1_r8,.96022104837560466155e-1_r8,&
       -.13486858218773139158_r8,.17360451334378319335_r8,&
       -.20913158018379167172_r8,.23907654949776246859_r8,&
       -.26192548161477102970_r8,.27699864052741948107_r8,&
       -.28432695853487095482_r8,.28448186552194680157_r8,&
       -.27839658337617157169_r8,.26720327422032545919_r8,&
       -.25209898969048886827_r8,.23424493493810997748_r8,&
       -.21469801828258118088_r8,.19437049669898026028_r8,&
       -.17401212913960327669_r8,.15420903436758014570_r8,&
       -.13539394342452291497_r8,.11786338661988688362_r8,&
       -.10179832039934106068_r8,.87285631081354269540e-1_r8,&
       -.74338768749844687828e-1_r8,.62916431192973831951e-1_r8,&
       -.52938730452196001459e-1_r8,.44300646423185523639e-1_r8,&
       -.36882824571960708148e-1_r8,.30559932289668292664e-1_r8,&
       -.25206873839963280882e-1_r8/
    DATA ast(4,0:30)/.49393729669138401882e-1_r8,-.16460159056974091096_r8,&
        .37671106920022966088_r8,-.70446154818843445310_r8,&
        1.1517696715633945107_r8,-1.7071311993240408422_r8,&
        2.3458575390237919962_r8,-3.0341470817455636226_r8,3.7338648123687769017_r8,&
       -4.4071025283752592817_r8,5.0199043219968235593_r8,-5.5448551987388746055_r8,&
        5.9624855820903817327_r8,-6.2616193802379769230_r8,6.4388904169360639081_r8,&
       -6.4976858367265413542_r8,6.4467642666502533981_r8,-6.2987589550594794269_r8,&
        6.0687264338988964841_r8,-5.7728499953831028737_r8,5.4273610811194061563_r8,&
       -5.0477040039578825218_r8,4.6479413723592056878_r8,-4.2403787420074181144_r8,&
        3.8353760898208164707_r8,-3.4413090294005663100_r8,3.0646425650883640354_r8,&
       -2.7100830787497833682_r8,2.3807788873442699130_r8,-2.0785451219757146248_r8,&
        1.8040941616937906026_r8/
    DATA ast(5,0:28)/-.36006972101648543486_r8,1.3987393559898456251_r8,&
       -3.6892362137755083799_r8,7.8737303343225028183_r8,&
       -14.569386493592020890_r8,24.261162048097651407_r8,&
       -37.213581231065030993_r8,53.417291724674799405_r8,&
       -72.576260862138309707_r8,94.133040133549441048_r8,&
       -117.32368949349208833_r8,141.25103106722442453_r8,&
       -164.96454465710299113_r8,187.53671898263579271_r8,-208.12824122199390069_r8,&
        226.03734163390257598_r8,-240.73138542763382112_r8,251.86109378733459125_r8,&
       -259.25942968790957949_r8,262.92818670198517035_r8,-263.01574432667643324_r8,&
        259.78942292465782198_r8,-253.60552167078665437_r8,244.87958311821722765_r8,&
       -234.05880699965059900_r8,221.59791652287623744_r8,-207.93921917262608293_r8,&
        193.49713382040369616_r8,-178.64709066614237649_r8/
    DATA ast(6,0:25)/4.0084065348046763094_r8,-17.642621241600336767_r8,&
        52.311391233887190108_r8,-124.66354447894727130_r8,256.05306717453359351_r8,&
       -470.80863230664034499_r8,793.63941447186478508_r8,-1246.6101080073623060_r8,&
        1846.1783018847890331_r8,-2600.7449029420965852_r8,3509.0494384975085126_r8,&
       -4559.5844595960940674_r8,5731.0427527316521420_r8,-6993.6748735645957022_r8,&
        8311.3386629781516967_r8,-9643.9722024686523983_r8,10950.214002409855203_r8,&
       -12189.920564766787596_r8,13326.380783686538323_r8,-14328.087748340125799_r8,&
        15169.991649192512856_r8,-15834.215245207554837_r8,16310.260903736659595_r8,&
       -16594.773183346838010_r8,16690.942899479998399_r8,-16607.648660184592203_r8/
    DATA ast(7,0:22)/-63.290222545472642845_r8,309.49063765889768949_r8,&
       -1013.7206553394096370_r8,2655.5227786862968090_r8,-5969.3059766708673407_r8,&
        11964.743360104188382_r8,-21906.959390128289116_r8,37252.392923816787966_r8,&
       -59544.188680047082730_r8,90278.109242533664710_r8,-130754.87888773039296_r8,&
        181936.88302880890521_r8,-244326.24121616523617_r8,317877.96339339844357_r8,&
       -401957.01966033293854_r8,495342.63338040766157_r8,-596277.79843283893916_r8,&
        702557.58142432196808_r8,-811646.60776554538511_r8,920814.40627216124150_r8,&
       -1027276.9472424641737_r8,1128333.5453824016662_r8,-1221490.0112516386311_r8/
    DATA bst(0,0:36)/-.32142857142857142857e-1_r8,.15555555555555555556e-1_r8,&
       -.10085343228200371058e-1_r8,.68923076923076923077e-2_r8,&
       -.47973770749280953363e-2_r8,.33672793089263677499e-2_r8,&
       -.23740671688092001398e-2_r8,.16782309959031527659e-2_r8,&
       -.11883320683635189568e-2_r8,.84238779139280588983e-3_r8,&
       -.59762438660929409956e-3_r8,.42422138636439454428e-3_r8,&
       -.30126023626926855833e-3_r8,.21400917293042855689e-3_r8,&
       -.15206634775554541912e-3_r8,.10807399278785359819e-3_r8,&
       -.76820938020122548526e-4_r8,.54612918802955620593e-4_r8,&
       -.38829207461698758135e-4_r8,.27609666542412784727e-4_r8,&
       -.19633469331004862369e-4_r8,.13962435063909612290e-4_r8,&
       -.99300043010631579079e-5_r8,.70625008660394052299e-5_r8,&
       -.50232596984969135714e-5_r8,.35729626884250873769e-5_r8,&
       -.25414708010873361644e-5_r8,.18078147381341295781e-5_r8,&
       -.12859778274667428677e-5_r8,.91479251956680900289e-6_r8,&
       -.65075912613628042224e-6_r8,.46294093040524422199e-6_r8,&
       -.32933491162045547170e-6_r8,.23429130240268590734e-6_r8,&
       -.16667872738882887922e-6_r8,.11857940849768171016e-6_r8,&
       -.84361252507627603914e-7_r8/
    DATA bst(1,0:30)/.11616363324175824176e-1_r8,-.10738690547767928720e-1_r8,&
        .11273908748814210999e-1_r8,-.11330200650226966016e-1_r8,&
        .10882662128717442695e-1_r8,-.10073431983586506837e-1_r8,&
        .90528041420056262576e-2_r8,-.79437161494535131346e-2_r8,&
        .68354205437865891183e-2_r8,-.57866877852097172840e-2_r8,&
        .48319188837933835377e-2_r8,-.39875048709666930997e-2_r8,&
        .32573797517277923555e-2_r8,-.26374472093943733107e-2_r8,&
        .21188963477292799694e-2_r8,-.16905597364664082581e-2_r8,&
        .13405067192519113713e-2_r8,-.10570579364055743550e-2_r8,&
        .82938060092290405418e-3_r8,-.64779233018712147103e-3_r8,&
        .50387108159099364466e-3_r8,-.39044280607171355606e-3_r8,&
        .30149757745508396450e-3_r8,-.23206892067503710671e-3_r8,&
        .17809916581662199089e-3_r8,-.13630510145467900797e-3_r8,&
        .10405223392871875512e-3_r8,-.79241928092168523010e-4_r8,&
        .60213081558437480541e-4_r8,-.45658357387345341779e-4_r8,&
        .34554059520366810843e-4_r8/
    DATA bst(2,0:30)/-.17619791271984698754e-1_r8,.22460738436828023930e-1_r8,&
       -.31095536623705537180e-1_r8,.39843610073496917522e-1_r8,&
       -.47521689852425072253e-1_r8,.53475980779992135034e-1_r8,&
       -.57414884418372632722e-1_r8,.59320052235418060306e-1_r8,&
       -.59362779310300406734e-1_r8,.57828321677111474388e-1_r8,&
       -.55053943676791592811e-1_r8,.51382942898551505119e-1_r8,&
       -.47133919046907858690e-1_r8,.42582911713303221181e-1_r8,&
       -.37955441784732551387e-1_r8,.33425554856704672102e-1_r8,&
       -.29119366448832191622e-1_r8,.25121136014826860211e-1_r8,&
       -.21480423732279927589e-1_r8,.18219346128000877782e-1_r8,&
       -.15339318528807313005e-1_r8,.12826952106907783895e-1_r8,&
       -.10658970903423837305e-1_r8,.88061445510104236352e-2_r8,&
       -.72363110956985624537e-2_r8,.59166055051532985471e-2_r8,&
       -.48150249286128209578e-2_r8,.39014607485083487893e-2_r8,&
       -.31483167600335778722e-2_r8,.25308172470515622861e-2_r8,&
       -.20270915098905906085e-2_r8/
    DATA bst(3,0:30)/.60909763139637582786e-1_r8,-.96541486771214866841e-1_r8,&
        .16264377675453827437_r8,-.24915885452321001414_r8,&
        .35009667669419877168_r8,-.45836730339191387087_r8,&
        .56648575100308131913_r8,-.66748754950937166485_r8,&
        .75561704712993985977_r8,-.82671466079887381067_r8,&
        .87832844401420206021_r8,-.90961364830930523024_r8,&
        .92109362061208216666_r8,-.91434938542624464396_r8,&
        .89169172213458344119_r8,-.85585378144126278784_r8,&
        .80972748910911579700_r8,-.75615478142621795241_r8,&
        .69777564969133675222_r8,-.63692892519825209993_r8,&
        .57559825199333248396_r8,-.51539418379669499028_r8,&
        .45756322033392356537_r8,-.40301535937866445856_r8,&
        .35236298156274064514_r8,-.30596531074994938002_r8,&
        .26397410487303626816_r8,-.22637751012529200117_r8,&
        .19304009378463816988_r8,-.16373793768410157069_r8,&
        .13818833224511412946_r8/
    DATA bst(4,0:29)/-.37829872479673768094_r8,.70699348356478890943_r8,&
        -1.3865797541213790568_r8,2.4460312321724494685_r8,&
        -3.9207588293947925720_r8,5.8080313668527459768_r8,&
        -8.0630182513775208787_r8,10.603669134030102425_r8,&
       -13.320620699455056397_r8,16.089517395703050953_r8,&
       -18.783542996117974847_r8,21.284502840495644149_r8,&
       -23.491416999286173127_r8,25.326169364855214008_r8,&
       -26.736231662339052310_r8,27.694813083767226802_r8,&
       -28.198977054536741489_r8,28.266337819667546672_r8,&
       -27.930931277990815536_r8,27.238778053930476484_r8,&
       -26.243549745521884151_r8,25.002633225307177664_r8,&
       -23.573777915301600688_r8,22.012416302016647166_r8,&
       -20.369672976431342515_r8,18.691023002223376311_r8,&
       -17.015524914036864532_r8,15.375534273694109775_r8,&
       -13.796797055459840802_r8,12.298824764252335256_r8/
    DATA bst(5,0:26)/3.7008098833796096974_r8,-7.8943179119715393505_r8,&
        17.523666638205992170_r8,-34.732110691221732984_r8,&
        62.145391921905206236_r8,-102.17014739446708717_r8,&
        156.60003201213783498_r8,-226.31321551629452367_r8,311.08916186164296289_r8,&
       -409.56133862117968200_r8,519.30074205116124335_r8,-637.00727007547042847_r8,&
        758.77505311063660762_r8,-880.39404139944726093_r8,997.65225418279181001_r8,&
       -1106.6093103403454609_r8,1203.8202784542850274_r8,-1286.4978451886305591_r8,&
        1352.6090567441321809_r8,-1400.9096495972083876_r8,1430.9238827471423665_r8,&
       -1442.8807756116574511_r8,1437.6189354564596975_r8,-1416.4720513211302631_r8,&
        1381.1460162370242338_r8,-1333.5968879195104838_r8,1275.9168367360627308_r8/
    DATA bst(6,0:22)/-52.440232872505911846_r8,124.91492555422010956_r8,&
       -307.80505874513857089_r8,673.65965906811191089_r8,-1324.7001416829727542_r8,&
        2383.2377235388771970_r8,-3981.7258172782825530_r8,6249.8655385949204720_r8,&
       -9300.3111310734301499_r8,13214.918162292874992_r8,-18033.350332471047341_r8,&
        23745.434126825075803_r8,-30288.059814249641659_r8,37546.796886266045900_r8,&
       -45361.823914328331689_r8,53537.336176016916081_r8,-61853.323266034541002_r8,&
        70078.506118453102599_r8,-77983.267923314792836_r8,85351.571541509186577_r8,&
       -91991.086631042572063_r8,97741.013229975714811_r8,-102477.35110602911554_r8/
    DATA fik/1.0_r8,-.1_r8,.42142857142857142857e-1_r8,-.21912698412698412698e-1_r8,&
         .12448755411255411255e-1_r8,-.74288549743906886764e-2_r8,&
         .45746350641201831678e-2_r8,-.28791700676332804184e-2_r8,&
         .18414150026089362744e-2_r8,-.11923004331977092900e-2_r8,&
         .77957025714134722714e-3_r8,-.51376250832318276629e-3_r8,&
         .34081277524518474479e-3_r8,-.22733442456475481892e-3_r8,&
         .15235591836762461067e-3_r8,-.10252250859554482318e-3_r8,&
         .69234005415816302727e-4_r8,-.46900261834713094097e-4_r8,&
         .31859086498603514106e-4_r8,-.21695264658176838875e-4_r8,&
         .14806797927678100534e-4_r8,-.10125785547514577093e-4_r8,&
         .69372508648302511184e-5_r8,-.47606638519588906843e-5_r8,&
         .32719619838504728166e-5_r8,-.22519366067543023254e-5_r8,&
         .15519011644445928413e-5_r8,-.10707548626041065912e-5_r8,&
         .73959982335703643779e-6_r8,-.51138836489547516290e-6_r8,&
         .35393386504561175831e-6_r8/
    DATA chik/-.1_r8,.74285714285714285714e-1_r8,-.54095238095238095238e-1_r8,&
         .39063615749330035044e-1_r8,-.28085509411223696938e-1_r8,&
         .20139984242514854760e-1_r8,-.14417859514984444956e-1_r8,&
         .10309442176430897606e-1_r8,-.73655122348323755815e-2_r8,&
         .52589148288898677708e-2_r8,-.37529918053932661301e-2_r8,&
         .26772689420573039019e-2_r8,-.19092899409669881562e-2_r8,&
         .13612619157861153785e-2_r8,-.97033149280715834007e-3_r8,&
         .69154709536062995295e-3_r8,-.49278576070368617280e-3_r8,&
         .35110626274417083999e-3_r8,-.25013277029410498759e-3_r8,&
         .17818059964612554284e-3_r8,-.12691507492553136865e-3_r8,&
         .90392689282224654705e-4_r8,-.64376056356811587852e-4_r8,&
         .45844739898716106234e-4_r8,-.32646108214340078356e-4_r8,&
         .23246222222742878232e-4_r8,-.16552150417025056676e-4_r8,&
         .11785263078514715153e-4_r8,-.83908989293805415224e-5_r8,&
         .59739749862631200831e-5_r8,-.42530962002148905661e-5_r8/
    eps=epss
    mu=sqrt(2.0_r8*a); t=abs(x)/(mu*sqrt2);  t2=t*t
    ierr=0
    IF (mode==0) THEN   
      IF (t <= 1) THEN
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
        t2p1=1.0_r8-t2
        IF (a*abs(asin(t)+t*sqrt(t2p1))>log(giant)) THEN
          ierr=1
          waxm=giant
          waxmp=giant
          wax=dwarf
          waxp=dwarf
        ENDIF
      ELSE
        ffa=pi*a*0.5_r8
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
        IF ((ffa < log(dwarf)).OR.(ffa>log(giant))) THEN
          ierr=1
          waxm=giant
          waxmp=giant
          wax=dwarf 
          waxp=dwarf
        ENDIF
      ENDIF
    ENDIF
    IF (ierr==0) THEN
      mu2=mu*mu
      mu4=mu2*mu2
      mu8=mu4*mu4
      mu13=xpowy(mu,onethird)
      mu23=mu13*mu13
      mu43=xpowy(mu4,onethird)
      mu83=xpowy(mu8,onethird)
      lk(0)=1.0_r8
      lk(1)=-.86805555555555555555555555555555555555555555555556e-3_r8;
      lk(2)=-.40496700585133744855967078189300411522633744855967e-3_r8;
      lk(3)=-.55970973301118646421903237746859146036100768611057e-3_r8;
      lk(4)=-16.943102322943963245388661733269423877173517985977_r8;
      lk(5)=-.94037146647019565951398813260001779563624393715061e-2_r8;
      lk(6)=-.84327434386741291387140061090400783753265461929755e-1_r8;
      lk(7)=-1.1162982841593720168762155424414027479385891257438_r8;   
      IF (abs(1.0_r8-t)<eps*100) THEN
        psi=0.0_r8
        zeta=0.0_r8
      ELSEIF (t > 1.1_r8) THEN
        psi=0.5_r8*(t*sqrt(t2-1.0_r8)-log(t+sqrt(t2-1.0_r8)))
        zeta=xpowy(3*psi*0.5_r8,twothird)
      ELSEIF (t < 0.9_r8) THEN
        sq=sqrt((1.0_r8-t)*(1.0_r8+t))
        etta=0.5_r8*(acos(t)-t*sq)
        zeta=-xpowy(3*etta*0.5_r8,twothird)
      ELSE
        ! Taylor series
        y=t-1.0_r8
        zmasf=cozmas(12)
        DO k=0,11
          j=11-k
          zmasf=cozmas(j)+y*zmasf
        ENDDO
        zmasf=0.5_r8*xpowy(abs(y),1.5_r8)*zmasf
        IF (t < 1) THEN
          zeta=-xpowy(1.5_r8*zmasf,twothird)
        ELSE
          zeta=xpowy(1.5_r8*zmasf,twothird)
        ENDIF
      ENDIF
      argu=-zeta*mu43
      ! -----------------------------
      ! Computation of l(mu) 
      ! -----------------------------
      DO k=0,7
        IF (k==0) THEN
          mus=1.0_r8
          lmus=0.0_r8
        ELSE
          mus=mus*mu4
        ENDIF
        lmus=lmus+lk(k)/mus
      ENDDO
      lmus=sqrt(sqrt2)/sqrt(mu)*lmus      
      IF ((t > 0.9_r8).AND.(t < 1.1_r8)) THEN
      ! Taylor series for the functions phi, chi
        phis=0.0_r8
        chi=0.0_r8
        twom13=xpowy(2.0_r8,-onethird)
        eta =zeta*twom13
        etal=eta
        DO k=0,20
          IF (k==0) THEN
            etal=1.0_r8
            phis=0.0_r8
            chi=0.0_r8
          ELSE
            etal=etal*eta
          ENDIF
          phis=phis+fik(k)*etal
          chi=chi+chik(k)*etal
        ENDDO
        phis=xpowy(2.0_r8,-onesix)*phis
        chi=twom13*chi
      ELSE
       ! Exact values
        phis=xpowy(zeta/((t-1.0_r8)*(t+1.0_r8)),0.25_r8)
        chi=0.25_r8*(1.0_r8-2.0_r8*t*xpowy(phis,6.0_r8))/zeta
      ENDIF
      sas=1.0_r8
      sbs=0.0_r8
      scs=0.0_r8
      sds=1.0_r8
      bs=0.0_r8
      bsp=0.0_r8
      twom13=xpowy(2.0_r8,-onethird)
      eta =zeta*twom13
      DO k=0,6
        bso=bs
        bspo=bsp
        as=0.0_r8
        bs=0.0_r8
        asp=0.0_r8
        bsp=0.0_r8
        IF (k > 0) THEN
!         mu2k=-1.0_r8*mu2k*mu2
          mu2k=-1.0_r8*mu2k*mu4
        ELSE
          mu2k=1.0_r8
        ENDIF
        f2=1.0_r8/mu2k
        IF (k <= 5) THEN
          DO l=0,40
            IF (l.EQ.0) THEN
              etal=1.0_r8
            ELSE
              etal=etal*eta
            ENDIF
            IF (k >0 ) THEN
              IF (l < inda(k)) THEN
                as=as+ast(k,l)*etal
              ENDIF
              IF ((l+1) < inda(k)) THEN
                asp=asp+(l+1)*ast(k,l+1)*etal
              ENDIF
            ENDIF
            IF (l < indb(k)) THEN
              bs=bs+bst(k,l)*etal
            ENDIF
            IF ((l+1) < indb(k)) THEN
              bsp=bsp+(l+1)*bst(k,l+1)*etal
            ENDIF
          ENDDO
          bs=bs/twom13
          asp=asp*twom13
          sas=sas+as*f2
          sbs=sbs+bs*f2    
        ENDIF   
        IF (k>=1) THEN
          sds=sds+(as+chi*bso+bspo)*f2
          scs=scs+(chi*as+asp+zeta*bs)*f2
        ELSE
          scs=scs+(chi+zeta*bs)*f2
        ENDIF
      ENDDO
      IF ((t<=1).AND.(mode==1)) THEN
        CALL airsca(argu,air,bir,dair,dbir)
        expia=1.0_r8
      ELSE
        CALL aibi(argu,air,bir)
        CALL aibip(argu,dair,dbir)
        IF (mode==0) THEN
          expia=exp(pi*a*0.5_r8)
        ELSE
          expia=1.0_r8
        ENDIF
      ENDIF
      mulm1=sqrtpi*mu13*lmus*phis
      mulm2=sqrtpi*mu23*lmus/phis
      dfacu=mulm1/(sqrt2*expia)
      dfacup=mulm2*0.5_r8/expia
      dfacv=mulm1*sqrt2*expia
      dfacvp=mulm2*expia
      wax=dfacu*(bir*sas+dbir*sbs/mu83)  
      waxp=-dfacup*(-bir*scs/mu43+dbir*sds)
      waxm=dfacv*(air*sas+dair*sbs/mu83)  
      waxmp=dfacvp*(-air*scs/mu43+dair*sds)
    ENDIF
    END SUBROUTINE expair


    SUBROUTINE wammx(a, x, mode, wax, waxp, waxm, waxmp, ierr)
    ! --------------------------------------------------
    ! Calculation of W(a,x), W'(a,x) by using asymptotic
    ! expansions in terms of elementary functions for 
    ! a positive and t=x/(2a**(1/2))>1
    ! --------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode,  mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   wax,  W(a,x)
    !   waxp, W'(a,x)
    !   ierr, error flag
    ! ---------------------------------------------------------
    USE someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
     INTEGER, INTENT(IN)  :: mode
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    REAL(r8), INTENT(OUT) :: waxm
    REAL(r8), INTENT(OUT) :: waxmp
    INTEGER, INTENT(INOUT) :: ierr
    REAL(r8) ka, kka, sqka, x2, x24, a2, a2xit, taut, w, &
             sr, si, tr, ti, sqa, p, q, r, s, c, chit
    REAL(r8) phik(0:8), psik(0:8)
    REAL(r8), DIMENSION(0:20) :: pjk1, qjk1, pjk2, qjk2, pjk3, qjk3
    REAL(r8), DIMENSION(0:20) :: pjk4, qjk4, pjk5, qjk5, pjk6, qjk6
    REAL(r8), DIMENSION(0:20) :: pjk7, qjk7, pjk8, qjk8
    INTEGER k, k2, k2m1
    DATA pjk1/-0.75_r8,-2.5_r8,-1.6666666666666666667_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk1/1.25_r8,3.5_r8,2.3333333333333333333_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk2/3.28125_r8,27.875_r8,67.375_r8,64.166666666666666667_r8,&
              21.388888888888888889_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8/
    DATA qjk2/-4.21875_r8,-33.625_r8,-79.958333333333333333_r8,&
              -75.833333333333333333_r8,-25.277777777777777778_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk3/-27.0703125_r8,-396.328125_r8,-1814.5125_r8,&
              -3861.2291666666666667_r8,-4254.25_r8,&
              -2363.4722222222222222_r8,-525.21604938271604938_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk3/31.9921875_r8,453.109375_r8,2047.1125_r8,&
              4330.5208333333333333_r8,4759.0277777777777778_r8,&
              2641.5277777777777778_r8,587.00617283950617284_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk4/329.91943359375_r8,7077.45703125_r8,49877.23984375_r8,&
              173584.125_r8,344694.42447916666667_r8,&
              411244.16666666666667_r8,292637.25231481481481_r8,&
              114759.70679012345679_r8,19126.617798353909465_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk4/-373.90869140625_r8,-7859.49609375_r8,-54833.22265625_r8,&
              -189793.22083333333333_r8,-375736.00052083333333_r8,&
              -447530.41666666666667_r8,-318189.01311728395062_r8,&
              -124738.81172839506173_r8,-20789.801954732510288_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk5/-5328.1988525390625_r8,-153266.177490234375_r8,&
              -1496065.7216099330357_r8,-7447882.350390625_r8,&
              -22100856.077647569444_r8,-42015571.25390625_r8,&
              -52708852.666377314815_r8,-43565653.690200617284_r8,&
              -22880216.541280864198_r8,-6933398.9519032921811_r8,&
              -924453.19358710562414_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk5/5889.0618896484375_r8,167081.378173828125_r8,&
              1618334.6748744419643_r8,8018883.473046875_r8,&
              23725519.159852430556_r8,45020597.459635416667_r8,&
              56413247.888020833333_r8,46595060.471836419753_r8,&
              24461987.833204732510_r8,7411564.3968621399177_r8,&
              988208.58624828532236_r8,0.0_r8,0.0_r8,0.0_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk6/107230.00190734863281_r8,3911743.4463500976562_r8,&
              49592558.348275320871_r8,327327256.24445452009_r8,&
              1320314667.6672144717_r8,3526512665.7394314236_r8,&
              6506075517.7748770255_r8,8449099135.7631655093_r8,&
              7726873009.8828125000_r8,4881459532.0875128601_r8,&
              2031139222.9600694444_r8,501515857.52100480110_r8,&
              55723984.169000533455_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA qjk6/-116554.34989929199219_r8,-4209790.5294799804688_r8,&
              -53046844.358616420201_r8,-348738575.26079101562_r8,&
              -1402928860.3276453993_r8,-3740348588.5545789931_r8,&
              -6891971112.2487340856_r8,-8942654096.4852912809_r8,&
              -8173643819.4593942901_r8,-5161879657.2836291152_r8,&
              -2147389212.0536479767_r8,-530173906.52220507545_r8,&
              -58908211.835800563938_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk7/-2585008.9745521545410_r8,-115109601.28864288330_r8,&
              -1814958496.4086914062_r8,-15118535470.703432792_r8,&
              -78228459666.499029134_r8,-273482436440.61774631_r8,&
              -677676841014.59605577_r8,-1223138659577.8502062_r8,&
              -1629074225394.3409650_r8,-1602934313802.8913484_r8,&
              -1152270884857.3288805_r8,-588917533590.47791281_r8,&
              -202879881562.49714220_r8,-42266641992.186904626_r8,&
              -4025394475.4463718691_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/       
    DATA qjk7/2776491.1208152770996_r8,122708978.79421234131_r8,&
              1925312041.3459734235_r8,15983788916.836100551_r8,&
              82508725486.563457380_r8,287952809783.26722765_r8,&
              712650296975.96442744_r8,1285107966580.4986798_r8,&
              1710504064563.1047815_r8,1682287248047.2539706_r8,&
              1208933528167.4145608_r8,617751113141.81053884_r8,&
              212787605947.74543705_r8,44328429406.439924364_r8,&
              4221755181.5657070823_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    DATA pjk8/72622595.878824591637_r8,3839296326.6825199127_r8,&
              72972869796.835744858_r8,741010050728.21005118_r8,&
              4728818470068.3977036_r8,20667683341938.663361_r8,&
              65097477007191.150894_r8,152479523376491.21695_r8,&
              270645546827389.17508_r8,367482071900019.28684_r8,&
              382220257297374.97871_r8,302446999586890.60270_r8,&
              179020321106247.51147_r8,76811168492401.261776_r8,&
              22575519201784.001194_r8,4067661117438.5587738_r8,&
              338971759786.54656448_r8,&
              0.0_r8,0.0_r8,0.0_r8,0.0_r8/        
    DATA qjk8/-77307924.645200371742_r8,-4062968417.3650646210_r8,&
              -76913846048.558694839_r8,-778781085218.40286473_r8,&
              -4959356633780.3453088_r8,-21641132068083.787434_r8,&
              -68082808844916.255468_r8,-159329371436569.03712_r8,&
              -282612499203202.34409_r8,-383535732852395.08975_r8,&
              -398767118978192.21973_r8,-315453595729638.25935_r8,&
              -186682859410003.59807_r8,-80088585766097.861952_r8,&
              -23536984672244.367115_r8,-4240753079882.7527641_r8,&
              -353396089990.22939701_r8,&
               0.0_r8,0.0_r8,0.0_r8,0.0_r8/
    ierr=0
    IF (mode==0) THEN
      IF (-0.5_r8*pi*a<log(dwarf)) THEN
        ierr=1
        wax=dwarf
        waxp=dwarf
        waxm=giant
        waxmp=giant
      ENDIF
    ENDIF
    IF (2.0_r8*log(x*0.5_r8)>log(giant)) THEN
      ierr=1
      wax=dwarf
      waxp=dwarf
      waxm=giant
      waxmp=giant
    ENDIF
    IF (ierr==0) THEN
      x2= x*0.5_r8; 
      x24= x2*x2; a2= 2._r8*a; sqa= sqrt(a);
      p= sqrt(x24-a);
      w= 4._r8*p*(x2+p); taut= a2/w; a2xit= x2*p;
      a2xit= a2xit-a*log((x2+p)/sqa);
      phik(0)= 1; psik(0)= 1; r= -1._r8/w;
      phik(1)= r*polysum(pjk1, taut, 2); psik(1)= r*polysum(qjk1, taut,2); r= -r/w;
      phik(2)= r*polysum(pjk2, taut, 4); psik(2)= r*polysum(qjk2, taut,4); r= -r/w;
      phik(3)= r*polysum(pjk3, taut, 6); psik(3)= r*polysum(qjk3, taut,6); r= -r/w;
      phik(4)= r*polysum(pjk4, taut, 8); psik(4)= r*polysum(qjk4, taut,8); r= -r/w;
      phik(5)= r*polysum(pjk5, taut,10); psik(5)= r*polysum(qjk5,taut,10); r= -r/w;
      phik(6)= r*polysum(pjk6, taut,12); psik(6)= r*polysum(qjk6,taut,12); r= -r/w;
      phik(7)= r*polysum(pjk7, taut,14); psik(7)= r*polysum(qjk7,taut,14); r= -r/w;
      phik(8)= r*polysum(pjk8, taut,16); psik(8)= r*polysum(qjk8, taut,16);
      r= 1.0_r8; sr= 1.0_r8; tr= 1.0_r8; si= 0.0_r8; ti= 0.0_r8;
      DO k= 1,4 
        r= -r;
        k2=2*k;
        k2m1=k2-1; 
        sr= sr+r*phik(k2); si= si+r*phik(k2m1);  
        tr= tr+r*psik(k2); ti= ti+r*psik(k2m1)
      ENDDO
      chit= piquart+rhostar(a)+a2xit; 
      s= sin(chit); c= cos(chit); p= sqrt(p);
      IF (mode==0) THEN
        q=exp(-pi*a*0.5_r8);
      ELSE
        q=1.0_r8 
      ENDIF 
      IF (-2.0_r8*pi*a<log(dwarf)) THEN
        kka=0.0_r8
      ELSE
        kka=exp(-2.0_r8*pi*a)
      ENDIF
      ka=sqrt(1.0_r8+sqrt(1.0_r8+kka))  
      sqka=q/ka;
      wax= sqka*(c*sr-s*si)/p;
      waxp=-sqka*p*(s*tr+c*ti);
      waxm= (s*sr+c*si)/(p*sqka);
      waxmp= -p*(c*tr-s*ti)/sqka;
    ENDIF
    END SUBROUTINE wammx

    FUNCTION waxatx0(a)
    ! ---------------------------------
    ! Computation of W(a,0), any a 
    ! ---------------------------------
    USE someconstants
    IMPLICIT NONE
    REAL(r8) a, waxatx0
    REAL(r8) b, g, r, s, s2, s4, s6, s8, c, w, ws, gn(0:5)
    REAL(r8) n075, n075d, n025, n025d
    INTEGER n, m
    IF (abs(a)<1.0e-18_r8) THEN
      waxatx0= 1.0227656721131686716_r8 !This is sqrt(pi/2)/gamma(0.75)
    ELSE
      b= a*a*0.25_r8;
      IF (abs(a) > 20) THEN 
        m= 0 
      ELSE
        m= 1+nint(sqrt(100-b));
      ENDIF
      w= 1.0_r8/(m*m+b); ws= sqrt(w);
      s= a*ws*0.5_r8; s2= s*s; s4= s2*s2; s6= s2*s4; s8= s4*s4;
      c= 2.0_r8*s2-1.0_r8;
      gn(0)= 1.0_r8;
      gn(1)= c/32.0_r8;
      gn(2)= (11.0_r8-84*s2+84*s4)/2048.0_r8;
      gn(3)= c*(2684.0_r8*s4-2684.0_r8*s2+173.0_r8)/65536.0_r8;
      gn(4)= (2885168.0_r8*s8-5770336.0_r8*s6+22931.0_r8-723976.0_r8*s2+&
              3609144.0_r8*s4)/8388608.0_r8;
      gn(5)= c*(334374768.0_r8*s8-668749536.0_r8*s6+397250360.0_r8*s4-&
              62875592.0_r8*s2+1319183.0_r8)/268435456.0_r8;
      g= gn(5); 
      DO n = 4,0,-1 
        g= gn(n)+w*g; 
      ENDDO
      g= sqrt(ws*g);
      r= 1; m= m-1;
      DO n= 0,m 
        n075=n+0.75_r8;
        n075d=n075*n075;
        n025=n+0.25_r8;
        n025d=n025*n025;
        r= r*(n075d+b)/(n025d+b);
      ENDDO
      g= sqrt(r)*g;
      waxatx0= sqrt(g)*exp(-0.75_r8*log(2.0_r8))
    ENDIF  
    END FUNCTION waxatx0

    SUBROUTINE waxsma(a, x, wax, waxp, waxm, waxmp)
    ! --------------------------------------------------
    ! Calculation of W(a,x), W'(a,x) by using McLaurin
    ! series, for small |a| and |x|.
    ! --------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    ! Outputs:
    !   wax,  W(a,x)
    !   waxp, W'(a,x)
    !   waxm, W(a,-x)
    !   waxmp, W'(a,-x)
    ! ---------------------------------------------------
    USE someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    REAL(r8), INTENT(OUT) :: waxm
    REAL(r8), INTENT(OUT) :: waxmp    
    REAL(r8) wa0, wa0p, w1, w1p, w2, w2p, &
             t, t1, t2, t3, t4, x2, delta 
    REAL(r8) an(0:100), bn(0:100)
    INTEGER n, n2
    wa0= waxatx0(a); wa0p= -0.5_r8/wa0;
    IF (abs(x)<1.0e-20_r8) THEN
      wax= wa0; waxp= wa0p; waxm= wa0; waxmp= wa0p;  
    ELSE
      an(0)= 1.0_r8; an(1)= a*0.5_r8; bn(0)= 1.0_r8; bn(1)= a/6.0_r8;
      x2= x*x; t= x2; n= 0; n2= 0;
      w1= 1.0_r8+an(1)*t; w2= 1.0_r8+bn(1)*t;
      w1p= a; w2p= 1.0_r8+a*t*0.5_r8;
      delta= 1.0_r8; t3= 1.0_r8; t4= 1.0_r8;
      DO WHILE ((delta > 1.0e-17_r8).AND.(n<50)) 
        an(n+2)= (a*an(n+1)-an(n)*0.25_r8)/((n2+3.0_r8)*(n2+4.0_r8));
        bn(n+2)= (a*bn(n+1)-bn(n)*0.25_r8)/((n2+4.0_r8)*(n2+5.0_r8));
        w1p= w1p+(n2+4)*an(n+2)*t;
        t= x2*t; t1= an(n+2)*t; t2= bn(n+2)*t;
        w1= w1 + t1; w2= w2 + t2;
        w2p= w2p+(n2+5)*t2;
        n= n+1; n2= n2+2;
        delta= (abs(t1)+abs(t2)+abs(t3)+abs(t4))/(abs(w1)+abs(w2));
        t3= t1; t4= t2;
      ENDDO
      w2= x*w2; w1p= x*w1p;
      wax= wa0*w1+wa0p*w2; waxp=wa0*w1p+wa0p*w2p;
      waxm= wa0*w1-wa0p*w2;  waxmp= -wa0*w1p+wa0p*w2p;
    ENDIF
    END SUBROUTINE waxsma

    SUBROUTINE taylor(a,xi,wax,waxp)
    ! --------------------------------------------------
    !  Integration of the ODE for the W functions using
    !  the Taylor method
    ! --------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   xi,    argument of the functions
    ! Outputs:
    !   wax,  W(a,x)
    !   waxp, W'(a,x)
    ! --------------------------------------------------
    USE someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: xi
    REAL(r8), INTENT(OUT) :: wax
    REAL(r8), INTENT(OUT) :: waxp
    REAL(r8) x, xx, wa0, wa0p, waxm, waxmp, y0, y1, h, s, sd
    INTEGER ierr,k,m,n,j
    k=-1 
    IF (xi>0) k=1
    x=abs(xi)
    IF ((k.eq.1).AND.(a>0.0_r8)) THEN
      m=-1
    ELSE
      m=1
    ENDIF
    IF (m.eq.1) THEN
      xx=0.0_r8
      wa0= waxatx0(a); wa0p= 0.5_r8/(wa0);
      y0=wa0; y1=wa0p;
      IF (k>0) y1=-y1
    ELSE
      xx=50.0_r8
      CALL wammx(a,xx,0,wa0,wa0p,waxm,waxmp,ierr)
      y0=wa0; y1=wa0p;
    ENDIF
    h=m*0.0625_r8
    n=int(abs((xx-x)/h));
    DO j=0,n-1
      CALL step(a,xx,h,y0,y1,s,sd) 
      xx=xx+h
      y0=s
      y1=sd
    ENDDO
    IF (m.eq.1) THEN
      IF (xx<x) THEN 
        h=x-xx;       
        CALL step(a,xx,h,y0,y1,s,sd)
      ENDIF
    ELSE
      IF (xx>x) THEN 
        h=x-xx;       
        CALL step(a,xx,h,y0,y1,s,sd)
      ENDIF
    ENDIF
    wax=s; waxp=sd
    IF (xi<0) waxp=-waxp
    END SUBROUTINE taylor

    SUBROUTINE step(a,x,h,y0,y1,s,sd)
    ! Auxiliary routine for taylor
    USE someconstants
    IMPLICIT NONE
    REAL(r8) y0, y1, y2, y3, y4, a, x, di, xin, xifac, &
             h2, h3, h2d2, lt, h, s, coun, cound, sd, ltd
    INTEGER i,imax
    imax=20        
    xin=x
    h2=h*h
    h2d2=h2*0.5_r8
    h3=h2*h
    xifac=xin*xin*0.25_r8-a
    y2=-xifac*y0
    y3=-xifac*y1-x*0.5_r8*y0
    s=y0+y1*h+y2*h2d2+y3*h3*onesix;
    sd=y1+y2*h+y3*h2d2
    lt=h3*onesix;
    ltd=h2*0.5_r8;
    DO i=0,imax
      di=float(i)
      coun=di+4.0_r8
      cound=di+3.0_r8
      y4=xifac*y2+(di+2.0_r8)*xin*0.5_r8*y1+&
             (di+2.0_r8)*(di+1.0_r8)*0.25_r8*y0
      y4=-y4
      lt=lt*h/coun
      ltd=ltd*h/cound
      s=s+lt*y4
      sd=sd+ltd*y4
      y0=y1
      y1=y2
      y2=y3
      y3=y4
    ENDDO
    END SUBROUTINE step
 END MODULE Parabolw 

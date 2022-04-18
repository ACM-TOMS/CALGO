    MODULE constants_639
      USE prec, ONLY : wp
      REAL (wp), PARAMETER :: one = 1.0_wp, two = 2.0_wp, zero = 0.0_wp, &
        four = 4.0_wp, sixteen = 16.0_wp, half = 0.5_wp, &
        three_eights = 0.375_wp, six = 6.0_wp, one_twenty = 120.0_wp, &
        ten_to_minus_4 = 10.0E-4_wp, ten_to_minus_8 = 10.0E-8_wp, &
        quarter = 0.25_wp, hundred = 100.0_wp
    END MODULE constants_639

    PROGRAM driver_639
! .. Use Statements ..
      USE calgo_639, ONLY : oscint
      USE prec, ONLY : wp
      USE oscsub, ONLY : alpha, const6, ebase, halfpi, nb, nix, numfun, pi
      USE constants_639, ONLY : four, half, one, two, zero
! ..
! .. Local Scalars ..
      REAL (wp) :: azero, cin, eps, period, result, rfirst
      INTEGER :: alloc, endv, i, jp, jpmk, jpr, kpr, ndim1, ndim2, nquad
      LOGICAL :: comments
      CHARACTER (50) :: form
      CHARACTER (80) :: gname
      CHARACTER (80) :: hname, notes, error_mess
! ..
! .. Local Arrays ..
      REAL (wp), ALLOCATABLE :: work(:,:)
      INTEGER :: istate(6)
      CHARACTER (80) :: com(3)
! ..
! .. External Functions ..
      REAL (wp), EXTERNAL :: exact, gper, hfun
! ..
! .. External Subroutines ..
      EXTERNAL g5and9
! ..
! .. Intrinsic Functions ..
      INTRINSIC atan, exp, log, sqrt, min



!      CONSTANTS, ETC.USED IN SUBROUTINES.
      ebase = exp(one)
      pi = four*atan(one)
      halfpi = pi*half
      cin = one + sqrt(two)
      const6 = one/log(cin)
      alpha = zero
      nb = 2
10    READ (*,*,end=20) numfun
      READ (*,'(a)') gname
      READ (*,'(a)') hname
      READ (*,'(a)') notes
      READ (*,*) azero, rfirst, eps, period
      READ (*,*) ndim1, ndim2, nquad
      READ (*,*) nix(1:3)
      READ (*,*) comments
      IF (comments) THEN
        READ (*,'(a)') com(1:3)
      END IF

      ALLOCATE (work(ndim1,ndim2),STAT=alloc)
      IF (alloc/=0) THEN
        WRITE (*,'("Allocation fails in driver_639")')
        STOP
      END IF

      WRITE (*,90000) numfun, gname, hname
      WRITE (*,90040) notes
      WRITE (*,90050) azero, period, rfirst, eps, nquad
      WRITE (*,90060) ndim1, ndim2
      WRITE (*,90070)
      IF (nix(1)>0 .OR. nix(3)>0) WRITE (*,90080)

      CALL oscint(azero,period,rfirst,eps,nquad,work,result,istate,g5and9, &
        hfun,gper)

      WRITE (*,90090)
      IF (istate(1)>=0) THEN
        WRITE (*,90110) result, exact(), exact() - result
      END IF

      WRITE (*,90120) (istate(i),i=1,5)
      WRITE (*,90140) istate(6)

      WRITE (*,90020)
      IF (istate(1)>=0) THEN
        WRITE (*,90100)
      ELSE
        WRITE (*,'(1x,''ERROR: '',a)') error_mess(istate(1))
      END IF

      IF (comments) THEN
        WRITE (*,90130) com(1), com(2), com(3)
      END IF

      WRITE (*,90030)

      IF (nix(2)>0) THEN
        WRITE (*,90010)
        jp = istate(3)
        jpmk = jp - istate(4)
        DO jpr = 1, min(istate(3),nix(2))
          endv = jp - jpr
          IF (jpr>jpmk) THEN
            endv = endv + 1
          END IF

          endv = min(endv,ndim2)
          WRITE (form,fmt='("(1x,i3,",i3,"(7d10.2,/4X))")') endv
          WRITE (*,fmt=form) jpr, (work(jpr,kpr),kpr=1,endv)
        END DO

      END IF

      DEALLOCATE (work)
      GO TO 10

90000 FORMAT (/' NUMFUN = ',1X,I1,10X,'F(X) = G(X).H(X).'/22X,'G(X) = ',A/, &
        22X,'H(X) = ',A)
90010 FORMAT (/25X,'Finite Average Table'/25X,'--------------------')
90020 FORMAT ( &
        ' ************************** Exit status **************************'/)
90030 FORMAT (/ &
        ' *****************************************************************')
90040 FORMAT (/' INPUT PARAMETERS-'//1X,A80/)
90050 FORMAT (4X,'AZERO',12X,'PERIOD',12X,'RFIRST',12X,'EPS',12X,'NQUAD', &
        12X/D13.6,5X,D12.6,6X,D12.6,6X,D12.6,6X,I3/)
90060 FORMAT (/4X,'NDIM1',10X,'NDIM2'/5X,I3,11X,I3/)
90070 FORMAT ('     ----------------------')
90080 FORMAT (/'    Code check output from HFUN, GPER, QRULE'/)
90090 FORMAT (/' PRINCIPAL OUTPUT PARAMETERS-'//)
90100 FORMAT (1X,' Successful exit from OSCINT'/)
90110 FORMAT (5X,'RESULT',20X,'Exact',21X,'Difference'/D17.8,10X,D16.8,10X, &
        D16.8,10X//)
90120 FORMAT (4X,'ISTATE(1)',5X,'ISTATE(2)',5X,'ISTATE(3)',5X,'ISTATE(4)',8X, &
        'ISTATE(5)'/5X,'(IFAIL)',7X,'(LSIGCH)',8X,'(N)',7X,'(Column No.)',6X, &
        '(Row No.)'/I9,11X,I3,11X,I3,11X,I3,13X,I4/)
90130 FORMAT (3(/1X,A))
90140 FORMAT (4X,' No. of function values (ISTATE(6)) = ',I6/)

20  END PROGRAM driver_639

    CHARACTER (len=80) FUNCTION error_mess(ifail)
      USE calgo_639, ONLY : err_gauss_routine, err_max_intervals, &
        err_sign_change, err_input_period, err_gper_function, err_input_params
      INTEGER, INTENT (IN) :: ifail

      SELECT CASE (ifail)
      CASE (err_gauss_routine)
        error_mess = 'Error reported from user supplied function GAUSS'
      CASE (err_max_intervals)
        error_mess = 'Maximum number of intervals exceeded'
      CASE (err_sign_change)
        error_mess = &
          'Normal sign change pattern is violated after the grace period'
      CASE (err_input_period)
        error_mess = 'Supplied value of PERIOD <= 10**(-5)'
      CASE (err_gper_function)
        error_mess = &
          'Routine has detected that GPER is not periodic with given period'
      CASE (err_input_params)
        error_mess = 'Input parameter error (NQUAD = 0,1 or WORK &
          &array not at least (1,1))'
      CASE DEFAULT
        error_mess = 'ILLEGAL ERROR RETURN -- PLEASE REPORT'
      END SELECT

    END FUNCTION error_mess


    FUNCTION exact()
      USE prec, ONLY : wp
      USE oscsub, ONLY : alpha, ebase, halfpi, numfun, pi
      USE constants_639, ONLY : four, hundred, one, quarter, sixteen, two, &
        zero

!      EXACT IS THE VALUE OF THE INTEGRAL OF F(X) OVER (0,INFINITY).
!      HERE, OF COURSE, F(X) = H(X).G(X) GIVEN BY HFUN AND GPER BELOW.

! .. Function Return Value ..
      REAL (wp) :: exact
! ..
! .. Intrinsic Functions ..
      INTRINSIC cosh, exp

      SELECT CASE (numfun)
      CASE (1)
        exact = halfpi
      CASE (2)
        exact = pi/(two*ebase)
      CASE (3)
        exact = exp(-(one+alpha)*halfpi)/cosh(halfpi)
      CASE (4)
        exact = one
      CASE (5)
        exact = one + hundred*(one-exp(-pi*pi*quarter))
      CASE (6)
        exact = pi/sixteen*(ebase**(-4)+four*ebase**(-2))
      CASE DEFAULT
        exact = zero
      END SELECT
    END FUNCTION exact

    FUNCTION hfun(x)
      USE prec, ONLY : wp
      USE oscsub, ONLY : alpha, const6, ebase, halfpi, nb, nix, numfun, pi
      USE constants_639, ONLY : hundred, one, one_twenty, six, ten_to_minus_4, &
        ten_to_minus_8, two, zero
! .. Function Return Value ..
      REAL (wp) :: hfun
! ..
! .. Scalar Arguments ..
      REAL (wp) :: x
! ..
! .. Local Scalars ..
      INTEGER :: ncalc
! ..
! .. Local Arrays ..
      REAL (wp) :: b(nb)
! ..
! .. External Subroutines ..
      EXTERNAL rjbesl
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, exp, sin


!      ALPHA - FRACTIONAL PART OF ORDER FOR WHICH BESSEL FUNCTION IS TO
!              BE CALCULATED, 0 <= ALPHA < 1
!      NB      NUMBER OF BESSEL FUNCTIONS IN SEQUENCE.


!      IF (NUMFUN .EQ. 1) F(X) = SIN(X)/X
!      IF (NUMFUN .EQ. 2) F(X) = COS(X)/(1+X**2)
!      IF (NUMFUN .EQ. 3) F(X) = EXP(-X*SINH(PI/2)*J1(X)
!      IF (NUMFUN .EQ. 4) F(X) = QUADPACK 5.2.3. PAGE118
!      IF (NUMFUN .EQ. 5) F(X) = (1 + 100*EXP(-X**2/PI**2))*J1(X)
!      IF (NUMFUN .EQ. 6) F(X) = (COS**4(X) - 3/8)/(1 + X**2)

!      WHEN NUMFUN IS 1,3,4 OR 5, H(X) = F(X).
!      WHEN NUMFUN IS 2 OR 6,     H(X) = 1/(1 + X**2)
      SELECT CASE (numfun)
      CASE (1)
        IF (abs(x)>ten_to_minus_4) THEN
          hfun = sin(x)/x
        ELSE
          hfun = one - x**2/six + x**4/one_twenty
        END IF
      CASE (2)
        hfun = one/(one+(x**two))
      CASE (3)
        CALL rjbesl(x,alpha,nb,b,ncalc)
        hfun = ebase**(-x*sinh(halfpi))*b(2)
      CASE (4)
        CALL rjbesl(x,alpha,nb,b,ncalc)
        hfun = const6*b(1)
        IF (x>=ten_to_minus_8) THEN
          hfun = ((one-exp(-x))/x)*hfun
        END IF
      CASE (5)
        CALL rjbesl(x,alpha,nb,b,ncalc)
        hfun = (one+hundred*exp(-(x**two)/(pi**two)))*b(2)
      CASE (6)
        hfun = one/(one+(x**two))
      CASE DEFAULT
        hfun = zero
      END SELECT

      IF (nix(3)>0) THEN
        nix(3) = nix(3) - 1
        WRITE (*,'(''   X = '',D16.8,''   HFUN = '',D16.8)') x, hfun
      END IF

    END FUNCTION hfun


    FUNCTION gper(x)
      USE prec, ONLY : wp
      USE oscsub, ONLY : nix, numfun
      USE constants_639, ONLY : four, one, three_eights, zero
! .. Function Return Value ..
      REAL (wp) :: gper
! ..
! .. Scalar Arguments ..
      REAL (wp) :: x

      SELECT CASE (numfun)
      CASE (1,3,4,5)
        gper = one
      CASE (2)
        gper = cos(x)
      CASE (6)
        gper = (cos(x))**four - three_eights
      CASE DEFAULT
        gper = zero
      END SELECT
      IF (nix(3)/=0) THEN
        WRITE (*,'(''   X = '',D16.8,''   GPER = '',D16.8)') x, gper
      END IF
      RETURN

    END FUNCTION gper



    SUBROUTINE g5and9(nquad,weight,abscis,ierr)
      USE prec, ONLY : wp

!      THIS IS AN EXTRACT FROM A QUADRATURE ROUTINE CONSTRUCTED ONLY FOR
!      USE IN A DRIVER WHICH ILLUSTRATES OSCINT. A NAG LIBRARY
!      SUBSCRIBER MAY REPLACE THIS BY D01BCF FOR GENERAL USE.

!      THE FOLLOWING DATA ARE WEIGHTS AND ABSCISSAS OF THE FIVE POINT
!      AND THE NINE POINT GAUSS-LEGENDRE QUADRATURE RULES RESPECTIVELY.
!      NORMALISED TO THE INTERVAL (-1,1).
! .. Scalar Arguments ..
      INTEGER :: ierr, nquad
! ..
! .. Array Arguments ..
      REAL (wp) :: abscis(nquad), weight(nquad)
! ..
! .. Local Arrays ..
      REAL (wp) :: absc5(5), absc9(9), wt5(5), wt9(9)
! ..
! .. Data Statements ..
      DATA wt5/.2369268850561891_wp, .4786286704993665_wp, &
        .5688888888888889_wp, .4786286704993665_wp, .2369268850561891_wp/
      DATA absc5/ -.9061798459386640_wp, -.5384693101056831_wp, 0.0_wp, &
        .5384693101056830_wp, .9061798459386640_wp/
      DATA wt9/.0812743883615745_wp, .1806481606948574_wp, &
        .2606106964029355_wp, .3123470770400029_wp, .3302393550012598_wp, &
        .3123470770400029_wp, .2606106964029356_wp, .1806481606948575_wp, &
        .8127438836157467D-01/
      DATA absc9/ -.9681602395076261_wp, -.8360311073266358_wp, &
        -.6133714327005904_wp, -.3242534234038090_wp, 0.0_wp, &
        .3242534234038087_wp, .6133714327005902_wp, .8360311073266357_wp, &
        .9681602395076260_wp/
! ..

      SELECT CASE (nquad)
      CASE (5)
        abscis(1:5) = absc5
        weight(1:5) = wt5
        ierr = 0
      CASE (9)
        abscis(1:9) = absc9
        weight(1:9) = wt9
        ierr = 0
      CASE DEFAULT
        ierr = 59
      END SELECT

    END SUBROUTINE g5and9

    MODULE calgo_639
!!!!!!!!!!!!!
!!! Documentation
!!
! The routine OSCINT attempts to compute an approximation, RESULT, to the 
! integral of a user defined overall integrand F(X) = HFUN(X)*GFUN(X) over
! a semi-infinite interval (AZERO,INFINITY). The overall integrand 
! function should be ultimately oscillating in sign and periodic
! with user-defined period of PERIOD.

! The routine considers a sequence of intervals. All the intervals,
! except the first, are of length 0.5*PERIOD. A quadrature rule is used
! to approximate the integral over each interval. These approximations
! are stored in the array QLIST(J), J=1,2,3,... and the values ultimately
! oscillate in sign. The routine then uses a series acceleration technique
! based on the euler transformation to sum this sequence.

! Input parameters to OSCINT

! AZERO  - REAL, INTENT(IN)
!         defines the lower limit of the integration

! PERIOD - REAL, INTENT(IN)
!         defines the least positive period or ultimate period of
!         the overall integrand function

! RFIRST - REAL, INTENT(IN)
!         defines the right-hand endpoint of the first interval.
!         This allows the user to locate his subdivision.  For cautious
!         running, arrange RFIRST to coincide with an "ultimate zero".  For
!         slightly more adventurous but less reliable running, arrange RFIRST
!         to coincide with an "ultimate" peak.  Otherwise, set RFIRST < AZERO,
!         in which case, OSCINT uses AZERO instead of RFIRST.

! EPS    - REAL, INTENT(IN)
!         defines the required accuracy.

! NQUAD  - REAL, INTENT(IN)
!         defines the number of abscissa points to be used by the quadrature
!         rule in each interval. If NQUAD is set > 1 the trapezoidal rule
!         is used otherwise an ABS(NQUAD) rule specified by the user
!         provided function GAUSS (see below) is used.

! WORK   - REAL(:,:), INTENT(OUT)
!         used to store the finite average table. This table may provide useful
!         information in the event of the routine failing to produce a result
!         to the desired accuracy.

!         The size of the array controls the maximum amount of work that the
!         routine will perform before abandoning the computation.

!         NDIM1 = SIZE(WORK,1) defines the maximum number of intervals to be


!         NDIM2 = SIZE(WORK,2) defines the number of finite averages to be
!         computed. A suggested value for NDIM2 is 15; generally 4<= NDIM2 <= 20.
!         There is no point in having NDIM2 > NDIM1.

! RESULT - REAL, INTENT(OUT)
!         used to return the final approximation to the user defined integral.
!         The user should check that the computation was successful by
!         ensuring that ISTATE(1) is also returned as zero.

! ISTATE - INTEGER(6), INTENT(OUT)
!         vector of integers that provides the status of the given result

!         ISTATE(1): failure indicator
!                    0 : successful 
!                    >0 : apparently successful but unsatisfactory features of
!                 possible interest to a sophisticated user have been
!                 detected. See notes on termination below for more details.

!                    <0 :
!            err_max_intervals : maximum number of intervals exceeded;
!                 could increase  leading dimension of WORK array and rerun.
!            err_sign_change : the normal sign change pattern is violated 
!                 after the grace period. See notes below for details.
!            err_input_period : supplied value of PERIOD <= 10**(-5).
!            err_gauss_routine : failure indicator set by user supplied 
!                 GAUSS routine.
!            err_gper_function : routine has detected that GPER is not
!                 periodic.
!            err_input_params : input parameter error (NQUAD = 0,1 or 
!                 NDIM1, NDIM2 defined above < 1).
!            err_alloc: ALLOCATE failed to provide requested array
!                       space.
!            err_dealloc: DEALLOCATE failed to free array space.
!                         In this case RESULT may contain useful
!                         information.

!          ISTATE(2): (LSIGCH) indicates last interval in which the 
!                     sign of the integral coincided with that of the
!                     integral over the previous interval

!          ISTATE(3): (NOW) the number of finite integrals, QLIST(Q),
!                     evaluated in the calculation.
!          ISTATE(4): (NCOL) the column of the finite average table,
!                                (WORK), on which the  result is based.
!          ISTATE(5): (NROW) the row of the finite average table,
!                                (WORK), on which the  result is based.
!          ISTATE(6): (NCOUNT) the number of calls to function HFUN.


! GAUSS  - SUBROUTINE (N, WEIGHT, ABSCIS, IFAIL)
!        name of a user defined subroutine which provides weights and abscissa
!        points. This routine is not called when NQUAD > 0. When NQUAD is negative
!        it is called with N = ABS(NQUAD) and should return a set of weights and 
!        abscissa points suitable for integration over the interval [-1,1] in
!        the real arrays WEIGHT(1:N) and ABSCIS(1:N) respectively. 

!        The argument IFAIL may be used as an error return. If GAUSS returns
!        with IFAIL/=0 the OSCINT aborts with ISTATE(1) = err_gauss_routine.

!        For users with the NAG library D01BCF may be used to compute the relevant
!        data.

! HFUN    - REAL FUNCTION (X)
! GFUN    - REAL FUNCTION (X)

!         the two user defined functions that define the overall integrand
!        to be of the form F(X)=HFUN(X)*GFUN(X). GPER must be periodic with period 
!        coinciding with the input parameter PERIOD defined above.

!        It is always possible to choose HFUN to be the integrand function F(X) and
!        code GPER to return the value 1.0. However, when F(X) has a periodic factor,
!        J(X), repetitive evaluation of J(X) in each interval may be avoided by
!        setting GPER(X)=J(X) and HFUN(X)=H(X). The routine checks that GPER(X) is
!        periodic by making some evaluations in the third and fourth intervals.
!        If it discovers that GPER(X) is not periodic, OSCINT terminates with
!        ISTATE(1)=err_gper_function.

! Termination: 

!    In what follows, the "normal sign change pattern" is one in which
! successive values of QLIST(Q), the value of the approximation to
! the integral over the Qth interval, oscillate in sign.  The "grace
! period" is Q<=10 when the normal sign change pattern is not
! insisted on.  LSIGCH is the highest value of Q for which
! QLIST(Q)*QLIST(Q-1) is positive.  Termination comes about
! after calculating QLIST(NOW), when either

! 1.  an approximation of sufficient accuracy is currently available
!     (the routine sets ISTATE(1) = MAX(0,4-(NOW-LSIGCH))
! or

! 2.  the normal sign change pattern is violated after the grace
!      period.   (The routine sets ISTATE(1) = err_sign_change)
! or

! 3.  the user set limit, NDIM1, of intervals have been calculated
!     i.e.  NOW = NMAX.  (The routine sets ISTATE(1) = err_max_intervals)

!    The routine checks 1., then 2., and then 3.  Termination
! under 1. may occur without any normal sign change pattern
! emerging.  This may be due to misuse of the routine, but the
! problem is small enough to be correctly handled.  In this case,
! ISTATE(1) = MAX(0,4-(NOW-LSIGCH)) may be a small positive  integer.


! .. Use Statements ..
      USE prec, ONLY : wp
! ..
! .. Local Arrays ..
      REAL (wp), ALLOCATABLE :: abscis(:), qlist(:), savper(:), weight(:)
! ..
! .. Parameters ..
      REAL (wp), PARAMETER, PRIVATE :: half = 0.5_wp
      REAL (wp), PARAMETER, PRIVATE :: one = 1.0_wp
      REAL (wp), PARAMETER, PRIVATE :: permin = 1.0E-5_wp
      REAL (wp), PARAMETER, PRIVATE :: two = 2.0_wp
      REAL (wp), PARAMETER, PRIVATE :: zero = 0.0_wp
      INTEGER, PARAMETER :: err_alloc = -7000, err_dealloc = -8000, &
        err_gauss_routine = -4000, err_gper_function = -5000, &
        err_input_params = -6000, err_input_period = -3000, &
        err_max_intervals = -100, err_sign_change = -200
! ..
    CONTAINS

      SUBROUTINE oscint(azero,period,rfirst,eps,nquad,work,result,istate, &
          gauss,hfun,gper)
! .. Scalar Arguments ..
        REAL (wp), INTENT (IN) :: azero, eps, period, rfirst
        REAL (wp), INTENT (OUT) :: result
        INTEGER, INTENT (IN) :: nquad
! ..
! .. Array Arguments ..
        REAL (wp), INTENT (OUT) :: work(:,:)
        INTEGER, INTENT (OUT) :: istate(6)
! ..
! .. Function Arguments ..
        REAL (wp), EXTERNAL :: gper, hfun
! ..
! .. Subroutine Arguments ..
        EXTERNAL gauss
! ..
! .. Local Scalars ..
        REAL (wp) :: hasper, wmin
        INTEGER :: jp, k, mj1, ndim1, ndim2, nml, now, nowjp, nrow, status
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, dabs, max, min, size, sum
! ..
        ndim1 = size(work,1)
        ndim2 = size(work,2)
        hasper = half*period
        wmin = one
        istate = 0
        result = zero
        IF (nquad==0 .OR. nquad==1 .OR. ndim1<1 .OR. ndim2<1) THEN
          istate(1) = err_input_params
          RETURN
        ELSE IF (period<permin) THEN
          istate(1) = err_input_period
          RETURN
        END IF
        ALLOCATE (qlist(ndim1),savper(2*abs(nquad)),weight(max(1, &
          -nquad)),abscis(max(1,-nquad)),STAT=status)
        IF (status/=0) THEN
          istate(1) = err_alloc
          RETURN
        END IF

! loop to construct table
        nrow = 0
        now = 0
JPLOOP: DO jp = 1, ndim1

          DO k = 1, min(jp,ndim2)

            IF (k==1) THEN
              work(jp,k) = qrule(jp-1,azero,hasper,rfirst,nquad,gauss,hfun, &
                gper,istate)
              IF (istate(1)/=0) THEN
                EXIT JPLOOP
              END IF
            ELSE
              mj1 = jp - k + 1
              work(mj1,k) = (work(mj1,k-1)+work(jp-k+2,k-1))*half
              IF (dabs(work(mj1,k))<wmin) THEN
                wmin = dabs(work(mj1,k))
                now = k
                nrow = mj1
                nowjp = jp
              END IF
              IF (jp/=k) THEN
                IF ((abs(work(mj1,k))<eps) .AND. (abs(work(jp-k,k))<eps)) THEN
                  now = k
                  nrow = mj1
                  nowjp = jp
                  EXIT JPLOOP
                END IF
              END IF
            END IF
          END DO
        END DO JPLOOP

        istate(4) = now
        istate(5) = nrow
        IF (istate(1)<0) THEN
          istate(3) = jp - 1
        ELSE
          IF (jp>ndim1) THEN
            istate(1) = err_max_intervals
            istate(3) = ndim1
          ELSE
            result = sum(qlist(1:nowjp-now)) + half*sum(work(nrow,1:now-1)) + &
              work(nrow,now)
            istate(3) = jp
            nml = istate(3) - istate(2)
            IF (nml<4 .AND. istate(1)==0) THEN
              istate(1) = max(0,4-nml)
            END IF
          END IF
        END IF

        DEALLOCATE (qlist,savper,abscis,weight,STAT=status)
        IF (status/=0) THEN
          istate(1) = err_dealloc
          RETURN
        END IF

        RETURN
      END SUBROUTINE oscint

      FUNCTION qrule(j,azero,hasper,rfirst,nquad,gauss,hfun,gper,istate)
! .. Function Return Value ..
        REAL (wp) :: qrule
! ..
! .. Scalar Arguments ..
        REAL (wp), INTENT (IN) :: azero, hasper, rfirst
        INTEGER, INTENT (IN) :: j, nquad
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: istate(6)
! ..
! .. Function Arguments ..
        REAL (wp), EXTERNAL :: gper, hfun
! ..
! .. Subroutine Arguments ..
        EXTERNAL gauss
! ..
! .. Local Scalars ..
        REAL (wp) :: a, diff, fun, wsum, wt, xi, y
        REAL (wp), SAVE :: b, gmax, tendpt
        INTEGER :: i, ielm, ifail, npts
        INTEGER, SAVE :: index
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, max, mod
! ..
!!!
!          This routine evaluates the integral of the function,
!     HFUN(X)*GPER(X) over the interval (A,B) where:
!          generally (J>0)     A = previous value of B
!                              B = A + HASPER
!                    (J=0)     A = AZERO
!                              B = RFIRST when RFIRST > AZERO
!                              B = AZERO + HASPER when RFIRST <= AZERO
! Input parameters:
!        J - INTEGER, INTENT(IN)
!            defines which term, QLIST(J), of the series is being evaluated.
!   HASPER - REAL, INTENT(IN)
!            defines half the period.
!      other input and output parameters are identical to those in  OSCINT
!      (described above).
!!!
        npts = abs(nquad)

        IF (j==0) THEN
          index = 0
! call GAUSS if NQUAD < 0
          IF (nquad<0) THEN
            ifail = 0
            CALL gauss(npts,weight,abscis,ifail)
            IF (ifail/=0) THEN
              istate(1) = err_gauss_routine
              qrule = zero
              RETURN
            END IF
          END IF
        END IF

        qrule = zero
!  set interval endpoints, A and B
        IF (j/=0) THEN
          a = b
        ELSE
          a = azero
          gmax = zero
        END IF

        IF (rfirst<=azero .OR. j/=0) THEN
          b = a + hasper
        ELSE
          b = rfirst
        END IF

! loop on abscissas for interval #J - starts to calculate QRULE
        DO i = 1, npts
          xi = i
! calculate abscissa Y
          IF (nquad<0) THEN
            y = (b-a)*abscis(i)*half + (b+a)*half
            wt = weight(i)
          ELSE
            y = ((xi-1)*b+a*(npts-i))/(npts-1)
            wt = one
          END IF

! IELM is location in SAVPER for GPER function values
! J=0: ignores. J=1,2: where to put value. J>2: where to get
! value from.
          ielm = npts*mod(j-1,2) + i

          IF (j>=3 .AND. i==1 .AND. nquad>0) THEN
            GO TO 10
          END IF
          SELECT CASE (j)
          CASE (0)
            fun = hfun(y)*gper(y)

          CASE (1:2)
            savper(ielm) = gper(y)
            IF (i==1 .AND. nquad>0) THEN
              gmax = abs(savper(ielm))
              GO TO 10
            END IF
! check for constant function
            IF (savper(ielm)==savper(1)) THEN
              index = index + 1
            END IF
            gmax = max(gmax,abs(savper(ielm)))
            fun = hfun(y)*savper(ielm)

          CASE (3:4)
!  check that GPER is periodic with PERIOD=PERIOD
            diff = abs(savper(ielm)-gper(y))

            IF (diff<gmax*permin) THEN
              fun = hfun(y)*savper(ielm)
            ELSE
              istate(1) = err_gper_function
            END IF

          CASE (5:)
            IF (index==2*npts) THEN
              fun = hfun(y)
            ELSE
              fun = hfun(y)*savper(ielm)
            END IF
          END SELECT

          istate(6) = istate(6) + 1

10        CONTINUE
          IF (nquad>0) THEN
            IF (i==1) THEN
              IF (j==0) THEN
                fun = fun*half
              END IF
              IF (j>0) THEN
                fun = tendpt
              END IF
            END IF

            IF (i==npts) THEN
              fun = fun*half
              tendpt = fun
            END IF
            qrule = qrule + fun
          ELSE
            qrule = qrule + wt*fun
          END IF

        END DO
! loop on abscissa for interval #J ends
        IF (index==2*npts) THEN
          qrule = qrule*savper(1)
        END IF

        IF (nquad>0) THEN
          wsum = npts - 1
        ELSE
          wsum = two
        END IF
        qrule = qrule*(b-a)/wsum
        qlist(j+1) = qrule

        IF (j>0) THEN
          IF (qrule*qlist(j)>0) THEN
            istate(2) = j
          END IF
        END IF
        IF (istate(2)>9) THEN
          istate(1) = err_sign_change
        END IF

      END FUNCTION qrule

    END MODULE calgo_639

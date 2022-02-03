      module MARCOMOD
      private
      public:: sl01f,sl02f,sl02fm,sl03f,sl03fm,sl04f,sl04fm,sl05f,
     +         sl05fm,sl06fm,sl07f,sl07fm

      contains

C SL01F & SL02F complete text (including NAG routines)
      SUBROUTINE SL01F(A,B,EMIN,EMAX,KMIN,KMAX,COEFFN,SETUP,AINFO,BINFO,
     &                 SYM,WK,IWK,KNOB,TOL,IFAIL)
C This routine counts the number of eigenvalues in a given range
C [emin,emax]. If emax lies in the continuous spectrum and emin
C does not, then an estimate will be made of the infimum of the
C continuous spectrum and this will be returned in lieu of emax.
C The eigenvalues are counted by returning kmin and kmax, the
C indices of the lowest and highest index eigenvalues lying in
C the interval.
C
C Brief specification:
C
C a,b: real. On entry, a and b specify the actual interval (a,b)on which
C            the (regular or singular) problem is posed. On exit, they
C            specify the (possibly shorter) interval on which it was solved.
C emin, emax: real. On entry, [emin,emax] is the interval in which the
C            code must count the eigenvalues. If this interval intersects
C            the continuous spectrum then provided emin does not lie in the
C            continuous spectrum the routine will locate a number elam such
C            that the infimum of the continuous spectrum is believed
C            to occur in (elam,elam+tol). The routine will then use the
C            interval [emin,elam] in place of [emin,emax]: it will count
C            the eigenvalues in this interval, and return elam in place
C            of emax. If emin is in the continuous spectrum or within a
C            distance tol of the continuous spectrum, the routine will
C            return with kmin = kmax = -10.
C coeffn: subroutine supplied by the user. Specification:
C            SUBROUTINE coeffn(x,p,q,w)
C            DOUBLE PRECISION x,p,q,w
C      C     x: On entry, a value of the independent variable in (a,b).
C               Must be unchanged on exit.
C      C     p,q,w: On exit these must specify the values of the
C      C        coefficients p(x), q(x), w(x).
C setup: subroutine supplied by the user. Specification:
C            SUBROUTINE setup(y,pdy,elam,iend)
C            INTEGER iend
C            DOUBLE PRECISION y,pdy,elam
C      C     iend: Specifies the end of the interval at which the b.c.'s
C      C           are to be applied. If iend = 0, then B.C.'s are required
C      C           at the left hand end, while iend = 1 requests b.c.'s at
C      C           the right-hand end. Must not be changed on exit.
C      C     y,pdy: On exit, these must specify values of y(.) and py'(.)
C      C           which are consistent with the B.C.'s at the appropriate
C      C           end.
C      C     elam: On entry, elam specifies the current value of the
C      C           eigenparameter. It is included since in many problems
C      C           the boundary conditions are dependent on this parameter.
C kmin, kmax: integers. On exit, these specify respectively the highest and
C            lowest indices of eigenvalues in [emin,elam], where elam is emax
C            if emax is not in the continuous spectrum, and is an estimate of
C            the infimum of the continuous spectrum otherwise. (see
C            description of emax above). If emin is within a distance tol of
C            the inf of the continuous spectrum then kmin and kmax will both
C            be set to -10 on exit. If the routine detects no eigenvalues in
C            the range [emin,emax] then it will return with kmin > kmax.
C knob: integer. The maximum number of endpoint shifts which the routine
C            is allowed to make in the course of the computation. Default
C            value of 60 should only be increased in special circumstances.
C            To select default value set knob < 10 on entry. The minimum
C            value of knob within the routine is therefore 10. On exit,
C            knob will have either the default value if selected, or the
C            value assigned by the user on entry.
C ainfo: character*1. On entry ainfo must be set as follows:
C            If x=a is a finite regular endpoint set ainfo = 'R' or ainfo = 'r'.
C            If x=a is a finite singular endpoint set ainfo ='S' or ainfo = 's'
C            If a=-infinity set ainfo = 'I' or ainfo = 'i'.
C binfo: character*1. On entry binfo must be set as follows:
C            If x=b is a finite regular endpoint set binfo = 'R' or binfo = 'r'.
C            If x=b is a finite singular endpoint set binfo ='S' or binfo = 's'
C            If b=+infinity set binfo = 'I' or ainfo = 'i'.
C sym: logical. If the problem is symmetric about the midpoint (a+b)/2
C            (taken as 0 if the endpoints are infinite) then set sym = .true.,
C             otherwise sym should be .false.
C wk: real array of dimension precisely (0:iwk,1:4). Used as workspace.
C iwk: integer. On entry, iwk must have a value which satisfies
C                        iwk > 10 + knob,
C            where knob has the value defined by the user if this is > 9,
C            and is otherwise equal to 60. Unchanged on exit.
C tol: real. This parameter controls the computation in various ways. It is
C            used in the location of the infimum of the continuous spectrum
C            if necessary (see description of emax given above) and generally
C            any eigenvalue which lies within the interval [emin,elam] (where
C            elam = min(emax,inf{Continuous Spectrum}) with at least tol
C            clearance from each end, should be counted. Unchanged on exit.
C ifail: integer. The error flag.

C     .. Parameters ..
      INTEGER IPARAM
      PARAMETER (IPARAM=1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,EMAX,EMIN,TOL
      INTEGER IFAIL,IWK,KMAX,KMIN,KNOB
      LOGICAL SYM
      CHARACTER AINFO*1,BINFO*1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WK(0:IWK,1:4)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL COEFFN,SETUP
C     ..
C     .. Local Scalars ..
      INTEGER ASTAT,BSTAT,ICOFUN,IFO,ISING,N,NREC
      CHARACTER SRNAME*6
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION PARAMS(1:IPARAM)
      CHARACTER REC(2)*80
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION X02AJF
cc      INTEGER P01ABF
cc      EXTERNAL X02AJF,P01ABF
C     ..
C     .. External Subroutines ..
cc      EXTERNAL CNTER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      SRNAME = 'SL01F'
      ICOFUN = 2
      IFO = IFAIL
      IFAIL = 0

      IF (AINFO.EQ.'R' .OR. AINFO.EQ.'r') THEN
          ASTAT = 0

      ELSE IF (AINFO.EQ.'S' .OR. AINFO.EQ.'s') THEN
          ASTAT = 1

      ELSE IF (AINFO.EQ.'I' .OR. AINFO.EQ.'i') THEN
          ASTAT = 2

      ELSE
          ASTAT = 3
      END IF

      IF (BINFO.EQ.'R' .OR. BINFO.EQ.'r') THEN
          BSTAT = 0

      ELSE IF (BINFO.EQ.'S' .OR. BINFO.EQ.'s') THEN
          BSTAT = 1

      ELSE IF (BINFO.EQ.'I' .OR. BINFO.EQ.'i') THEN
          BSTAT = 2

      ELSE
          BSTAT = 3
      END IF

C Check input for errors.

      IF (IWK.LT.10 .OR. (ASTAT.EQ.3) .OR. (BSTAT.EQ.3) .OR.
     &    (SYM.AND.ASTAT.NE.BSTAT) .OR. EMAX.LT.EMIN+TOL .OR.
     &    TOL.LT.0.D0 .OR. (ASTAT.NE.2.AND.BSTAT.NE.2.AND.B.LE.A) .OR.
     &    ABS(IFO).GT.1) THEN

          IFAIL = 1
          NREC = 1

          IF (ASTAT.EQ.3) THEN
              WRITE (REC,FMT=9010)
              GO TO 100

          END IF

          IF (BSTAT.EQ.3) THEN
              WRITE (REC,FMT=9020)
              GO TO 100

          END IF

          IF (SYM .AND. ASTAT.NE.BSTAT) THEN
              WRITE (REC,FMT=9030)
              GO TO 100

          END IF

          IF (IWK.LT.10) THEN
              WRITE (REC,FMT=9000) IWK
              GO TO 100

          END IF

          NREC = 2
          IF (EMAX.LT.EMIN+TOL) THEN
              WRITE (REC,FMT=9040) EMIN,EMAX
              GO TO 100

          END IF

          NREC = 1
          IF (TOL.LT.0.D0) THEN
              WRITE (REC,FMT=9050) TOL
              GO TO 100

          END IF

          NREC = 2
          IF (ASTAT.NE.2 .AND. BSTAT.NE.2 .AND. B.LE.A) THEN
              WRITE (REC,FMT=9060) A,B
              GO TO 100

          END IF

          IF (ABS(IFO).GT.1) THEN
              WRITE (REC,FMT=9070) IFO
              IFO = -1
              GO TO 100

          END IF

      END IF

9000  FORMAT ('** IWK must be at least 10, IWK =',i8)
9010  FORMAT ('** Invalid input value of AINFO')
9020  FORMAT ('** Invalid input value of BINFO')
9030  FORMAT ('** For symmetric problems AINFO must equal BINFO')
9040  FORMAT ('** EMAX-EMIN-TOL must be non-negative,',/,'** EMIN = ',
     &       G18.8,' EMAX = ',G18.8)
9050  FORMAT ('** TOL must be strictly positive, TOL =',g18.8)
9060  FORMAT ('** For regular problems A < B is required ',/,'** A = ',
     &       g18.8,' B = ',g18.8)
9070  FORMAT ('** On entry IFAIL must be -1, 0 or +1 ',/,
     &       '** Input value of IFAIL = ',i8)

C END of error trapping

      IF (ASTAT.NE.0 .AND. BSTAT.EQ.0) THEN
          ISING = 1

      ELSE IF (BSTAT.NE.0 .AND. ASTAT.EQ.0) THEN
          ISING = 2

      ELSE IF (ASTAT.NE.0 .AND. BSTAT.NE.0) THEN
          ISING = 3

      ELSE
          ISING = 0
      END IF

      IF (ASTAT.EQ.2 .AND. BSTAT.EQ.2) THEN
          A = -1.D0
          B = 1.D0
          ICOFUN = 1
          GO TO 10

      END IF

      IF (ASTAT.EQ.2 .AND. BSTAT.NE.2) THEN
          A = -1.D0
          B = B/ (1.D0+ABS(B))
          ICOFUN = 1
      END IF

      IF (BSTAT.EQ.2 .AND. ASTAT.NE.2) THEN
          A = A/ (1.D0+ABS(A))
          B = 1.D0
          ICOFUN = 1
      END IF

10    N = IWK
      CALL CNTER(A,B,EMIN,EMAX,KMIN,KMAX,COEFFN,SETUP,ISING,SYM,WK(0,1),
     &           WK(1,2),WK(1,3),WK(1,4),N,KNOB,TOL,PARAMS,IPARAM,
     &           ICOFUN,IFAIL)
      IF (IFAIL.NE.0) THEN
          GO TO (20,30,40,50,60,70,80,90) IFAIL

20        WRITE (REC,FMT=9080)

9080      FORMAT ('** Input parameter error')

          GO TO 100

30        WRITE (REC,FMT=9090)

9090      FORMAT ('** Integration halted: step-size too small')

          GO TO 100

40        WRITE (REC,FMT=9100)

9100      FORMAT ('** Invalid boundary conditions')

          GO TO 100

50        WRITE (REC,FMT=9110)

9110      FORMAT ('** IWK was too small for given TOL')

          GO TO 100

60        WRITE (REC,FMT=9120)

9120      FORMAT ('** Uncertainty about start of continuous spectrum')

          GO TO 100

70        WRITE (REC,FMT=9130)

9130      FORMAT (
     &       '** Too many attempts to find start of continuous spectrum'
     &           )

          GO TO 100

80        WRITE (REC,FMT=9140)

9140      FORMAT ('Tighter TOL required for computation of KMAX')

          GO TO 100

90        WRITE (REC,FMT=9150)

9150      FORMAT ('May be eigenvalues close to EMIN or EMAX')

          GO TO 100

      END IF

      IF (ICOFUN.EQ.1) THEN
          A = A/MAX((1.D0-ABS(A)),X02AJF(1.D0))
          B = B/MAX((1.D0-ABS(B)),X02AJF(1.D0))
      END IF

      RETURN

100   IFAIL = P01ABF(IFO,IFAIL,SRNAME,NREC,REC)
      RETURN

      END SUBROUTINE SL01F

C --------------------------------------------------------------------
      SUBROUTINE CNTER(A,B,EMIN,EMAX,KMIN,KMAX,COEFFN,SETUP,ISING,ISYMM,
     &                 XMESH,PP,QP,WP,N,MAXITS,TOL,PARAMS,IPARAM,ICOFUN,
     &                 IFAIL)
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,EMAX,EMIN,TOL
      INTEGER ICOFUN,IFAIL,IPARAM,ISING,KMAX,KMIN,MAXITS,N
      LOGICAL ISYMM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PARAMS(1:IPARAM),PP(1:N),QP(1:N),WP(1:N),
     &                 XMESH(0:N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL COEFFN,SETUP
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AT,ATRUNC,BT,BTRUNC,CSPEC,DLAM,ELAM,FX,H,
     &                 PI,RMAX,RTRY,TRUNC,XS,YS
      INTEGER I,IFLAG,IND,IP1,IR,MAXEVS,NEVS,NP
      LOGICAL IRELOC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(17)
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION DACC,X01AAF,X02ALF
cc      EXTERNAL DACC,X01AAF,X02ALF
C     ..
C     .. External Subroutines ..
cc      EXTERNAL C05AZF,CEKSPK,COARSE,INITIA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,INT,LOG,MAX,SQRT
C     ..
C Initial error traps:
      IF (ABS(IFAIL).GT.1) THEN
          IFAIL = -ABS(IFAIL)
          GO TO 60

      END IF

      IF (N.LT.10 .OR. TOL.LE.0.D0 .OR. TOL.GE.1.D0 .OR. B.LE.A .OR.
     &    EMIN.GE. (EMAX-TOL)) THEN
          IFAIL = 1
          GO TO 60

      END IF

C End of input checks.
C Set some constants:
      PI = X01AAF(1.D0)
      IF (MAXITS.LT.10) MAXITS = 60
      MAXEVS = (LOG(EMAX-EMIN)-LOG(TOL))/LOG(2.D0) + 20
      NP = 10
C Stage 1: set up an initial discretisation of the problem.
      IF (ISING.EQ.0) THEN
C       Omit the check for continuous spectrum, go straight to the
C       counting stage.
          ATRUNC = A
          BTRUNC = B
          AT = A
          BT = B
          CALL COARSE(NP,ATRUNC,BTRUNC,XMESH,0)
          CALL INITIA(COEFFN,PP,QP,WP,XMESH,NP,PARAMS,IPARAM,ICOFUN)
          GO TO 30

      ELSE IF (ISING.EQ.1) THEN
          ATRUNC = A + (B-A)/6.D0
          BTRUNC = B

      ELSE IF (ISING.EQ.2) THEN
          BTRUNC = B - (B-A)/6.D0
          ATRUNC = A

      ELSE
          TRUNC = (B-A)/6.D0
          ATRUNC = A + TRUNC
          BTRUNC = B - TRUNC
      END IF
C  Get the first mesh and coefficient functions:
      CALL COARSE(NP,ATRUNC,BTRUNC,XMESH,0)
      CALL INITIA(COEFFN,PP,QP,WP,XMESH,NP,PARAMS,IPARAM,ICOFUN)
C Check to see if emax is in continuous spectrum:
      IRELOC = .true.
      AT = ATRUNC
      BT = BTRUNC
      CALL CEKSPK(EMAX,A,B,AT,BT,XMESH,PP,QP,WP,NP,N,MAXITS,TOL,COEFFN,
     &            SETUP,IRELOC,ISYMM,CSPEC,DLAM,PARAMS,IPARAM,ICOFUN,
     &            IFAIL)
      IF (IFAIL.NE.0) GO TO 60
C If emax is not in the continuous spectrum then set the truncated
C endpoints and go to the eigenvalue counting stage:
      IF (CSPEC.LT.0.D0) THEN
          ATRUNC = AT
          BTRUNC = BT
          GO TO 30

      END IF
C If we are here then it means that emax is in the continuous spectrum.
C In such a case, we must also check emin.
      IRELOC = .true.
      AT = ATRUNC
      BT = BTRUNC
      CALL CEKSPK(EMIN,A,B,AT,BT,XMESH,PP,QP,WP,NP,N,MAXITS,TOL,COEFFN,
     &            SETUP,IRELOC,ISYMM,CSPEC,DLAM,PARAMS,IPARAM,ICOFUN,
     &            IFAIL)
      IF (IFAIL.NE.0) GO TO 60
      IF (CSPEC.GT.0.D0) THEN
C In this case emin is also in the continuous spectrum. There are therefore
C no eigenvalues in [emin,emax]. Denote this by setting kmin = kmax = -10.
          KMIN = -10
          KMAX = KMIN
          RETURN

      END IF
C If we are here then it means that the infimum of the continous spectrum
C lies somewhere in (emin,emax]. We must find this infimum. We will do so
C by using a rootfinding routine (C05AZF). Because we are rootfinding on a
C discontinuous function, C05AZF will be using bisection.

      XS = EMIN
      YS = EMAX
      IND = -1
      IR = 1
      NEVS = 0
      FX = -1.D0
      C(1) = 1.D0
      IFLAG = 1

10    CALL C05AZF(XS,YS,FX,TOL/20.D0,IR,C,IND,IFLAG)

      IF (IND.EQ.0) THEN
          IF (IFLAG.NE.0 .AND. IFLAG.NE.5) THEN
              IFAIL = 5
              GO TO 60

          END IF

          GO TO 20

      END IF

      IF (IND.LT.2 .OR. IND.GT.4) THEN
          IFAIL = 5
          GO TO 60

      END IF

      NEVS = NEVS + 1
      IF (NEVS.GT.MAXEVS) THEN
          IFAIL = 6
          GO TO 60

      END IF

      IRELOC = .true.
      AT = ATRUNC
      BT = BTRUNC
      CALL CEKSPK(XS,A,B,AT,BT,XMESH,PP,QP,WP,NP,N,MAXITS,TOL,COEFFN,
     &            SETUP,IRELOC,ISYMM,FX,DLAM,PARAMS,IPARAM,ICOFUN,IFAIL)
      IF (IFAIL.NE.0) GO TO 60
      GO TO 10

20    ELAM = XS

C We have now found elam, an approximation to the infimum of the continuous
C spectrum. We will proceed to the counting stage with
      EMAX = MAX(EMIN,ELAM-TOL/10.D0)
      IF (EMAX.LE.EMIN) THEN
          KMIN = -10
          KMAX = KMIN
          RETURN

      END IF

C Counting stage: we now count the number of eigenvalues in [emin,emax].
C We do this by finding kmin and kmax, respectively the lowest and highest
C indices of eigenvalues in this interval. We need to be able to compute
C D(emin) and D(emax) to fairly high accuracy.
C Stage 1 of counting: find truncated ends.
30    IRELOC = .false.
      AT = XMESH(0)
      BT = XMESH(NP)
      CALL CEKSPK(EMAX,A,B,AT,BT,XMESH,PP,QP,WP,NP,N,MAXITS,TOL,COEFFN,
     &            SETUP,IRELOC,ISYMM,CSPEC,DLAM,PARAMS,IPARAM,ICOFUN,
     &            IFAIL)
      IF (IFAIL.NE.0) GO TO 60
      IF (CSPEC.GT.0.D0) THEN
          IFAIL = 7
          GO TO 60

      END IF
C We now have truncated ends given by at,bt
C and so we can make two extra-accurate evaluations of D(elam). This
C will require the use of our stepsize-controlling version of the
C Prufer theta-integrator.

C Evaluation of D(emin):
      ELAM = EMIN
40    IF (ISYMM) THEN
          RMAX = 0.5D0* (A+B)
          IP1 = NP/2

      ELSE
          RMAX = (ELAM*WP(1)-QP(1))/PP(1)
          IP1 = 1
          DO 50 I = 1,NP - 1
              RTRY = (ELAM*WP(I)-QP(I))/PP(I)
              IF (RMAX.LT.RTRY) THEN
                  IP1 = I
                  RMAX = RTRY
              END IF

50        CONTINUE
          RMAX = XMESH(IP1)
      END IF

C      H = MIN((XMESH(1)-XMESH(0)), (XMESH(NP)-XMESH(NP-1)))
      H = (XMESH(NP)-XMESH(0))/DBLE(NP)
      DLAM = DACC(ELAM,AT,BT,RMAX,COEFFN,SETUP,H,TOL/1.D1,PARAMS,IPARAM,
     &       ICOFUN,IFAIL)
      IF (IFAIL.NE.0) GO TO 60
C      WRITE (6,*) 'ELAM,D/PI:',ELAM,DLAM/(4.D0*ATAN(1.D0))

      IF (ELAM.LT.EMAX) THEN
          IF (DLAM.LE.0.D0) THEN
              KMIN = 0
              IF (DLAM.GT. (-TOL)) IFAIL = 8

          ELSE
              DLAM = DLAM/PI
              KMIN = INT(DLAM)
              RMAX = ABS(DLAM-KMIN)
              IF (RMAX.GT. (1.D0-TOL) .OR. RMAX.LT.TOL) IFAIL = 8
              IF ((DLAM-DBLE(KMIN)*PI).NE.0.D0) KMIN = KMIN + 1
          END IF

          ELAM = EMAX
          GO TO 40

      ELSE
          DLAM = DLAM/PI
          KMAX = INT(DLAM)
          RMAX = ABS(DLAM-KMAX)
          IF (RMAX.GT. (1.D0-TOL) .OR. RMAX.LT.TOL) IFAIL = 8
      END IF

      IF (KMIN.GT.KMAX) THEN
          KMIN = -10
          KMAX = -10
      END IF

C Ifail = 8 means that there may be an eigenvalue lurking very close to
C one of the endpoints, which would be difficult to detect.

C Successful termination:
      A = AT
      B = BT
      RETURN

60    RETURN

      END SUBROUTINE CNTER

C -------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION DACC(ELAM,A,B,C,COEFFN,SETUP,H,TOL,
     &                               PARAMS,IPARAM,ICOFUN,IFAIL)
C     .. Parameters ..
      DOUBLE PRECISION ONE,HALF,TRD
      PARAMETER (ONE=1.D0,HALF=5.D-1,TRD=ONE/3.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,ELAM,H,TOL
      INTEGER ICOFUN,IFAIL,IPARAM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PARAMS(1:IPARAM)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL COEFFN,SETUP
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ABSERR,DI,DTH1,ERR,FAC,HMAX,HO,P,PDY,PI,Q,SCALE,
     &                 SNEW,THETA,THETA1,THETA2,W,X1,XFIN,XMID,XO,Y
      INTEGER ICNT,KNTRY
      LOGICAL STEP1,TRANS
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION SCL,X01AAF,X02AJF
cc      EXTERNAL SCL,X01AAF,X02AJF
C     ..
C     .. External Subroutines ..
cc      EXTERNAL BCONS,COEFUN,ONESTP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,MAX,MIN,SIGN,SQRT
C     ..
      PI = X01AAF(1.D0)
      IFAIL = 0
      STEP1 = .TRUE.
      HMAX = ABS(H)
      HO = SIGN(1.D0,H)*MIN(SQRT(TOL),ABS(H))
      H = HO
      ICNT = 0
C Step from a to c
      XO = A
      XFIN = C
      DI = ONE
      KNTRY = 0
      X1 = XO + H
      CALL COEFUN((XO+X1)*HALF,P,Q,W,COEFFN,PARAMS,IPARAM,ICOFUN)
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      TRANS = ICOFUN .EQ. 1 .OR. ICOFUN .EQ. 3
      CALL BCONS(SETUP,Y,PDY,ELAM,XO,P,Q,W,0,TRANS,IFAIL)
      IF (IFAIL.NE.0) RETURN
      THETA = ATAN2(Y,PDY)
C ------------------------- Stepping stage -------------------------
10    THETA1 = THETA
      THETA2 = THETA
      SCALE = ONE
      SNEW = ONE
      X1 = XO + H
      IF ((X1-XFIN)*DI.GT.0.D0) X1 = XFIN
C Do one-step + two half-steps to integrate with error estimate.
      XMID = HALF* (XO+X1)
      CALL COEFUN(XMID,P,Q,W,COEFFN,PARAMS,IPARAM,ICOFUN)
      CALL ONESTP(XO,X1,DI,SCALE,THETA1,P,ELAM*W-Q)
      CALL COEFUN(HALF* (XO+XMID),P,Q,W,COEFFN,PARAMS,IPARAM,ICOFUN)
      CALL ONESTP(XO,XMID,DI,SNEW,THETA2,P,ELAM*W-Q)
      CALL COEFUN(HALF* (XMID+X1),P,Q,W,COEFFN,PARAMS,IPARAM,ICOFUN)
      CALL ONESTP(XMID,X1,DI,SNEW,THETA2,P,ELAM*W-Q)
      THETA1 = SCL(THETA1,ONE/SCALE)
      THETA2 = SCL(THETA2,ONE/SNEW)
C Error estimate:
      ERR = (THETA1-THETA2)*TRD
      ABSERR = ABS(ERR)

C      WRITE (6,6000) XO,XMID,X1,H,ABSERR
C6000  FORMAT ('XO,XMID,X1,H,ABSERR:',5G16.6)

      IF (ABSERR.GT.TOL) THEN
C          WRITE (6,*) 'Error too large, repeating step'
          KNTRY = KNTRY + 1
          IF (ABS(H).LT.SQRT(X02AJF(H)) .AND. KNTRY.GT.16) THEN
              IFAIL = 2
              RETURN

          END IF

          H = H*MIN(8.D-1, (TOL/ABSERR)**TRD)
          GO TO 10

      END IF

      IF (ABSERR.LT.0.8D0*TOL .AND. KNTRY.EQ.0) THEN
C         WRITE (6,*) 'Error too small'
C         IF (step1) WRITE (6,*) 'Repeating step'
          H = H/MAX(2.D-1, (ABSERR/ (0.8D0*TOL))**TRD)
          FAC = HMAX*MAX(ONE,XMID**2)
          H = SIGN(1.D0,H)*MIN(FAC,ABS(H))
          IF (STEP1 .AND. ABS(H).LT.FAC) GO TO 10
      END IF

      STEP1 = .FALSE.

      KNTRY = 0
      THETA = THETA2 - ERR
      XO = X1
      ICNT = ICNT + 1
      IF ((XFIN-X1)*DI.GT.0.D0) GO TO 10
C --------------------- END of stepping stage ------------------------

C If we have got here then it means that we have finished an integration.

      IF (DI.GT.0.D0) THEN
          DTH1 = THETA
          XO = B
          XFIN = C
          DI = -ONE
          H = HO*DI
          X1 = XO + H
          CALL COEFUN((XO+X1)*HALF,P,Q,W,COEFFN,PARAMS,IPARAM,ICOFUN)
          CALL BCONS(SETUP,Y,PDY,ELAM,XO,P,Q,W,1,TRANS,IFAIL)
          IF (IFAIL.NE.0) RETURN
          THETA = PI - ATAN2(Y,PDY)
          STEP1 = .TRUE.
          GO TO 10

      END IF

      DACC = DTH1 - THETA
      RETURN

      END FUNCTION DACC

C -------------------------------------------------------------------

      SUBROUTINE CEKSPK(ELAM,A,B,AT,BT,XMESH,PP,QP,WP,N,NMAX,MAXITS,TOL,
     &                  COEFFN,SETUP,IRELOC,ISYMM,CSPEC,DLAM,PARAMS,
     &                  IPARAM,ICOFUN,IFAIL)

C Routine for spectral analysis. Tells when a value of elam is in the
C continuous spectrum. Used as part of a rootfinding process to find
C the infimum of the continuous spectrum.
C
C Brief Specification:
C
C elam: real. On entry elam must specify the value of lambda to be
C             tested. Unchanged on exit.
C a,b:  real. On entry, the endpoints of the whole interval on which
C             the problem is posed. Unchanged on exit.
C at,bt: real. On entry, proposed initial truncation points from which
C             stepping towards the singular endpoint(s) will take place.
C             On exit, the points at which the algorithm stopped the
C             truncation procedure.
C xmesh: real array of dimension at least (0:nmax). On entry, the
C             entries xmesh(0),..,xmesh(n) specify a mesh covering the
C             interval (at,bt), and must be strictly increasing with
C             xmesh(0) = at, xmesh(n) = bt. On exit, if ireloc was
C             .true. on entry then xmesh will be unchanged; otherwise,
C             xmesh will contain the mesh which was in use when the
C             interval truncation algorithm stopped, and n will be the
C             number of intervals in this mesh.
C pp,qp,wp: real arrays of dimension at least (1:nmax). On entry, the
C             entries in pp(1),..,pp(n), qp(1),..,qp(n), wp(1),..,wp(n)
C             specify coefficient values. On exit they will be unchanged
C             if ireloc was .true. on entry, otherwise they will be
C             altered to agree with the new xmesh array (see above).
C n: integer. On entry, n specifies the number of intervals in the
C             current mesh. If ireloc is .true. on entry then n will be
C             unchanged on exit; otherwise n will be the number of
C             intervals in the new mesh (see xmesh above).
C nmax: integer. On entry, nmax specifies the maximum amount of storage
C             which is available for meshpoints and coefficient values.
C             nmax must be strictly greater than n, and must satisfy
C                             nmax > n + maxits
C             if there is to be no risk failure. Unchanged on exit.
C maxits: integer. The maximum number of moves of truncated endpoints
C             which will be made before declaring that no suitable
C             truncated endpoints can be found and that elam must lie
C             in the continuous spectrum. In general maxits = 60 will be
C             adequate, but at tight tolerances tol it may be advisable
C             to increase maxit. Unchanged on exit.
C tol: real.  On entry, the tolerance to which two successive approximations
C             to D(elam) with different truncated endpoints must agree before
C             the endpoints are deemed to be satisfactory and elam is
C             declared not to lie in the continuous spectrum. Unchanged on
C             exit.
C COEFFN, SETUP: subroutines supplied by the user for the evaluation of
C             coefficients and boundary conditions.
C ireloc: logical. Controls whether or not the arrays xmesh, pp, qp, wp,
C             and the integer n are changed on exit. See above. Unchanged
C             on exit.
C isymm: logical. On entry, isymm = .true. tells the routine that the
C             problem is symmetric about (a+b)/2. In this case the truncated
C             interval (at,bt) should be symmetric in the same way.
C             Unchanged on exit.
C cspec: real. This parameter is for output only and serves to tell the user
C             whether or not elam was found to be in the continuous spectrum.
C             cspec = -1 on exit ==> elam is not in continuous spectrum.
C             cspec = +1 on exit ==> elam is in the continuous spectrum.
C dlam: real. On exit, this parameter gives the final approximation obtained
C             for D(elam). Only meaningful if cspec = -1 on exit.
C ifail: integer. The error flag.
C
C END of brief specification.
C Declarations:
C     .. Parameters ..
      DOUBLE PRECISION HALF,TQTR
      PARAMETER (HALF=5.D-1,TQTR=7.5d-1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,AT,B,BT,CSPEC,DLAM,ELAM,TOL
      INTEGER ICOFUN,IFAIL,IPARAM,MAXITS,N,NMAX
      LOGICAL IRELOC,ISYMM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PARAMS(1:IPARAM),PP(1:NMAX),QP(1:NMAX),
     &                 WP(1:NMAX),XMESH(0:NMAX)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL COEFFN,SETUP
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ATRUNC,BTRUNC,DL1,DL2,DTHA,DTHB,FYNYT,OVFLOW,PDY,
     &                 SC,SHIFT,TH,TOLRED,Y,H
      INTEGER I,I1,I2,IP1,ISING,ISINGO,ITIME,KNTR
      LOGICAL TRANS
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION D,DTHETA,X01AAF,X02AHF,X02ALF
cc      EXTERNAL D,DTHETA,X01AAF,X02AHF,X02ALF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,MIN,SQRT
C     ..
C     .. External Subroutines ..
cc      EXTERNAL BCONS,COEFUN
C     ..
      OVFLOW = SQRT(X02ALF(1.D0))
      IFAIL = 0
      TOLRED = TOL*1.D-1
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      TRANS = ICOFUN .EQ. 1 .OR. ICOFUN .EQ. 3
C Stage 1: Find out which ends are singular.
      IF (A.EQ.AT .AND. B.EQ.BT) THEN
          CSPEC = -1
          RETURN

      ELSE IF (A.LT.AT .AND. B.EQ.BT) THEN
          ISING = 1

      ELSE IF (A.EQ.AT .AND. B.GT.BT) THEN
          ISING = 2

      ELSE
          ISING = 3
      END IF

      ISINGO = ISING
C Stage 2: Initialise the stepping procedure.

C  Set size of coarse mesh:
      SHIFT = XMESH(1) - XMESH(0)
C  Mark ends of input arrays:
      I1 = 1
      I2 = N
C  END of marking; special action if two singular points:
      ITIME = 0
      KNTR = 0
      IF (ISING.EQ.3) THEN
          ITIME = 1
          ISING = 1
      END IF
C First approximation to D(elam):
      DL1 = D(ELAM,N,N/2,PP,QP,WP,XMESH,SETUP,TRANS,IFAIL)
      IF (IFAIL.NE.0) RETURN

10    IF (ISING.EQ.1) THEN
C Add a point at the left hand end.
          IF ((N+1).GT.NMAX) THEN
              IFAIL = 4
              RETURN

          END IF

          DO 20 I = N,1,-1
              IP1 = I + 1
              XMESH(IP1) = XMESH(I)
              PP(IP1) = PP(I)
              QP(IP1) = QP(I)
              WP(IP1) = WP(I)
20        CONTINUE
          XMESH(1) = XMESH(0)
          XMESH(0) = XMESH(0) - MIN(SHIFT,TQTR* (XMESH(0)-A))
          CALL COEFUN(HALF* (XMESH(1)+XMESH(0)),PP(1),QP(1),WP(1),
     &                COEFFN,PARAMS,IPARAM,ICOFUN)
          N = N + 1
          I1 = I1 + 1
          I2 = I2 + 1
      END IF

C END of adding point at left-hand side.

      IF (ISING.EQ.2 .OR. ISYMM) THEN
C Add point at right hand side.
          N = N + 1
          IF (N.GT.NMAX) THEN
              IFAIL = 4
              RETURN

          END IF

          XMESH(N) = XMESH(N-1) + MIN(SHIFT,TQTR* (B-XMESH(N-1)))
          CALL COEFUN(HALF* (XMESH(N)+XMESH(N-1)),PP(N),QP(N),WP(N),
     &                COEFFN,PARAMS,IPARAM,ICOFUN)
      END IF

      DL2 = D(ELAM,N,N/2,PP,QP,WP,XMESH,SETUP,TRANS,IFAIL)
      IF (IFAIL.NE.0) RETURN

C Check for convergence
      IF (ABS(DL2-DL1).LT.TOLRED) GO TO 30

C If we got past the last block then we have to make another step
      KNTR = KNTR + 1
      IF (KNTR.GE.MAXITS) GO TO 50
      H = ABS(XMESH(1) - XMESH(0))
      FYNYT = (ELAM*WP(1)-QP(1))/PP(1)
      IF ((FYNYT.GE.0.0D0.AND.H*SQRT(ABS(FYNYT)).GE.1.0D-2*X02AHF(0.D0))
     &    .OR.ABS(FYNYT).GE.OVFLOW) GO TO 50
      H = ABS(XMESH(N)-XMESH(N-1))
      FYNYT = (ELAM*WP(N)-QP(N))/PP(N)
      IF ((FYNYT.GE.0.0D0.AND.H*SQRT(ABS(FYNYT)).GE.1.0D-2*X02AHF(0.D0)) 
     &    .OR.ABS(FYNYT).GE.OVFLOW) GO TO 50
C If we got past the last line then we have not yet decided whether or
C not we are in the continuous spectrum. Do another endpoint shift:
      DL1 = DL2
      GO TO 10
C Various exits now follow

C First the exit in case where dl2-dl1 has converged - we are not in the
C continuous spectrum.

30    IF (ITIME.EQ.1 .AND. .NOT.ISYMM) THEN
C  We have two endpoints to do and only one has been done -- go back & do
C  the other.
          ITIME = 2
          ISING = 2
          DL1 = DL2
          GO TO 10

      END IF

      AT = XMESH(0)
      BT = XMESH(N)
C One last check for continuous spectrum:

      IF (ISINGO.NE.0) THEN
          DTHA = 0.D0
          DTHB = 0.D0
          IF (ISINGO.GT.1) THEN
              CALL BCONS(SETUP,Y,PDY,ELAM,BT,PP(N),QP(N),WP(N),1,TRANS,
     &                   IFAIL)
              IF (IFAIL.NE.0) RETURN
              TH = ATAN2(Y,PDY)
              SC = 1.D0
              BTRUNC = BT
              DTHB = DTHETA(TH,SC,ELAM,BTRUNC,B,COEFFN,PARAMS,IPARAM,
     &               ICOFUN)
          END IF

          IF (ISINGO.NE.2) THEN
              CALL BCONS(SETUP,Y,PDY,ELAM,AT,PP(1),QP(1),WP(1),0,TRANS,
     &                   IFAIL)
              IF (IFAIL.NE.0) RETURN
              TH = ATAN2(Y,PDY)
              SC = 1.D0
              ATRUNC = XMESH(0)
              DTHA = DTHETA(TH,SC,ELAM,ATRUNC,A,COEFFN,PARAMS,IPARAM,
     &               ICOFUN)
          END IF

          IF ((DTHA+DTHB).GE.1.D1*X01AAF(1.D0)) GO TO 50
      END IF

      IF (IRELOC) THEN
          N = I2 - I1 + 1
          XMESH(0) = XMESH(I1-1)
          DO 40 I = 1,N
              IP1 = I1 + I - 1
              XMESH(I) = XMESH(IP1)
              PP(I) = PP(IP1)
              QP(I) = QP(IP1)
              WP(I) = WP(IP1)
40        CONTINUE
      END IF

      DLAM = DL2
      CSPEC = -1.D0
      RETURN

C Next the case where dl2-dl1 did not converge - we are in the continuous
C spectrum
50    AT = XMESH(0)
      BT = XMESH(N)
      IF (IRELOC) THEN
          N = I2 - I1 + 1
          XMESH(0) = XMESH(I1-1)
          DO 60 I = 1,N
              IP1 = I1 + I - 1
              XMESH(I) = XMESH(IP1)
              PP(I) = PP(IP1)
              QP(I) = QP(IP1)
              WP(I) = WP(IP1)
60        CONTINUE
      END IF

      DLAM = DL2
      CSPEC = 1.D0
      RETURN
C END of routine
      END SUBROUTINE CEKSPK
C ---------------------------------------------------------------------
C --------------------- Source from SL02F -----------------------------
C ---------------------------------------------------------------------

      SUBROUTINE sl02f(elam,a,b,k,ainfo,binfo,sym,tol,coeffn,setup,n,wk,
     &                 iwk,wksmal,ismal,knobs,ifail)
C This routine finds the kth eigenvalue elam of a regular or singular Sturm-
C Liouville problem
C
C               -(p(x)y'(x))'+q(x)y(x) = elam.w(x).y(x)  (a<x<b)
C
C with general boundary conditions at regular end points.
C
C VARIABLES:
C ELAM: REAL ARRAY of DIMENSION (1:2). On entry elam(1) must be set to an
C       approximation to the eigenvalue and elam(2) to a strictly non-zero
C       estimate of the error in the approximation. Neither of these numbers
C       need be at all accurate; the code will always find an approximation to
C       the kth eigenvalue rather than the eigenvalue closest to the initial
C       approximation. On exit, elam(1) contains the
C       approximation found and elam(2) an estimate of the error in this
C       approximation.
C A,B: REAL. The end-points of the interval. These must be finite. If the
C       problem is posed on an infinite or semi-infinite interval then it must
C       be transformed onto a finite interval by means of a suitable
C       transformation. If either A or B is singular then on exit it will be set
C       to the value of the corrseponding truncated endpoint used by the
C       algorithm.
C   K: INTEGER. The index of the eigenvalue to be found when these are arranged
C             in ascending order starting with index 0. Unchanged on exit.
C TOL: REAL. The error control parameter. EIGEN will attempt to ensure the its
C        approximation to the eigenvalue satisfies the inequality
C              |elam-elamtrue| < TOL.max(1,|elam|).
C       Unchanged on exit.
C AINFO: CHARACTER*1, input. On entry AINFO should be set as follows:
C        If A is a finite regular endpoint set AINFO = 'R' or AINFO = 'r'.
C        If A is a finite singular endpoint set AINFO = 'S' or AINFO = 's'.
C        If A is -infinity set AINFO = 'I' or AINFO = 'i'
C        AINFO is unchanged on exit.
C BINFO: CHARACTER*1, input. On entry BINFO should be set as follows:
C        If B is a finite regular endpoint set BINFO = 'R' or BINFO = 'r'.
C        If B is a finite singular endpoint set BINFO = 'S' or BINFO = 's'.
C        If B is +infinity set BINFO = 'I' or BINFO = 'i'
C        BINFO is unchanged on exit.
C SYM:  LOGICAL. If the problem is symmetric about the midpoint (A+B)/2
C        set SYM = .TRUE., otherwise set SYM = .FALSE. If the problem is
C        posed on (-oo,+oo) and SYM = .TRUE. then the point x=0 will be
C        taken as the midpoint and the problem must be symmetric about
C        x=0.
C        SYM is unchanged on exit.
C COEFFN: SUBROUTINE supplied by the user. Its specification is
C              SUBROUTINE coeffn(x,p,q,w)
C              real x,p,q,w
C        This subroutine serves to supply the coefficient functions which define
C       the problem.
C SETUP: SUBROUTINE supplied by the user. Its specification is
C              SUBROUTINE setup(y,pdy,eig,x,iend,ising)
C              INTEGER iend
C              real y,pdy
C               LOGICAL ising
C       This routine supplies values of y(x) and p(x)y'(x) which satisfy the
C       boundary conditions at the (iend+1)th end of the interval. That is,
C       with iend=0 y and pdy must satisfy the boundary conditions at x=a,
C       while for iend=1 they must satisfy the boundary conditions at x=b.
C       If the end is regular then the user should set y and pdy as described,
C       then set ising = .false. and return. If the end is singular the
C       user need not set y or pdy: it suffices to set ising = .true. then
C       return. Boundary conditions will then be generated automatically.
C       x is the point at which boundary conditions are actually required
C       and is equal to the endpoint unless the problem is singular. Its
C       value should not be changed. It is included to allow the user to
C       monitor the interval truncation process for singular problems.
C   N : INTEGER. In order to achieve computational economy the routine
C       stores the values of the coefficient functions and meshpoints in an
C       array,rather than re-evaluate repeatedly. There is no way of knowing in
C       advance how large this array might need to be since this will depend
C       on the problem concerned, on the eigenvalue requested, and on the
C       accuracy required, although on exit N will be reset to allow the user
C       to see how much storage was actually required. N is the maximum number
C       of mesh intervals which the user is prepared to allow the program to
C       use. In order to minimise the probability of the required accuracy not
C       being achieved, the user is advised to set N to a large number
C       (e.g. 5000; this copes with all the problems we have tested so far,
C       to reasonably high accuracies). If the required accuracy cannot be
C       achieved on a mesh with less than N intervals then the routine will
C       either stop (if IFAIL=0 on entry) or continue to compute an
C       approximation on a mesh of min(N,800) intervals (if IFAIL=1 on entry).
C       In the later case IFAIL will be set to an appropriate value on exit.
C       On exit, N contains the number of intervals used in the computation.
C    WK: real ARRAY of DIMENSION at least (0:N,4). Used as storage space.
C        If the eigenfunction is also required it may be obtained by passing
C        this array, after the call to SL02F, through the routines NORONE
C        and SNGLFN.
C   IWK: the first dimension of WK as declared in the calling (sub)program.
C         IWK must be bigger than N.
C WKSMAL: real ARRAY of DIMENSION (0:ISMAL,1:7). Used as workspace.
C ISMAL: On entry, ISMAL must specify a first dimension for the array
C       WKSMAL as declared in the calling (sub)program. ISMAL must satisfy
C       the following criteria:
C    1) in order to guarantee that no failures will occur, it must not
C       be less than 2*knobs(1)+knobs(2), where knobs(1) and knobs(2) are as
C       below. If the user opts for the default values of knobs(1) and knobs(2)
C       then the appropriate value of ISMAL is 130;
C    2) even if the user reduces knobs(1), say, ISMAL must never be set to
C       a value which is less than 100.
C KNOBS: INTEGER ARRAY of DIMENSION at least 2. This array may be used to
C       change the values of certain parameters which control computation.
C       The user who wishes to use all the default opttions should include
C       a statement
C                             DATA KNOBS/2*0/
C       in his (sub)program before the call to SL02F. The user who wishes
C       change the defaults controlled by KNOBS will require the following
C       information:
C      knobs(1) is only used for singular problems. It is the maximum number
C        of steps which will be performed in attemting to locate suitable
C        truncated endpoints for a singular problem.  At each step the
C        tentative truncated endpoint is moved towards the singular endpoint,
C        and the maximum shift is controlled by knobs(2) as indicated below.
C        The default value of knobs(1) is 60, and it is NOT RECOMMENDED that
C        this be increased, unless knobs(2) has been set greater than the
C        default value.
C      knobs(2) is the number of mesh-points in an initial coarse mesh which
C        obtains a first estimate of the eigenvalue and which, in the case of a
C        singular problem, is the first of a sequence of meshes which will be
C        used to find suitable truncated end points. The default value of
C        knobs(2) is 10. The minimum value which can be set is 6. Increasing
C        knobs(2) may occasionally result in improved accuracy, but will always
C        cause an increase in run-times. The value of knobs(2) also controls
C        the maximum stepsize which can be used in approaching a singular
C        endpoint, which is 2.(B-A)/(3.knobs(2)) if there are 2 singular
C        endpoints in the problem and  5.(B-A)/(6.knobs(2)) if there is 1
C        singular point. Thus if knobs(2) has been set particularly large, it
C        may be advisable to increase knobs(1) since a large number of steps
C        may be necessary before the enpoint is sufficiently close. Finally
C        in both regular and singular problems, knobs(2) also controls the
C        maximum stepsize during the final (high accuracy) integration: this
C        is (B-A)/knobs(2), where A and B denote truncated or input endpoints,
C        depending on whether or not the latter are singular.
C        The array KNOBS is unchanged on exit.
C IFAIL: INTEGER. On entry this parameter must be set to 0 or 1 in the usual
C        way. On exit, IFAIL will be set as follows:
C       IFAIL=1: Parameter error. TOL<=0, ELAM(2)=0,N<=0, B<=A all
C        produce this value of IFAIL on exit. ELAM will not be changed in this
C        case.
C       IFAIL=2: p(x) <= 0 for some x in (a,b)
C       IFAIL=3: Error in the routine SETUP: both y and pdy are zero.
C       IFAIL=4: The eigenvalue appears not to exist due to onset of continuous
C                 spectrum.
C       IFAIL=5: The truncated endpoints are unsatisfactory but cannot be
C                adjusted further without risking overflow. elam(1) contains
C                the eigenvalue approimation obtained on the last coarse mesh
C                used before failure; elam(2) is unchanged.
C       IFAIL=6: Too many evaluations of the Prufer miss-distance while
C       rootfinding. This is extremely unlikely since the internal value set
C       by the code is 100.
C       IFAIL=7: Too many attempts to bracket the eigenvalue. This error exit
C                may occur if elam(2) was too small on entry.
C       IFAIL=8: More than knobs(1) attempts have been made to find (a)
C       satisfactory truncated endpoint(s) without success. This error exit
C       can only be taken if the problem is declared to be singular by the
C       user.
C       IFAIL=9: The Prufer miss distance function appears to be non-monotone.
C        This is a serious failure and Marco Central Office should be informed.
C        John Pryce will then be subjected to severe persecution since he wrote
C        the miss-distance function routine. Actually I have messed it about
C        slightly since then, but any bugs should still be all his fault.
C       IFAIL=10: The meshing routine cannot produce a mesh fine enough to meet
C        the user's tolerance, but not because of storage problems. This may
C        happen if the coefficient functions are of a particularly pathological
C        nature. Check the routine COEFFN for errors and there are none increase
C        TOL.
C       IFAIL=11: The required accuracy cannot be met within the user-supplied
C        storage. Increase TOL or increase N.
C       IFAIL=12: If IFAIL was 1 on entry then the failure indicated by IFAIL
C        =11 will in this case be indicated by IFAIL = 12, and the
C        routine will have produced an estimate of the eigenvalue on an
C        emergency mesh, together with an estimate of the error.
C       IFAIL = 13. A serious rootfinding problem; contact God. (God did not
C        write C05AZF but he is probably the only person who knows how it
C        works).
C       IFAIL = 14. The eigenfunctions are too ill-conditioned or the problem
C        too pathological for an acceptable eigenvalue estimate to be
C        generated with such a coarse initial mesh. Try increasing knobs(2).
C     .. Parameters ..
      INTEGER iparam
      DOUBLE PRECISION one,p6
      PARAMETER (iparam=1,one=1.D0,p6=6.D-1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,tol
      INTEGER ifail,ismal,iwk,k,n
      LOGICAL sym
      CHARACTER ainfo*1,binfo*1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION elam(1:2),wk(0:iwk,1:4),wksmal(0:ismal,1:7)
      INTEGER knobs(1:2)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn,setup
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION dummy1,dummy2,elam0,eps
      INTEGER astat,bstat,icofun,ifo,imatch,ising,isymm,lknt,maxit,np,
     &        nrec
      LOGICAL noxtrp
      CHARACTER srname*6
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION params(1:iparam)
      CHARACTER rec(2)*80
C     ..
C     .. External Subroutines ..
cc      EXTERNAL eigen
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf
cc      INTEGER p01abf
cc      EXTERNAL x02ajf,p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max
C     ..
      ifo = ifail
      eps = x02ajf(one)**p6
      srname = 'SL02F'
C Check for input parameter errors.

      IF (ainfo.EQ.'R' .OR. ainfo.EQ.'r') THEN
          astat = 0

      ELSE IF (ainfo.EQ.'S' .OR. ainfo.EQ.'s') THEN
          astat = 1

      ELSE IF (ainfo.EQ.'I' .OR. ainfo.EQ.'i') THEN
          astat = 2

      ELSE
          astat = 3
      END IF

      IF (binfo.EQ.'R' .OR. binfo.EQ.'r') THEN
          bstat = 0

      ELSE IF (binfo.EQ.'S' .OR. binfo.EQ.'s') THEN
          bstat = 1

      ELSE IF (binfo.EQ.'I' .OR. binfo.EQ.'i') THEN
          bstat = 2

      ELSE
          bstat = 3
      END IF
      IF (k.LT.0 .OR. n.GT.iwk .OR. n.LT.10 .OR. ismal.LT.10 .OR.
     &    astat.EQ.3 .OR. bstat.EQ.3 .OR. (sym.AND.astat.NE.bstat) .OR.
     &    elam(2).LE.0.D0 .OR. tol.LT.eps .OR.
     &    (astat.NE.2.AND.bstat.NE.2.AND.b.LE.a) .OR.
     &    abs(ifo).GT.1) THEN

          ifail = 1
          nrec = 1
          IF (k.LT.0) THEN
              WRITE (rec,FMT=9000) k
              GO TO 160

          END IF

          IF (astat.EQ.3) THEN
              WRITE (rec,FMT=9040)
              GO TO 160

          END IF

          IF (bstat.EQ.3) THEN
              WRITE (rec,FMT=9050)
              GO TO 160

          END IF

          IF (astat.NE.bstat .AND. sym) THEN
              WRITE (rec,FMT=9060)
              GO TO 160

          END IF
          IF (n.GT.iwk) THEN
              WRITE (rec,FMT=9010) n,iwk
              GO TO 160

          END IF
          IF (n.LT.10) THEN
              WRITE (rec,FMT=9020) n
              GO TO 160

          END IF
          IF (ismal.LT.10) THEN
              WRITE (rec,FMT=9030) ismal
              GO TO 160

          END IF

          IF (elam(2).LE.0.0D0) THEN
              WRITE (rec,FMT=9070) elam(2)
              GO TO 160

          END IF

          nrec = 2
          IF (tol.LT.eps) THEN
              WRITE (rec,FMT=9080) tol,eps
              GO TO 160

          END IF

          IF (astat.NE.2 .AND. bstat.NE.2 .AND. b.LE.a) THEN
              WRITE (rec,FMT=9090) a,b
              GO TO 160

          END IF

          IF (abs(ifo).GT.1) THEN
              WRITE (rec,FMT=9100) ifo
              ifo = -1
              GO TO 160

          END IF

      END IF
C END of input error trapping

      IF (astat.EQ.0 .AND. bstat.EQ.0) ising = 0
      IF (astat.NE.0 .AND. bstat.EQ.0) ising = 1
      IF (bstat.NE.0 .AND. astat.EQ.0) ising = 2
      IF (astat.NE.0 .AND. bstat.NE.0) ising = 3

      icofun = 2

      IF (astat.EQ.2 .AND. bstat.EQ.2) THEN
          a = -1.D0
          b = 1.D0
          icofun = 1
          GO TO 10

      END IF

      IF (astat.EQ.2) THEN
          a = -1.D0
          b = b/ (1.D0+abs(b))
          icofun = 1
      END IF

      IF (bstat.EQ.2) THEN
          b = 1.D0
          a = a/ (1.D0+abs(a))
          icofun = 1
      END IF
C Sneaky storage of the un-truncated endpoints:
      wk(0,3) = a
      wk(0,4) = b
C END of sneaky storage

10    isymm = 0
      IF (sym) isymm = 1
      noxtrp = .false.
      np = knobs(2)
      lknt = knobs(1)
      maxit = 100
      imatch = 0
      dummy1 = elam(1)
      dummy2 = elam(2)
      CALL eigen(dummy1,elam0,dummy2,a,b,k,n,ismal,imatch,np,ising,
     &           isymm,lknt,maxit,coeffn,setup,wk(1,2),wk(1,3),wk(1,4),
     &           wk(0,1),wksmal(1,2),wksmal(1,3),wksmal(1,4),
     &           wksmal(0,1),wksmal(0,5),wksmal(0,6),wksmal(0,7),tol,
     &           noxtrp,params,iparam,icofun,ifail)
      elam(1) = dummy1
      elam(2) = dummy2
C Sneaky storage:
      wk(0,2) = imatch
C End of sneaky storage.
C *********************************************************************
C WARNING: for infinite interval problems, the mesh is in the user's  *
C coordinates but the endpoints are in finite transformed coordinates *
C given by xfinite = xinfinite/(1+abs(xinfinite))                     *
C *********************************************************************
      IF (icofun.EQ.1) THEN
C Transform endpoints to user's coordinates.
          a = a/max((1.D0-abs(a)),x02ajf(1.D0))
          b = b/max((1.D0-abs(b)),x02ajf(1.D0))
      END IF

      IF (ifail.NE.0) THEN
C Output error trapping
          nrec = 1
          GO TO (20,30,40,50,60,70,80,90,100,
     &           110,120,130,140,150) ifail

20        WRITE (rec,FMT=9110)
          GO TO 160

30        WRITE (rec,FMT=9120)
          GO TO 160

40        WRITE (rec,FMT=9130)
          GO TO 160

50        WRITE (rec,FMT=9140)
          GO TO 160

60        WRITE (rec,FMT=9150)
          GO TO 160

70        WRITE (rec,FMT=9160) maxit
          GO TO 160

80        WRITE (rec,FMT=9170) maxit
          GO TO 160

90        WRITE (rec,FMT=9180) knobs(1)
          GO TO 160

100       WRITE (rec,FMT=9190)
          GO TO 160

110       WRITE (rec,FMT=9200)
          GO TO 160

120       WRITE (rec,FMT=9210)
          GO TO 160

130       WRITE (rec,FMT=9220)
          GO TO 160

140       WRITE (rec,FMT=9230)
          GO TO 160

150       WRITE (rec,FMT=9240)
          GO TO 160

      END IF
C END of error trapping
C Elam(1) is our good approximation to the eigenvalue.
      RETURN

160   ifail = p01abf(ifo,ifail,srname,nrec,rec)
      RETURN

9000  FORMAT ('** K must not be negative, K = ',i8)
9010  FORMAT ('** N must not exceed IWK, N =',i8,'IWK =',i8)
9020  FORMAT ('** N must be at least 10, N =',i2)
9030  FORMAT ('** ISMAL must be at least 10, ISMAL=',i2)
9040  FORMAT ('** AINFO has an invalid value')
9050  FORMAT ('** BINFO has an invalid value')
9060  FORMAT ('** AINFO and BINFO must be the same if SYM = .TRUE.')
9070  FORMAT ('** ELAM(2) must be strictly positive, ELAM(2) =',G18.8)
9080  FORMAT ('** TOL is too small for optimal accuracy ',/,'** TOL =',
     &       g18.8,' smallest allowed = ',g18.8)
9090  FORMAT ('** For finite interval problems A < B is required ',/,
     &       '** A = ',g18.8,' B = ',g18.8)
9100  FORMAT ('** On entry IFAIL must be -1, 0 or +1 **',/,
     &       '** Input value of IFAIL = ',i8)
9110  FORMAT (' ** Parameter error')
9120  FORMAT (' ** P(X) is not strictly positive')
9130  FORMAT (' ** Invalid boundary conditions')
9140  FORMAT (' ** Cannot find eigenvalue with this index')
9150  FORMAT (' ** Danger of floating overflow near endpoints')
9160  FORMAT ('**',i8,' iterations were not enough to locate eigenvalue'
     &       )
9170  FORMAT ('** More than ',i8,' attempts to bracket eigenvalue')
9180  FORMAT ('** More than ',i8,' attempts to truncate interval')
9190  FORMAT (' ** Miss-distance non-monotone over a wide range')
9200  FORMAT (' ** Meshing routine could make no further progress')
9210  FORMAT (' ** Meshing routine ran out of storage')
9220  FORMAT (' ** Meshing routine ran out of storage: used ad-hoc mesh'
     &       )
9230  FORMAT (' ** Serious internal error: try different TOL')
9240  FORMAT (' ** Cannot obtain required accuracy')

      END SUBROUTINE sl02f
C ---------------------------------------------------------------------
C --------------------- Source from SL02FM-----------------------------
C ---------------------------------------------------------------------
      SUBROUTINE sl02fm(elam,a,b,m,k,ainfo,binfo,sym,tol,coeffn,setup,n,
     &                  wk,iwk,wksmal,ismall,knobs,ifail)
C     .. Parameters ..
      DOUBLE PRECISION smatch,zero
      INTEGER iparam
      DOUBLE PRECISION p6,two,half
      PARAMETER (smatch=1.D0,zero=0.D0,iparam=1,p6=6.D-1,two=2.D0,
     &          half=5.D-1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,tol
      INTEGER ifail,ismall,iwk,m,n
      LOGICAL sym
      CHARACTER ainfo*1,binfo*1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION elam(1:2,1:m),wk(0:iwk,1:4),wksmal(0:ismall,1:7)
      INTEGER k(1:m),knobs(1:2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION compar,dummy,elam0,elam2,eps,toloc
      INTEGER astat,bstat,i,icfn,icofun,iend,iflag,ifo,ii,imatch,iq,
     &        ising,istart,isymm,lknt,maxit,newmsh,no,np,nrec
      LOGICAL noxtrp,shift,trans
      CHARACTER srname*6
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION params(1:iparam)
      CHARACTER rec(2)*80
C     ..
C     .. External Subroutines ..
cc      EXTERNAL admesh,check,eigen,fillin,initia,norfun,solve
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02agf,x02ajf
cc      INTEGER p01abf
cc      EXTERNAL x02agf,x02ajf,p01abf
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn,setup
C     ..
C     .. Scalars in Common ..
      INTEGER icall,igt0,ilt0,int,ip,ismal,ival
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,max,min
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival,igt0,ilt0,ismal
C     ..
      ifo = ifail
      srname = 'SL02FM'
      eps = x02ajf(zero)**p6

      IF (ainfo.EQ.'R' .OR. ainfo.EQ.'r') THEN
          astat = 0

      ELSE IF (ainfo.EQ.'S' .OR. ainfo.EQ.'s') THEN
          astat = 1

      ELSE IF (ainfo.EQ.'I' .OR. ainfo.EQ.'i') THEN
          astat = 2

      ELSE
          astat = 3
      END IF

      IF (binfo.EQ.'R' .OR. binfo.EQ.'r') THEN
          bstat = 0

      ELSE IF (binfo.EQ.'S' .OR. binfo.EQ.'s') THEN
          bstat = 1

      ELSE IF (binfo.EQ.'I' .OR. binfo.EQ.'i') THEN
          bstat = 2

      ELSE
          bstat = 3
      END IF

      elam2 = 1.D0
      DO 10 i = 1,m
          elam2 = min(elam2,elam(2,i))
10    CONTINUE

      IF (n.GT.iwk .OR. n.LT.10 .OR. ismall.LT.10 .OR. (astat.EQ.3) .OR.
     &    (bstat.EQ.3) .OR. tol.LT.eps .OR. elam2.LE.zero .OR.
     &    (astat.NE.2.AND.bstat.NE.2.AND.b.LE.a) .OR.
     &    abs(ifo).GT.1) THEN

          ifail = 1
          nrec = 1

          IF (n.LT.10) THEN
              WRITE (rec,FMT=9040) n
              GO TO 250

          END IF
          IF (ismall.LT.10) THEN
              WRITE (rec,FMT=9050) ismall
              GO TO 250

          END IF
          IF (astat.EQ.3) THEN
              WRITE (rec,FMT=9010)
              GO TO 250

          END IF

          IF (bstat.EQ.3) THEN
              WRITE (rec,FMT=9020)
              GO TO 250

          END IF

          IF (sym .AND. astat.NE.bstat) THEN
              WRITE (rec,FMT=9030)
              GO TO 250

          END IF

          IF (n.GT.iwk) THEN
              WRITE (rec,FMT=9000) n,iwk
              GO TO 250

          END IF

          nrec = 2
          IF (tol.LT.eps) THEN
              WRITE (rec,FMT=9060) tol,eps
              GO TO 250

          END IF

          IF (elam2.LE.zero) THEN
              WRITE (rec,FMT=9070)
              GO TO 250

          END IF

          IF (astat.NE.2 .AND. bstat.NE.2 .AND. b.LE.a) THEN
              WRITE (rec,FMT=9080) a,b
              GO TO 250

          END IF

          IF (abs(ifo).GT.1) THEN
              WRITE (rec,FMT=9090) ifo
              ifo = -1
              GO TO 250

          END IF

      END IF

9000  FORMAT ('** N must not exceed IWK, N =',i8,'IWK =',i8)
9010  FORMAT ('** Invalid value of AINFO on entry')
9020  FORMAT ('** Invalid value of BINFO on entry')
9030  FORMAT ('** If SYM = .TRUE. then AINFO and BINO must be the same')
9040  FORMAT ('** N must be at least 10, N =',I2)
9050  FORMAT ('** ISMALL must be at least 10, ISMALL=',I2)
9060  FORMAT ('** TOL is too small for optimal accuracy ',/,'** TOL =',
     &       g18.8,' smallest allowed = ',g18.8)
9070  FORMAT ('** All the ELAM(2,i) must be strictly positive')
9080  FORMAT ('** For regular problems A < B is required ',/,'** A = ',
     &       g18.8,' B = ',g18.8)
9090  FORMAT ('** On entry IFAIL must be -1, 0 or +1 **',/,
     &       '** Input value of IFAIL = ',i8)

C Store input value of n, which will be changed on call to EIGEN:
      no = n

C Translate the control parameters for SL02FM into control parameters for
C EIGEN

      icofun = 2

      IF (astat.EQ.0 .AND. bstat.EQ.0) ising = 0
      IF (astat.NE.0 .AND. bstat.EQ.0) ising = 1
      IF (bstat.NE.0 .AND. astat.EQ.0) ising = 2
      IF (astat.NE.0 .AND. bstat.NE.0) ising = 3
      isymm = 0
      IF (sym) isymm = 1

      IF (astat.EQ.2 .AND. bstat.EQ.2) THEN
          a = -1.D0
          b = 1.D0
          icofun = 1
          GO TO 20

      END IF

      IF (astat.EQ.2) THEN
          a = -1.D0
          b = b/ (1.D0+abs(b))
          icofun = 1
      END IF

      IF (bstat.EQ.2) THEN
          b = 1.D0
          a = a/ (1.D0+abs(a))
          icofun = 1
      END IF

C Check the input to make sure that the eigenvalues are being requested in
C a suitable order
20    IF (k(1).LT.0) THEN
          ifail = 1
          nrec = 1
          WRITE (rec,FMT=9100)

9100      FORMAT ('** The indices K(i) must be non-negative')

          GO TO 250

      END IF

      DO 30 i = 1,m - 1
          IF (k(i).GE.k(i+1)) THEN
              ifail = 1
              nrec = 1
              WRITE (rec,FMT=9110)

9110          FORMAT ('** The indices K(i) must be in ascending order')

              GO TO 250

          END IF

30    CONTINUE

C Get the first approximation to each of the eigenvalues.
C 1: Get an approximation to the eigenvalue with the highest index:

      elam0 = elam(1,m)
      eps = elam(2,m)
      iflag = ifo
      noxtrp = .true.
      np = knobs(2)
      lknt = knobs(1)
      maxit = 100
C The call to EIGEN gives the truncated endpoints a and b and a mesh
C but the eigenvalue approximation returned is computed on the coarse
C mesh which gives the initial approximations to the eigenfunctions.
      CALL eigen(elam(1,m),dummy,eps,a,b,k(m),n,ismall,imatch,np,ising,
     &           isymm,lknt,maxit,coeffn,setup,wk(1,2),wk(1,3),wk(1,4),
     &           wk(0,1),wksmal(1,2),wksmal(1,3),wksmal(1,4),
     &           wksmal(0,1),wksmal(0,5),wksmal(0,6),wksmal(0,7),tol,
     &           noxtrp,params,iparam,icofun,iflag)
      ifail = iflag
C *********************************************************************
C WARNING: for infinite interval problems, the mesh is in the user's  *
C coordinates but the endpoints are in finite transformed coordinates *
C given by xfinite = xinfinite/(1+abs(xinfinite))                     *
C *********************************************************************
      IF (icofun.EQ.1) THEN
C Transform endpoints to user's coordinates
          a = a/max((1.D0-abs(a)),x02agf(1.D0))
          b = b/max((1.D0-abs(b)),x02agf(1.D0))
      END IF
C ERROR TRAPPING:
C   If (ifail>0) and (ifail.ne.12.OR.14) the calculation must be stopped:
      IF (ifail.EQ.12 .OR. ifail.EQ.14) GO TO 40
      IF (ifail.NE.0) GO TO 100
C END OF ERROR TRAPPING.

C This has generated a mesh which is appropriate for the highest index
C eigenvalue. We now test the mesh on all the other eigenfunctions to make
C sure that it is appropriate for the lower index eigenvalues. In the
C process we generate coarse-mesh approximations to the lower index
C eigenpairs.

      istart = 0
      iend = no

      icfn = icofun
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      IF (icofun.EQ.1 .OR. icofun.EQ.3) THEN
          icfn = icofun + 1
      END IF

40    IF (m.GT.1) THEN
C IFAIL has one of the values 0, 12, 14.
          IF (ifail.NE.0) GO TO 60
          DO 50 i = 1,m
C Adapt the mesh appropriately for the K(i)th eigenfunction.
C First compute K(i)th coarse-mesh eigenvalue.
              IF (i.EQ.1) THEN
                  IF (elam(1,i).GT.elam(1,m)) THEN
                      eps = max(eps,half* (elam(1,i)-elam(1,m)))
                      elam(1,i) = elam(1,m) - two*eps
                  END IF

              END IF

              IF (i.GT.1 .AND. i.LT.m) THEN
                  IF (elam(1,i).LT.elam(1,i-1) .OR.
     &                elam(1,i).GT.elam(1,m)) THEN
                      elam(1,i) = dble((k(i)-k(i-1))/ (k(m)-k(i-1)))**2
                      elam(1,i) = elam(1,i-1) +
     &                            elam(1,i)* (elam(1,m)-elam(1,i-1))
                      eps = half* (elam(1,m)-elam(1,i-1))
                  END IF

              END IF

              iq = np/2
              trans = .false.
              toloc = 0.D0
              CALL solve(elam(1,i),eps,k(i),np,iq,maxit,setup,
     &                   wksmal(1,2),wksmal(1,3),wksmal(1,4),
     &                   wksmal(0,1),trans,toloc,iflag)
C Error trap:
              IF (iflag.NE.0) THEN
                  ifail = iflag
                  GO TO 100

              END IF
C END of error trap.

C Next compute coarse mesh eigenfunction:
              elam0 = elam(1,i)
              CALL norfun(wksmal(0,5),wksmal(0,6),wksmal(0,7),smatch,
     &                    elam0,zero,k(i),np,iq,wksmal(1,2),wksmal(1,3),
     &                    wksmal(1,4),wksmal(0,1),setup,trans,iflag)
              elam(1,i) = elam0
C Error trap (only possibility is IFLAG = 3)
              IF (iflag.NE.0) THEN
                  ifail = iflag
                  GO TO 100

              END IF
C END of error trap.

C Now use coarse mesh eigenfunction to adapt the mesh.
              shift = .true.
              IF (i.EQ.1) shift = .false.
              CALL admesh(elam(1,i),wksmal(1,2),wksmal(1,3),wksmal(1,4),
     &                    wksmal(0,1),wksmal(0,5),wksmal(0,6),
     &                    wksmal(0,7),wk(1,2),wk(1,3),wk(1,4),wk(0,1),n,
     &                    np,tol,iq,shift,istart,iend,params,iparam,
     &                    coeffn,icfn,iflag)
              IF (iflag.NE.0) THEN
C Only possibilities are IFAIL = 10 and IFAIL = 11
                  ifail = iflag
                  GO TO 100

              END IF

              n = iend
              iend = no
50        CONTINUE
C Put the correct function values back into the arrays pp,qp,wp:
          CALL initia(coeffn,wk(1,2),wk(1,3),wk(1,4),wk(0,1),n,params,
     &                iparam,icfn)
C We now have an adapted mesh, and it is stored in the places 0,1,..,iend
C on the mesh. This is our fine mesh with which we continue the computation.

          CALL check(n,wk(1,2),iflag)
          IF (iflag.NE.0) THEN
              ifail = iflag
              GO TO 100

          END IF

60        DO 80 i = 1,m
              eps = elam(2,i)/10.D0
              imatch = 1
              dummy = (elam(1,i)*wk(imatch,4)-wk(imatch,3))/wk(imatch,2)
              DO 70 ii = 1,n
                  compar = (elam(1,i)*wk(ii,4)-wk(ii,3))/wk(ii,2)
                  IF (compar.GT.dummy) THEN
                      imatch = ii
                      dummy = compar
                  END IF

70            CONTINUE
              IF (sym) imatch = n/2
              IF (i.EQ.1) THEN
                  IF (elam(1,i).GT.elam(1,m)) THEN
                      eps = max(eps,half* (elam(1,i)-elam(1,m)))
                      elam(1,i) = elam(1,m) - two*eps
                  END IF

              END IF

              IF (i.GT.1 .AND. i.LT.m) THEN
                  IF (elam(1,i).LT.elam(1,i-1) .OR.
     &                elam(1,i).GT.elam(1,m)) THEN
                      elam(1,i) = dble((k(i)-k(i-1))/ (k(m)-k(i-1)))**2
                      elam(1,i) = elam(1,i-1) +
     &                            elam(1,i)* (elam(1,m)-elam(1,i-1))
                      eps = half* (elam(1,m)-elam(1,i-1))
                  END IF

              END IF

              trans = .false.
              toloc = tol/10.D0
              CALL solve(elam(1,i),eps,k(i),n,imatch,maxit,setup,
     &                   wk(1,2),wk(1,3),wk(1,4),wk(0,1),trans,toloc,
     &                   iflag)
              IF (iflag.NE.0) THEN
                  ifail = iflag
                  GO TO 100

              END IF

80        CONTINUE
      END IF
C Now halve the mesh and compute improved approximations:

      newmsh = 2*n
      icfn = 2
      CALL fillin(n,newmsh,wk(0,1),wk(1,2),wk(1,3),wk(1,4),coeffn,
     &            params,iparam,icfn)
      imatch = 2*imatch
      IF (sym) imatch = newmsh/2
      n = newmsh

C Recompute the eigenvalues:

      DO 90 i = 1,m
          eps = elam(2,i)/20.D0
          elam0 = elam(1,i)
          trans = .false.
          toloc = 0.D0
          CALL solve(elam0,eps,k(i),n,imatch,maxit,setup,wk(1,2),
     &               wk(1,3),wk(1,4),wk(0,1),trans,toloc,iflag)
          IF (iflag.NE.0) THEN
              ifail = iflag
              GO TO 100

          END IF
C Aitken extrapolation:
          elam(2,i) = (elam(1,i)-elam0)/3.D0
          elam(1,i) = elam0 - elam(2,i)
90    CONTINUE
      a = wk(0,1)
      b = wk(n,1)
      wk(0,2) = imatch
      IF (ifail.NE.12 .AND. ifail.NE.14) ifail = 0
      RETURN

100   IF (ifail.NE.0) THEN
C Output error trapping
          nrec = 1
          GO TO (110,120,130,140,150,160,170,180,190,
     &           200,210,220,230,240) ifail

110       WRITE (rec,FMT=9120)

9120      FORMAT ('** Parameter error')

          GO TO 250

120       WRITE (rec,FMT=9130)

9130      FORMAT ('** P(X) is not strictly positive')

          GO TO 250

130       WRITE (rec,FMT=9140)

9140      FORMAT ('** Invalid boundary conditions')

          GO TO 250

140       WRITE (rec,FMT=9150) k(m)

9150      FORMAT ('** Cannot find eigenvalue with index',I6)

          GO TO 250

150       WRITE (rec,FMT=9160)

9160      FORMAT ('** Danger of floating overflow near endpoints')

          GO TO 250

160       WRITE (rec,FMT=9170) maxit

9170      FORMAT ('**',i8,
     &           ' iterations were not enough to locate an eigenvalue')

          GO TO 250

170       WRITE (rec,FMT=9180) maxit

9180      FORMAT ('** More than ',i8,
     &           ' attempts to bracket an eigenvalue')

          GO TO 250

180       WRITE (rec,FMT=9190) knobs(1)

9190      FORMAT ('** More than ',i8,' attempts to truncate interval')

          GO TO 250

190       WRITE (rec,FMT=9200)

9200      FORMAT (' ** Miss-distance non-monotone over a wide range')

          GO TO 250

200       WRITE (rec,FMT=9210)

9210      FORMAT ('** Meshing routine could make no further progress')

          GO TO 250

210       WRITE (rec,FMT=9220)

9220      FORMAT ('** Meshing routine ran out of storage')

          GO TO 250

220       WRITE (rec,FMT=9230)

9230      FORMAT (
     &      '** Meshing routine ran out of storage and used ad-hoc mesh'
     &           )

          GO TO 250

230       WRITE (rec,FMT=9240)

9240      FORMAT ('** Serious internal error: try different TOL')

          GO TO 250

240       WRITE (rec,FMT=9250)

9250      FORMAT ('** Cannot obtain required accuracy ')

          GO TO 250

      END IF

250   ifail = p01abf(ifo,ifail,srname,nrec,rec)
      RETURN

      END SUBROUTINE sl02fm

C ---------------------------------------------------------------------

      SUBROUTINE admesh(elam,ptemp,qtemp,wtemp,xtemp,rlog,theta,scale,
     &                  pp,qp,wp,xmesh,n,ntemp,tol,imatch,shift,istart,
     &                  iend,params,iparam,coeffn,icofun,ifail)
C The mesh is stored in the places xmesh(n),..,xmesh(istart) and the rest of
C the array contains no useful information. RLOG(0:ntemp), SCALE(0:ntemp),
C THETA(0:ntemp) and ELAM are available to represent the eigenfunction with
C which the adaptation of the mesh will be carried out.

C Initialising Stage:
C   Ifail:
C     .. Parameters ..
      DOUBLE PRECISION half,safe1,safe2,one,sixth,trd,twelf,zero,two,
     &                 tenth,fifth
      PARAMETER (half=5.D-1,safe1=9.D-1,safe2=8.D-1,one=1.D0,
     &          sixth=one/6.D0,trd=one/3.D0,twelf=half*sixth,zero=0.D0,
     &          two=2.D0,tenth=1.D-1,fifth=two*tenth)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,tol
      INTEGER icofun,iend,ifail,imatch,iparam,istart,n,ntemp
      LOGICAL shift
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam),pp(0:iend),ptemp(1:ntemp),
     &                 qp(1:iend),qtemp(1:ntemp),rlog(0:ntemp),
     &                 scale(0:ntemp),theta(0:ntemp),wp(1:iend),
     &                 wtemp(1:ntemp),xmesh(0:iend),xtemp(0:ntemp)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION errest,errpm,errpp,errqm,errqp,errw1,errw2,errwm,
     &                 errwp,gamma,h,omega1,omega2,omega3,p1,p2,p3,q1,
     &                 q2,q3,toloc,tolow,w1,w2,w3,weight,x1,x2,x3
      INTEGER i,idi,j,k,kntr
      LOGICAL phase1,repeat
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION pdyval(0:1),xvals(0:1),yvals(0:1)
C     ..
C     .. External Subroutines ..
cc      EXTERNAL coefun,efun
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min,sign,sqrt
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION hloc
cc      EXTERNAL hloc
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
      ifail = 0
      IF (iend.EQ.istart) RETURN
C   Local Tolerance:
      toloc = tol
      tolow = toloc*safe2
      phase1 = .true.
      repeat = .false.
      k = 0
      kntr = 0
      i = istart
      idi = 1
      IF (iend.LT.istart) idi = -1
C We use the array PP as storage for the new mesh being formed.
      x1 = xmesh(i)
      x3 = xmesh(i+idi)
      pp(i) = x1

10    h = x3 - x1
      x2 = half* (x1+x3)
      IF (repeat) THEN
          p1 = p3
          q1 = q3
          w1 = w3

      ELSE
          CALL coefun(x1,p1,q1,w1,coeffn,params,iparam,icofun)
      END IF

      CALL coefun(x2,p2,q2,w2,coeffn,params,iparam,icofun)
      CALL coefun(x3,p3,q3,w3,coeffn,params,iparam,icofun)

C Do error estimation:  ******************************************

      IF (x1.LE.x3) THEN
          errw1 = w1 - w2
          errw2 = w3 - w2
          errwm = elam*errw1
          errwp = elam*errw2
          errpp = (one/p3-one/p2)
          errpm = (one/p1-one/p2)
          errqp = (q3-q2)
          errqm = (q1-q2)
          xvals(0) = x1
          xvals(1) = x3

      ELSE
          errw1 = w3 - w2
          errw2 = w1 - w2
          errwp = elam*errw2
          errwm = elam*errw1
          errpm = (one/p3-one/p2)
          errpp = (one/p1-one/p2)
          errqm = (q3-q2)
          errqp = (q1-q2)
          xvals(0) = x3
          xvals(1) = x1
      END IF

      IF (phase1) THEN
          CALL efun(xvals,yvals,pdyval,elam,rlog,theta,scale,xtemp,
     &              ptemp,qtemp,wtemp,1,ntemp,imatch)

      ELSE

          IF (x1.LE.x3) THEN
              IF (.NOT.repeat) THEN
                  yvals(0) = yvals(1)
                  pdyval(0) = pdyval(1)
              END IF

              CALL efun(xvals(1),yvals(1),pdyval(1),elam,rlog,theta,
     &                  scale,xtemp,ptemp,qtemp,wtemp,0,ntemp,imatch)

          ELSE
              IF (.NOT.repeat) THEN
                  yvals(1) = yvals(0)
                  pdyval(1) = pdyval(0)
              END IF

              CALL efun(xvals(0),yvals(0),pdyval(0),elam,rlog,theta,
     &                  scale,xtemp,ptemp,qtemp,wtemp,0,ntemp,imatch)
          END IF

      END IF
C END of eigenfunction evaluation.

      errest = twelf*abs(h* (errw1*min(one,yvals(0)**2)+errw2*min(one,
     &         yvals(1)**2)))
      errwm = errwm - errqm
      errwp = errwp - errqp
C There are a number of other different error monitors to be used, depending
C on whether the behaviour of solutions of the differential equation.
C Compute Prufer radius for use:
      weight = half*sqrt(yvals(0)**2+pdyval(0)**2+yvals(1)**2+
     &         pdyval(1)**2)
C Measure rate of oscillation:
      omega2 = (elam*w2-q2)/p2
C WEIGHT may be quite small near the endpoints, where other eigenfunctions
C do not decay just as fast. Since weight is supposed to take account of these
C other eigenfunctions we adjust its value where appropriate:
      IF (omega2.LT.zero) THEN
          IF (q2.LT.two*max(one,abs(elam))*w2) THEN
C It is not clear that there is any representative rate of decay for nearby
C eigenfunctions: a small increase in elam could change everything.
              weight = max(weight,tenth)

          ELSE
C We can be reasonably confident that the rate of decay of nearby
C eigenfunctions will be at least half as fast as that of the present
C eigenfunction.
              gamma = sqrt((two*max(one,elam)*w2-q2)/ (elam*w2-q2))
              weight = max(weight,weight**gamma)
          END IF

      ELSE
          weight = max(weight,min(tenth,two*weight))
      END IF
C First error monitor, valid in regions which are not highly oscillatory.
C This controls the components of the eigenfunction which are not linearly
C dependent on the exact eigenfunction:
      errest = max(errest,sixth*weight*
     &         abs(h* (abs(yvals(0)*errwm+yvals(1)*errwp)/sqrt(max(one,
     &         omega2))+abs(errpm*pdyval(0)+errpp*pdyval(1)))))
C Second error monitor, uses Simpson's rule to measure a local contribution
C to the Green's formula for eigenvalue error:
      errest = max(errest,sixth*abs(h* ((yvals(0)**2)*errwm+
     &         (yvals(1)**2)*errwp+errpm* (pdyval(0)**2)+
     &         errpp* (pdyval(1)**2))))/max(one,abs(elam))
      IF ((h**2)*omega2.GT.half) THEN
C Safeguard in highly oscillatory regions: in such regions we can control
C eigenvalue error by controlling error in Prufer theta:
          omega1 = sqrt(abs(elam*w1-q1)/p1)
          omega3 = sqrt(abs(elam*w3-q3)/p3)
          omega2 = sqrt(omega2)
          errest = max(errest,abs(two* (p2/w2)*omega2* (omega1-
     &             two*omega2+omega3)*h)/max(one,elam))
      END IF
C All the error monitors here are O(h**3) for small h.

C END of error estimation. ********************************************

C Decide what to do next, depending on size of error:
      IF (errest.LE.toloc) THEN
C Check that the error is not TOO small
          IF (errest.LT.tolow .AND. k.EQ.0) THEN
              h = h/max(fifth, (errest/tolow)**trd)
              IF (phase1) THEN
                  repeat = .true.
                  kntr = kntr + 1
                  IF (kntr.LT.8) GO TO 10
              END IF

          END IF
C        the error test has been passed; move to next step if there is one.
          repeat = .false.
          phase1 = .false.
          i = i + idi
          pp(i) = x3
          IF ((xmesh(n)-x3)*idi.LE.zero) THEN
C The end of the range has been reached
              GO TO 20

          END IF
C If we are here then it means that the meshing must continue.
          x1 = x3
          x3 = x3 + sign(one,h)*min(abs(h),hloc(x3,xmesh,n))
          IF ((xmesh(n)-x3)*idi.LE.zero) THEN
              x3 = xmesh(n)
          END IF

          k = 0
          GO TO 10

      ELSE
          repeat = .true.
C        the error test has been failed; reduce the stepsize and repeat step.
C        Error traps:
          IF ((iend-i)*idi.LE.zero) THEN
C         Run out of space
              ifail = 11
              GO TO 40

          END IF

          k = k + 1
          IF (k.GE.8) THEN
C         Coefficients are too nasty
              ifail = 10
              GO TO 40

          END IF
C        END of error traps
          h = h*min(safe1, (toloc/errest)**trd)
          x3 = x1 + h
          GO TO 10

      END IF

20    CONTINUE

C We now have an adapted mesh, which is stored in the array spaces
C pp(istart),..,pp(i). We copy it into the corresponding slots in the array
C xmesh.
      iend = i
      DO 30 j = istart,iend,idi
          xmesh(j) = pp(j)
30    CONTINUE
40    RETURN

      END SUBROUTINE admesh

C ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION HLOC(X,XMESH,N)
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE
      PARAMETER (HALF=5.D-1,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XMESH(0:N)
C     ..
C     .. Local Scalars ..
      INTEGER J
      DOUBLE PRECISION XHI,XLO,HHI,HLO,S
C     ..

      IF (X.LE.HALF*(XMESH(0)+XMESH(1))) THEN
          HLOC = XMESH(1)-XMESH(0)
      ELSE IF (X.GE.HALF*(XMESH(N-1)+XMESH(N))) THEN
          HLOC = XMESH(N)-XMESH(N-1)
      ELSE
          DO 10 J = 1,N-1
             XLO = HALF*(XMESH(J-1)+XMESH(J))
             XHI = HALF*(XMESH(J)+XMESH(J+1))
             IF (X.GE.XLO.AND.X.LT.XHI) THEN
                HHI = XMESH(J+1)-XMESH(J)
                HLO = XMESH(J)-XMESH(J-1)
                S = (X-XLO)/(XHI-XLO)
                HLOC = S*HHI + (ONE-S)*HLO 
                RETURN
             END IF
10        CONTINUE
      END IF
      RETURN

      END FUNCTION HLOC
C ---------------------------------------------------------------------
C --------------------- Source from SNEW  -----------------------------
C ---------------------------------------------------------------------

      SUBROUTINE eigen(elam0,elam,eps,a,b,k,n,ntemp,imatch,np,ising,
     &                 isymm,lknt,maxit,coeffn,setup,pp,qp,wp,xmesh,
     &                 ptemp,qtemp,wtemp,xtemp,rlog,theta,scale,tol,
     &                 noxtrp,params,iparam,icofun,ifail)
C     .. Parameters ..
      DOUBLE PRECISION three,six,smatch,one,zero,two,half,tqtr,p9
      PARAMETER (three=3.D0,six=6.D0,smatch=1.D0,one=1.D0,zero=0.D0,
     &          two=2.D0,half=5.D-1,tqtr=three/4.D0,p9=9.D-1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,elam,elam0,eps,tol
      INTEGER icofun,ifail,imatch,iparam,ising,isymm,k,lknt,maxit,n,np,
     &        ntemp
      LOGICAL noxtrp
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam),pp(1:n),ptemp(1:ntemp),qp(1:n),
     &                 qtemp(1:ntemp),rlog(0:ntemp),scale(0:ntemp),
     &                 theta(0:ntemp),wp(1:n),wtemp(1:ntemp),xmesh(0:n),
     &                 xtemp(0:ntemp)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn,setup
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION atrunc,btrunc,diff,dtha,dthb,elamo,epso,fynyt,
     &                 ovflow,pc,pi,qc,rat1,rat2,rat3,rmatch,sc,shift,
     &                 store,th,tolmsh,toloc,wc,xc
      INTEGER i,icfn,idmmy,ifc,ifo,imesh,iq,iref,isingo,itime,itw,
     &        itwm1,knt,lknto,newmsh,nmesh,npo
      LOGICAL trans
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION dtheta,x01aaf,x02ahf,x02alf
cc      EXTERNAL dtheta,x01aaf,x02ahf,x02alf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL check,coarse,coefun,fillin,initia,mmesh,norfun,solve
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,max,min,nint,sqrt
C     ..
C Set a big number to approximate infinity:
      ovflow = sqrt(x02alf(1.D0))
C     ..
C Check input for errors (and change in certain cases)
C
      IF (abs(ifail).GT.1) THEN
          ifail = -abs(ifail)
          GO TO 150

      END IF
C
      IF (n.LT.10 .OR. tol.LE.0.0D0 .OR. eps.LE.0.0D0 .OR. b.LE.a) THEN
          ifail = 1
          GO TO 150

      END IF
C
C End of input checks.
C
      pi = x01aaf(pi)
C
C Store the input values of certain parameters which may be changed by the
C code:
C
      ifo = ifail
      lknto = lknt
      isingo = ising
      elamo = elam0
      epso = eps
      ifail = 0
C trans determines when we are to perform transformation of independent
C variables:
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      trans = icofun .EQ. 1 .OR. icofun .EQ. 3
C
C Set default values on certain parameters:
C
      IF (maxit.EQ.0) maxit = 100
      IF (lknt.EQ.0) lknt = 60
      IF (np.LE.5) np = 10
C np is stored so that it may be used to control the maximum mesh size:
      npo = np
C The paramter itime is needed to count how many ends have been truncated if
C the problem is singular:
      itime = 0
C The parameter iref is needed to count how many times we have tried to
C find a mesh to give decent accuracy, without success:
      iref = 0
C
C Set the match-point integer for the computation of the initial approximation
C to the eigenvalue.
C
      iq = np/2
C
C Set up artificial end-points according to the value of the ISING parameter
C
      IF (ising.EQ.0) THEN
C The problem is regular
          atrunc = a
          btrunc = b

      ELSE IF (ising.EQ.1) THEN
C Singular point at x=a; truncate:
          atrunc = a + (b-a)/six
          btrunc = b

      ELSE IF (ising.EQ.2) THEN
C Singular point at x=b; truncate:
          atrunc = a
          btrunc = b - (b-a)/six

      ELSE
C Both ends are singular; truncate:
          atrunc = a + (b-a)/six
          btrunc = b - (b-a)/six
      END IF
C Set the initial step size for approaching a singular endpoint:
      shift = (btrunc-atrunc)/np

C Compute initial approximation with truncated end-points:
C 1: Set up an initial mesh of np points:
      CALL coarse(np,atrunc,btrunc,xmesh,0)
C 2: Set up the arrays of coefficient function values at the mesh centres:
      CALL initia(coeffn,pp,qp,wp,xmesh,np,params,iparam,icofun)
C 3: Check the array PP to see that p is an admissible coefficient function (>0)
      ifc = 0
      CALL check(np,pp,ifc)
      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF
C
C 4:  Solve the discrete S-L problem on the small mesh of np points to get
C   initial approximation to the eigenvalue:
C
      toloc = tol/10.D0
      CALL solve(elam0,eps,k,np,iq,maxit,setup,pp,qp,wp,xmesh,trans,
     &           toloc,ifail)
C Special treatment of regular problems: do not need to change the
C end-points inview of the estimate obtained. Note also that no failures
C are tolerated for regular problems.
      IF (ising.EQ.0) THEN
          nmesh = n/2
          IF (ifail.NE.0) GO TO 150
C         Go straight to the final two sweeps on adapted meshes, with
C         Richardson extrapolation.
          GO TO 50

      END IF

C ********************************************************************
C **Only singular problems from here until the statement labelled 50**
C ********************************************************************

C
      IF (ifail.EQ.7) THEN
C This can happen if we truncate too far in and so we reset the ifail and
C elam0 to their old values and choose new truncated end-points.
          ifail = 0
          elam0 = elamo
      END IF
C Other types of failure are not to be tolerated!
      IF (ifail.NE.0) GO TO 150
C We now have an initial estimate of the eigenvalue, from which
C we can decide how to truncate the interval.

      IF (ising.EQ.3) THEN
C We pretend that the problem is only singular at x=a and deal with truncation
C at that end first. We set itime=1 to denote the fact that we are dealing with
C the first of two singular points in this case.
          itime = 1
          ising = 1
      END IF
C We can now start the sequence of computations which leads to the
C choice of the truncated endpoint(s). knt counts the number of attempts
C to find each suitable endpoint, and must not exceed lknt.
      knt = 1
C Add an extra point closer to the singular point under consideration:
20    np = np + 1
      iq = np/2
      IF (ising.EQ.1) THEN
C Add the point between atrunc and a:
          DO 30 i = np,2,-1
              pp(i) = pp(i-1)
              qp(i) = qp(i-1)
              wp(i) = wp(i-1)
              xmesh(i) = xmesh(i-1)
30        CONTINUE
          xmesh(1) = xmesh(0)
          xmesh(0) = xmesh(0) - min(shift,tqtr* (xmesh(0)-a))
          atrunc = xmesh(0)
          xc = xmesh(0) + half* (xmesh(1)-xmesh(0))
          CALL coefun(xc,pc,qc,wc,coeffn,params,iparam,icofun)
          pp(1) = pc
          qp(1) = qc
          wp(1) = wc
      END IF

      IF (isymm.EQ.1) THEN
C This is a special case: the interval has two singular endpoints and must be
C truncated symmetrically, so we must add an extra point at the other end too:
          np = np + 1
          iq = np/2
      END IF

      IF (ising.EQ.2 .OR. isymm.EQ.1) THEN
C Add in an extra mesh-point between btrunc and b. In the case isymm = 1 this
C is done purely to preserve symmetry.
          xmesh(np) = xmesh(np-1) + min(shift,tqtr* (b-xmesh(np-1)))
          btrunc = xmesh(np)
          xc = xmesh(np-1) + half* (xmesh(np)-xmesh(np-1))
          CALL coefun(xc,pc,qc,wc,coeffn,params,iparam,icofun)
          pp(np) = pc
          qp(np) = qc
          wp(np) = wc
      END IF
C Check that p(x) never vanishes or becomes negative:
      ifc = 0
      CALL check(np,pp,ifc)
      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF
C Extra point or points have been added; compute the eigenvalue and see if
C the the extra point(s) have made it change by much (store the old value in
C elam0 for comparison)
      elam = elam0
      CALL solve(elam,eps,k,np,iq,maxit,setup,pp,qp,wp,xmesh,trans,
     &           toloc,ifail)

C      WRITE (6,FMT=9000) xmesh(0),xmesh(np),elam,ifail

      IF (ifail.EQ.7) THEN
C This failure can be caused by truncating too much when certain asymptotic
C forms are used for the boundary conditions. We must see at which end of the
C interval the eigenfunction is oscillating most rapidly, and concentrate on
C that end provided this will not risk overflow.
          rat1 = (elam*wp(1)-qp(1))/pp(1)
          rat2 = (elam*wp(np)-qp(np))/pp(np)
C Check for danger of overflow:
          rat3 = max(rat1,rat2)
          IF ((rat3.GE.0.D0.AND.sqrt(abs(rat3)).GE.1.0D-2*x02ahf(0.D0)) 
     &        .OR.max(abs(rat1),abs(rat2)).GE.ovflow) THEN
              a = xmesh(0)
              b = xmesh(np)
              GO TO 150

          END IF
C End of overflow check. Now select the end of the interval with most rapid
C oscillations. Note that if the problem has only one singular end then
C itime = 0 so there is no danger of treating a regular end as singular.
          IF (itime.EQ.1 .AND. isymm.NE.1) THEN
              IF (rat2.GT.rat1) THEN
                  ising = 2

              ELSE
                  ising = 1
              END IF

          END IF

          ifail = 0
          elam = elam0
          diff = 10.0D0*tol

          GO TO 40
C End of special treatment for IFAIL = 7 failures.
      END IF

C Other failures (ifail \neq 7) are not tolerated!
      IF (ifail.NE.0) GO TO 150

      diff = abs(elam-elam0)/max(one,abs(elam))
C There is no danger in setting eps = diff; for if diff=0 then we will
C not be calling SOLVE without resetting eps to its original value epso.
      eps = diff
40    knt = knt + 1
      elam0 = elam
      IF (knt.GE.lknt) THEN
C Cannot find a suitable end-point in lknt steps
          ifail = 8
          a = xmesh(0)
          b = xmesh(np)
          GO TO 150

      END IF
C Less than lknt steps taken, but endpoints not yet satisfactory:
      IF (diff.GT.tol) THEN
C Check to see if endpoint(s) can be moved without danger of overflow.
          fynyt = (elam*wp(1)-qp(1))/pp(1)
          IF ((fynyt.GE.0.D0.AND.sqrt(abs(fynyt)).GE.
     &        1.0D-2*x02ahf(0.D0)) .OR. abs(fynyt).GE.ovflow) THEN
              ifail = 5
              a = xmesh(0)
              b = xmesh(np)
              GO TO 150

          END IF

          fynyt = (elam*wp(np)-qp(np))/pp(np)
          IF ((fynyt.GE.0.D0.AND.sqrt(abs(fynyt)).GE.
     &        1.0D-2*x02ahf(0.D0)) .OR. abs(fynyt).GE.ovflow) THEN
              ifail = 5
              a = xmesh(0)
              b = xmesh(np)
              GO TO 150

          END IF

C If we get through the last loop then we can indeed move an endpoint. Loop
C back to 20 to do this.
          GO TO 20

      END IF
C
C If we are doing a problem with TWO SINGULAR POINTS and itime=1 then only the
C end of the interval where the eigenfunction oscillates most rapidly has been
C dealt with so far.
      IF (itime.EQ.1 .AND. isymm.NE.1) THEN
          IF (ising.EQ.1) THEN
              ising = 2

          ELSE
              ising = 1
          END IF

          elam0 = elam
C See, I told you we would reset eps!
          eps = epso
          knt = 1
          itime = 2
C Loop back to start truncation procedure for the other end.
          GO TO 20

      END IF

C *******************************************************************
C ****   End of interval truncation for singular problems ***********
C *******************************************************************

C We now have some approximations to the eigenvalue, as well as
C suitable truncated endpoints.
      atrunc = xmesh(0)
      btrunc = xmesh(np)
C We are now ready to solve our problem on the truncated interval. We will
C compute two approximations in order to get an error estimate.
      elam0 = elam
      eps = epso
C Our first fine mesh will have only half as many points as our eventual mesh
      nmesh = n/2
C Find a suitable matching point (it is at this stage that regular problems
C come back under consideration -- see the GO TO 50 above).
50    IF (isymm.EQ.1) THEN
          rmatch = half* (atrunc+btrunc)

      ELSE
          rmatch = xmesh(iq)
          store = elam0*wp(iq) - qp(iq)
          DO 60 idmmy = 1,np - 1
              IF ((elam0*wp(idmmy)-qp(idmmy)).GT.store) THEN
                  store = elam0*wp(idmmy) - qp(idmmy)
                  rmatch = xmesh(idmmy)
              END IF

60        CONTINUE
      END IF
C
C In order to form the mesh for the final two sweeps, we need to know a
C decent approximation to the eigenfunction, which we get by using a coarse
C mesh of at most ntemp or 100 points.
70    IF (np.LT.min(n/2,ntemp/2,50)) THEN
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
          IF (icofun.EQ.1 .OR. icofun.EQ.3) THEN
              icfn = icofun + 1
              DO 80 i = np,1,-1
                  itw = 2*i
                  itwm1 = itw - 1
                  xmesh(itw) = xmesh(i)
                  xmesh(itwm1) = half* (xmesh(i)+xmesh(i-1))
80            CONTINUE
              np = 2*np
              GO TO 70

          ELSE
              CALL fillin(np,2*np,xmesh,pp,qp,wp,coeffn,params,iparam,
     &                    icofun)
              np = 2*np
              GO TO 70

          END IF

      END IF

C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      IF (icofun.EQ.1 .OR. icofun.EQ.3) THEN
          xmesh(0) = xmesh(0)/ (one-abs(xmesh(0)))
          icfn = icofun + 1
C From now on we abandon transformed coordinates:
          trans = .false.
          DO 90 i = 1,np
              xmesh(i) = xmesh(i)/ (one-abs(xmesh(i)))
              CALL coefun(half* (xmesh(i-1)+xmesh(i)),pp(i),qp(i),wp(i),
     &                    coeffn,params,iparam,icfn)
90        CONTINUE
      END IF

      xtemp(0) = xmesh(0)
      DO 100 i = 1,np
          xtemp(i) = xmesh(i)
          ptemp(i) = pp(i)
          qtemp(i) = qp(i)
          wtemp(i) = wp(i)
100   CONTINUE

      ifc = 0
      CALL check(np,ptemp,ifc)
      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF

      ising = isingo
      eps = epso
      iq = np/2
      toloc = 0.D0
      CALL solve(elam0,eps,k,np,iq,maxit,setup,ptemp,qtemp,wtemp,xtemp,
     &           trans,toloc,ifc)
C      WRITE (6,FMT=*) elam0,np,ifc

      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF

      CALL norfun(rlog,theta,scale,smatch,elam0,zero,k,np,iq,ptemp,
     &            qtemp,wtemp,xtemp,setup,trans,ifc)
      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF

C Check all singular problems to make sure we have not strayed into the
C continuous spectrum.
      IF (isingo.NE.0) THEN
          dtha = zero
          dthb = zero
          IF (isingo.GT.1) THEN
              th = theta(np)
              sc = scale(np)
              dthb = dtheta(th,sc,elam0,btrunc,b,coeffn,params,iparam,
     &               icofun)
          END IF

          IF (isingo.NE.2) THEN
              th = theta(0)
              sc = scale(0)
              dtha = dtheta(th,sc,elam0,atrunc,a,coeffn,params,iparam,
     &               icofun)
          END IF

          IF ((dtha+dthb).GE.two*pi) THEN
              ifail = 4
              ising = isingo
              lknt = lknto
              a = atrunc
              b = btrunc
              RETURN

          END IF

      END IF
C End of check for continuous spectrum.

C We can now choose the new mesh using an appropriate equidistribution
C process. This is done by calling the routine MMESH.
C Note the special status of the IFC parameter for this routine. On
C entry, it equals the value which IFAIL had on entry to EIGEN, and
C if that value was 1, then an emergency mesh will be used if the
C routine runs out of space. If this happens then ifc will be set to
C 12 on exit.
C
      imesh = nmesh
      IF (isymm.EQ.1) imesh = (imesh-1)/2
      ifc = ifo
C Store the eigenvalue corresponding to the eigenfunction which is
C used for meshing; this will be useful later.
      elamo = elam0
C
      tolmsh = tol
110   CALL mmesh(rmatch,atrunc,0,elam0,coeffn,rlog,theta,scale,xtemp,
     &          ptemp,qtemp,wtemp,np,iq,xmesh,pp,qp,wp,imesh,tolmsh,npo,
     &          params,iparam,icofun,ifc)

C
C Error handling:
C
      IF (ifc.NE.0) THEN
C The routine cannot produce a mesh to meet the users tolerance
C within the allocated storage, or has ground to a halt because the
C coefficient functions are too pathological
C
          ifail = ifc
          IF (ifail.NE.12) GO TO 150
      END IF
C
C End of error handling; produce second half of mesh;
C in symmetric case we just need to reflect the mesh through the midpoint:
C
      IF (isymm.EQ.1) THEN
C Just reflect the mesh through its right-hand endpoint if the problem
C is symmetric.
          nmesh = 2*imesh
          xmesh(nmesh) = xmesh(0)
          DO 120 i = 1,imesh
              xmesh(imesh+i) = two*rmatch - xmesh(imesh-i)
              pp(imesh+i) = pp(imesh-i+1)
              qp(imesh+i) = qp(imesh-i+1)
              wp(imesh+i) = wp(imesh-i+1)
120       CONTINUE

      ELSE
C End of symmetric case.
C
          ifc = ifo
          CALL mmesh(rmatch,btrunc,imesh,elam0,coeffn,rlog,theta,scale,
     &              xtemp,ptemp,qtemp,wtemp,np,iq,xmesh,pp,qp,wp,nmesh,
     &              tolmsh,npo,params,iparam,icofun,ifc)
C
C Error handling:
C
          IF (ifc.NE.0) THEN
C The routine cannot produce a mesh to meet the users tolerance
C within the allocated storage, or has ground to a halt because the
C coefficient functions are too pathological
C
              ifail = ifc
              IF (ifail.NE.12) GO TO 150
          END IF
C
C End of error handling
C if no error has occured then nmesh is now set to a new value which
C is equal to the actual number of intervals in the computed mesh
C and is less than or equal to its input value.
C
      END IF
C We now have a mesh which can be stored within the user-specified
C storage, and which is adapted in a suitable manner to help achieve
C the user-specified tolerance. we compute the first approxmation to
C the eigenvalue by using this mesh.
C
C Check that the array pp has no zero or negative entries:
      ifc = 0
      CALL check(nmesh,pp,ifc)
      IF (ifc.NE.0) THEN
          ifail = ifc
          GO TO 150

      END IF

C
C We are now ready for the first of the two final computations. IFAIL now
C has one of two values 0 or 12, the second (if it arises) coming from the
C use of an emergency mesh.
      ifc = 0
      toloc = tol/100.D0
      CALL solve(elam0,eps,k,nmesh,imesh,maxit,setup,pp,qp,wp,xmesh,
     &           trans,toloc,ifc)
      eps = epso/1.D2
C
C Store the eigenvalue estimate for extrapolation later
C
      elam = elam0
C
C Error handling:
      IF (ifc.NE.0) THEN
          IF (ifail.EQ.12) GO TO 130
          ifail = ifc
          GO TO 150

      END IF

      IF (noxtrp) THEN
          a = atrunc
          b = btrunc
          n = nmesh
          imatch = imesh
          ising = isingo
          lknt = lknto
          RETURN

      END IF

C Now we halve the mesh-size and repeat the process to obtain a
C better eigenvalue approximation. the mesh-size is halved by
C adding in all the midpoints of intervals in the current mesh.
C the arrays pp,qp,wp also must be updated, and all this is done
C by the routine fillin. storage is always adequate since, if it
C were not, this would have been detected by the mesh routine above.
130   newmsh = 2*nmesh
C Reset match-point index:
      imesh = 2*imesh
      icfn = icofun
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      IF (icofun.EQ.1) THEN
          icfn = 2

      ELSE IF (icofun.EQ.3) THEN
          icfn = 4
      END IF

      CALL fillin(nmesh,newmsh,xmesh,pp,qp,wp,coeffn,params,iparam,icfn)
C      eps = epso/100.d0
      ifc = 0
      toloc = 0.D0
      CALL solve(elam0,eps,k,newmsh,imesh,maxit,setup,pp,qp,wp,xmesh,
     &           trans,toloc,ifc)
C
C Error handling:
      IF (ifc.NE.0) THEN
C          ic = imesh
          n = newmsh
          a = atrunc
          b = btrunc
          IF (ifail.GT.1) THEN
              GO TO 150

          ELSE
              ifail = ifc
              GO TO 150

          END IF

      END IF
C
C Error estimate:
      eps = (elam-elam0)/three
      elam0 = elam0 - eps
      IF (abs(eps)/max(one,abs(elam0)).GT.tol** (two/three) .AND.
     &    ifail.EQ.0) THEN
C The accuracy achieved is unacceptable. Go back to the meshing stage with
C a reduced minimum mesh size.
          IF (iref.GT.1) THEN
C           We have already done this too many times. The problem is simply
C           too nasty. Fail with appropriate error flag
              ifail = 14
              GO TO 140

          END IF

          iref = iref + 1
          rat1 = sqrt(abs(max(one,abs(elam0))* (tol** (two/three))/eps))
          npo = nint(dble(npo)/ (tqtr*rat1))
          tolmsh = tolmsh*min(p9,rat1** (three/two))
          elam0 = elamo
          nmesh = n/2
          imesh = nmesh
          IF (isymm.EQ.1) imesh = (imesh-1)/2
          GO TO 110

      END IF

C Use elam0 to store the eigenvalue corresponding to the coarse-mesh
C eigenfunction which was used for automatic meshing:
140   elam = elamo
C
C Prepare to return to the user; make sure a and b have the same values
C as on entry, not to mention ising. also we set n = newmesh so that the
C user can see how many mesh-points were used.

      lknt = lknto
      ising = isingo
      n = newmsh
      imatch = imesh
      a = atrunc
      b = btrunc

C     write(6,*) 'np before exiting EIGEN: ',np

C Either an emergency mesh has been used (IFAIL = 12), or else the accuracy
C achieved is inadequate (IFAIL = 14), or else the routine has completed
C its run successfully.
      IF (ifail.NE.12 .AND. ifail.NE.14) ifail = 0
150   RETURN
C9000  FORMAT (' ALFA: ',G18.12,' BETA: ',G18.12,' ELAM: ',D10.4,
C     &       ' IFAIL ',I2)

      END SUBROUTINE eigen
C--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION dtheta(theta,scale,elam,bt,b,coeffn,
     &                 params,iparam,icofun)
C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=5.D-1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION b,bt,elam,scale,theta
      INTEGER icofun,iparam
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION di,ovflow,p,q,qbig,s,th,w,x,xend,xo
      INTEGER icount
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02alf
cc      EXTERNAL x02alf
C     ..
C     .. External Subroutines ..
c      EXTERNAL coefun,onestp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sign,sqrt
C     ..
      ovflow = sqrt(x02alf(1.D0))
      th = theta
      s = scale
      xo = bt
      xend = half* (b+bt)
      dtheta = 0.D0
      di = sign(1.D0,b-bt)
      icount = 0

10    x = half* (xo+xend)
      CALL coefun(x,p,q,w,coeffn,params,iparam,icofun)
      qbig = elam*w - q
      CALL onestp(xo,xend,di,s,th,p,qbig)
      dtheta = (th-theta)*di
      icount = icount + 1
      IF (icount.GE.5) RETURN
      IF (abs(qbig).GE.ovflow*p) RETURN
      xo = xend
      xend = half* (xend+b)
      GO TO 10

      END FUNCTION dtheta

C -------------------------------------------------------------------------

      SUBROUTINE solve(elam,eps,k,n,ic,maxit,setup,pp,qp,wp,xmesh,trans,
     &                 tol,ifail)
C     .. Parameters ..
      DOUBLE PRECISION ten
      PARAMETER (ten=1.D1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,eps,tol
      INTEGER ic,ifail,k,maxit,n
      LOGICAL trans
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(n),qp(n),wp(n),xmesh(0:n)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. Scalars in Common ..
      INTEGER icall,int,ip,ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION da,db,fx,pi,pik,rl1,rl2,tight,xs,ys
      INTEGER ifd,iflag,ind,ir,nevs
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION c(17)
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION d,x01aaf,x02ajf
cc      EXTERNAL d,x01aaf,x02ajf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL c05azf,intval
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC max
C     ..
      pi = x01aaf(pi)
      pik = k*pi
      tight = max(x02ajf(0.D0)*ten,tol)
      ifail = 0
      iflag = 0

C     WE NOW SEEK TO FIND AN INTERVAL [RL1,RL2] CONTAINING THE
C     ZERO OF THE FUNCTION D(RL)-K*PI

      int = ip
      CALL intval(elam,rl1,rl2,eps,da,db,k,n,ic,pp,qp,wp,xmesh,setup,
     &            trans,iflag)
      int = ip - int

C     CHECK THAT THE ERROR FLAG HAS NOT BEEN CHANGED

      IF (iflag.NE.0) THEN
          ifail = iflag
          RETURN

      END IF

C     WE NOW HAVE AN INTERVAL [A,B] WITH THE PROPERTY THAT
C     D(A) <= 0 AND D(B) >= 0

      xs = rl1
      ys = rl2

      ind = -1
      ir = 0
      nevs = 0
      fx = da
      c(1) = db
      iflag = 1

20    CALL c05azf(xs,ys,fx,tight,ir,c,ind,iflag)

      IF (ind.EQ.0) THEN
          IF (iflag.NE.0 .AND. iflag.NE.5) ifail = 13
          GO TO 40

      END IF

      IF (ind.LT.2 .OR. ind.GT.4) THEN
          ifail = 13
          RETURN

      END IF

      ifd = 0
      fx = d(xs,n,ic,pp,qp,wp,xmesh,setup,trans,ifd) - pik
      IF (ifd.NE.0) THEN
          ifail = ifd
          RETURN

      END IF

      nevs = nevs + 1

      IF (nevs.GT.maxit) GO TO 30

      GO TO 20

30    ifail = 6
      RETURN

40    elam = xs

      RETURN

      END SUBROUTINE solve

C---------------------------------------------------------------------

      SUBROUTINE fillin(nmesh,newmsh,xmesh,pp,qp,wp,coeffn,params,
     &                  iparam,icofun)
C     HALVES THE SIZE OF A MESH.
C     .. Scalar Arguments ..
      INTEGER icofun,iparam,newmsh,nmesh
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam),pp(newmsh),qp(newmsh),
     &                 wp(newmsh),xmesh(0:newmsh)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION pc,qc,wc,xc
      INTEGER i,itw,itwm1
C     ..
C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=5.D-1)
C     ..
C     .. External Subroutines ..
cc      EXTERNAL coefun
C     ..
      DO 10 i = nmesh,1,-1
          itw = 2*i
          itwm1 = itw - 1
          xmesh(itw) = xmesh(i)
          xmesh(itwm1) = half* (xmesh(itw)+xmesh(i-1))
          xc = half* (xmesh(itw)+xmesh(itwm1))
          CALL coefun(xc,pc,qc,wc,coeffn,params,iparam,icofun)
          pp(itw) = pc
          qp(itw) = qc
          wp(itw) = wc
          xc = half* (xmesh(itwm1)+xmesh(i-1))
          CALL coefun(xc,pc,qc,wc,coeffn,params,iparam,icofun)
          pp(itwm1) = pc
          qp(itwm1) = qc
          wp(itwm1) = wc
10    CONTINUE

      RETURN

      END SUBROUTINE fillin
C---------------------------------------------------------------------

      SUBROUTINE coarse(n,a,b,xmesh,ising)

C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
      INTEGER ising,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION xmesh(0:n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION h
      INTEGER i
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble
C     ..
      xmesh(0) = a
      xmesh(n) = b

      IF (ising.EQ.0) THEN
C     THE PROBLEM IS REGULAR, USE A UNIFORM MESH

          h = (b-a)/dble(n)
          DO 10 i = 1,n - 1
              xmesh(i) = xmesh(i-1) + h
10        CONTINUE

      ELSE IF (ising.EQ.1) THEN
C     X = A IS SINGULAR BUT X = B IS NOT

          h = (b-a)/dble(n-4)
          xmesh(1) = xmesh(0) + h/16.0D0
          xmesh(2) = xmesh(1) + h/16.0D0
          xmesh(3) = xmesh(2) + h/8.0D0
          xmesh(4) = xmesh(3) + h/4.0D0
          xmesh(5) = xmesh(4) + h/2.0D0
          DO 20 i = 6,n - 1
              xmesh(i) = xmesh(i-1) + h
20        CONTINUE

      ELSE IF (ising.EQ.2) THEN
C     X=A IS REGULAR BUT X = B IS SINGULAR

          h = (b-a)/dble(n-4)

          xmesh(n-1) = xmesh(n) - h/16.0D0
          xmesh(n-2) = xmesh(n-1) - h/16.0D0
          xmesh(n-3) = xmesh(n-2) - h/8.0D0
          xmesh(n-4) = xmesh(n-3) - h/4.0D0
          xmesh(n-5) = xmesh(n-4) - h/2.0D0
          DO 30 i = n - 6,1,-1
              xmesh(i) = xmesh(i+1) - h
30        CONTINUE

      ELSE IF (ising.EQ.3) THEN
C     BOTH ENDS OF THE INTERVAL ARE SINGULAR POINTS.

          h = (b-a)/dble(n-8)

          xmesh(1) = xmesh(0) + h/16.0D0
          xmesh(2) = xmesh(1) + h/16.0D0
          xmesh(3) = xmesh(2) + h/8.0D0
          xmesh(4) = xmesh(3) + h/4.0D0
          xmesh(5) = xmesh(4) + h/2.0D0

          xmesh(n-1) = xmesh(n) - h/16.0D0
          xmesh(n-2) = xmesh(n-1) - h/16.0D0
          xmesh(n-3) = xmesh(n-2) - h/8.0D0
          xmesh(n-4) = xmesh(n-3) - h/4.0D0
          xmesh(n-5) = xmesh(n-4) - h/2.0D0

          DO 40 i = 6,n - 6
              xmesh(i) = xmesh(i-1) + h
40        CONTINUE

      ELSE
C     THERE IS AN ERROR IN THE VALUE OF ISING
      END IF

      RETURN

      END SUBROUTINE coarse
C---------------------------------------------------------------------
      SUBROUTINE initia(coeffn,pp,qp,wp,xmesh,n,params,iparam,icofun)
C     .. Scalar Arguments ..
      INTEGER icofun,iparam,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam),pp(n),qp(n),wp(n),xmesh(0:n)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
C     .. Scalars in Common ..
      INTEGER icall,int,ip,ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION pc,qc,wc,xc
      INTEGER i
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival
C     ..
C     .. External Subroutines ..
cc      EXTERNAL coefun
C     ..
      DO 10 i = 1,n
          xc = 0.5D0* (xmesh(i)+xmesh(i-1))
          CALL coefun(xc,pc,qc,wc,coeffn,params,iparam,icofun)
          pp(i) = pc
          qp(i) = qc
          wp(i) = wc
10    CONTINUE

      RETURN

      END SUBROUTINE initia
C---------------------------------------------------------------------

      SUBROUTINE intval(elc,ela,elb,eps,fa,fb,k,n,ic,pp,qp,wp,xmesh,
     &                  setup,trans,iflag)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ela,elb,elc,eps,fa,fb
      INTEGER ic,iflag,k,n
      LOGICAL trans
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(n),qp(n),wp(n),xmesh(0:n)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION h,pi,sub
      INTEGER ifail,inits,nits
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION d,x01aaf
cc      EXTERNAL d,x01aaf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..

      pi = x01aaf(pi)
      sub = k*pi

      ela = elc - abs(eps)
      elb = elc + abs(eps)
      nits = 0
      inits = 0
      ifail = 0
      iflag = 0

      fa = d(ela,n,ic,pp,qp,wp,xmesh,setup,trans,ifail) - sub
      IF (ifail.NE.0) GO TO 50
      fb = d(elb,n,ic,pp,qp,wp,xmesh,setup,trans,ifail) - sub
      IF (ifail.NE.0) GO TO 50

10    h = elb - ela
      IF (nits.GT.50) THEN
          GO TO 40

      END IF

      IF (fa.GT.0.0D0 .AND. fb.LT.0.0D0) THEN
          ela = elb
          elb = ela + 2.D0*h
          fa = fb
          fb = d(elb,n,ic,pp,qp,wp,xmesh,setup,trans,ifail) - sub
          IF (ifail.NE.0) GO TO 50
          inits = inits + 1
          IF (inits.GT.4) GO TO 30
          nits = nits + 1
          GO TO 10

      END IF

      IF (fa*fb.LE.0.0D0) GO TO 20

      IF (fb.LT.0.0D0) THEN
          ela = elb
          fa = fb
          elb = elb + 2.0D0*h
          fb = d(elb,n,ic,pp,qp,wp,xmesh,setup,trans,ifail) - sub
          IF (ifail.NE.0) GO TO 50
          nits = nits + 1
          GO TO 10

      END IF

      IF (fa.GT.0.0D0) THEN
          elb = ela
          fb = fa
          ela = ela - 2.0D0*h
          fa = d(ela,n,ic,pp,qp,wp,xmesh,setup,trans,ifail) - sub
          IF (ifail.NE.0) GO TO 50
          nits = nits + 1
          GO TO 10

      END IF

20    iflag = 0
      RETURN

30    iflag = 9
      RETURN

40    iflag = 7
      RETURN

50    iflag = ifail
      RETURN

      END SUBROUTINE intval

C---------------------------------------------------------------

      SUBROUTINE nrmlse(yl,pdyl,ifail)
C     .. Scalar Arguments ..
      DOUBLE PRECISION pdyl,yl
      INTEGER ifail
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION s
C     ..
      s = yl**2 + pdyl**2
      IF (s.EQ.0.0D0) THEN
          ifail = 3
          RETURN

      END IF

      IF (yl.LT.0.0D0) THEN
          yl = -yl
          pdyl = -pdyl
      END IF

      IF (yl.EQ.0.0D0 .AND. pdyl.LT.0.0D0) THEN
          pdyl = -pdyl
      END IF

      ifail = 0
      RETURN

      END SUBROUTINE nrmlse

C---------------------------------------------------------------

      SUBROUTINE check(n,pp,ifail)
C THIS ROUTINE CHECKS THE INPUT TO ROUTINE D TO ENSURE THAT IT MAKES SENSE.
C     .. Scalar Arguments ..
      INTEGER ifail,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION rmini
      INTEGER i
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC min
C     ..
      ifail = 0
      rmini = 1.0D0
      DO 10 i = 1,n
          rmini = min(rmini,pp(i))
10    CONTINUE
      IF (rmini.LE.0.D0) ifail = 2
      RETURN

      END SUBROUTINE check

C -----------------------------------------------------------

      SUBROUTINE coefun(x,p,q,w,coeffn,params,iparam,icofun)
C     .. Scalar Arguments ..
      DOUBLE PRECISION p,q,w,x
      INTEGER icofun,iparam
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ex,fac
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..

C Params is an unused array introduced to enable the code to be
C easily adapted to handle problems with parameter-dependent
C coefficients.

      IF (icofun.EQ.1) THEN
          fac = 1.D0 - abs(x)
          ex = x/fac
          fac = fac**2
          CALL coeffn(ex,p,q,w)
          p = p*fac
          q = q/fac
          w = w/fac
          RETURN

      ELSE

          CALL coeffn(x,p,q,w)
          RETURN

      END IF

      END SUBROUTINE coefun
C ---------------------------------------------------------------------
C --------------------- Source from NEWCOMP ---------------------------
C ---------------------------------------------------------------------

C Routine for single eigenvalue case:
      SUBROUTINE sl04f(xval,yvals,pdyval,ival,elam,eigfns,wk,iwk,nval,n,
     &                 ifail)
C Brief Specification:
C This routine takes in appropriate information about the normalised
C eigenfunction from the routine SL03F and returns the value of the
C normalised eigenfunction at the points xval(i),i=1,nval.
C VARIABLES:
C      XVAL: REAL ARRAY of DIMENSION (0:IVAL). A set of points supplied
C            by the user. The eigenfunction will be evaluated at the
C            points XVAL(i) for 0 <= i <= NVAL. (n.b. IVAL is the
C            DIMENSION of XVAL as declared in the calling (sub)program,
C            NVAL the number of these points at which the eigenfunction
C            is to be evaluated). The points in XVAL MUST BE ARRANGED
C            in ASCENDING ORDER,
C                XVAL(0) <= XVAL(1) <= XVAL(2) <= ... <= XVAL(NVALS).
C            otherwise SL04F will rearrange them in this order.
C     YVALS: REAL ARRAY of DIMENSION (0:IVAL). On exit, YVALS(i) is, for each
C            0 <= i <= NVALS,  the value of the approximate eigenfunction at
C            the point XVAL(i).
C    PDYVAL: REAL ARRAY of DIMENSION (0:IVAL). On exit, PDYVAL(i) is, for each
C            0 <= i <= NVALS, the value of the approximation to p(x)y'(x) at
C            x = XVAL(i). Here p(x) is the coefficient from the defining
C            equation of the Sturm-Liouville problem,
C
C                     -(p(x)y')' + q(x)y = Lambda.w(x)y
C      ELAM: REAL ARRAY of DIMENSION (1:2). Supplied by the preceding call
C            to the routine SL02F. Unchanged on exit.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4). Supplied by the preceding
C            call to the routine SL02F. The preceding call to SL02F need
C            not occur in the same subprogram, but if it does not then the
C            array wk must be declared to have the same dimensions in the
C            subprogram from which SL02F is called as in the subprogrm from
C            which SL04F is called.
C            Unchanged on exit.
C    EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3). Supplied by the preceding
C            call to the routine SL03F.
C            Unchanged on exit.
C       IWK: The first dimension of wk as declared in the calling
C            subprogram. Unchanged on exit.
C     NVALS: One less than the number of points at which the eigenfunction
C            is required. Unchanged on exit.
C         N: The number of mesh intervals used by SL02F in computing the
C            eigenvalue approximation. Supplied by the preceding call to
C            SL02F. Unchanged on exit.
C END of brief specification.
C
C Additional Information (for me):
C                eigfns(i,1) = rlog(i) for i=0,1
C                eigfns(i,2) = theta(i) for i=0,1
C                eigfns(i,3) = scale(i) for i=0,1.
C     wk: *real* ARRAY of DIMENSION (0:iwk,1:4). The entries in wk are as
C            follows:
C                wk(i,1) = xmesh(i) for i=0,n
C                wk(i,2) = pp(i)    for i=1,n
C                wk(i,3) = qp(i)    for i=1,n
C                wk(i,4) = wp(i)    for i=1,n.
C
C     .. Scalar Arguments ..
      INTEGER ifail,ival,iwk,n,nval
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3),elam(1:2),pdyval(0:ival),
     &                 wk(0:iwk,1:4),xval(0:ival),yvals(0:ival)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam0
      INTEGER ic
C     ..
C     .. External Subroutines ..
cc      EXTERNAL efn1
C     ..
C Retrieve matchpoint index from sneaky storage:
      ic = wk(0,2)
C Un-do the Richardson extrapolation:
      elam0 = elam(1) + elam(2)
C Compute the eigenfunction and store in arrays:
      CALL efn1(xval,yvals,pdyval,elam0,eigfns(0,1),eigfns(0,2),
     &          eigfns(0,3),wk(0,1),wk(1,2),wk(1,3),wk(1,4),nval,n,ic)
C IFAIL is included just for compatibility.
      ifail = 0

      END SUBROUTINE sl04f
C ----------------------------------------------------------------------------

      SUBROUTINE sl03f(elam,k,eigfns,wk,iwk,n,setup,ifail)
C This routine computes information required by the routine SL04F in order
C to evaluate the normalised eigenfunction.
C
C Brief Specification:
C VARIABLES
C     ELAM: REAL ARRAY of DIMENSION (1:2). Supplied by the preceding call to
C           SL02F. Unchanged on exit.
C        K: The index of the eigenvalue and eigenfunction in question. As
C           supplied to SL02F for the computation of ELAM. Unchanged on exit.
C   EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3) This array must be supplied
C           by the user and will be used to store the information about the
C           normalised eigenfunction required by SL03F.
C       WK: REAL ARRAY of DIMENSION (0:IWK,1:4). Supplied by the preceding
C            call to the routine SL02F. The preceding call to SL02F need
C            not occur in the same subprogram, but if it does not then the
C            array WK must be declared to have the same dimensions in the
C            subprogram from which SL02F is called as in the subprogrm from
C         N: INTEGER. The number of mesh intervals used by the routine SL02F
C            in the computation of the eigenvalue approximation. Supplied by
C            the preceding call to SL02F and unchanged on exit.
C     SETUP: User-supplied SUBROUTINE. The same as supplied to SL02F at the
C            preceding call to SL02F.
C   IFAIL: INTEGER. This parameter should be zero on exit. If it is not then
C          it will have the value 3 indicating an error in the routine SETUP.
C END of brief specification

C     .. Parameters ..
      DOUBLE PRECISION smatch
      PARAMETER (smatch=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER ifail,iwk,k,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3),elam(1:2),wk(0:iwk,1:4)
C     ..
C     .. Local Scalars ..
      INTEGER ic,iflag,nrec
      CHARACTER srname*6
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. External Functions
C     .. External Subroutines ..
cc      EXTERNAL norfn1
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      srname = ' SL03F'
      nrec = 1
      iflag = 0
      IF (abs(ifail).GT.1) THEN
          WRITE (rec,FMT=9000)
          ifail = -1
          iflag = 1

      ELSE

C Retrive matchpoint index from sneaky storage:
          ic = wk(0,2)
C This routine just acts as an interface for NORFUN
          CALL norfn1(eigfns(0,1),eigfns(0,2),eigfns(0,3),smatch,
     &                elam(1),elam(2),k,n,ic,wk(1,2),wk(1,3),wk(1,4),
     &                wk(0,1),setup,iflag)
          IF (iflag.NE.0) THEN
              WRITE (rec,FMT=9010)

          ELSE

C IF iflag is nonzero iflag will have the value 3
C     Eigenfunction has been normalised as required
              ifail = 0
              RETURN

          END IF

      END IF

      ifail = p01abf(ifail,iflag,srname,nrec,rec)
9000  FORMAT (' ** Parameter error: ifail out of range')
9010  FORMAT (' ** Invalid boundary conditions')
C$st$ Unreachable comments ...

      END SUBROUTINE sl03f

C------------------------------------------------------------------------------
C Routine for many eigenvalue case:

      SUBROUTINE sl04fm(xval,yvals,pdyval,ival,elam,k,kvals,m,mvals,
     &                  eigfns,wk,iwk,nval,n,ifail)
C Brief Specification:
C This routine takes in appropriate information about the normalised
C eigenfunctions from the routine NRMANY and returns the values of the
C normalised eigenfunctions at the points xval(i),i=1,nval.
C VARIABLES:
C      XVAL: REAL ARRAY of DIMENSION (0:IVAL). A set of points supplied
C            by the user. The eigenfunctions with indices KVALS(i),
C            i = 1,..,MVALS, will be evaluated at the
C            points XVAL(i) for 0 <= i <= NVAL. (n.b. IVAL is the
C            DIMENSION of XVAL as declared in the calling (sub)program,
C            NVAL the number of these points at which the eigenfunction
C            is to be evaluated). The points in XVAL MUST BE ARRANGED
C            in ASCENDING ORDER,
C                XVAL(0) <= XVAL(1) <= XVAL(2) <= ... <= XVAL(NVALS).
C            Unchanged on exit.
C         K: INTEGER array of DIMENSION at least (1:M). This array must
C            contain the same values as the array by the same name supplied to
C            the routine NRMANY to generate the eigenfunction information
C            array EIGFNS. Unchanged on exit.
C     MVALS: INTEGER. There are M eigenfunctions stored in the array EIGFNS.
C            MVALS specifies how many of these are to be evaluated at each point
C            in the array XVALS. MVALS <= M. Unchanged on exit.
C     KVALS: INTEGER array of DIMENSION at least (1:MVALS). On entry, KVALS
C            must be set by the user to specify the indices KVALS(i), i =
C            1,..,MVALS, of the eigenfunctions to be evaluated. Unchanged
C            on exit.
C     YVALS: REAL ARRAY of DIMENSION (0:IVAL,1:p) where p>= MVALS.
C            On exit, YVALS(i,j) is, for each 0 <= i <= NVALS, 1<=j<=MVALS,
C            the value of the approximate KVALS(j)th eigenfunction at the
C            point XVAL(i).
C    PDYVAL: REAL ARRAY of DIMENSION (0:IVAL,1:p) where p>=MVALS
C            On exit, PDYVAL(i,j) is, for each 0 <= i <= NVALS, 1<=j<=MVALS
C            the value of the approximation to p(x)y_{KVALS(j)}'(x) at
C            x = XVAL(i). Here p(x) is the coefficient from the defining
C            equation of the Sturm-Liouville problem,
C
C                     -(p(x)y')' + q(x)y = Lambda.w(x)y
C            and y_{KVALS(j)} is the KVALS(j)th eigenfunction of this problem.
C      ELAM: REAL ARRAY of DIMENSION (1:2,1:p) where p >= m.
C            Supplied by the preceding call to the routine SL02FM. Unchanged
C            on exit.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4)
C            Supplied by the preceding call to the routine SL02FM. The
C            preceding call to SL02FM need not occur in the same subprogram,
C            but if it does not then the array wk must be declared to have the
C            same dimensions in the subprogram from which SL02FM is called as
C            in the subprogrm from which SL04F is called.
C            Unchanged on exit.
C    EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3,1:p) where p>=m.
C            Supplied by the preceding call to the routine NRMANY. Note that
C            IWK is the dimension of this routine as declared in the calling
C            (sub)program, and must be the same as first dimension of wk as
C            declared in the calling subprogram. The preceding call to NRMANY
C            need not occur in the same subprogram, but if it does not then the
C            array EIGFNS must be declared to have the same dimensions in the
C            subprogram from which NRMANY is called as in the subprogrm from
C            which SL04FM is called.
C            Unchanged on exit.
C       IWK: INTEGER. The first dimension of wk  as declared in the
C            calling subprogram. Unchanged on exit.
C     NVALS: INTEGER. The number of points in XVAL at which the eigenfunction
C            approximation is to be computed. Unchanged on exit.
C         N: INTEGER. The number of mesh intervals used by SL02FM in computing
C            the eigenvalue approximation. Supplied by the preceding call to
C            SL02FM. Unchanged on exit.
C         M: INTEGER. The number of eigenvalues computed at the preceding call
C            to SL02FM. As supplied to SL02FM by the user. Unchanged on exit.
C END of brief spec.
C
C Additional Information (for me):
C                eigfns(i,1,k) = rlog_{k}(i) for i=0,1,k=1,m,
C                eigfns(i,2,k) = theta_{k}(i) for i=0,1,k=1,m,
C                eigfns(i,3,k) = scale_{k}(i) for i=0,1,k=1,m.
C     wk: *real* ARRAY of DIMENSION (0:n,1:4). The entries in wk are as
C            follows:
C                wk(i,1) = xmesh(i) for i=0,n
C                wk(i,2) = pp(i)    for i=1,n
C                wk(i,3) = qp(i)    for i=1,n
C                wk(i,4) = wp(i)    for i=1,n.
C
C     .. Scalar Arguments ..
      INTEGER ifail,ival,iwk,m,mvals,n,nval
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3,1:m),elam(1:2,1:m),
     &                 pdyval(0:ival,1:m),wk(0:iwk,1:4),xval(0:ival),
     &                 yvals(0:ival,1:m)
      INTEGER k(1:m),kvals(1:mvals)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam0
      INTEGER i,ic,iflag,j,jc,nrec
      CHARACTER srname*6
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL efn1
C     ..
C Retrieve matchpoint index from sneaky storage:
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      ic = wk(0,2)
      iflag = 0
      nrec = 2
      IF (mvals.GT.m .OR. mvals.LE.0 .OR. abs(ifail).GT.1) THEN
          iflag = 1
          IF (abs(ifail).GT.1) ifail = -1
          WRITE (rec,FMT=9000)

      ELSE

          DO 30 i = 1,mvals
              DO 10 j = 1,m
                  IF (kvals(i).EQ.k(j)) GO TO 20
10            CONTINUE
              GO TO 40

20            jc = j
C Undo each Richarson extrapolation before computing the eigenfunctions:
              elam0 = elam(1,jc) + elam(2,jc)
              CALL efn1(xval,yvals(0,i),pdyval(0,i),elam0,
     &                  eigfns(0,1,jc),eigfns(0,2,jc),eigfns(0,3,jc),
     &                  wk(0,1),wk(1,2),wk(1,3),wk(1,4),nval,n,ic)
30        CONTINUE

          ifail = 0
          RETURN

40        CONTINUE

          iflag = 2
          WRITE (rec,FMT=9010)
      END IF

      ifail = p01abf(ifail,iflag,srname,nrec,rec)
9000  FORMAT (' ** Parameter error: MVALS or IFAIL out of range')
9010  FORMAT (' ** Requested eigenfunctions not all available')

      END SUBROUTINE sl04fm

C ------------------------------------------------------------------------

      SUBROUTINE sl03fm(elam,k,m,eigfns,wk,iwk,n,setup,ifail)
C This routine computes information required by the routine SL04FM in order
C to evaluate the normalised eigenfunction.
C
C Brief Specification:
C VARIABLES
C     ELAM: REAL ARRAY of DIMENSION (1:2,1:p) where p >= m.  Supplied by the
C           preceding call to SL02FM. Unchanged on exit.
C        K: INTEGER ARRAY of DIMENSION at least (1:M). The array of indices
C           of the eigenvalues and eigenfunctions in question. As supplied to
C           SL02FM for the computation of ELAM. Unchanged on exit.
C   EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3,1:p) where p >= m. This array
C           must be supplied by the user and will be used to store the
C           information about the normalised eigenfunctions required by SL04FM.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4)
C            Supplied by the preceding call to the routine SL02FM. The
C            preceding call to SL02FM need not occur in the same subprogram,
C            but if it does not then the array wk must be declared to have the
C            same dimensions in the subprogram from which SL02FM is called as
C            in the subprogrm from which SL04F is called.
C            Unchanged on exit.
C       IWK: INTEGER. The first dimension of WK as declared in the calling
C            (sub)program. Unchanged on exit.
C         N: INTEGER. The number of mesh intervals used by the routine SL02FM
C            in the computation of the eigenvalue approximation. Supplied by
C            the preceding call to SL02FM and unchanged on exit.
C         M: INTEGER. The number of eigenvalues computed at the preceding call
C            to SL02FM. As supplied to SL02FM by the user. Unchanged on exit.
C     SETUP: User-supplied SUBROUTINE. The same as supplied to SL02FM at the
C            preceding call to SL02FM.
C   IFAIL: INTEGER. This parameter should be zero on exit. If it is not then
C          it will have the value 3 indicating an error in the routine SETUP.
C END of brief specification
C
C Additional information (for me):
C This routine interfaces to NORFN1 and allows the normalisation of many
C eigenfunctions at once.
C     .. Parameters ..
      DOUBLE PRECISION smatch
      PARAMETER (smatch=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER ifail,iwk,m,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3,1:m),elam(1:2,1:m),wk(0:iwk,1:4)
      INTEGER k(1:m)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. Local Scalars ..
      INTEGER i,ic,iflag,nrec
      CHARACTER srname*6
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL norfn1
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      iflag = 0
      srname = 'SL03FM'
      nrec = 2
      IF (abs(ifail).GT.1) THEN
          WRITE (rec,FMT=9000)
          ifail = -1
          iflag = 1

      ELSE
C     Retrieve matchpoint index from sneaky storage:
          ic = wk(0,2)
C This routine is just an interface for NORFUN
          DO 10 i = 1,m
              CALL norfn1(eigfns(0,1,i),eigfns(0,2,i),eigfns(0,3,i),
     &                    smatch,elam(1,i),elam(2,i),k(i),n,ic,wk(1,2),
     &                    wk(1,3),wk(1,4),wk(0,1),setup,iflag)

              IF (iflag.NE.0) GO TO 20
10        CONTINUE
C If ifail is nonzero ifail will have the value 3.
C     Eigenfunctions have been normalised as required
          ifail = 0
          RETURN

20        WRITE (rec,FMT=9010)
          iflag = 2
      END IF

      ifail = p01abf(ifail,iflag,srname,nrec,rec)
9000  FORMAT (' ** Invalid input value of IFAIL')
9010  FORMAT (' ** Invalid boundary conditions')
C$st$ Unreachable comments ...

      END SUBROUTINE sl03fm

C------------------------------------------------------------------------------

      SUBROUTINE efun(xval,yvals,pdyval,elam,rlog,theta,scale,xmesh,pp,
     &                qp,wp,nval,n,ic)
C     .. Parameters ..
      DOUBLE PRECISION zero,one
      PARAMETER (zero=0.0D0,one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam
      INTEGER ic,n,nval
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pdyval(0:nval),pp(1:n),qp(1:n),rlog(0:n),
     &                 scale(0:n),theta(0:n),wp(1:n),xmesh(0:n),
     &                 xval(0:nval),yvals(0:nval)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION di,dummy,erl,rlgtmp,sqrts,stemp,thtmp,xend,
     &                 xo
      INTEGER i,j,jtemp,kdummy
C     ..
C     .. External Subroutines ..
cc      EXTERNAL fulstp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos,sin,sqrt
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION spexp
cc      EXTERNAL spexp
C     ..
C      rnorm = zero
      i = 0
      j = 0
      di = one
10    CONTINUE
      IF (xval(i).LT.xmesh(0)) di = -one
      jtemp = j
      DO 20 kdummy = j,n - 1
          IF (xval(i).GT.xmesh(kdummy)) THEN
              jtemp = kdummy

          ELSE
              GO TO 30

          END IF

20    CONTINUE
      GO TO 40

30    CONTINUE

40    j = jtemp
      xend = xval(i)
      IF (xend.GT.xmesh(ic)) THEN
          xo = xmesh(j+1)
          rlgtmp = rlog(j+1)
          thtmp = theta(j+1)
          stemp = scale(j+1)
          di = -one

      ELSE
          rlgtmp = rlog(j)
          thtmp = theta(j)
          stemp = scale(j)
          xo = xmesh(j)
      END IF

      CALL fulstp(xo,xend,di,stemp,thtmp,rlgtmp,pp(j+1),
     &            elam*wp(j+1)-qp(j+1),wp(j+1),dummy,1)
      erl = spexp(rlgtmp)
      sqrts = sqrt(stemp)
      yvals(i) = sin(thtmp)*erl/sqrts
      pdyval(i) = cos(thtmp)*erl*sqrts
      i = i + 1
      IF (i.LE.nval) GO TO 10

      END SUBROUTINE efun
C------------------------------------------------------------------------------

      SUBROUTINE efn1(xval,yvals,pdyval,elam,rlog,theta,scale,xmesh,pp,
     &                qp,wp,nval,n,ic)
C     .. Parameters ..
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam
      INTEGER ic,n,nval
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pdyval(0:nval),pp(1:n),qp(1:n),rlog(0:1),
     &                 scale(0:1),theta(0:1),wp(1:n),xmesh(0:n),
     &                 xval(0:nval),yvals(0:nval)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION di,dummy,p,q,rlgloc,scloc,sgnh,thtloc,w,xend,xo
      INTEGER i,idi,ifail,ii,j,jsplit
C     ..
C     .. External Subroutines ..
cc      EXTERNAL fulstp,m01caf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos,dble,max,min,sign,sin,sqrt
C     ..

C First order points in ascending order
C     .. External Functions ..
cc      DOUBLE PRECISION spexp
cc      EXTERNAL spexp
C     ..
      CALL m01caf(xval(0),1,nval+1,'A',ifail)

C Next find which points lie to the right of the match-point and which to
C the left.
      jsplit = -1
      DO 10 i = 0,nval
C Modify a point if it lies outwith the permitted range:
          xval(i) = min(xmesh(n),max(xmesh(0),xval(i)))
          IF (xval(i).LE.xmesh(ic)) jsplit = i
10    CONTINUE
      i = 1
      j = 0

      IF (jsplit.EQ.-1) THEN
          jsplit = 0
          GO TO 30

      END IF

      idi = 1
      ii = i
C Initial conditions
      rlgloc = rlog(0)
      thtloc = theta(0)
      scloc = scale(0)
      xo = xmesh(0)
      di = dble(idi)
      sgnh = sign(one,di)

20    p = pp(ii)
      q = qp(ii)
      w = wp(ii)
      IF (xmesh(i)*sgnh.LT.xval(j)*sgnh) THEN
          xend = xmesh(i)
          CALL fulstp(xo,xend,sign(one,xend-xo),scloc,thtloc,rlgloc,p,
     &                elam*w-q,w,dummy,1)
          i = i + idi
          ii = ii + idi

      ELSE
          xend = xval(j)
          CALL fulstp(xo,xend,sign(one,xend-xo),scloc,thtloc,rlgloc,p,
     &                elam*w-q,w,dummy,1)
          dummy = spexp(rlgloc)
          yvals(j) = sin(thtloc)*dummy/sqrt(scloc)
          pdyval(j) = cos(thtloc)*dummy*sqrt(scloc)
          IF (xmesh(i).EQ.xval(j)) THEN
              i = i + idi
              ii = ii + idi
              i = min(i,n)
              i = max(i,0)
              ii = min(ii,n)
              ii = max(ii,0)
          END IF

          j = j + idi
      END IF

C If we have evaluated at all required points then bail out:
      IF (idi.GT.0 .AND. j.GT.jsplit) GO TO 30
      IF (idi.LT.0 .AND. j.LT.jsplit) GO TO 40
      xo = xend
      GO TO 20

30    IF (j.LE.nval) THEN
C We still have more points to do.
          rlgloc = rlog(1)
          thtloc = theta(1)
          scloc = scale(1)
          xo = xmesh(n)
          i = n - 1
          idi = -1
          ii = i + 1
          j = nval
          di = dble(idi)
          sgnh = sign(one,di)
          GO TO 20

      END IF

40    RETURN

C If we are here then all the required values have been computed.

      END SUBROUTINE efn1

C------------------------------------------------------------------------------
      SUBROUTINE norfun(rlog,theta,scale,smatch,elam0,ep,k,n,ic,pp,qp,
     &                  wp,xmesh,setup,trans,ifail)
C
C This routine computes arrays containing nodal values of the scale-factor,
C Prufer radius and Prufer angle, all normalised to agree with each other and
C yield an eigenfunction with L2(w)-norm equal to unity. The mesh may be in
C the user's coordinates (trans = .false.) or transformed coordinates
C (trans = .true.)
C VARIABLES:
C          rlog:  *real* ARRAY of DIMENSION at least (0:n) On exit, this contain
C                the required nodal values of log(r).
C      theta: *real* ARRAY of DIMENSION at least (0:n). On exit, this contains
C                the required nodal values of theta.
C      scale: *real* ARRAY of DIMENSION at least (0:n). On exit, this contains
C                the required nodal values of the scale-factor.
C     ELAM0: *real*. The eigenvalue corresponding to the required eigenfunction.
C             Unchanged on exit.
C        N,IC: INTEGERs. N is the number of mesh-intervals used by  EIGEN in its
C             computation of the eigenfunction. IC is the integer such that the
C             point XMESH(IC) was used as the matching point for the computation
C             of the Prufer miss-distance function.
C   PP,QP,WP: *real* ARRAYs of DIMENSIONs at least N. These store the values
C              of the coefficient functions at the midpoints XMESH(i+1/2).
C      XMESH: *real* ARRAY of DIMENSION at least N. Stores the mesh-points used
C             by EIGEN in the computation of the eigenfunction.
C     smatch: the value of the scale-factor at the matching point xmesh(ic).
C   IFAIL: INTEGER. This parameter should be zero on exit. If it is not then
C          it will have the value 3 indicating an error in the routine SETUP.
C     .. Parameters ..
C Set up left boundary conditions:
      DOUBLE PRECISION zero,one,two,half
      PARAMETER (zero=0.0D0,one=1.0D0,two=2.0D0,half=one/two)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam0,ep,smatch
      INTEGER ic,ifail,k,n
      LOGICAL trans
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(n),qp(n),rlog(0:n),scale(0:n),theta(0:n),
     &                 wp(n),xmesh(0:n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam,eps,facl,facr,pdy,pi,rgtlog,rlftlg,rnorml,
     &                 rnormr,y
      INTEGER i
C     ..
C     .. External Subroutines ..
cc      EXTERNAL bcons,fulprf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC atan2,log
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION spexp,x01aaf
cc      EXTERNAL spexp,x01aaf
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
      ifail = 0
      pi = x01aaf(0.0D0)
      elam = elam0 + ep
      CALL bcons(setup,y,pdy,elam,xmesh(0),pp(1),qp(1),wp(1),0,trans,
     &           ifail)
      IF (ifail.EQ.0) THEN
C If ifail is nonzero ifail = 3
          theta(0) = atan2(y,pdy)
          rlog(0) = zero
C Shoot from left-hand end to matching point:
          rnorml = zero
          scale(0) = one
          i = 0
          CALL fulprf(i,theta,rlog,scale,rnorml,ic,smatch,elam,pp,qp,wp,
     &                xmesh,n,n,-1)
          rlftlg = rlog(ic)
C Set up right boundary conditions:
          CALL bcons(setup,y,pdy,elam,xmesh(n),pp(n),qp(n),wp(n),1,
     &               trans,ifail)
          IF (ifail.EQ.0) THEN
C If ifail is nonzero ifail = 3
              theta(n) = (k+1)*pi - atan2(y,pdy)
              rlog(n) = zero
C Shoot from right-hand end to matching point:
              rnormr = zero
              scale(n) = one
              i = n
              CALL fulprf(i,theta,rlog,scale,rnormr,ic,smatch,elam,pp,
     &                    qp,wp,xmesh,n,n,-1)
              rgtlog = rlog(ic)
C Choose the normalising factors facl and facr to match up the two halves
C of the solution; actually facl and facr here denote the logarithms of the
C normalising factors.
C To reduce the probability of overflow we consider two cases.
C
              IF (rgtlog.GT.rlftlg) THEN
                  eps = spexp(rlftlg-rgtlog)
                  facl = -log((rnormr*eps)**2+ (rnorml**2))*half
                  facr = facl + rlftlg - rgtlog

              ELSE
                  eps = spexp(rgtlog-rlftlg)
                  facr = -log((rnorml*eps)**2+ (rnormr**2))*half
                  facl = facr + rgtlog - rlftlg
              END IF
C
C The normalising factors have now been computed.
C
              DO 10 i = 0,ic - 1
                  rlog(i) = rlog(i) + facl
10            CONTINUE
              DO 20 i = ic,n
                  rlog(i) = rlog(i) + facr
20            CONTINUE
C
C We now have the following information:
C
C The scale factors scale(i), the values of rlog(i), and the values of
C theta(i), for i=1,..,n, with the theta values scaled in a manner which is
C appropriate for the scale factors and the rlog values scaled to agree with
C scale-factors and yield an eigenfunction witth L2(w)-norm equal to unity.
C
          END IF

      END IF

      END SUBROUTINE norfun
C-------------------------------------------------------------------------------
      SUBROUTINE norfn1(rlog,theta,scale,smatch,elam0,ep,k,n,ic,pp,qp,
     &                  wp,xmesh,setup,ifail)
C
C This routine computes arrays containing end-point values of the scale-factor,
C Prufer radius and Prufer angle, all normalised to agree with each other and
C yield an eigenfunction with L2(w)-norm equal to unity.
C The mesh is assumed to be supplied in the user's original coordinates
C and not in transformed coordinates.
C VARIABLES:
C    rlog:  *real* ARRAY of DIMENSION at least (0:1) On exit, this contains
C             the required end-point values of log(r).
C   theta: *real* ARRAY of DIMENSION at least (0:1). On exit, this contains
C             the required end-point values of theta.
C   scale: *real* ARRAY of DIMENSION at least (0:1). On exit, this contains
C             the required end-point values of the scale-factor.
C   ELAM0: *real*. The eigenvalue corresponding to the required eigenfunction.
C          Unchanged on exit.
C    N,IC: INTEGERs. N is the number of mesh-intervals used by  EIGEN in its
C          computation of the eigenfunction. IC is the integer such that the
C          point XMESH(IC) was used as the matching point for the computation
C          of the Prufer miss-distance function.
CPP,QP,WP: *real* ARRAYs of DIMENSIONs at least N. These store the values
C          of the coefficient functions at the midpoints XMESH(i+1/2).
C   XMESH: *real* ARRAY of DIMENSION at least N. Stores the mesh-points used
C          by EIGEN in the computation of the eigenfunction.
C  smatch: the value of the scale-factor at the matching point xmesh(ic).
C   IFAIL: INTEGER. This parameter should be zero on exit. If it is not then
C          it will have the value 3 indicating an error in the routine SETUP.
C     .. Parameters ..
C Set up left boundary conditions:
      DOUBLE PRECISION zero,one,two,half
      PARAMETER (zero=0.0D0,one=1.0D0,two=2.0D0,half=one/two)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam0,ep,smatch
      INTEGER ic,ifail,k,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(n),qp(n),rlog(0:1),scale(0:1),theta(0:1),
     &                 wp(n),xmesh(0:n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam,eps,facl,facr,pdy,pi,rgtlog,rlftlg,rlo,
     &                 rnorml,rnormr,rro,sclo,scro,thlo,thro,y
      INTEGER i
      LOGICAL trans
C     ..
C     .. External Subroutines ..
cc      EXTERNAL bcons,fulprf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC atan2,log
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION spexp,x01aaf
cc      EXTERNAL spexp,x01aaf
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
      ifail = 0
      pi = x01aaf(0.0D0)
      elam = elam0 + ep
C This routine never works with transformed coordinates.
      trans = .false.
      CALL bcons(setup,y,pdy,elam,xmesh(0),pp(1),qp(1),wp(1),0,trans,
     &           ifail)
      IF (ifail.EQ.0) THEN
C If ifail is nonzero ifail = 3
          theta(0) = atan2(y,pdy)
          rlog(0) = zero
C Shoot from left-hand end to matching point:
          rnorml = zero
          scale(0) = one
          i = 0
          CALL fulprf(i,theta,rlog,scale,rnorml,ic,smatch,elam,pp,qp,wp,
     &                xmesh,1,n,-2)
          rlftlg = rlog(1)
          rlo = rlog(0)
          thlo = theta(0)
          sclo = scale(0)
C Set up right boundary conditions:
          CALL bcons(setup,y,pdy,elam,xmesh(n),pp(n),qp(n),wp(n),1,
     &               trans,ifail)
          IF (ifail.EQ.0) THEN
C If ifail is nonzero ifail = 3
              theta(0) = (k+1)*pi - atan2(y,pdy)
              rlog(0) = zero
C Shoot from right-hand end to matching point:
              rnormr = zero
              scale(0) = one
              i = n
              CALL fulprf(i,theta,rlog,scale,rnormr,ic,smatch,elam,pp,
     &                    qp,wp,xmesh,1,n,-2)
              rgtlog = rlog(1)
              rro = rlog(0)
              thro = theta(0)
              scro = scale(0)
C Choose the normalising factors facl and facr to match up the two halves
C of the solution; actually facl and facr here denote the logarithms of the
C normalising factors.
C To reduce the probability of overflow we consider two cases.
C
              IF (rgtlog.GT.rlftlg) THEN
                  eps = spexp(rlftlg-rgtlog)
                  facl = -log((rnormr*eps)**2+ (rnorml**2))*half
                  facr = facl + rlftlg - rgtlog

              ELSE
                  eps = spexp(rgtlog-rlftlg)
                  facr = -log((rnorml*eps)**2+ (rnormr**2))*half
                  facl = facr + rgtlog - rlftlg
              END IF
C
C The normalising factors have now been computed.
C

              rlog(0) = rlo + facl
              rlog(1) = rro + facr
              theta(0) = thlo
              theta(1) = thro
              scale(0) = sclo
              scale(1) = scro
C
          END IF

      END IF

      END SUBROUTINE norfn1
C-------------------------------------------------------------------------------
      SUBROUTINE fulprf(i,thetrr,rlogrr,scale,rnorm,iend,send,elam,pp,
     &                  qp,wp,xmesh,idim,ndim,iswtch)
C
C Integrates the Prufer phase and radius (forward or backward depending on
C whether iend islarger or smaller than i) advancing the values of i,theta and s
C and accumulating an approximation to the L2(w) integral of the eigenfunction
C at each step.
C iswtch is a switch parameter. It works as follows:
C
C iswtch = 1 or 2: both the Prufer radius and the Prufer angle are advanced from
C               i to iend.
C iswtch= -1 or -2: all the information returned by iswtch = 1 is returned,
C and in addition the norm of the corresponding eigenfunction is returned.
C The arrays thetrr and rlogrr serve as follows:
C IF iswtch=(+/-)1, then on entry rlogrr(i), thetrr(i) and scale(i) must
C                   contain the initial values of the log of Prufer radius,
C                   the Prufer angle and scale-factor. On exit, rlogrr(j),
C                   thetrr(j) and scale(j) will contain these Prufer
C                   quantities at each of the points xmesh(j) for j=i,iend.
C                   (n.b.we refer here to the input i; i is changed on exit).
C                   If a rescaling has been necessary, then the initial
C                   log Prufer radius rlog(i) will have been changed to
C                   correspond to the new scaling.
C IF iswtch=(+/-)2, then on entry rlogrr(0), thetrr(0) and scale(0) must
C                   contain the initial values of the log of the Prufer radius,
C                   the Prufer angle and the scale-factor. On exit, they
C                   will be unchanged unless rlog(0) has been changed to
C                   correspond to a new scaling of the eigenfunction. The
C                   elements rlogrr(1), thetrr(1) and scale(1) contain the
C                   appropriate Prufer quantities at xmesh(iend).
C
C In both cases, the scale-factor s at the final point will be send, specified
C by the user.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,rnorm,send
      INTEGER i,idim,iend,iswtch,ndim
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(ndim),qp(ndim),rlogrr(0:idim),scale(0:idim),
     &                 thetrr(0:idim),wp(ndim),xmesh(0:ndim)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION delta,di,hlfsub,p,q,rlog,s,safe,sbtrct,theta,
     &                 twosaf,w
      INTEGER idi,idummy,ii,iold,ioldo
      LOGICAL ltf
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION rescal,scl,spexp,x02amf
cc      EXTERNAL rescal,scl,spexp,x02amf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL fulstp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,log,max,sign,sqrt
C     ..
C     .. External Functions ..
C     ..
C Set a parameter which can be safely exponentiated:
C     .. Parameters ..
      DOUBLE PRECISION one,zero,half
      PARAMETER (one=1.0D0,zero=0.0D0,half=5.0D-1)
C     ..
      twosaf = -log(x02amf(one))*half
      safe = half*twosaf
      ltf = abs(iswtch) .EQ. 1
C Store input values of certain parameters:
      ioldo = i
C Set the L2(w) norm to zero before starting the integration
      rnorm = zero
C Start stepping
      idi = sign(1,iend-i)
      di = dble(idi)
      ii = i
      IF (idi.LT.0) ii = i + 1
      IF (ltf) THEN
          s = scale(i)
          theta = thetrr(i)
          rlog = rlogrr(i)

      ELSE
          s = scale(0)
          theta = thetrr(0)
          rlog = rlogrr(0)
C Set scale(1), thetrr(1) and rlogrr(1) in case the range of integration
C is null:
          scale(1) = scale(0)
          thetrr(1) = thetrr(0)
          rlogrr(1) = rlogrr(0)

      END IF

10    CONTINUE

C
      IF (i.NE.iend) THEN
          iold = i
          i = i + idi
          ii = ii + idi
          p = pp(ii)
          w = wp(ii)
C  Q here has same sign as in writeup of 12/1/89
          q = wp(ii)*elam - qp(ii)
C  Integrate over the interval [xmesh(iold),xmesh(i)]
C          h = xmesh(i) - xmesh(iold)
          CALL fulstp(xmesh(iold),xmesh(i),di,s,theta,rlog,p,q,w,delta,
     &                sign(1,iswtch))
C Update the logr and theta arrays, and the scale-factor array:
          IF (ltf) THEN
              rlogrr(i) = rlog
              scale(i) = s
              thetrr(i) = theta

          ELSE
              rlogrr(1) = rlog
              scale(1) = s
              thetrr(1) = theta
          END IF
C Do some scaling to avoid overflow problems:
          IF (rlog.GT.safe) THEN
C Subtract a suitable amount from every value of logr, and divide the
C norm accordingly:
              IF (sign(1,iswtch).EQ.-1) THEN
C The contribution to the norm from a single step cannot exceed 1 for the
C normalised eigenfunction:
                  sbtrct = max(delta,twosaf)
                  hlfsub = half*sbtrct
                  rnorm = spexp(-sbtrct)*rnorm
                  delta = delta - sbtrct
                  rlog = rlog - hlfsub
                  IF (ltf) THEN
                      DO 20 idummy = ioldo,i,idi
                          rlogrr(idummy) = rlogrr(idummy) - hlfsub
20                    CONTINUE

                  ELSE
                      rlogrr(0) = rlogrr(0) - hlfsub
                      rlogrr(1) = rlogrr(1) - hlfsub
                  END IF

              ELSE
                  rlog = rlog - safe
                  IF (ltf) THEN
                      DO 30 idummy = ioldo,i,idi
                          rlogrr(idummy) = rlogrr(idummy) - safe
30                    CONTINUE

                  ELSE
                      rlogrr(0) = rlogrr(0) - safe
                      rlogrr(1) = rlogrr(1) - safe
                  END IF

              END IF

          END IF
          IF (sign(1,iswtch).EQ.-1) rnorm = rnorm + spexp(delta)

          GO TO 10

      END IF
C Integration of the Prufer equations is now complete
      rnorm = sqrt(rnorm)

      IF (ltf) THEN
          thetrr(i) = scl(thetrr(i),send/s)
          scale(i) = send
          rlog = rlogrr(i)
          rlogrr(i) = rescal(s,send,thetrr(i),rlog)

      ELSE
          thetrr(1) = scl(thetrr(1),send/s)
          scale(1) = send
          rlog = rlogrr(1)
          rlogrr(1) = rescal(s,send,thetrr(1),rlog)
      END IF
      END SUBROUTINE fulprf

C-------------------------------------------------------------------------------

      SUBROUTINE fulstp(xo,xend,di,s,theta,rlog,p,qbig,w,sqrnrm,iswtch)
C This routine advances the solutions of both the Prufer r and theta equations
C over an interval [xo,xend] on which the coefficient functions are constant.
C It also stores in  sqrnrm the logarithm of the contribution to the square of
C the eigenfunction norm arising from integration across the interval (xo,xend).
C VARIABLES:
C
C xo, xend: the endpoints of the interval over which integration must proceed;
C           unchanged on exit.
C        s: the scale-factor in use when the routine is called; this is the
C           scale-factor appropriate to the input value of theta. On exit,
C           s will be reset to the output value of the scale-factor if
C           this is different from the input value
C    theta: the Prufer angle. On entry, theta must specify the value of the
C           Prufer angle at the point xo, normalised according to the input
C           value of the scale-factor s. On exit, theta contains the value of
C           the Prufer angle at the point xend, normalised according to the
C           output value of the scale-factor s.
C     rlog: the logarithm of the Prufer radius, normalised to be appropriate
C           to the scale-factor s. On input, rlog contains the logarithm of the
C           Prufer radius at the point xo. On output, rlog contains the
C           logarithm of the Prufer radius at the point xend.
C    p,w,qbig: the values of the coefficient functions p and w, and of
C           elam*w-q = qbig on the interval [xo,xend]. Unchanged on exit.
C       di: indicates whether xo is greater than or less than xend;
C           di = 1.0 if xo<xend, di = -1.0 if xo>xend.
C  iswtch: a switch parameter. With iswtch = 1 the Prufer theta-equation
C           and the Prufer r-equation are integrated, and with iswtch = -1
C           the norm is also accumulated.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION di,p,qbig,rlog,s,sqrnrm,theta,w,xend,xo
      INTEGER iswtch
C     ..
C     .. Scalars in Common ..
      INTEGER icall,igt0,ilt0,int,ip,ismall,ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alfa,alfah,c2th,c2thn,cth,h,hlog,oldtrm,pdy,
     &                 pdynew,phi2t,phit,pi,pi4,psit,rlgnew,rmess,snew,
     &                 spth,sqrts,sth,t,tah,term,thtnew,y,ynew
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION chi,phi,psi,rescal,scl,spexp,sptanh,x01aaf,x02amf
cc      EXTERNAL chi,phi,psi,rescal,scl,spexp,sptanh,x01aaf,x02amf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,atan2,cos,log,sign,sin,sqrt,tan
C     ..
C     .. External Functions ..
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival,igt0,ilt0,ismall
C     ..
C
C     .. Parameters ..
      DOUBLE PRECISION one,two,half,four,qtr,smal,zero
      PARAMETER (one=1.0D0,two=2.0D0,half=one/two,four=4.0D0,
     &          qtr=one/four,smal=1.0D-1,zero=0.0D0)
C     ..
      pi = x01aaf(zero)
      pi4 = pi*qtr

      h = xend - xo

      IF (h.NE.zero) THEN
          hlog = log(abs(h))
          t = h*h*qbig/p
C
          IF (abs(t).GT.smal) THEN
              snew = sqrt(abs(qbig*p))
              alfa = snew/p
              alfah = h*snew/p
              thtnew = scl(theta,snew/s)
              rlgnew = rescal(s,snew,thtnew,rlog)
              s = snew
              IF (t.GT.zero) THEN
                  igt0 = igt0 + 1
                  IF (iswtch.EQ.-1) THEN
C Compute the contribution to the square of the L2-norm also. Note the
C special precaution which has to be taken to deal with indefinite
C problems where w may be zero.
                      rmess = half*w*abs(one- (sin((thtnew+alfah)*two)-
     &                        sin(two*thtnew))* (half/alfah))/s
                      IF (rmess.GT.0.0D0) THEN
                          sqrnrm = log(rmess) + two*rlgnew + hlog

                      ELSE
C (n.b. x02amf(*) is a very small number, O(1e-39)).
                          sqrnrm = log(x02amf(one))
                      END IF

                  END IF
C Update theta and rlog:
                  theta = thtnew + alfah
                  rlog = rlgnew

              ELSE
                  ilt0 = ilt0 + 1
                  c2thn = cos(two*thtnew)
                  theta = scl(thtnew-pi4*di,spexp(-abs(alfah)*two)) +
     &                    pi4*di
                  c2th = cos(two*theta)
C Update rlog:
                  tah = cos(thtnew-pi4*di)
                  IF (abs(tah).GT.smal) THEN
                      tah = (tan(thtnew-pi4*di))**2
                      tah = (one+tah)/ (spexp(-abs(alfah)*four)*tah+one)
                      rlog = abs(alfah) + rlgnew - log(tah)*half

                  ELSE
                      oldtrm = log(abs(c2thn))
                      term = log(abs(c2th))
                      rlog = (oldtrm-term)*half + rlgnew
                  END IF

                  IF (iswtch.EQ.-1) THEN

                      cth = spexp(-abs(alfah))
                      cth = two*cth/ (one-cth**2)
                      spth = spexp(-two*abs(alfah))
                      sth = (sign(one,alfah)/alfa)* (one+spth)/
     &                      (one-spth) - cth**2*h
                      IF (rlgnew.LE.rlog) THEN
                          sth = half*sth* (spexp((rlgnew-rlog)*two)*
     &                          sin(thtnew)**2+sin(theta)**2)/s
                          sth = spexp(rlgnew-rlog)*sign(one,alfah)*
     &                          sin(theta)* (sin(thtnew)/s)*abs(cth)*
     &                          (h/sptanh(alfah)-one/alfa) + sth

                      ELSE
                          sth = half*sth* (sin(thtnew)**2+
     &                          spexp((rlog-rlgnew)*two)*sin(theta)**2)/
     &                          s
                          sth = spexp(rlog-rlgnew)*sign(one,alfah)*
     &                          sin(theta)* (sin(thtnew)/s)*abs(cth)*
     &                          (h/sptanh(alfah)-one/alfa) + sth
                      END IF

                      IF (abs(w*sth).GT.zero) THEN
                          sqrnrm = log(abs(w*sth)) + two*rlog
                          IF (rlgnew.GT.rlog) sqrnrm = sqrnrm +
     &                        two* (rlgnew-rlog)

                      ELSE
                          sqrnrm = log(x02amf(one))
                      END IF

                  END IF

              END IF

          ELSE
              ismall = ismall + 1
              sth = sin(theta)
              cth = cos(theta)
              phit = phi(-t)
              psit = psi(-t)
              sqrts = sqrt(s)
              IF (iswtch.EQ.-1) THEN
C Compute the contribution to the square of the L2-norm also.
                  phi2t = phit*psit
                  y = sth/sqrts
                  pdy = cth*sqrts
                  rmess = half*w*abs((pdy**2)* (one-phi2t)/ (p*qbig)+
     &                    (y**2)* (one+phi2t)+two*y*pdy* (phit**2)*h/p)
                  IF (rmess.GT.0.0D0) THEN
                      sqrnrm = log(rmess) + rlog*two + hlog

                  ELSE
                      sqrnrm = log(x02amf(one))
                  END IF

              END IF

              theta = atan2((s/p*cth**2+qbig/s*sth**2)*h,
     &                (s/p-qbig/s)*h*sth*cth+chi(-t)) + theta
C Update rlog:
              pdy = sqrts*cth
              y = sth/sqrts
              pdynew = pdy*psit - y*h*qbig*phit
              ynew = h*phit*pdy/p + y*psit
              rlog = log((pdynew**2)/s+ (ynew**2)*s)*half + rlog
          END IF

      END IF

      END SUBROUTINE fulstp

C----------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION rescal(s,snew,thtnew,rlog)
C     .. Scalar Arguments ..
      DOUBLE PRECISION rlog,s,snew,thtnew
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sinsq
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,log,sin
C     ..
      sinsq = sin(thtnew)**2
      rescal = rlog - log(abs(1.0D0-sinsq)*snew/s+sinsq*s/snew)*0.5D0

      END FUNCTION rescal
C----------------------------------------------------------------------------

      SUBROUTINE onestp(xo,xend,di,s,theta,p,qbig)
C This routine advances the solution of the scaled Prufer equation
C across an interval [xo,xend] on which the coefficient functions p, q,
C and w are constants. It will be called from the meshing routine, from the
C prufer routine below, and also from the routine which computes the
C eigenfunction.
C
C VARIABLES:
C
C xo, xend: the endpoints of the interval over which integration must proceed;
C           unchanged on exit.
C        s: the scale-factor in use when the routine is called; this is the
C            scale-factor appropriate to the input value of theta. On exit,
C           s will be reset to the output value of the scale-factor if
C           this is different from the input value
C    theta: the Prufer angle. On entry, theta must specify the value of the
C            Prufer angle at the point xo, normalised according to the input
C           value of the scale-factor s. On exit, theta contains the value of
C            the Prufer angle at the point xend, normalised according to the
C             output value of the scale-factor s.
C      p,Q: the values of the coefficient function p and of elam*w-q = Q
C           on the interval [xo,xend]. Unchanged on exit.
C       di: indicates whether xo is greater than or less than xend;
C           di = 1.0 if xo<xend, di = -1.0 if xo>xend.
C     .. Scalar Arguments ..
      DOUBLE PRECISION di,p,qbig,s,theta,xend,xo
C     ..
C     .. Scalars in Common ..
      INTEGER icall,igt0,ilt0,int,ip,ismall,ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alfah,cth,h,pi,pi4,snew,sqrtp,sqrtq,sth,t
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,atan2,cos,sin,sqrt
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION chi,scl,spexp,x01aaf
cc      EXTERNAL chi,scl,spexp,x01aaf
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival,igt0,ilt0,ismall
C     ..

      pi = x01aaf(0.0D0)
      pi4 = pi/4.0D0

      h = xend - xo
      t = h*h*qbig/p
C
      IF (abs(t).GT.0.1D0) THEN
          sqrtq = sqrt(abs(qbig))
          sqrtp = sqrt(p)
          snew = sqrtq*sqrtp
          alfah = (h*sqrtq)/sqrtp
          theta = scl(theta,snew/s)
          s = snew
          IF (t.GT.0.0D0) THEN
              igt0 = igt0 + 1
              theta = theta + alfah

          ELSE
              ilt0 = ilt0 + 1
              theta = scl(theta-pi4*di,spexp(-abs(alfah)*2.D0)) + pi4*di
          END IF

      ELSE
          ismall = ismall + 1
          sth = sin(theta)
          cth = cos(theta)
          theta = atan2((s/p*cth**2+qbig/s*sth**2)*h,
     &            (s/p-qbig/s)*h*sth*cth+chi(-t)) + theta
      END IF
      END SUBROUTINE onestp

C-----------------------------------------------------------------------

      SUBROUTINE prufer(i,theta,s,iend,send,elam,pp,qp,wp,xmesh)
C
C Integrates the Prufer phase (forward or backward depending on whether iend is
C larger or smaller than i) advancing the values of i,theta and s
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,s,send,theta
      INTEGER i,iend
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(*),qp(*),wp(*),xmesh(0:*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION di,p,q
      INTEGER ii,iold
C     ..
C     .. External Subroutines ..
cc      EXTERNAL onestp
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION scl
cc      EXTERNAL scl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sign
C     ..
      di = sign(1,iend-i)
      ii = i
      IF (di.LT.0.0D0) ii = i + 1
10    CONTINUE
C
      IF (i.NE.iend) THEN
          iold = i
          i = i + di
          ii = ii + di
          p = pp(ii)
C  Q here has same sign as in writeup of 12/1/89
          q = wp(ii)*elam - qp(ii)

          CALL onestp(xmesh(iold),xmesh(i),di,s,theta,p,q)

          GO TO 10

      END IF

      theta = scl(theta,send/s)
      s = send

      END SUBROUTINE prufer

C-----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION d(elam,nmesh,imatch,pp,qp,wp,xmesh,
     &                            setup,trans,iflag)
C        DOUBLE PRECISION FUNCTION d(elam,pp,qp,wp,xmesh,nmesh,
C     *      imatch,smatch)
C
C Computes the miss-distance with xmesh(imatch) taken as the matching point
C and smatch as the scale-factor there
C
C     .. Parameters ..
      DOUBLE PRECISION smatch
      PARAMETER (smatch=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam
      INTEGER iflag,imatch,nmesh
      LOGICAL trans
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(nmesh),qp(nmesh),wp(nmesh),xmesh(0:nmesh)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. Scalars in Common ..
      INTEGER icall,igt0,ilt0,int,ip,ismall,ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION pdyl,pdyr,pi,s,thetal,thetar,yl,yr
      INTEGER i,ifail
C     ..
C     .. External Subroutines ..
cc      EXTERNAL bcons,prufer
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC atan2
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x01aaf
cc      EXTERNAL x01aaf
C     ..
C     .. Common blocks ..
      COMMON icall,int,ip,ival,igt0,ilt0,ismall
C     ..
      pi = x01aaf(0.0D0)
      IF (nmesh.GT.99) ip = ip + 1
      iflag = 0
      ifail = 0
C Left leg:
      CALL bcons(setup,yl,pdyl,elam,xmesh(0),pp(1),qp(1),wp(1),0,trans,
     &           ifail)
      IF (ifail.NE.3) THEN
          i = 0
          thetal = atan2(yl,pdyl)
          s = 1.0D0
          CALL prufer(i,thetal,s,imatch,smatch,elam,pp,qp,wp,xmesh)
C Right leg:
          CALL bcons(setup,yr,pdyr,elam,xmesh(nmesh),pp(nmesh),
     &               qp(nmesh),wp(nmesh),1,trans,ifail)
          IF (ifail.NE.3) THEN
              i = nmesh
              thetar = pi - atan2(yr,pdyr)
              s = 1.0D0
              CALL prufer(i,thetar,s,imatch,smatch,elam,pp,qp,wp,xmesh)
C
              d = thetal - thetar
          END IF

      END IF

      IF (ifail.NE.0) iflag = ifail

      END FUNCTION d

C -------------------------------------------------------------------

      SUBROUTINE bcons(setup,y,pdy,elam,x,p,q,w,iend,trans,ifail)
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,p,pdy,q,w,x,y
      INTEGER iend,ifail
      LOGICAL trans
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL setup
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION qbig,s,xend
      LOGICAL ising
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
C     ..
C     .. Parameters ..
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
      IF (trans) THEN
          xend = x/ (one-abs(x))

      ELSE
          xend = x
      END IF

      CALL setup(y,pdy,elam,xend,iend,ising)
      IF (ising) THEN
C      This end of the interval is singular from the point of view of
C      boundary conditions.
          ifail = 0
          qbig = elam*w - q
          IF (qbig.LT.0.0D0) THEN
              y = 1.0D0
              s = -p*qbig
              pdy = sqrt(s)
              s = sqrt(y**2+s)
              y = y/s
              pdy = pdy/s

          ELSE
              s = 1.0D0/p
              IF (s.GT.qbig) THEN
                  y = 1.0D0
                  pdy = 0.0D0

              ELSE
                  y = 0.0D0
                  pdy = 1.0D0
              END IF

          END IF

      ELSE
C The end is regular.
          s = sqrt(y**2+pdy**2)
          IF (s.EQ.0.0D0) THEN
              ifail = 3

          ELSE
              y = y/s
              pdy = pdy/s
              IF (iend.EQ.1) pdy = -pdy
              IF (y.LT.0.0D0) THEN
                  y = -y
                  pdy = -pdy
              END IF

              IF (y.EQ.0.0D0 .AND. pdy.LT.0.0D0) pdy = -pdy

              ifail = 0

          END IF

      END IF

      END SUBROUTINE bcons
      SUBROUTINE mmesh(xo,xend,ic,elam,coeffn,rlog,theta,scale,xtemp,
     &                ptemp,qtemp,wtemp,ntemp,imatch,xmesh,pp,qp,wp,n,
     &                tol,npo,params,iparam,icofun,ifail)
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,tol,xend,xo
      INTEGER ic,icofun,ifail,imatch,iparam,n,npo,ntemp
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION params(1:iparam),pp(1:n),ptemp(1:ntemp),qp(1:n),
     &                 qtemp(1:ntemp),rlog(0:ntemp),scale(0:ntemp),
     &                 theta(0:ntemp),wp(1:n),wtemp(1:ntemp),xmesh(0:n),
     &                 xtemp(0:ntemp)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL coeffn
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION big,di,errest,errpm,errpp,errqm,errqp,errw1,
     &                 errw2,errwm,errwp,gamma,h,hmax,omega1,omega2,
     &                 omega3,p1,p2,p3,q1,q2,q3,ratio,toloc,tolow,w1,w2,
     &                 w3,weight,x1,x2,x3,xende,xoo
      INTEGER i,icfn,icnt,idi,idmy,ifin,ii,irf,kntr,npts
      LOGICAL infty,phase1,repeat
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION pdyval(0:1),xvals(0:1),yvals(0:1)
C     ..
C     .. External Subroutines ..
cc      EXTERNAL coarse,coefun,efun,initia
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf,x02amf
cc      EXTERNAL x02ajf,x02amf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,exp,int,log,max,min,sin,sqrt
C     ..
C     .. Parameters ..
      DOUBLE PRECISION hlf,two,one,safe,trd,five,ten,safe3,p008,sixth,
     &                 zero,tenth,twelf
      PARAMETER (hlf=5.d-1,two=2.d0,one=1.d0,safe=0.8D0,trd=one/3.d0,
     &          five=5.d0,ten=10.d0,safe3=7.29d-1,p008=8.d-3,
     &          sixth=one/6.d0,zero=0.d0,tenth=one/ten,twelf=hlf*sixth)
C     ..
C This subroutine covers [xo,xend] with a suitable mesh for eigenvalue
C computation.
      toloc = tol
      tolow = safe*toloc
      ifin = ifail
      ifail = 0
      i = ic
      ii = ic + 1
      di = one
      idi = 1
      icfn = icofun
      IF (icofun.EQ.1) THEN
          icfn = 2
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      ELSE IF (icofun.EQ.3) THEN
          icfn = 4
      END IF

      infty = .false.
C REMARK:RWB: In the NAG version of the code icofun NEVER has the value 3.
      IF (icofun.EQ.1 .OR. icofun.EQ.3) infty = .true.
      IF (xo.EQ.xend) RETURN
      xoo = xo
      xende = xend
      IF (infty) THEN
          xoo = xoo/ (one-abs(xoo))
          xende = xende/ (one-abs(xende))
      END IF

      IF (xoo.GT.xende) THEN
          i = n
          ii = n
          di = -one
          idi = -1
      END IF

      icnt = 0
      xmesh(i) = xoo
C Set the initial mesh-size; first setting:
      IF (imatch.GE.1) THEN
          h = abs(xtemp(imatch)-xtemp(imatch-1))

      ELSE
          h = abs(xtemp(imatch)-xtemp(imatch+1))
      END IF

C Second setting:
      big = abs((elam*wtemp(imatch)-qtemp(imatch))*
     &      (sin(theta(imatch))**2)/scale(imatch)+
     &      (one/ptemp(imatch))* (cos(theta(imatch))**2)*scale(imatch))*
     &      exp(two*rlog(imatch))

      hmax = abs(xend-xo)/npo
      h = min(hmax,h,tol/max(x02amf(one),big))
C Initial mesh-size has now been set.
C maximum no. of attempts to get right stepsize (set to avoid stepsize
C becoming less than machine precision x02ajf(*) in one step):
      irf = int(log(x02ajf(one))*hlf/log(safe))
      kntr = 0
      phase1 = .true.
      repeat = .false.
C
C This purpose of this routine is to provide a suitably adapted mesh for
C the solution of the Sturm-Liouville eigenvalue problem. Mesh adaptation
C is based on keeping an estimate of the error per interval below a specified
C tolerance.
C
C General Step (from i to i+idi):
10    x3 = xmesh(i) + h*di
      x1 = xmesh(i)
      IF ((x3-xende)*di.GT.zero) x3 = xende
      xmesh(i+idi) = x3
      x2 = hlf* (xmesh(i)+xmesh(i+idi))
      IF (.NOT.repeat .AND. .NOT.phase1) THEN
C If the step is a repeat after a failure then we dont update at the end
C from which we are shooting, otherwise we do.
          p1 = p3
          q1 = q3
          w1 = w3
      END IF
C We only need to evaluate the coefficients at the end from which we are
C shooting on the first step, otherwise they are available from the
C previous step and have been set above.
      IF (phase1) CALL coefun(x1,p1,q1,w1,coeffn,params,iparam,icfn)
      CALL coefun(x2,p2,q2,w2,coeffn,params,iparam,icfn)
      CALL coefun(x3,p3,q3,w3,coeffn,params,iparam,icfn)
      IF (x1.LE.x3) THEN
          errw1 = w1 - w2
          errw2 = w3 - w2
          errwm = elam*errw1
          errwp = elam*errw2
          errpp = (one/p3-one/p2)
          errpm = (one/p1-one/p2)
          errqp = (q3-q2)
          errqm = (q1-q2)
          xvals(0) = x1
          xvals(1) = x3

      ELSE
          errw2 = w1 - w2
          errw1 = w3 - w2
          errwp = elam*errw1
          errwm = elam*errw2
          errpm = (one/p3-one/p2)
          errpp = (one/p1-one/p2)
          errqm = (q3-q2)
          errqp = (q1-q2)
          xvals(0) = x3
          xvals(1) = x1
      END IF

C We only need to evaluate the eigenfunction at the end from which we are
C shooting on the first step, otherwise it is  available from the
C previous step:
      IF (phase1) THEN
          CALL efun(xvals,yvals,pdyval,elam,rlog,theta,scale,xtemp,
     &              ptemp,qtemp,wtemp,1,ntemp,imatch)

      ELSE
C Two cases to consider depending on stepping direction:
          IF (idi.GT.0) THEN
C Case of repeated step is special; the eigenfunction is not updated.
              IF (.NOT.repeat) THEN
                  yvals(0) = yvals(1)
                  pdyval(0) = pdyval(1)
              END IF

              CALL efun(xvals(1),yvals(1),pdyval(1),elam,rlog,theta,
     &                  scale,xtemp,ptemp,qtemp,wtemp,0,ntemp,imatch)

          ELSE
C Case of repeated step is special; the eigenfunction is not updated.
              IF (.NOT.repeat) THEN
                  yvals(1) = yvals(0)
                  pdyval(1) = pdyval(0)
              END IF

              CALL efun(xvals(0),yvals(0),pdyval(0),elam,rlog,theta,
     &                  scale,xtemp,ptemp,qtemp,wtemp,0,ntemp,imatch)
          END IF

      END IF
C END of eigenfunction evaluation.
C The most basic error monitor ensures that we get the right normalisation:
      errest = twelf*abs(h* (errw1*min(one,yvals(0)**2)+errw2*min(one,
     &         yvals(1)**2)))
      errwm = errwm - errqm
      errwp = errwp - errqp
C There are a number of other different error monitors to be used, depending
C on whether the behaviour of solutions of the differential equation.
C Compute Prufer radius for use:
      weight = hlf*sqrt(yvals(0)**2+pdyval(0)**2+yvals(1)**2+
     &         pdyval(1)**2)
C Measure rate of oscillation:
      omega2 = (elam*w2-q2)/p2
C WEIGHT may be quite small near the endpoints, where other eigenfunctions
C do not decay just as fast. Since weight is supposed to take account of these
C other eigenfunctions we adjust its value where appropriate:
      IF (omega2.LT.zero) THEN
          IF (q2.LT.two*max(one,abs(elam))*w2) THEN
C It is not clear that there is any representative rate of decay for nearby
C eigenfunctions: a small increase in elam could change everything.
              weight = max(weight,tenth)

          ELSE
C We can be reasonably confident that the rate of decay of nearby
C eigenfunctions will be at least half as fast as that of the present
C eigenfunction.
              gamma = sqrt((two*max(one,elam)*w2-q2)/ (elam*w2-q2))
              weight = max(weight,weight**gamma)
          END IF

      ELSE
          weight = max(weight,min(tenth,two*weight))
      END IF
C First error monitor, valid in regions which are not highly oscillatory.
C This controls the components of the eigenfunction which are not linearly
C dependent on the exact eigenfunction:
      errest = max(errest,sixth*weight*
     &         abs(h* (abs(yvals(0)*errwm+yvals(1)*errwp)/sqrt(max(one,
     &         omega2))+abs(errpm*pdyval(0)+errpp*pdyval(1)))))
C Second error monitor, uses Simpson's rule to measure a local contribution
C to the Green's formula for eigenvalue error:
      errest = max(errest,sixth*abs(h* ((yvals(0)**2)*errwm+
     &         (yvals(1)**2)*errwp+errpm* (pdyval(0)**2)+
     &         errpp* (pdyval(1)**2))))/max(one,abs(elam))
      IF ((h**2)*omega2.GT.hlf) THEN
C Safeguard in highly oscillatory regions: in such regions we can control
C eigenvalue error by controlling error in Prufer theta:
          omega1 = sqrt(abs(elam*w1-q1)/p1)
          omega3 = sqrt(abs(elam*w3-q3)/p3)
          omega2 = sqrt(omega2)
          errest = max(errest,abs(two* (p2/w2)*omega2* (omega1-
     &             two*omega2+omega3)*h)/max(one,elam))
      END IF
C All the error monitors here are O(h**3) for small h.

      repeat = .false.

      IF (errest.GT.toloc) THEN
          repeat = .true.
          ratio = safe
          IF (toloc.LT.errest*safe3) ratio = (toloc/errest)**trd
          h = h*ratio
          icnt = icnt + 1
          IF (icnt.LT.irf) GO TO 10
          ifail = 10
          RETURN

      END IF

      IF (errest.LT.tolow .AND. icnt.EQ.0) THEN
          ratio = five
          IF (errest.GT.tolow*p008) ratio = (tolow/errest)**trd
          h = h*ratio
          IF (infty) THEN
              h = min(h,hmax* (one+x2**2))

          ELSE
              h = min(h,hmax)
          END IF
C The next piece of code gets the meshing routine 'on-scale' at the first step.
          IF (phase1) THEN
              repeat = .true.
              kntr = kntr + 1
              IF (kntr.LT.5) GO TO 10
          END IF

      END IF
C General step has been completed successfully, prepare for next step.
      phase1 = .false.
      icnt = 0
      pp(ii) = p2
      qp(ii) = q2
      wp(ii) = w2
      i = i + idi
      ii = ii + idi
C CHECK THAT we have not reached the end:
      IF (xmesh(i)*di.GE.xende*di) GO TO 20
C If we have not reached the end then we must make sure that there is at least
C one more space in each array for us to continue.
      IF (ii.GE.n .OR. i.LE.ic) THEN
          ifail = 11
          IF (ifin.EQ.1) GO TO 50
          GO TO 40

      END IF
C Everything checked.
      GO TO 10
C We have now successfully completed the meshing. If we were meshing backwards
C then we must ensure that the mesh nodes obtained are indexed correctly.
20    IF (idi.LT.0) THEN
          xmesh(ic) = xende
          DO 30 idmy = 1,n - i
              xmesh(ic+idmy) = xmesh(i+idmy)
              pp(ic+idmy) = pp(i+idmy)
              qp(ic+idmy) = qp(i+idmy)
              wp(ic+idmy) = wp(i+idmy)
30        CONTINUE
          x1 = hlf* (xende+xmesh(ic+1))
          CALL coefun(x1,p1,q1,w1,coeffn,params,iparam,icfn)
          pp(ic+1) = p1
          qp(ic+1) = q1
          wp(ic+1) = w1
          n = ic + n - i

      ELSE
          n = i
          xmesh(n) = xende
          x1 = hlf* (xende+xmesh(n-1))
          CALL coefun(x1,p1,q1,w1,coeffn,params,iparam,icfn)
          pp(n) = p1
          qp(n) = q1
          wp(n) = w1
      END IF

40    RETURN
C Eventually we hope to have an emergency mesh installed here

50    IF (ic.EQ.0) THEN
          npts = min(n/2,400)
          IF (npts.LT.5) RETURN
          CALL coarse(npts,xende,xoo,xmesh,1)
          n = npts
          CALL initia(coeffn,pp,qp,wp,xmesh,n,params,iparam,icfn)

      ELSE IF (ic.GT.0) THEN
          npts = min((n-ic)/2,400)
          IF (npts.LT.5) RETURN
          CALL coarse(npts,xoo,xende,xmesh(ic),2)
          CALL initia(coeffn,pp(ic+1),qp(ic+1),wp(ic+1),xmesh(ic),npts,
     &                params,iparam,icfn)
          n = npts + ic

      ELSE
          RETURN

      END IF

      ifail = 12
      RETURN

      END SUBROUTINE mmesh
C ---------------------------------------------------------------------
C --------------------- Source from EXPECT ----------------------------
C ---------------------------------------------------------------------

      SUBROUTINE eigint(fun,tol,iopt,k,kvals,elam,expval,m,mvals,wk,
     &                  eigfns,iwk,n,ifail)
C August 1990: This routine has been superceded by SL05FM
C Brief specification:
C This routine computes mvals integrals of Fourier type if IOPT = 1,
C of expectation type if IOPT = 2; that is, we compute the integrals
C                 \int_{a}^{b}fun(x)(y_{Kvals(j)}(x))^{IOPT}dx   (*)
C VARIABLES:
C       FUN: REAL FUNCTION supplied by the user. Its specification is
C            FUNCTION FUN(x)
C            REAL x
C
C            and it must define FUN for all a <= x <=b, where [a,b] is the
C            interval on which the Sturm-Liouville problem is posed.
C       TOL: REAL. The routine will attempt to ensure that the approximations
C            to the integrals which are obtained are within TOL of the actual
C            value of the integral when the exact eigenfunctions are replaced
C            by the approximate eigenfunctions. However it is impossible to
C            guarantee that the answer will be within TOL of the actual
C            integral with the exact eigenfunction because the eigenfunctions
C            may be very ill-conditioned.
C            Unchanged on exit.
C      IOPT: INTEGER. The value of IOPT decides the type of integrals to be
C            computed (see (*) above). Unchanged on exit.
C    EXPVAL: REAL ARRAY of DIMENSION at least (0:1,1:m). On exit, EXPVAL(0,j)
C            (1 <= j <= m) contains an approximation to the integral in
C            (*) where the eigenfunction concerned is the K(j)th eigenfunction.
C            EXPVAL(1,j) contains an estimate of the error due to numerical
C            integration, to allow the user to assess whether or not TOL
C            was sufficiently small.
C     ELAM: REAL ARRAY of DIMENSION (0:1,1:p) where p >= m.  Supplied by the
C           preceding call to sl02fm. Unchanged on exit.
C        K: INTEGER ARRAY of DIMENSION at least (1:M). The array of indices
C           of the eigenvalues and eigenfunctions stored in the arrays EIGFNS
C           and ELAM. As supplied to  sl02fm for the computation of ELAM.
C           Unchanged on exit.
C    MVALS: INTEGER. The number of integrals to be computed. Unchanged on
C            exit.
C    KVALS: INTEGER ARRAY of DIMENSION at least (1:MVALS). The indices of the
C           eigenfunctions for which the integrals are to be computed. Every
C           integer which appears as an element of this array must also
C           appear as an element of the array K. Unchanged on exit.
C   EIGFNS: REAL ARRAY of DIMENSION (0:IWK,1:3,1:p) where p >= m. This array
C           must be supplied by the user and will be used to store the
C           information about the normalised eigenfunctions required by EIGINT.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4,1:p) where p >= m.
C            Supplied by the preceding call to the routine sl02fm. The
C            preceding call to sl02fm need not occur in the same subprogram,
C            but if it does not then the array wk must be declared to have the
C            same dimensions in the subprogram from which sl02fm is called as
C            in the subprogrm from which EIGINT is called.
C            Unchanged on exit.
C       IWK: INTEGER. The first dimension of WK as declared in the calling
C            (sub)program. Unchanged on exit.
C         N: INTEGER. The number of mesh intervals used by the routine sl02fm
C            in the computation of the eigenvalue approximation. Supplied by
C            the preceding call to sl02fm and unchanged on exit.
C         M: INTEGER. The number of eigenvalues computed at the preceding call
C            to sl02fm. As supplied to sl02fm by the user. Unchanged on exit.
C     IFAIL: INTEGER. On entry, IFAIL should be assigned the value 0 or 1. On
C            successful exit, IFAIL will have the value 0. The only failure
C            possible is that the routine may not be able to achieve the
C            accuracy TOL requested by the user, and then IFAIL will have the
C            value 1 on exit.
C END of brief specification.
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER ifail,iopt,iwk,m,mvals,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:iwk,1:3,1:m),elam(0:1,1:m),
     &                 expval(0:1,1:m),wk(0:iwk,1:4)
      INTEGER k(1:m),kvals(1:mvals)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam0
      INTEGER i,ic,j,jc
C     ..
C     .. External Subroutines ..
cc      EXTERNAL entgrl
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
      ifail = 0
      IF (iopt.LT.1 .OR. iopt.GT.2) THEN
          ifail = 1
          RETURN

      END IF

C Retrieve the matchpoint index from the sneaky storage:
      ic = wk(0,2)

      DO 30 i = 1,mvals
C Search for KVALS(i) in the array K
          DO 10 j = 1,m
              IF (kvals(i).EQ.k(j)) THEN
                  jc = j
                  GO TO 20

              END IF

10        CONTINUE

          ifail = 2
          RETURN

20        elam0 = elam(0,jc) + elam(1,jc)
          CALL entgrl(expval(0,i),expval(1,i),elam0,fun,tol,wk(0,1),
     &                wk(1,2),wk(1,3),wk(1,4),eigfns(0,1,jc),
     &                eigfns(0,2,jc),eigfns(0,3,jc),n,ic,iopt,ifail)
C ifail is only non-zero on exit if an error has been detected.
          IF (ifail.NE.0) RETURN

30    CONTINUE
C
C All the required integrals have now been computed
      RETURN

      END SUBROUTINE eigint

C --------------------------------------------------------------------

      SUBROUTINE sl05fm(fun,tol,iopt,k,kvals,elam,expval,m,mvals,wk,
     &                  eigfns,iwk,n,ifail)
C Brief specification:
C This routine computes mvals integrals of Fourier type if IOPT = 1,
C of expectation type if IOPT = 2; that is, we compute the integrals
C                 \int_{a}^{b}fun(x)(y_{Kvals(j)}(x))^{IOPT}dx   (*)
C VARIABLES:
C       FUN: REAL FUNCTION supplied by the user. Its specification is
C            FUNCTION FUN(x)
C            REAL x
C
C            and it must define FUN for all a <= x <=b, where [a,b] is the
C            interval on which the Sturm-Liouville problem is posed.
C       TOL: REAL. The routine will attempt to ensure that the approximations
C            to the integrals which are obtained are within TOL of the actual
C            value of the integral when the exact eigenfunctions are replaced
C            by the approximate eigenfunctions. However it is impossible to
C            guarantee that the answer will be within TOL of the actual
C            integral with the exact eigenfunction because the eigenfunctions
C            may be very ill-conditioned.
C            Unchanged on exit.
C      IOPT: INTEGER. The value of IOPT decides the type of integrals to be
C            computed (see (*) above). Unchanged on exit.
C    EXPVAL: REAL ARRAY of DIMENSION at least (0:1,1:m). On exit, EXPVAL(0,j)
C            (1 <= j <= m) contains an approximation to the integral in
C            (*) where the eigenfunction concerned is the K(j)th
C            eigenfunction, while EXPVAL(1,j) contains an estimate of the
C            error in EXPVAL(0,j) due to numerical integration, allowing
C            the user to decide whether or not TOL was sufficiently small.
C     ELAM: REAL ARRAY of DIMENSION (0:1,1:p) where p >= m.  Supplied by the
C           preceding call to sl02fm. Unchanged on exit.
C        K: INTEGER ARRAY of DIMENSION at least (1:M). The array of indices
C           of the eigenvalues and eigenfunctions stored in the arrays EIGFNS
C           and ELAM. As supplied to  sl02fm for the computation of ELAM.
C           Unchanged on exit.
C    MVALS: INTEGER. The number of integrals to be computed. Unchanged on
C            exit.
C    KVALS: INTEGER ARRAY of DIMENSION at least (1:MVALS). The indices of the
C           eigenfunctions for which the integrals are to be computed. Every
C           integer which appears as an element of this array must also
C           appear as an element of the array K. Unchanged on exit.
C   EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3,1:p) where p >= m. This array
C           must be supplied by the user and will be used to store the
C           information about the normalised eigenfunctions required by SL05FM.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4)
C            Supplied by the preceding call to the routine sl02fm. The
C            preceding call to sl02fm need not occur in the same subprogram,
C            but if it does not then the array wk must be declared to have the
C            same dimensions in the subprogram from which sl02fm is called as
C            in the subprogrm from which EIGINT is called.
C            Unchanged on exit.
C       IWK: INTEGER. The first dimension of WK as declared in the calling
C            (sub)program. Unchanged on exit.
C         N: INTEGER. The number of mesh intervals used by the routine sl02fm
C            in the computation of the eigenvalue approximation. Supplied by
C            the preceding call to sl02fm and unchanged on exit.
C         M: INTEGER. The number of eigenvalues computed at the preceding call
C            to SL02FM. As supplied to sl02fm by the user. Unchanged on exit.
C     IFAIL: INTEGER. On entry, IFAIL should be assigned the value 0 or 1. On
C            successful exit, IFAIL will have the value 0.
C END of brief specification.
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER ifail,iopt,iwk,m,mvals,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3,1:m),elam(0:1,1:m),
     &                 expval(0:1,1:m),wk(0:iwk,1:4)
      INTEGER k(1:m),kvals(1:mvals)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam0
      INTEGER i,ic,iflag,j,jc,nrec
      CHARACTER srname*6
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. External Subroutines ..
cc      EXTERNAL nntgrl
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      srname = 'SL05FM'
      IF (iopt.LT.1 .OR. iopt.GT.2 .OR. abs(ifail).GT.1 .OR.
     &    tol.LE.0.0D0 .OR. mvals.GT.m .OR. mvals.LT.0) THEN
          iflag = 1
          IF (abs(ifail).GT.1) ifail = -1
          WRITE (rec,FMT=9000)

9000      FORMAT (
     &  '** Parameter error: IOPT or IFAIL or TOL or MVALS out of range'
     &           )

          nrec = 1
          GO TO 40

      END IF

C Retrieve the matchpoint index from the sneaky storage:
      ic = wk(0,2)

      DO 30 i = 1,mvals
C Search for KVALS(i) in the array K
          DO 10 j = 1,m
              IF (kvals(i).EQ.k(j)) THEN
                  jc = j
                  GO TO 20

              END IF

10        CONTINUE

          iflag = 2
          WRITE (rec,FMT=9010)

9010      FORMAT ('** Required eigenfunctions not available')

          nrec = 1
          GO TO 40

20        elam0 = elam(0,jc) + elam(1,jc)
          CALL nntgrl(expval(0,i),expval(1,i),elam0,fun,tol,wk(0,1),
     &                wk(1,2),wk(1,3),wk(1,4),eigfns(0,1,jc),
     &                eigfns(0,2,jc),eigfns(0,3,jc),n,ic,iopt,iflag)
C ifail is only non-zero on exit if an error has been detected.
          IF (iflag.NE.0) THEN
              iflag = 3
              WRITE (rec,FMT=9020)

9020          FORMAT ('** Automatic step-size control has failed',/,
     &               '** Try increasing TOL')

              nrec = 2
              GO TO 40

          END IF

30    CONTINUE
C
C All the required integrals have now been computed
      ifail = 0
      RETURN

40    ifail = p01abf(ifail,iflag,srname,nrec,rec)
      RETURN

      END SUBROUTINE sl05fm

C --------------------------------------------------------------------
      SUBROUTINE sl05f(fun,tol,iopt,k,elam,expval,wk,eigfns,iwk,n,ifail)
C Brief specification:
C This routine computes m integrals of Fourier type if IOPT = 1,
C of expectation type if IOPT = 2; that is, we compute the integrals
C                 \int_{a}^{b}fun(x)(y_{K(j)}(x))^{IOPT}dx   (*)
C VARIABLES:
C       FUN: REAL FUNCTION supplied by the user. Its specification is
C            FUNCTION FUN(x)
C            REAL x
C
C            and it must define FUN for all a <= x <=b, where [a,b] is the
C            interval on which the Sturm-Liouville problem is posed.
C       TOL: REAL. The routine will attempt to ensure that the approximation
C            to the integral which is obtained is within TOL of the actual
C            value of the integral when the exact eigenfunction is  replaced
C            by the approximate eigenfunction. However it is impossible to
C            guarantee that the answer will be within TOL of the actual
C            integral with the exact eigenfunction because the eigenfunction
C            may be very ill-conditioned.
C            Unchanged on exit.
C      IOPT: INTEGER. The value of IOPT decides the type of integral to be
C            computed (see (*) above). Unchanged on exit.
C    EXPVAL: REAL.  On exit, EXPVAL(0) contains an approximation to the
C            integral in (*) where the eigenfunction concerned is the Kth
C            eigenfunction. EXPVAL(1) contains an estimate of the error in
C            EXPVAL(0) due to numerical integration, allowing the user to
C            decide whether or not TOL was sufficiently small.
C     ELAM: REAL ARRAY of DIMENSION (0:1).  Supplied by the
C           preceding call to sl02f. Unchanged on exit.
C        K: INTEGER. The index of the eigenvalue and eigenfunction in question.
C           As supplied to sl02f for the computation of ELAM. Unchanged on
C           exit.
C   EIGFNS: REAL ARRAY of DIMENSION (0:1,1:3). This array
C           must be supplied by the user and will be used to store the
C           information about the normalised eigenfunctions required by SL05F.
C       wk: REAL ARRAY of DIMENSION (0:IWK,1:4).
C            Supplied by the preceding call to the routine sl02f. The
C            preceding call to sl02f need not occur in the same subprogram,
C            but if it does not then the array wk must be declared to have the
C            same dimensions in the subprogram from which sl02f is called as
C            in the subprogram from which SL05F is called.
C            Unchanged on exit.
C       IWK: INTEGER. The first dimension of WK as declared in the calling
C            (sub)program. Unchanged on exit.
C         N: INTEGER. The number of mesh intervals used by the routine sl02f
C            in the computation of the eigenvalue approximation. Supplied by
C            the preceding call to sl02f and unchanged on exit.
C     IFAIL: INTEGER. On entry, IFAIL should be assigned the value 0 or 1. On
C            successful exit, IFAIL will have the value 0.
C END of brief specification.
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER ifail,iopt,iwk,k,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3),elam(0:1),expval(0:1),
     &                 wk(0:iwk,1:4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam0
      INTEGER ic,iflag,nrec
      CHARACTER srname*6
C     ..
C     .. External Subroutines ..
cc      EXTERNAL nntgrl
C     ..
C     ..  External Functions ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      srname = ' SL05F'

      IF (iopt.LT.1 .OR. iopt.GT.2 .OR. abs(ifail).GT.1 .OR.
     &    tol.LE.0.0D0) THEN
          iflag = 1
          IF (abs(ifail).GT.1) ifail = -1
          WRITE (rec,FMT=9000)

9000      FORMAT (
     &           '** Parameter error: IOPT or IFAIL or TOL out of range'
     &           )

          nrec = 1
          GO TO 10

      END IF

C Retrieve the matchpoint index from the sneaky storage
      ic = wk(0,2)

      elam0 = elam(0) + elam(1)
      CALL nntgrl(expval(0),expval(1),elam0,fun,tol,wk(0,1),wk(1,2),
     &            wk(1,3),wk(1,4),eigfns(0,1),eigfns(0,2),eigfns(0,3),n,
     &            ic,iopt,iflag)
C ifail is only non-zero on exit if an error has been detected.
      IF (iflag.NE.0) THEN
          iflag = 2
          WRITE (rec,FMT=9010)

9010      FORMAT ('** Automatic step-size control has failed',/,
     &           '** Try increasing TOL')

          nrec = 2
          GO TO 10

      END IF

      ifail = 0
C
C All the required integrals have now been computed
      RETURN

10    ifail = p01abf(ifail,iflag,srname,nrec,rec)
      RETURN

      END SUBROUTINE sl05f
C --------------------------------------------------------------------
      SUBROUTINE entgrl(expval,errest,elam,fun,tol,xmesh,pp,qp,wp,rlog,
     &                  theta,scale,n,ic,iopt,ifail)
C     .. Parameters ..
      DOUBLE PRECISION one,two,half,qtr,safe,tosafe,bound,zero,
     &                 ffteen
      INTEGER nvec
      PARAMETER (one=1.D0,two=2.D0,half=one/two,
     &          qtr=half*half,safe=0.95D0,tosafe=0.8D0,bound=2.D-1,
     &          zero=0.D0,ffteen=1.5D1,nvec=1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,errest,expval,tol
      INTEGER ic,ifail,iopt,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(1:n),qp(1:n),rlog(0:n),scale(0:n),theta(0:n),
     &                 wp(1:n),xmesh(0:n)
C     ..
C     .. Subroutine Arguments ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION abserr,ans4,di,epsilo,err,floc,h,hmax,p,q,ratio,
     &                 rl,rln,sc,scn,tend,th,thn,tmid,to,toloc,value,w,
     &                 xend,xfin,xmid,xo
      INTEGER i,icount,idi,ifin,ii,init,ip1,itime,nstep
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ans1(1:nvec),ans2(1:nvec),ans3(1:nvec),
     &                 fend(1:nvec),fmid(1:nvec),fo(1:nvec)
C     ..
C     .. External Subroutines ..
cc      EXTERNAL quad
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf
cc      EXTERNAL x02ajf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
      epsilo = x02ajf(1.D0)
      toloc = safe*tol
      value = zero
      errest = zero
      nstep = 0

      itime = 1
      init = 0
      ifin = ic - 1
      idi = 1
      di = one
      icount = 0
      h = xmesh(1) - xmesh(0)

10    DO 40 i = init,ifin,idi
          xo = xmesh(i)
          ip1 = i + idi
          hmax = xmesh(ip1) - xmesh(i)
          h = di*min(abs(h),abs(hmax))
          xfin = xmesh(ip1)
          xend = xo + h
          ii = i
          IF (idi.GT.0) ii = ip1
          p = pp(ii)
          q = qp(ii)
          w = wp(ii)
          rl = rlog(i)
          th = theta(i)
          sc = scale(i)

20        rln = rl
          thn = th
          scn = sc
C Do one step
          xmid = half* (xo+xend)
          to = xo
          tmid = xmid
          tend = xend
          fo(1) = fun(to)
          fmid(1) = fun(tmid)
          fend(1) = fun(tend)
          CALL quad(xo,xend,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,iopt,
     &              ans1)
C Do two half-steps
          rl = rln
          th = thn
          sc = scn
          floc = fend(1)
          fend(1) = fmid(1)
          to = half* (xo+xmid)
          fmid(1) = fun(to)
          CALL quad(xo,xmid,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,iopt,
     &              ans2)
          fo(1) = fend(1)
          fend(1) = floc
          to = half* (xo+xmid)
          fmid(1) = fun(to)
          CALL quad(xmid,xend,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,
     &              iopt,ans3)
          ans4 = ans2(1) + ans3(1)
C Estimate the error
          err = (ans4-ans1(1))/ffteen
          errest = errest + err
          abserr = abs(err)

          IF (abserr.GT.toloc) THEN
              ratio = min(half, (toloc/abserr)**qtr)
C It is conceivable that a failed step may occur when we are stepping exactly
C to the end xfin, in which case (xend-xo) \neq h. In such a case we must
C set h = ratio*(xend-xo) for the repeat step, rather than h = ratio*h.
C In other cases xend-xo=h and there is nothing to worry about.
              h = ratio* (xend-xo)
              xend = xo + h
              icount = icount + 1
C At most ten attempts to get the stepsize right:
              IF (icount.LT.10) THEN
                  rl = rln
                  th = thn
                  sc = scn
                  GO TO 20

              END IF

              ifail = 1
              GO TO 50

          ELSE
              value = value + ans4 + err
              IF (((xend-xfin)*di+epsilo).GE.zero) GO TO 30
              IF (abserr.LT.toloc*tosafe .AND. icount.EQ.0) THEN
C The error is very small and the step before last was not a failure so
C we may increase the stepsize.
                  ratio = one/max(bound, (abserr/toloc)**qtr)
                  h = ratio*h
              END IF
C We have completed the last step successfully, and provided the step before
C was not a failure we may have increased the step length.
              icount = 0
              xo = xend
              xend = xo + h
              IF (((xend-xfin)*di+epsilo).GE.zero) xend = xfin
              GO TO 20

          END IF

C Have successfully completed the integration from xmesh(i) to xmesh(i+idi)

30        nstep = nstep + 1
          icount = 0
40    CONTINUE

C If itime = 1 then we have only done the shooting from xmesh(0) to
C xmesh(ic), and we must therefore do the shooting from xmesh(n) to
C xmesh(ic).

      IF (itime.EQ.1) THEN
          itime = 2
          idi = -1
          di = -one
          init = n
          ifin = ic + 1
          h = xmesh(n) - xmesh(n-1)
          GO TO 10

      END IF

C Successful completion:
      ifail = 0
      expval = value
C      PRINT *,'ENTGRL error estimate:',errest
      RETURN
C Unsuccessful completion:
50    expval = zero
      RETURN

      END SUBROUTINE entgrl
C ---------------------------------------------------------------------

      SUBROUTINE nntgrl(expval,errest,elam,fun,tol,xmesh,pp,qp,wp,rlog,
     &                  theta,scale,n,ic,iopt,ifail)

C     .. Parameters ..
      DOUBLE PRECISION one,two,half,qtr,safe,tosafe,bound,zero,
     &                 ffteen
      INTEGER nvec
      PARAMETER (one=1.D0,two=2.D0,half=one/two,
     &          qtr=half*half,safe=0.95D0,tosafe=0.8D0,bound=2.D-1,
     &          zero=0.D0,ffteen=1.5D1,nvec=1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,errest,expval,tol
      INTEGER ic,ifail,iopt,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp(1:n),qp(1:n),rlog(0:1),scale(0:1),theta(0:1),
     &                 wp(1:n),xmesh(0:n)
C     ..
C     .. Subroutine Arguments ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION abserr,ans4,di,epsilo,err,floc,h,hmax,p,q,ratio,
     &                 rl,rln,sc,scn,tend,th,thn,tmid,to,toloc,value,w,
     &                 xend,xfin,xmid,xo
      INTEGER i,icount,idi,ifin,ii,init,ip1,itime,nstep
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ans1(1:nvec),ans2(1:nvec),ans3(1:nvec),
     &                 fend(1:nvec),fmid(1:nvec),fo(1:nvec)
C     ..
C     .. External Subroutines ..
cc      EXTERNAL quad
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf
cc      EXTERNAL x02ajf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
      epsilo = x02ajf(1.D0)
      toloc = safe*tol
      value = zero
      errest = zero
      nstep = 0

      itime = 1
      init = 0
      ifin = ic - 1
      idi = 1
      di = one
      icount = 0
      h = xmesh(1) - xmesh(0)

      rl = rlog(0)
      th = theta(0)
      sc = scale(0)

10    DO 40 i = init,ifin,idi
          xo = xmesh(i)
          ip1 = i + idi
          hmax = xmesh(ip1) - xmesh(i)
          h = di*min(abs(h),abs(hmax))
          xfin = xmesh(ip1)
          xend = xo + h
          ii = i
          IF (idi.GT.0) ii = ip1
          p = pp(ii)
          q = qp(ii)
          w = wp(ii)
C These three lines removed to comply with reduced storage requirements
C          rl = rlog(i)
C          th = theta(i)
C          sc = scale(i)
C
20        rln = rl
          thn = th
          scn = sc
C Do one step
          xmid = half* (xo+xend)
          to = xo
          tmid = xmid
          tend = xend
          fo(1) = fun(to)
          fmid(1) = fun(tmid)
          fend(1) = fun(tend)
          CALL quad(xo,xend,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,iopt,
     &              ans1)
C Do two half-steps
          rl = rln
          th = thn
          sc = scn
          floc = fend(1)
          fend(1) = fmid(1)
          to = half* (xo+xmid)
          fmid(1) = fun(to)
          CALL quad(xo,xmid,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,iopt,
     &              ans2)
          fo(1) = fend(1)
          fend(1) = floc
          to = half* (xend+xmid)
          fmid(1) = fun(to)
          CALL quad(xmid,xend,p,q,w,elam,rl,th,sc,fo,fmid,fend,nvec,
     &              iopt,ans3)
          ans4 = ans2(1) + ans3(1)
C Estimate the error
          err = (ans4-ans1(1))/ffteen
          errest = errest + err
          abserr = abs(err)

          IF (abserr.GT.toloc) THEN
              ratio = min(half, (toloc/abserr)**qtr)
C It is conceivable that a failed step may occur when we are stepping exactly
C to the end xfin, in which case (xend-xo) \neq h. In such a case we must
C set h = ratio*(xend-xo) for the repeat step, rather than h = ratio*h.
C In other cases xend-xo=h and there is nothing to worry about.
              h = ratio* (xend-xo)
              xend = xo + h
              icount = icount + 1
C At most ten attempts to get the stepsize right:
              IF (icount.LT.10) THEN
                  rl = rln
                  th = thn
                  sc = scn
                  GO TO 20

              END IF

              ifail = 1
              GO TO 50

          ELSE
              value = value + ans4 + err
              IF (((xend-xfin)*di+epsilo).GE.zero) GO TO 30
              IF (abserr.LT.toloc*tosafe .AND. icount.EQ.0) THEN
C The error is very small and the step before last was not a failure so
C we may increase the stepsize.
                  ratio = one/max(bound, (abserr/toloc)**qtr)
                  h = ratio*h
              END IF
C We have completed the last step successfully, and provided the step before
C was not a failure we may have increased the step length.
              icount = 0
              xo = xend
              xend = xo + h
              IF (((xend-xfin)*di+epsilo).GE.zero) xend = xfin
              GO TO 20

          END IF

C Have successfully completed the integration from xmesh(i) to xmesh(i+idi)

30        nstep = nstep + 1
          icount = 0
40    CONTINUE

C If itime = 1 then we have only done the shooting from xmesh(0) to
C xmesh(ic), and we must therefore do the shooting from xmesh(n) to
C xmesh(ic).

      IF (itime.EQ.1) THEN
          itime = 2
          idi = -1
          di = -one
          init = n
          ifin = ic + 1
          h = xmesh(n) - xmesh(n-1)
          rl = rlog(1)
          th = theta(1)
          sc = scale(1)
          GO TO 10

      END IF

C Successful completion:
      ifail = 0
      expval = value
C      PRINT *,'ENTGRL error estimate:',errest
      RETURN
C Unsuccessful completion:
50    expval = zero
      RETURN

      END SUBROUTINE nntgrl

C -------------------------------------------------------------------

      SUBROUTINE quad(xo,xend,p,q,w,elam,rlog,theta,scale,fo,fmid,fend,
     &                nvec,iopt,ans)
C     .. Parameters ..
      DOUBLE PRECISION one,two,half,qtr,three,twtrds,trd,sxteen
      PARAMETER (one=1.D0,two=2.D0,half=one/two,qtr=half*half,
     &          three=3.D0,twtrds=two/three,trd=one/three,sxteen=16.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam,p,q,rlog,scale,theta,w,xend,xo
      INTEGER iopt,nvec
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ans(1:nvec),fend(1:nvec),fmid(1:nvec),fo(1:nvec)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alfa,alfa2,alfa3,alfah,calfah,chqtrt,cntrb1,
     &                 cntrb2,cntrb3,cth,cthn,cthnp1,dummy,erl,f0,f1,f2,
     &                 gqtrt,h,onovsh,pdy,pdynew,phit,psit,qbig,
     &                 qtrt,rln,s2thn,s2thp1,salfah,snew,sqrts,sth,sthn,
     &                 sthnp1,t,talfah,term,thn,thnp1,x1,x2,xit,y,
     &                 yend,ynew,yo
      INTEGER i
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION chi,g,phi,psi,rescal,scl,x02ajf,xi
cc      EXTERNAL chi,g,phi,psi,rescal,scl,x02ajf,xi
C     ..
C     .. External Subroutines ..
cc      EXTERNAL fulstp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,atan2,cos,exp,log,sign,sin,sqrt,tanh
C     ..
      DO 10 i = 1,nvec
          ans(i) = 0.D0
10    CONTINUE

      h = xend - xo
      IF (abs(h).LE.x02ajf(1.D0)) RETURN
C      xmid = half* (xo+xend)
      qbig = elam*w - q
      t = (h**2)*qbig/p
      qtrt = t*qtr
C Interpolate to f using a quadratic expression
C f(x) = f0 + f1*(x-xmid) + f2*(x-xmid)**2
C where the coefficients are given by

C END of interpolation to f.

      IF (iopt.EQ.2) THEN

C  We require an expectation integral. There are three terms in the integral
C  arising from the three parts which approximate to f, and there are three
C  cases depending on the size of the parameter t.

          IF (abs(t).GT.0.1D0) THEN
              alfa2 = abs(qbig/p)
              alfa = sqrt(alfa2)
              alfah = alfa*h
              IF (t.GT.0.1D0) THEN
                  snew = sqrt(abs(qbig*p))
                  thn = scl(theta,snew/scale)
                  rln = rescal(scale,snew,thn,rlog)
C  In this case the eigenfunction has a very simple expression and all
C  three terms in the integral can be computed without divide by zero problems.
C 1: The constant contribution
                  thnp1 = thn + alfah
                  s2thn = sin(two*thn)
                  s2thp1 = sin(two*thnp1)
                  cntrb1 = half*h - qtr* (s2thp1-s2thn)/alfa
C 2: The linear contribution
                  salfah = sin(alfah)
                  calfah = cos(alfah)
                  cntrb2 = h*sin(two*thn+alfah)* (salfah/alfah-calfah)*
     &                     qtr/alfa
C 3: The quadratic contribution
                  cntrb3 = ((h**3)/sxteen)* (twtrds-
     &                     two*cos(two+thn+alfah)*
     &                     (salfah/alfah+two*calfah/abs(t)-
     &                     two*salfah/ (alfah**3)))
C 4: The entire contribution
                  DO 20 i = 1,nvec
                      f0 = fmid(i)
                      f1 = (fend(i)-fo(i))/h
                      f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                      ans(i) = w*exp(two*rln)*
     &                         (f0*cntrb1+f1*cntrb2+f2*cntrb3)/snew
20                CONTINUE
                  theta = thnp1
                  scale = snew
                  rlog = rln

              ELSE IF (t.LT. (-0.1D0)) THEN
C  In this case the integral is a pain to calculate because if it is not done
C  properly then there is the danger of floating overflow problems.
C  We need to know the value of the eigenfunction at each end of the interval
                  yo = exp(rlog)*sin(theta)/sqrt(scale)
                  x1 = xo
                  x2 = xend
C This call also does the updating of rlog, theta and scale:
                  CALL fulstp(x1,x2,sign(one,h),scale,theta,rlog,p,qbig,
     &                        w,dummy,-1)
                  yend = exp(rlog)*sin(theta)/sqrt(scale)
C  We can now proceed to compute each contribution in turn.
C 1: The constant contribution
C
                  cntrb1 = sign(one,h)*exp(dummy)
C 2: The linear contribution
                  talfah = tanh(alfah)
                  cntrb2 = qtr* (h**2)* (yend-yo)* (yend+yo)*
     &                     (one/talfah-one/alfah)/alfah
C 3: The quadratic contribution
                  onovsh = sign(one,alfah)*sqrt((one/talfah**2)-one)
                  term = one/alfah - two/ ((alfah**2)*talfah) +
     &                   two/ (alfah**3)
                  cntrb3 = ((half*h)**3)* ((yo**2+yend**2)*
     &                     (-trd* (onovsh**2)+ (one/talfah)*term)+
     &                     two*yo*yend* (trd/ (talfah)-term)*onovsh)
C 4: The full contribution
                  DO 30 i = 1,nvec
                      f0 = fmid(i)
                      f1 = (fend(i)-fo(i))/h
                      f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                      ans(i) = f0*cntrb1 + w* (f1*cntrb2+f2*cntrb3)
30                CONTINUE
C END of case t < -0.1
              END IF
C END of case ABS(t) > 0.1
          ELSE
C  Case ABS(t) <=0.1. This is by far the worst case since series expansions
C  of various terms are required for small alfah. Compute the contributions
C  from each of the terms involved.
              yo = exp(rlog)*sin(theta)/sqrt(scale)
              x1 = xo
              x2 = xend
              CALL fulstp(x1,x2,sign(one,h),scale,theta,rlog,p,qbig,w,
     &                    dummy,-1)
              yend = exp(rlog)*sin(theta)/sqrt(scale)
C 1: The constant term. We do this, and compute the values at the ends of the
C                    interval at the same time, by calling the fulstp routine.
              cntrb1 = sign(one,h)*exp(dummy)
C 2: The linear contribution
              cntrb2 = ((half*h)**2)* (yend-yo)* (yend+yo)*g(-t)
C 3: The quadratic contribution (Uugh!)
              xit = xi(t)
              cntrb3 = ((half*h)**3)* (((yend-yo)**2)*xit/ (phi(-t)**2)+
     &                  ((yend+yo)**2)/ (three* (one+psi(-t)))-
     &                 half* (yo**2+yend**2)*t*xit)
C 4: The full contribution:
              DO 40 i = 1,nvec
                  f0 = fmid(i)
                  f1 = (fend(i)-fo(i))/h
                  f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                  ans(i) = f0*cntrb1 + w* (f1*cntrb2+f2*cntrb3)
40            CONTINUE
C END of cases for iopt = 2
          END IF

      ELSE IF (iopt.EQ.1) THEN
C This case is not so easy -- we require a Fourier integral.
          IF (abs(t).GT.0.1D0) THEN
              snew = sqrt(abs(qbig*p))
              thn = scl(theta,snew/scale)
              rln = rescal(scale,snew,thn,rlog)
              alfa2 = abs(qbig/p)
              alfa = sqrt(alfa2)
              alfa3 = alfa*alfa2
              alfah = sqrt(abs(t))
              IF (t.GT.0.1D0) THEN
                  thnp1 = thn + alfah
C Compute the contributions from each of the terms approximating f
                  cthn = cos(thn)
                  cthnp1 = cos(thnp1)
                  sthn = sin(thn)
                  sthnp1 = sin(thnp1)
C 1: The constant term
                  cntrb1 = (cthn-cthnp1)/alfa
C 2: The linear term
                  cntrb2 = h* (cthnp1+cthn)*half/alfa -
     &                     (sthnp1-sthn)/ (alfa2)
C 3: The quadratic term
                  cntrb3 = (h**2)*qtr* (cthn-cthnp1)/alfa +
     &                     h* (sthnp1+sthn)/ (alfa2) -
     &                     two* (cthn-cthnp1)/ (alfa3)
C 4: The full contribution
                  DO 50 i = 1,nvec
                      f0 = fmid(i)
                      f1 = (fend(i)-fo(i))/h
                      f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                      ans(i) = w*exp(rln)* (f0*cntrb1+f1*cntrb2+
     &                         f2*cntrb3)/sqrt(snew)
50                CONTINUE
                  thn = thnp1

              ELSE
C             ** t < -0.1 **
C First compute the values of the eigenfunction at x = xo and xend.
                  yo = exp(rln)*sin(thn)/sqrt(snew)
                  x1 = xo
                  x2 = xend
C Advance the values of rln,thn, and the scalefactor snew.
                  CALL fulstp(x1,x2,sign(one,h),snew,thn,rln,p,qbig,w,
     &                        dummy,1)
C Compute the eigenfunction at the point xend
                  yend = exp(rln)*sin(thn)/sqrt(snew)

C Compute the contribution from each of these approximating terms
C 1: The Constant Term
                  talfah = tanh(half*alfa*h)
                  cntrb1 = (yend+yo)*talfah/alfa
C 2: The Linear Term
                  cntrb2 = (yend-yo)* (half*h/talfah-one/alfa)/alfa
C 3: The Quadratic Term
                  cntrb3 = two* (yo+yend)* ((((half*h)**2)+two/alfa2)*
     &                     talfah-h/alfa)/alfa
C 4: The full contribution
                  DO 60 i = 1,nvec
                      f0 = fmid(i)
                      f1 = (fend(i)-fo(i))/h
                      f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                      ans(i) = w* (f0*cntrb1+f1*cntrb2+f2*cntrb3)
60                CONTINUE
              END IF

              rlog = rln
              scale = snew
              theta = thn
C END of case where ABS(t) > 0.1
          ELSE
C The case where ABS(t) <= 0.1
C              hlog = log(abs(h))
              sth = sin(theta)
              cth = cos(theta)
              sqrts = sqrt(scale)
              phit = phi(-t)
              psit = psi(-t)
              erl = exp(rlog)
              y = sth/sqrts
              pdy = cth*sqrts
C Update rlog and theta:
              theta = theta + atan2(h* (scale/p*cth**2+
     &                qbig/scale*sth**2),chi(-t)+
     &                (scale/p-qbig/scale)*h*sth*cth)
              ynew = h*phit*pdy/p + y*psit
              pdynew = pdy*psit - y*h*qbig*phit
              rlog = rlog + half*log((pdynew**2)/scale+ (ynew**2)*scale)
C
C Compute the contribution from each of these approximating terms
C 1: The Constant Term
              chqtrt = chi(-qtrt)
              cntrb1 = half*h*erl* (y+ynew)/chqtrt
C 2: The Linear Term
              gqtrt = g(-qtrt)
              cntrb2 = (h**2)*qtr*erl* (ynew-y)*gqtrt
C 3: The Quadratic Term ???
              cntrb3 = ((h*half)**3)*erl* ((ynew+y)/chqtrt)*
     &                 (one-two*gqtrt)
C 4: The whole contribution:
              DO 70 i = 1,nvec
                  f0 = fmid(i)
                  f1 = (fend(i)-fo(i))/h
                  f2 = two* (fo(i)+fend(i)-two*fmid(i))/ (h**2)
                  ans(i) = w* (f0*cntrb1+f1*cntrb2+f2*cntrb3)

70            CONTINUE
          END IF

      END IF
C If we were integrating backwards then the next line will get the correct
C sign for the integral which we actually want.

      DO 80 i = 1,nvec
          ans(i) = ans(i)*sign(one,h)
80    CONTINUE

      RETURN

      END SUBROUTINE quad
C -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION g(t)
C     .. Parameters ..
      DOUBLE PRECISION fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10,
     &                 fac11,fac12,fac13,fac14,fac15
      PARAMETER (fac2=5.D-1,fac3=fac2/3.D0,fac4=fac3/4.D0,
     &          fac5=fac4/5.D0,fac6=fac5/6.D0,fac7=fac6/7.D0,
     &          fac8=fac7/8.D0,fac9=fac8/9.D0,fac10=fac9/10.D0,
     &          fac11=fac10/11.D0,fac12=fac11/12.D0,fac13=fac12/13.D0,
     &          fac14=fac13/14.D0,fac15=fac14/15.D0)
      DOUBLE PRECISION c0,c1,c2,c3,c4,c5,c6
      PARAMETER (c0=fac2-fac3,c1=fac4-fac5,c2=fac6-fac7,c3=fac8-fac9,
     &          c4=fac10-fac11,c5=fac12-fac13,c6=fac14-fac15)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION t
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION phi
c      EXTERNAL phi
C     ..
      g = c0 + t* (c1+t* (c2+t* (c3+t* (c4+t* (c5+t*c6)))))/phi(t)
      RETURN

      END FUNCTION g

C --------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION xi(t)
C     .. Parameters ..
      DOUBLE PRECISION fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10,
     &                 fac11,fac12,fac13,fac14
      PARAMETER (fac2=5.D-1,fac3=fac2/3.D0,fac4=fac3/4.D0,
     &          fac5=fac4/5.D0,fac6=fac5/6.D0,fac7=fac6/7.D0,
     &          fac8=fac7/8.D0,fac9=fac8/9.D0,fac10=fac9/10.D0,
     &          fac11=fac10/11.D0,fac12=fac11/12.D0,fac13=fac12/13.D0,
     &          fac14=fac13/14.D0)
      DOUBLE PRECISION c1,c2,c3,c4,c5,c6,c7
      PARAMETER (c1=fac2/5.D0,c2=-fac4/7.D0,c3=fac6/9.D0,c4=-fac8/11.D0,
     &          c5=fac10/13.D0,c6=-fac12/15.D0,c7=fac14/17.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION t
C     ..
      xi = c1 + t* (c2+t* (c3+t* (c4+t* (c5+t* (c6+t*c7)))))
      RETURN

      END FUNCTION xi
C ---------------------------------------------------------------------
C --------------------- Source from EXTRA -----------------------------
C ---------------------------------------------------------------------
      SUBROUTINE sl07fm(elam1,elam2,wk1,wk2,iwk1,iwk2,m1,m2,k1,k2,
     &                  eigfn1,eigfn2,n1,n2,i1,i2,fun,tol,ans,ifail)
C This routine provides an interface for the routine FCF for computing
C Franck-Condon factors. It is designed primarily for the multiple
C eigenvalue case where the eigenvalues of each problem have all been
C computed on a single mesh, but may also be used for the case where the
C eigenvalues have been computed singly. A single Franck-Condon factor,
C specified by the integers i1 and i2, is computed at each call.
C Failure flags:
C IFAIL = 1: Parameter out of range on entry.
C IFAIL = 3: An error in SUBROUTINE fcf. Cannot achieve required accuracy.
C IFAIL = 2: A parameter error. We do not have the eigenfunctions for the
C            integrals requested.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER i1,i2,ifail,iwk1,iwk2,m1,m2,n1,n2
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ans(0:1),eigfn1(0:1,1:3,1:m1),
     &                 eigfn2(0:1,1:3,1:m2),elam1(0:1,1:m1),
     &                 elam2(0:1,1:m2),wk1(0:iwk1,1:4),wk2(0:iwk2,1:4)
      INTEGER k1(1:m1),k2(1:m2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION eig1,eig2
      INTEGER i,ic1,ic2,iflag,j1,j2,nrec
      CHARACTER srname*6
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. External Subroutines ..
cc      EXTERNAL nfcf
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      iflag = 0
      srname = 'SL07FM'
      IF (abs(ifail).GT.1 .OR. tol.LE.0.D0) THEN
          WRITE (rec,FMT=9000)

9000      FORMAT ('** Parameter error: IFAIL or TOL out of range')

          IF (abs(ifail).GT.1) ifail = -1
          iflag = 1
          nrec = 1
          GO TO 50

      END IF
C Check to make sure we have available the eigenfunctions for the integral
C requested:
      DO 10 i = 1,m1
          IF (k1(i).EQ.i1) THEN
              j1 = i
              GO TO 20

          END IF

10    CONTINUE

      iflag = 2
      WRITE (rec,FMT=9010)

9010  FORMAT ('** Necessary eigenfunctions not both available')

      nrec = 1
      GO TO 50

20    DO 30 i = 1,m2
          IF (k2(i).EQ.i2) THEN
              j2 = i
              GO TO 40

          END IF

30    CONTINUE

      iflag = 2
      WRITE (rec,FMT=9020)

9020  FORMAT ('** Necessary eigenfunctions not both available')

      nrec = 1
      GO TO 50

40    ic1 = wk1(0,2)
      ic2 = wk2(0,2)

      iflag = 0
      eig1 = elam1(0,j1) + elam1(1,j1)
      eig2 = elam2(0,j2) + elam2(1,j2)
      CALL nfcf(eig1,wk1(0,1),wk1(1,2),wk1(1,3),wk1(1,4),eigfn1(0,1,j1),
     &          eigfn1(0,2,j1),eigfn1(0,3,j1),eig2,wk2(0,1),wk2(1,2),
     &          wk2(1,3),wk2(1,4),eigfn2(0,1,j2),eigfn2(0,2,j2),
     &          eigfn2(0,3,j2),n1,n2,ic1,ic2,fun,tol,ans(0),ans(1),
     &          iflag)
      IF (iflag.NE.0) THEN
          iflag = 3
          WRITE (rec,FMT=9030)

9030      FORMAT ('** Automatic step-size control has failed',/,
     &           '** Try a larger value of TOL')

          nrec = 2
          GO TO 50

      END IF

      ifail = 0
      RETURN

50    ifail = p01abf(ifail,iflag,srname,nrec,rec)
      END SUBROUTINE sl07fm
C ---------------------------------------------------------------------
      SUBROUTINE sl07f(elam1,elam2,wk1,wk2,iwk1,iwk2,eigfn1,eigfn2,n1,
     &                 n2,fun,tol,ans,ifail)
C This routine provides an interface for the routine FCF for computing
C Franck-Condon factors. It is designed  for the single
C eigenvalue case. A single Franck-Condon factor, is computed at each call.
C Failure flags:
C IFAIL = 1: Parameter out of range on entry.
C IFAIL = 2: An error in SUBROUTINE fcf. Cannot achieve required accuracy.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER ifail,iwk1,iwk2,n1,n2
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ans(0:1),eigfn1(0:1,1:3),eigfn2(0:1,1:3),
     &                 elam1(0:1),elam2(0:1),wk1(0:iwk1,1:4),
     &                 wk2(0:iwk2,1:4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION eig1,eig2
      INTEGER ic1,ic2,iflag,nrec
      CHARACTER srname*6
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
c      EXTERNAL p01abf
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. External Subroutines ..
cc      EXTERNAL nfcf
C     ..
C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      iflag = 0
      srname = ' SL07F'
      IF (abs(ifail).GT.1 .OR. tol.LE.0.D0) THEN
          WRITE (rec,FMT=9000)

9000      FORMAT ('** Parameter error: IFAIL or TOL out of range')

          IF (abs(ifail).GT.1) ifail = -1
          iflag = 1
          nrec = 1
          GO TO 10

      END IF

      ic1 = wk1(0,2)
      ic2 = wk2(0,2)

      iflag = 0
      eig1 = elam1(0) + elam1(1)
      eig2 = elam2(0) + elam2(1)
      CALL nfcf(eig1,wk1(0,1),wk1(1,2),wk1(1,3),wk1(1,4),eigfn1(0,1),
     &          eigfn1(0,2),eigfn1(0,3),eig2,wk2(0,1),wk2(1,2),wk2(1,3),
     &          wk2(1,4),eigfn2(0,1),eigfn2(0,2),eigfn2(0,3),n1,n2,ic1,
     &          ic2,fun,tol,ans(0),ans(1),iflag)
      IF (iflag.NE.0) THEN
          iflag = 2
          WRITE (rec,FMT=9010)

9010      FORMAT ('** Automatic step-size control has failed',/,
     &           '** Try a larger value of TOL')

          nrec = 2
          GO TO 10

      END IF

      ifail = 0
      RETURN

10    ifail = p01abf(ifail,iflag,srname,nrec,rec)
      END SUBROUTINE sl07f
C ---------------------------------------------------------------------

      SUBROUTINE nfcf(elam1,xmesh1,pp1,qp1,wp1,rlog1,theta1,scale1,
     &                elam2,xmesh2,pp2,qp2,wp2,rlog2,theta2,scale2,n1,
     &                n2,ic1,ic2,fun,tol,ans,errest,ifail)
C Declarations

C     .. Parameters ..
      DOUBLE PRECISION one,two,three,half,trd,qtr,safe,tosafe,bound,zero
      PARAMETER (one=1.D0,two=2.D0,three=3.D0,half=one/two,
     &          trd=one/three,qtr=half*half,safe=0.95D0,tosafe=0.8D0,
     &          bound=2.D-1,zero=0.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ans,elam1,elam2,errest,tol
      INTEGER ic1,ic2,ifail,n1,n2
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION pp1(1:n1),pp2(1:n2),qp1(1:n1),qp2(1:n2),
     &                 rlog1(0:1),rlog2(0:1),scale1(0:1),scale2(0:1),
     &                 theta1(0:1),theta2(0:1),wp1(1:n1),wp2(1:n2),
     &                 xmesh1(0:n1),xmesh2(0:n2)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION abserr,ans1,ans2,ans3,ans4,di,diff,dummy,err,h,
     &                 hmax,p1,p2,q1,q2,ratio,rl1,rl2,rln1,
     &                 rln2,sc1,sc2,scn1,scn2,th1,th2,
     &                 thn1,thn2,toloc,value,w1,w2,x1,x2,x3,xend,
     &                 xfin,xhigh,xlow,xmid,xo,xsplit,xstart,xstop
      INTEGER i,i1,i1high,i1low,i1m1,i2,i2high,i2low,i2m1,icount,idi,
     &        ii1,ii2,ipp,isplit,isplt1,isplt2,itest,iupdat,nstep
      LOGICAL frstme
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION rloc(0:1),scloc(0:1),thloc(0:1)
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf
cc      EXTERNAL x02ajf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL fulprf,fulstp,inrprd
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
C This subroutine computes Franck-Condon factors, given the appropriate
C eigenfunctions.

      errest = zero
      nstep = 0
C Very first bit: locate a suitable splitting point for the integration.
C Locate this point on the mesh XMESH1 (well it has to be somewhere!)

      toloc = tol*safe
      icount = 0

      i1low = ic1
      i1high = ic1
      i2low = ic2
      i2high = ic2
      DO 10 i = ic1,1,-1
          IF ((elam1*wp1(i)-qp1(i)).LT.zero) GO TO 20
          i1low = i - 1
10    CONTINUE
20    DO 30 i = ic1,n1
          IF ((elam1*wp1(i)-qp1(i)).LT.zero) GO TO 40
          i1high = i
30    CONTINUE
40    DO 50 i = ic2,1,-1
          IF ((elam2*wp2(i)-qp2(i)).LT.zero) GO TO 60
          i2low = i - 1
50    CONTINUE
60    DO 70 i = ic2,n2
          IF ((elam2*wp2(i)-qp2(i)).LT.zero) GO TO 80
          i2high = i
70    CONTINUE

80    IF ((xmesh1(i1high)-xmesh2(i2low))*
     &    (xmesh2(i2high)-xmesh1(i1low)).GE.zero) THEN
C There is an overlap of the classically allowed regions. We choose a
C matchpoint in this region.
          xhigh = min(xmesh2(i2high),xmesh1(i1high))
          xlow = max(xmesh2(i2low),xmesh1(i1low))
          xsplit = half* (xhigh+xlow)
          diff = abs(xsplit-xmesh1(i1high))
          isplt1 = i1high
          DO 90 i = i1low,i1high - 1
              IF (xmesh1(i).GE.xlow .AND. xmesh1(i).LE.xhigh) THEN
                  IF (abs(xmesh1(i)-xsplit).LE.diff) THEN
                      isplt1 = i
                      diff = abs(xmesh1(i)-xsplit)
                  END IF

              END IF

90        CONTINUE
          diff = abs(xsplit-xmesh2(i2high))
          isplt2 = i2high
          DO 100 i = i2low,i2high - 1
              IF (xmesh2(i).GE.xlow .AND. xmesh2(i).LE.xhigh) THEN
                  IF (abs(xmesh2(i)-xsplit).LE.diff) THEN
                      isplt2 = i
                      diff = abs(xmesh2(i)-xsplit)
                  END IF

              END IF

100       CONTINUE

      ELSE
C There is no overlap in the classically allowed regions. Use a
C heuristic.
          xsplit = half* (xmesh1(ic1)+xmesh2(ic2))
          diff = abs(xsplit-xmesh1(1))
          isplt1 = 1
          DO 110 i = 2,n1 - 1
              IF (abs(xmesh1(i)-xsplit).LE.diff) THEN
                  isplt1 = i
                  diff = abs(xmesh1(i)-xsplit)
              END IF

110       CONTINUE
          diff = abs(xsplit-xmesh2(i))
          isplt2 = 1
          DO 120 i = 2,n2 - 1
              IF (abs(xmesh2(i)-xsplit).LE.diff) THEN
                  isplt2 = i
                  diff = abs(xmesh2(i)-xsplit)
              END IF

120       CONTINUE

      END IF

      IF (abs(xmesh1(isplt1)-xsplit).LT.abs(xmesh2(isplt2)-xsplit)) THEN
          itest = 1
          isplit = isplt1

      ELSE IF (abs(xmesh1(isplt1)-xsplit).GT.
     &         abs(xmesh2(isplt2)-xsplit)) THEN
          itest = 2
          isplit = isplt2

      ELSE IF (xmesh1(isplt1).LE.xmesh2(isplt2)) THEN
          itest = 1
          isplit = isplt1

      ELSE
          itest = 2
          isplit = isplt2
      END IF
C What we just did is by no means foolproof. It is always possible that
C there is a part of the range of integration where shooting in one
C direction is unstable for one eigenfunction while shooting in a different
C direction is unstable for the other eigenfunction. This subroutine assumes
C that such an interval does not exist.

      rl1 = rlog1(0)
      th1 = theta1(0)
      sc1 = scale1(0)
      rl2 = rlog2(0)
      th2 = theta2(0)
      sc2 = scale2(0)
      idi = 1
      di = one
      i1 = 1
      i2 = 1
      i1m1 = i1 - idi
      i2m1 = i2 - idi

      IF (xmesh1(0).EQ.xmesh2(0)) THEN
          xstart = xmesh1(0)
          iupdat = 0

      ELSE IF (xmesh1(0).LT.xmesh2(0)) THEN
          xstart = xmesh2(0)
          iupdat = 1

      ELSE
          xstart = xmesh1(0)
          iupdat = 2
      END IF

      IF (iupdat.EQ.1) THEN

130       IF ((xmesh1(i1)-xstart)*di.LE.zero) THEN
              rloc(0) = rl1
              thloc(0) = th1
              scloc(0) = sc1
C WARNING: this changes i1m1
              CALL fulprf(i1m1,thloc,rloc,scloc,dummy,i1,one,elam1,pp1,
     &                    qp1,wp1,xmesh1,1,n1,2)
              rl1 = rloc(1)
              th1 = thloc(1)
              sc1 = scloc(1)
              i1m1 = i1
              i1 = i1 + idi
              GO TO 130

          END IF

C Integrate RL1, TH1 and SC1 from XMESH1(I1M1) to XSTART.

          ii1 = max(i1,i1m1)
          CALL fulstp(xmesh1(i1m1),xstart,di,sc1,th1,rl1,pp1(ii1),
     &                elam1*wp1(ii1)-qp1(ii1),wp1(ii1),dummy,1)

      ELSE IF (iupdat.EQ.2) THEN

140       IF ((xmesh2(i2)-xstart)*di.LE.zero) THEN
              rloc(0) = rl2
              thloc(0) = th2
              scloc(0) = sc2
C WARNING: this changes i2m1
              CALL fulprf(i2m1,thloc,rloc,scloc,dummy,i2,one,elam2,pp2,
     &                    qp2,wp2,xmesh2,1,n2,2)
              rl2 = rloc(1)
              th2 = thloc(1)
              sc2 = scloc(1)
              i2m1 = i2
              i2 = i2 + idi
              GO TO 140

          END IF

C Integrate RL2, TH2 and SC2 from XMESH2(I2M1) to XSTART.

          ii2 = max(i2,i2m1)
          CALL fulstp(xmesh2(i2m1),xstart,di,sc2,th2,rl2,pp2(ii2),
     &                elam2*wp2(ii2)-qp2(ii2),wp2(ii2),dummy,1)

      END IF

      frstme = .true.

      value = zero

      h = min(xmesh1(i1)-xstart,xmesh2(i2)-xstart)
C We have now located a suitable starting point and initial conditions.

C Integrate forward by a suitable step.
C Select correct coefficient function values:
150   ipp = max(i1,i1m1)
      p1 = pp1(ipp)
      q1 = qp1(ipp)
      w1 = wp1(ipp)
      ipp = max(i2,i2m1)
      p2 = pp2(ipp)
      q2 = qp2(ipp)
      w2 = wp2(ipp)
      IF ((xmesh1(i1)-xmesh2(i2))*di.LT.zero) THEN
          xstop = xmesh1(i1)
          i1m1 = i1
          i1 = i1 + idi

      ELSE
          xstop = xmesh2(i2)
          IF (xmesh1(i1).EQ.xmesh2(i2)) THEN
              i1m1 = i1
              i1 = i1 + idi
          END IF

          i2m1 = i2
          i2 = i2 + idi
      END IF

      hmax = min(abs(xmesh1(i1)-xmesh1(i1m1)),
     &       abs(xmesh2(i2)-xmesh2(i2m1)))
C Now perform integration from xtart to xstop:

      xo = xstart
      xfin = xstop
      IF (abs(xfin-xo).LT.x02ajf(one)) THEN
C This interval is trivial. Do the next one.
          xstart = xstop
          GO TO 150

      END IF

      xend = xstart + di*min(abs(xstop-xstart),abs(h),hmax)

160   rln1 = rl1
      thn1 = th1
      scn1 = sc1
      rln2 = rl2
      thn2 = th2
      scn2 = sc2

      xmid = half* (xo+xend)
      x1 = xo
      x2 = xmid

      CALL inrprd(x1,x2,p1,q1,w1,rl1,th1,sc1,elam1,p2,q2,w2,rl2,th2,sc2,
     &            elam2,one,ans2)

      ans1 = ans2
      ans2 = ans2*fun(half* (xo+xmid))
      x2 = xmid
      x3 = xend
      CALL inrprd(x2,x3,p1,q1,w1,rl1,th1,sc1,elam1,p2,q2,w2,rl2,th2,sc2,
     &            elam2,one,ans3)

      ans1 = (ans1+ans3)*fun(half* (xo+xend))
      ans4 = ans2 + ans3*fun(half* (xmid+xend))
C Estimate the error
      err = (ans4-ans1)/three
      errest = errest + err
      abserr = abs(err)
      IF (abserr.GT.toloc) THEN
          ratio = min(half, (toloc/abserr)**trd)
C It is conceivable that a failed step may occur when we are stepping exactly
C to the end xfin, in which case (xend-xo) neq h. In such a case we must
C set h = ratio*(xend-xo) for the repeat step, rather than h = ratio*h.
C In other cases xend-xo=h and there is nothing to worry about.
          xend = xo + ratio* (xend-xo)
          xend = di*min(di* (xend),di*xfin)
          icount = icount + 1
C At most three attempts to get the stepsize right:
          IF (icount.LT.4) THEN
              rl1 = rln1
              th1 = thn1
              sc1 = scn1
              rl2 = rln2
              th2 = thn2
              sc2 = scn2
              GO TO 160

          END IF

          ifail = 1
          GO TO 200

      ELSE
          value = value + (ans4+err)
          IF (xend*di.GE. (xfin-di*x02ajf(1.D0))*di) GO TO 170
          IF (abserr.LT.toloc*tosafe .AND. icount.EQ.0) THEN
C The error is very small and the step before last was not a failure so
C we may increase the stepsize.
              ratio = one/max(bound, (abserr/toloc)**qtr)
              h = di*min(ratio*abs(h),hmax)
          END IF
C We have completed the last step successfully, and provided the step before
C was not a failure we may have increased the step length.
          icount = 0
          xo = xend
          xend = xo + h
C Step exactly to the end of the interval if xend >= xfin - x02ajf. In
C this eventuality we do not set h = xend-xo since this might result in
C excessively small stepsizes for some time thereafter.
          IF (xend*di.GE. (xfin-di*x02ajf(1.D0))*di) THEN
              xend = xfin
          END IF

          GO TO 160

      END IF

C Have successfully completed the integration from xo to xfin.

170   nstep = nstep + 1
      icount = 0
C Integration from xstart to xstop has now been performed. What next?

      IF ((itest.EQ.1.AND.idi* (i1-isplit).LE.0) .OR.
     &    (itest.EQ.2.AND.idi* (i2-isplit).LE.0)) THEN
C Still another step to do (n.b. itest always has one of the values 1,2)
          xstart = xstop
          GO TO 150

      END IF

      IF (frstme) THEN
C We have only done one half of the integration. Let us do the other.

          rl1 = rlog1(1)
          th1 = theta1(1)
          sc1 = scale1(1)
          rl2 = rlog2(1)
          th2 = theta2(1)
          sc2 = scale2(1)
          idi = -1
          di = -one
          i1 = n1 - 1
          i2 = n2 - 1
          i1m1 = i1 - idi
          i2m1 = i2 - idi
          xstart = min(xmesh1(n1),xmesh2(n2))

          IF (xmesh1(n1).EQ.xmesh2(n2)) THEN
              xstart = xmesh1(n1)
              iupdat = 0

          ELSE IF (xmesh1(n1).GT.xmesh2(n2)) THEN
              xstart = xmesh2(n2)
              iupdat = 1

          ELSE
              xstart = xmesh1(n1)
              iupdat = 2
          END IF

          IF (iupdat.EQ.1) THEN

180           IF ((xmesh1(i1)-xstart)*di.LT.zero) THEN
                  rloc(0) = rl1
                  thloc(0) = th1
                  scloc(0) = sc1
C WARNING: this changes i1m1
                  CALL fulprf(i1m1,thloc,rloc,scloc,dummy,i1,one,elam1,
     &                        pp1,qp1,wp1,xmesh1,1,n1,2)
                  rl1 = rloc(1)
                  th1 = thloc(1)
                  sc1 = scloc(1)
                  i1m1 = i1
                  i1 = i1 + idi
                  GO TO 180

              END IF
C Integrate RL1, TH1 and SC1 from XMESH1(I1M1) to XSTART.
              ii1 = max(i1,i1m1)
              CALL fulstp(xmesh1(i1m1),xstart,di,sc1,th1,rl1,pp1(ii1),
     &                    elam1*wp1(ii1)-qp1(ii1),wp1(ii1),dummy,1)

          ELSE IF (iupdat.EQ.2) THEN

190           IF ((xmesh2(i2)-xstart)*di.LT.zero) THEN
                  rloc(0) = rl2
                  thloc(0) = th2
                  scloc(0) = sc2
C WARNING: this changes i2m1
                  CALL fulprf(i2m1,thloc,rloc,scloc,dummy,i2,one,elam2,
     &                        pp2,qp2,wp2,xmesh2,1,n2,2)
                  rl2 = rloc(1)
                  th2 = thloc(1)
                  sc2 = scloc(1)
                  i2m1 = i2
                  i2 = i2 + idi
                  GO TO 190

              END IF

C Integrate RL2, TH2 and SC2 from XMESH2(I2M1) to XSTART.
              ii2 = max(i2,i2m1)
              CALL fulstp(xmesh2(i2m1),xstart,di,sc2,th2,rl2,pp2(ii2),
     &                    elam2*wp2(ii2)-qp2(ii2),wp2(ii2),dummy,1)

          END IF

          frstme = .false.
          h = max(xmesh1(i1),xmesh2(i2)) - xstart
          GO TO 150
      END IF

C We have now completed the integration.

      ans = value

200   RETURN

      END SUBROUTINE nfcf
C --------------------------------------------------------------------

      SUBROUTINE sl06fm(elam,eigfns,wk,iwk,n,k,m,i1,i2,fun,tol,elem,
     &                  ifail)
C This routine computes the so-called matrix elements which are of interest
C in problems arising in quantum mechanics. It acts as an interface for the
C routine matel.
C     .. Scalar Arguments ..
      DOUBLE PRECISION tol
      INTEGER i1,i2,ifail,iwk,m,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION eigfns(0:1,1:3,1:m),elam(1:2,1:m),elem(0:1),
     &                 wk(0:iwk,1:4)
      INTEGER k(1:m)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION elam1,elam2
      INTEGER i,ic,iflag,j1,j2,nrec
      LOGICAL check1,check2
      CHARACTER srname*6
C     ..
C     .. External Subroutines ..
cc      EXTERNAL nmatel
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,nint
C     ..

C     .. Local Arrays ..
      CHARACTER rec(2)*80
C     ..
      srname = 'SL06FM'
      IF (abs(ifail).GT.1 .OR. tol.LE.0.D0) THEN
          IF (abs(ifail).GT.1) ifail = -1
          iflag = 1
          WRITE (rec,FMT=9000)

9000      FORMAT ('** Parameter error: IFAIL or TOL out of range')

          nrec = 1
          GO TO 30

      END IF

      ic = nint(wk(0,2))

C Check that the requested eigenfunctions are available:
      check1 = .false.
      check2 = .false.
      DO 10 i = 1,m
          IF (k(i).EQ.i1) THEN
              j1 = i
              check1 = .true.
          END IF

          IF (k(i).EQ.i2) THEN
              j2 = i
              check2 = .true.
          END IF

          IF (check1 .AND. check2) GO TO 20
10    CONTINUE

      iflag = 2
      WRITE (rec,FMT=9010)

9010  FORMAT ('** Necessary eigenfunctions not all available')

      nrec = 1
      GO TO 30
C End of check for eigenfunction availability

20    elam1 = elam(1,j1) + elam(2,j1)
      elam2 = elam(1,j2) + elam(2,j2)
      iflag = 0
      CALL nmatel(elam1,elam2,eigfns(0,1,j1),eigfns(0,1,j2),
     &            eigfns(0,2,j1),eigfns(0,2,j2),eigfns(0,3,j1),
     &            eigfns(0,3,j2),wk(1,2),wk(1,3),wk(1,4),wk(0,1),n,ic,
     &            tol,fun,elem,iflag)
      IF (iflag.NE.0) THEN
          iflag = 3
          WRITE (rec,FMT=9020)

9020      FORMAT ('** Autmatic step-size control has failed',/,
     &           '** Try using a larger value of TOL')

          nrec = 2
          GO TO 30

      END IF

      ifail = 0
      RETURN

30    ifail = p01abf(ifail,iflag,srname,nrec,rec)
      RETURN

      END SUBROUTINE sl06fm
C --------------------------------------------------------------------------
      SUBROUTINE nmatel(elam1,elam2,rlog1,rlog2,theta1,theta2,scale1,
     &                  scale2,pp,qp,wp,xmesh,n,ic,tol,fun,ans,ifail)
C This routine computes the value of the matrix element
C
C            Integral(a,b)w(x)Y(1)(x)fun(x)Y(2)(x)dx
C
C     .. Parameters ..
      DOUBLE PRECISION one,two,three,half,trd,qtr,safe,tosafe,bound,zero
      PARAMETER (one=1.D0,two=2.D0,three=3.D0,half=one/two,
     &          trd=one/three,qtr=half*half,safe=0.95D0,tosafe=0.8D0,
     &          bound=2.D-1,zero=0.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION elam1,elam2,tol
      INTEGER ic,ifail,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ans(0:1),pp(1:n),qp(1:n),rlog1(0:1),rlog2(0:1),
     &                 scale1(0:1),scale2(0:1),theta1(0:1),theta2(0:1),
     &                 wp(1:n),xmesh(0:n)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION fun
      EXTERNAL fun
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION abserr,ans1,ans2,ans3,ans4,di,epsilo,err,errest,
     &                 h,hmax,p,q,ratio,rl1,rl2,rln1,rln2,sc1,sc2,scn1,
     &                 scn2,th1,th2,thn1,thn2,toloc,value,w,xend,xfin,
     &                 xmid,xo
      INTEGER i,icount,idi,ifin,ii,init,ip1,itime,nstep
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf
cc      EXTERNAL x02ajf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL inrprd
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
      toloc = safe*tol
      value = zero
      errest = zero
      nstep = 0

      itime = 1
      init = 0
      ifin = ic - 1
      idi = 1
      di = one
      icount = 0
      h = xmesh(1) - xmesh(0)
      rl1 = rlog1(0)
      th1 = theta1(0)
      sc1 = scale1(0)
      rl2 = rlog2(0)
      th2 = theta2(0)
      sc2 = scale2(0)
10    DO 40 i = init,ifin,idi
          xo = xmesh(i)
          ip1 = i + idi
          hmax = xmesh(ip1) - xmesh(i)
          h = di*min(abs(h),abs(hmax))
          xfin = xmesh(ip1)
C          write(9,*) xo,xfin
          xend = xo + h
          ii = i
          IF (idi.GT.0) ii = ip1
          p = pp(ii)
          q = qp(ii)
          w = wp(ii)

20        rln1 = rl1
          thn1 = th1
          scn1 = sc1
          rln2 = rl2
          thn2 = th2
          scn2 = sc2

          xmid = half* (xo+xend)
          CALL inrprd(xo,xmid,p,q,w,rl1,th1,sc1,elam1,p,q,w,rl2,th2,sc2,
     &                elam2,one,ans2)
          ans1 = ans2
          ans2 = ans2*fun(half* (xo+xmid))

          CALL inrprd(xmid,xend,p,q,w,rl1,th1,sc1,elam1,p,q,w,rl2,th2,
     &                sc2,elam2,one,ans3)
          ans1 = (ans1+ans3)*fun(half* (xo+xend))
          ans4 = ans2 + ans3*fun(half* (xmid+xend))
C Estimate the error
          err = (ans4-ans1)/three
          errest = errest + err
          abserr = abs(err)

          IF (abserr.GT.toloc) THEN
              ratio = min(half, (toloc/abserr)**trd)
C It is conceivable that a failed step may occur when we are stepping exactly
C to the end xfin, in which case (xend-xo) neq h. In such a case we must
C set h = ratio*(xend-xo) for the repeat step, rather than h = ratio*h.
C In other cases xend-xo=h and there is nothing to worry about.
              h = ratio* (xend-xo)
              xend = xo + h
              icount = icount + 1
C At most three attempts to get the stepsize right:
              IF (icount.LT.4) THEN
                  rl1 = rln1
                  th1 = thn1
                  sc1 = scn1
                  rl2 = rln2
                  th2 = thn2
                  sc2 = scn2
                  GO TO 20

              END IF

              ifail = 1
              GO TO 50

          ELSE
              value = value + ans4 + err
              IF ((xend-xfin)*di.GE.0.0D0) GO TO 30
              IF (abserr.LT.toloc*tosafe .AND. icount.EQ.0) THEN
C The error is very small and the step before last was not a failure so
C we may increase the stepsize.
                  ratio = one/max(bound, (abserr/toloc)**qtr)
                  h = ratio*h
              END IF
C We have completed the last step successfully, and provided the step before
C was not a failure we may have increased the step length.
              icount = 0
              xo = xend
              xend = xo + h
C Step exactly to the end of the interval if xend >= xfin - x02ajf. In
C this eventuality we do not set h = xend-xo since this might result in
C excessively small stepsizes for some time thereafter.
              IF (xend*di.GE. (xfin*di-x02ajf(epsilo))) THEN
                  xend = xfin
              END IF

              GO TO 20

          END IF

C Have successfully completed the integration from xmesh(i) to xmesh(i+idi)

30        nstep = nstep + 1
          icount = 0
40    CONTINUE

C If itime = 1 then we have only done the shooting from xmesh(0) to
C xmesh(ic), and we must therefore do the shooting from xmesh(n) to
C xmesh(ic).

      IF (itime.EQ.1) THEN
          itime = 2
          idi = -1
          di = -one
          init = n
          ifin = ic + 1
          h = xmesh(n) - xmesh(n-1)
          rl1 = rlog1(1)
          th1 = theta1(1)
          sc1 = scale1(1)
          rl2 = rlog2(1)
          th2 = theta2(1)
          sc2 = scale2(1)
          GO TO 10

      END IF

C Successful completion:
      ifail = 0
      ans(0) = value
      ans(1) = errest
C      WRITE (6,FMT=*) 'MATEL error estimate:',errest
      RETURN
C Unsuccessful completion:
50    ans(0) = zero
      ans(1) = errest
      RETURN

      END SUBROUTINE nmatel
C --------------------------------------------------------------------------
      SUBROUTINE inrprd(xo,xend,p1,q1,w1,rlog1,theta1,scale1,elam1,p2,
     &                  q2,w2,rlog2,theta2,scale2,elam2,f,ans)
C This routine computes the inner product of two eigenfunctions over
C an interval (xo,xend), advancing the Prufer radius, angle and scalefactor
C for each eigenfunction as the integration proceeds.
C
C     .. Parameters ..
      DOUBLE PRECISION half,qtr,two
      PARAMETER (half=5.D-1,qtr=half*half,two=2.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ans,elam1,elam2,f,p1,p2,q1,q2,rlog1,rlog2,scale1,
     &                 scale2,theta1,theta2,w1,w2,xend,xo
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alfa1h,alfa2h,cth1,cth2,dalfa,dalfa2,diff1,diff2,
     &                 dtheta,erl1,erl2,h,h2,pdy1nd,pdy1o,pdy2nd,pdy2o,
     &                 r01,r11,salfa,signum,snew,sqrnrm,sqrt1,sqrt2,
     &                 sth1,sth2,stheta,t1,t1o,t2,t2o,tdiff,y1end,y1o,
     &                 y2end,y2o
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION phi,rescal,scl,x02ajf
cc      EXTERNAL phi,rescal,scl,x02ajf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL chidif,fulstp,phidif
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,exp,sign,sin,sqrt
C     ..
      h = xend - xo
      h2 = h**2
      signum = sign(1.D0,h)
      t1o = (elam1*w1-q1)/p1
      t2o = (elam2*w2-q2)/p2
      tdiff = t1o - t2o
      t1 = t1o*h2
      t2 = t2o*h2

      IF (abs(tdiff).GT.sqrt(x02ajf(1.D0))) THEN
C In this case there is a very simple formula for the integral involving
C only the values of the eigenfunctions and their derivatives at each end.

C STEP ONE: Compute the values of each eigenfunction and its derivative
C at xo and xend.
C At xo:
          erl1 = exp(rlog1)
          erl2 = exp(rlog2)
          sth1 = sin(theta1)
          sth2 = sin(theta2)
          cth1 = cos(theta1)
          cth2 = cos(theta2)
          sqrt1 = sqrt(scale1)
          sqrt2 = sqrt(scale2)

          y1o = erl1*sth1/sqrt1
          y2o = erl2*sth2/sqrt2
          pdy1o = erl1*cth1*sqrt1
          pdy2o = erl2*cth2*sqrt2

C At xend: first need to advance the eigenfunctions using the FULSTP routine.

          CALL fulstp(xo,xend,signum,scale1,theta1,rlog1,p1,elam1*w1-q1,
     &                w1,sqrnrm,1)
          CALL fulstp(xo,xend,signum,scale2,theta2,rlog2,p2,elam2*w2-q2,
     &                w2,sqrnrm,1)

          erl1 = exp(rlog1)
          erl2 = exp(rlog2)
          sth1 = sin(theta1)
          sth2 = sin(theta2)
          cth1 = cos(theta1)
          cth2 = cos(theta2)
          sqrt1 = sqrt(scale1)
          sqrt2 = sqrt(scale2)

          y1end = erl1*sth1/sqrt1
          y2end = erl2*sth2/sqrt2
          pdy1nd = erl1*cth1*sqrt1
          pdy2nd = erl2*cth2*sqrt2
C We can now write down the required formula for the integral:

          ans = - ((pdy1nd*y2end-pdy1o*y2o)/p1-
     &          (pdy2nd*y1end-pdy2o*y1o)/p2)*f*w1/tdiff

      ELSE
C         (if (ABS(tdiff).LT.0.0001) THEN
C This is the case where t1 and t2 are close, and we must use the routines
C designed earlier for this case. There are two cases to consider: the first
C occurs when both t1 and t2 lie in (-infty,0.1), the second when one of
C these quantities lies in (0.1,infty).

C The case where t1 and t2 are both less than 0.1 is what the
C routines phidif and chidif were designed to handle.

          IF (t1.LE.0.1D0 .AND. t2.LE.0.1D0) THEN
              CALL chidif(-t2,-t1,diff1)
              r11 = diff1*h
              CALL phidif(-t2,-t1,diff2)
              r01 = -diff2*h
              erl1 = exp(rlog1)
              erl2 = exp(rlog2)
              sth1 = sin(theta1)
              sth2 = sin(theta2)
              cth1 = cos(theta1)
              cth2 = cos(theta2)
              sqrt1 = sqrt(scale1)
              sqrt2 = sqrt(scale2)

              y1o = erl1*sth1/sqrt1
              y2o = erl2*sth2/sqrt2
              pdy1o = erl1*cth1*sqrt1
              pdy2o = erl2*cth2*sqrt2

C At xend: first need to advance the eigenfunctions using the FULSTP routine.

              CALL fulstp(xo,xend,signum,scale1,theta1,rlog1,p1,
     &                    elam1*w1-q1,w1,sqrnrm,1)
              CALL fulstp(xo,xend,signum,scale2,theta2,rlog2,p2,
     &                    elam2*w2-q2,w2,sqrnrm,1)

              erl1 = exp(rlog1)
              erl2 = exp(rlog2)
              sth1 = sin(theta1)
              sth2 = sin(theta2)
              cth1 = cos(theta1)
              cth2 = cos(theta2)
              sqrt1 = sqrt(scale1)
              sqrt2 = sqrt(scale2)

              y1end = erl1*sth1/sqrt1
              y2end = erl2*sth2/sqrt2
              pdy1nd = erl1*cth1*sqrt1
              pdy2nd = erl2*cth2*sqrt2

              ans = ((y1o*y2o+y1end*y2end)*r11+
     &              (y1o*y2end+y1end*y2o)*r01)*w1*f

          ELSE
C In this case we use a completely different formula, which requires the
C updating of the scalefactor. Because in this case at least one of
C t1, t2 is greater than 0.1 and the other is at most 0.0001*h**2 different,
C it follows that both are positive.
              snew = sqrt((elam1*w1-q1)*p1)
              theta1 = scl(theta1,snew/scale1)
              rlog1 = rescal(scale1,snew,theta1,rlog1)
              scale1 = snew
              alfa1h = sqrt(t1)

              snew = sqrt((elam2*w2-q2)*p2)
              diff1 = sqrt(t2o) - sqrt(t1o)
              theta2 = scl(theta2,snew/scale2)
              rlog2 = rescal(scale2,snew,theta2,rlog2)
              scale2 = snew
              alfa2h = sqrt(t2)

              dtheta = theta1 - theta2
              stheta = theta1 + theta2
              dalfa = alfa1h - alfa2h
              dalfa2 = dalfa**2
              salfa = alfa1h + alfa2h
              r11 = half*h* (cos(dtheta)*phi(-dalfa2)-
     &              sin(dtheta)*half*diff1* (phi(-qtr*dalfa2)**2)-
     &              cos(stheta)*sin(salfa)/salfa+
     &              signum*two*sin(stheta)* (sin(salfa*half)**2)/salfa)
              ans = w1*f*exp(rlog1+rlog2)*r11/sqrt(scale1*scale2)
C Rlog1, rlog2, scale1 and scale2 do not need updating, but theta1 and theta2
C require updating:

              theta1 = theta1 + h*sqrt(t1o)
              theta2 = theta2 + h*sqrt(t2o)
C
C END of case where ABS(tdiff).LT.sqrt(x02ajf(1.d0))
          END IF
C END of all cases
      END IF
C Compensate for the possible effect of backwards integration:

      ans = ans*sign(1.D0,xend-xo)
      RETURN

      END SUBROUTINE inrprd

C -------------------------------------------------------------------------

      SUBROUTINE phidrv(t,phi0,phi1,phi2,phi3,phi4,phi5)
C     .. Parameters ..
      DOUBLE PRECISION p15,p3,p6
      PARAMETER (p15=15.D0,p3=3.D0,p6=6.D0)
      DOUBLE PRECISION half,qtr,eighth
      PARAMETER (half=5.D-1,qtr=half*half,
C     &          sxtnth=qtr*qtr,
     &          eighth=half*qtr)
      DOUBLE PRECISION fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10,fac11,
     &                 fac12,fac13,fac14,fac15,fac16,fac17,fac18,fac19
      PARAMETER (fac3=1.D0/6.D0,fac4=fac3/4.D0,fac5=fac4/5.D0,
     &          fac6=fac5/6.D0,fac7=fac6/7.D0,fac8=fac7/8.D0,
     &          fac9=fac8/9.D0,fac10=fac9/1.D1,fac11=fac10/11.D0,
     &          fac12=fac11/12.D0,fac13=fac12/13.D0,fac14=fac13/14.D0,
     &          fac15=fac14/15.D0,fac16=fac15/16.D0,fac17=fac16/17.D0,
     &          fac18=fac17/18.D0,fac19=fac18/19.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION phi0,phi1,phi2,phi3,phi4,phi5,t
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION epos,psi0,sqrto,to,to2,to3
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION phi,psi,x02amf
cc      EXTERNAL phi,psi,x02amf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,cosh,log,min,sin,sinh,sqrt
C     ..
      epos = (log(x02amf(half))**2)*half
      to = t

      IF (abs(t).GT.0.1D0) THEN

          IF (t.GT.0.1D0) THEN
              to = min(t,epos)
              sqrto = sqrt(to)
              phi0 = sinh(sqrto)/sqrto
              psi0 = cosh(sqrto)

          ELSE IF (t.LT. (-0.1D0)) THEN
              sqrto = sqrt(-to)
              phi0 = sin(sqrto)/sqrto
              psi0 = cos(sqrto)

          ELSE
              phi0 = phi(to)
              psi0 = psi(to)
          END IF

          to2 = to**2
          to3 = to2*to
C          to4 = to3*to
C          to5 = to4*to

          phi1 = (psi0-phi0)*half/to
          phi2 = ((p3+to)*phi0-p3*psi0)*qtr/to2
          phi3 = ((p15+t)*psi0- (p6*t+p15)*phi0)*eighth/to3
C          phi4 = ((to2+p45*to+p105)*phi0- (p10*to+p105)*psi0)*sxtnth/to4
C          phi5 = ((to2+p105*to+p945)*psi0- (p15*to2+p420*to+p945)*phi0)
C     &           *trtsnd/to5
      ELSE
C          IF (ABS(t).LE.0.1) THEN use Taylor expansions
          phi0 = phi(t)

C Compute the higher derivatives

          phi1 = fac3 + t* (2.D0*fac5+t*
     &           (3.D0*fac7+t* (4.D0*fac9+t* (5.D0*fac11+
     &           t*6.D0*fac13))))
          phi2 = 2.D0*fac5 + t* (6.D0*fac7+
     &           t* (12.D0*fac9+t* (20.D0*fac11+t* (30.D0*fac13+
     &           t*42.D0*fac15))))
          phi3 = 6.D0*fac7 + t* (24.D0*fac9+
     &           t* (60.D0*fac11+t* (120.D0*fac13+t* (210.D0*fac15+
     &           t*336.D0*fac17))))
C          phi4 = 24.*fac9 + t* (120.*fac11+
C     &           t* (360.*fac13+t* (840.*fac15+t* (1680.*fac17+
C     &           t*3024.*fac19))))
C          phi5 = 120.*fac11 + t* (720.*fac13+
C     &           t* (2520.*fac15+t* (6720.*fac17+t* (15120.*fac19+
C     &           t*30240.*fac21))))

      END IF

      RETURN

      END SUBROUTINE phidrv

C ----------------------------------------------------------------------------

      SUBROUTINE chidrv(t,chi0,chi1,chi2,chi3,chi4,chi5)

C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION chi0,chi1,chi2,chi3,chi4,chi5,t
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION phi0,phi1,phi2,phi3,phi4,phi5,psi0,psi1,psi2,
     &                 psi3,rat0,rat1,rat2,rat3,rat4,rat5,rat6,sqrtt,to
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION chi,psi,x02amf
cc      EXTERNAL chi,psi,x02amf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL phidrv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,cosh,log,min,sqrt,tan,tanh
C     ..
      CALL phidrv(t,phi0,phi1,phi2,phi3,phi4,phi5)

      sqrtt = sqrt(abs(t))
      to = min(t, (log(x02amf(half))**2)*half)

      IF (t.GT.0.1D0) THEN
          chi0 = sqrtt/tanh(sqrtt)
          psi0 = cosh(sqrt(to))

      ELSE IF (t.LT. (-0.1D0)) THEN
          chi0 = sqrtt/tan(sqrtt)
          psi0 = cos(sqrtt)

      ELSE
          chi0 = chi(t)
          psi0 = psi(t)
      END IF

      psi1 = half*phi0
      psi2 = half*phi1
      psi3 = half*phi2
C      psi4 = half*phi3
C      psi5 = half*phi4

C We now simply type in the expressions which were obtained by MAPLE when
C working out the derivative. We must be careful to do this in such a way
C that overflow is avoided. This is achieved by defining the following ratios:

      rat0 = psi1/phi0
      rat1 = psi0/phi0
      rat2 = phi1/phi0
      rat3 = psi2/phi0
      rat4 = phi2/phi0
      rat5 = psi3/phi0
      rat6 = phi3/phi0
C      rat7 = psi4/phi0
C      rat8 = phi4/phi0
C      rat9 = psi5/phi0
C      rat10 = phi5/phi0

      chi1 = rat0 - rat1*rat2
      chi2 = rat3 - rat1*rat4 - 2.D0* (rat0-rat1*rat2)*rat2
      chi3 = rat5 - 3.D0*rat3*rat2 + 6.D0*rat0* (rat2**2) -
     &       3.D0*rat0*rat4 - 6.D0*rat1* (rat2**3) +
     &       6.D0*rat1*rat2*rat4 - rat1*rat6
C      chi4 = rat7 - 4.d0*rat6*rat2 + 12.d0*rat3* (rat2**2) -
C     &       6.d0*rat3*rat4 - 24.d0*rat0* (rat2**3) +
C     &       24.d0*rat0*rat2*rat4 - 4.d0*rat0*rat6 +
C     &       24.d0*rat1* (rat2**4) - 36.d0*rat1* (rat2**2)*rat4 +
C     &       6.d0*rat1* (rat4**2) + 8.d0*rat1*rat2*rat6 - rat1*rat8
C      chi5 = rat9 - 5.d0*rat8*rat2 + 20.d0* (rat2**2)*rat6 -
C     &       60.d0* (rat2**3)*rat3 - 10.d0*rat3*rat6 - 5.d0*rat0*rat8 -
C     &       rat1*rat10 + 60.d0*rat3*rat0*rat4 +
C     &       120.d0*rat0* (rat2**4) - 18.d1*rat0* (rat2**2)*rat4 +
C     &       30.d0*rat0* (rat4**2) + 40.d0*rat0*rat2*rat6 -
C     &       12.d1*rat1* (rat2**5) + 24.d1*rat1* (rat2**3)*rat4 -
C     &       90.d0*rat1*rat2* (rat4**2) - 60.d0*rat1* (rat2**2)*rat6 +
C     &       20.d0*rat1*rat4*rat6 + 10.d0*rat1*rat2*rat8

      RETURN

      END SUBROUTINE chidrv

C --------------------------------------------------------------------------

      SUBROUTINE chidif(t1,t2,diff)
C This routine computes the ratio
C
C            chi(t2)-chi(t1)
C            ---------------
C                t2 - t1
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION diff,t1,t2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION chi0,chi01,chi02,chi1,chi2,chi3,chi4,chi5,sqrt1,
     &                 sqrt2,tdiff,tmid
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION chi
cc      EXTERNAL chi
C     ..
C     .. External Subroutines ..
cc      EXTERNAL chidrv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt,tan,tanh
C     ..
      sqrt1 = sqrt(abs(t1))
      sqrt2 = sqrt(abs(t2))

      IF (t1.GT.0.1D0) THEN
          chi01 = sqrt1/tanh(sqrt1)

      ELSE IF (t1.LT. (-0.1D0)) THEN
          chi01 = sqrt1/tan(sqrt1)

      ELSE
          chi01 = chi(t1)
      END IF

      IF (t2.GT.0.1D0) THEN
          chi02 = sqrt1/tanh(sqrt2)

      ELSE IF (t2.LT. (-0.1D0)) THEN
          chi02 = sqrt2/tan(sqrt2)

      ELSE
          chi02 = chi(t2)
      END IF

      IF (abs(t1-t2).GT.0.01D0) THEN
          diff = (chi01-chi02)/ (t1-t2)

      ELSE
C          Expand in a Taylor series about (t1+t2)/2

          tmid = 0.5D0* (t1+t2)
          CALL chidrv(tmid,chi0,chi1,chi2,chi3,chi4,chi5)
          tdiff = t2 - t1
          diff = chi1 + chi3* (tdiff**2)/24.D0
C     &         +   chi5* (tdiff**4)/ (1920.d0)

      END IF

      RETURN

      END SUBROUTINE chidif

C ---------------------------------------------------------------------

      SUBROUTINE phidif(t1,t2,diff)
C     .. Scalar Arguments ..
      DOUBLE PRECISION diff,t1,t2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION drv1,drv3,phi0,phi01,phi02,phi1,phi2,phi3,phi4,
     &                 phi5,rat1,rat2,sqrto,tdiff,tmid,to1,to2
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION phi,x02amf
cc      EXTERNAL phi,x02amf
C     ..
C     .. External Subroutines ..
cc      EXTERNAL phidrv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,log,min,sin,sinh,sqrt
C     ..
      to1 = min(t1, (log(x02amf(t1))**2)*5.D-1)
      to2 = min(t2, (log(x02amf(t2))**2)*5.D-1)

      IF (abs(t1-t2).GT.0.01D0) THEN

          IF (t1.GT.0.1D0) THEN
              sqrto = sqrt(to1)
              phi01 = sinh(sqrto)/sqrto

          ELSE IF (t1.LT. (-0.1D0)) THEN
              sqrto = sqrt(-to1)
              phi01 = sin(sqrto)/sqrto

          ELSE
              phi01 = phi(to1)
          END IF

          IF (t2.GT.0.1D0) THEN
              sqrto = sqrt(to2)
              phi02 = sinh(sqrto)/sqrto

          ELSE IF (t2.LT. (-0.1D0)) THEN
              sqrto = sqrt(-to2)
              phi02 = sin(sqrto)/sqrto

          ELSE
              phi02 = phi(to2)
          END IF

          diff = (1.D0/phi01-1.D0/phi02)/ (t1-t2)

      ELSE
C          ABS(t1-t2) is very small
          tmid = 0.5D0* (t1+t2)
          CALL phidrv(tmid,phi0,phi1,phi2,phi3,phi4,phi5)
          tdiff = t2 - t1
          rat1 = phi1/phi0
          rat2 = phi2/phi0
C          rat3 = phi3/phi0
C          rat4 = phi4/phi0
C          rat5 = phi5/phi0

          drv1 = -rat1
          drv3 = -6.D0* ((rat1**2)-rat2)*rat1
C          drv5 = -12.d1* (rat1**5) + 24.d1* (rat1**3)*rat2 -
C     &           9.d1*rat1* (rat2**2) - 6.d1* (rat1**2)*rat3 +
C     &           2.d1*rat2*rat3 + 1.d1*rat1*rat4 - rat5

          diff = (drv1+ (tdiff**2)*drv3/24.D0)/phi0
C    &           + (tdiff**4)*drv5/192.d1)/phi0

      END IF

      RETURN

      END SUBROUTINE phidif
C ---------------------------------------------------------------------
C --------------------- Source from USEFUN ----------------------------
C ---------------------------------------------------------------------

C -------------------- Useful Functions -----------------------------

      DOUBLE PRECISION FUNCTION spexp(x)
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION eneg
      LOGICAL first
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02amf
cc      EXTERNAL x02amf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,log
C     ..
C     .. Save statement ..
      SAVE first,eneg
C     ..
C     .. Data statements ..
      DATA first/.true./
C     ..
      IF (first) THEN
          eneg = log(x02amf(x))
          first = .false.
      END IF

      spexp = x02amf(x)
      IF (x.GT.eneg) spexp = exp(x)

      END FUNCTION spexp

C -------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION sptanh(x)

C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION loc
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION spexp
cc      EXTERNAL spexp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sign
C     ..
      loc = spexp(-2.0D0*abs(x))
      sptanh = sign(1.0D0,x)* (1.0D0-loc)/ (1.0D0+loc)
      END FUNCTION sptanh

C ---------------- Other Useful Functions --------------------------

      DOUBLE PRECISION FUNCTION phi(v)
C     .. Parameters ..
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5
      PARAMETER (a0=1.D0,a1=a0/6.D0,a2=a1/20.D0,a3=a2/42.D0,a4=a3/72.D0,
     &          a5=a4/110.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION v
C     ..
      phi = a0 + v* (a1+v* (a2+v* (a3+v* (a4+v*a5))))

      END FUNCTION phi

      DOUBLE PRECISION FUNCTION psi(v)
C     .. Parameters ..
      DOUBLE PRECISION b0,b1,b2,b3,b4,b5
      PARAMETER (b0=1.D0,b1=b0/2.D0,b2=b1/12.D0,b3=b2/30.D0,b4=b3/56.D0,
     &          b5=b4/90.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION v
C     ..
      psi = b0 + v* (b1+v* (b2+v* (b3+v* (b4+v*b5))))

      END FUNCTION psi

      DOUBLE PRECISION FUNCTION chi(v)
C     .. Scalar Arguments ..
      DOUBLE PRECISION v
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION phi,psi
cc      EXTERNAL phi,psi
C     ..
      chi = psi(v)/phi(v)

      END FUNCTION chi

      DOUBLE PRECISION FUNCTION scl(th,u)
C     .. Parameters ..
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION th,u
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC atan2,cos,sin
C     ..

      scl = th + atan2((u-one)*sin(th)*cos(th),one+ (u-one)*sin(th)**2)

      END FUNCTION scl

C ---------------------------------------------------------------------
C --------------------- Source from USEFUN ----------------------------
C ---------------------------------------------------------------------

      SUBROUTINE c05azf(x,y,fx,tolx,ir,c,ind,ifail)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-496 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     .. Parameters ..
      CHARACTER*6 srname
      PARAMETER (srname='C05AZF')
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION fx,tolx,x,y
      INTEGER ifail,ind,ir
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION c(17)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ab,diff,diff1,diff2,rel,rmax,tol,tol1
      INTEGER i
      LOGICAL t
C     ..
C     .. Local Arrays ..
      CHARACTER p01rec(1)*1
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02ajf,x02akf
cc      INTEGER p01abf
cc      EXTERNAL x02ajf,x02akf,p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,int,max,sign,sqrt
C     ..
      i = 0
      IF ((ind.LE.0.OR.ind.GT.4) .AND. ind.NE.-1) THEN
C     USER NOT CHECKED IND OR CHANGED IT
          i = 2
          ind = 0

      ELSE IF (tolx.GT.0.0D0 .AND. (ir.EQ.0.OR.ir.EQ.1.OR.ir.EQ.2)) THEN
          rel = 1.0D0
          ab = 1.0D0
          IF (ir.EQ.1) rel = 0.0D0
          IF (ir.EQ.2) ab = 0.0D0
          IF (ind.EQ.-1) THEN
              c(3) = x

          ELSE
              GO TO (10,30,40,20) ind

10            c(3) = x
              ind = 2
              RETURN

20            IF (fx.EQ.0.0D0) THEN
                  GO TO 80

              ELSE
                  c(4) = fx
                  rmax = abs(fx)
                  IF (c(13)*rmax.LE.c(15)) THEN
                      c(16) = 0.0D0

                  ELSE
                      IF (c(16).EQ.1.0D0) c(16) = -1.0D0
                      IF (c(16).EQ.0.0D0) c(16) = 1.0D0
                  END IF

                  IF (c(2).GE.0.0D0) THEN
                      t = c(4) .GE. 0.0D0

                  ELSE
                      t = c(4) .LE. 0.0D0
                  END IF

                  IF (t) THEN
                      GO TO 50

                  ELSE
                      i = int(c(17)+0.1D0)
                      i = i + 1
                      IF (c(11).EQ.c(12)) i = 0
                      c(17) = dble(i)
                      GO TO 60

                  END IF

              END IF

          END IF

30        IF (fx.NE.0.0D0) THEN
              c(4) = fx
              c(15) = abs(fx)
              c(16) = 0.0D0
              x = y
              y = c(3)
              c(2) = c(4)
              c(5) = x
              IF (ind.EQ.-1) THEN
                  fx = c(1)
                  ind = 3

              ELSE
                  ind = 3
                  RETURN

              END IF

          ELSE
              GO TO 80

          END IF

40        IF (fx.EQ.0.0D0) THEN
              GO TO 80

          ELSE IF (sign(1.0D0,fx).NE.sign(1.0D0,c(2))) THEN
              c(6) = fx
              c(13) = sqrt(x02ajf(0.0D0))
              c(15) = max(c(15),abs(fx))
              c(14) = x02akf(0.0D0)
              c(16) = 0.0D0

          ELSE
              ind = 0
              i = 1
              GO TO 90

          END IF

50        c(1) = c(5)
          c(2) = c(6)
          c(17) = 0.0D0
60        IF (abs(c(2)).LT.abs(c(4))) THEN
              IF (c(1).NE.c(5)) THEN
                  c(7) = c(5)
                  c(8) = c(6)
              END IF

              c(5) = c(3)
              c(6) = c(4)
              x = c(1)
              c(3) = x
              c(4) = c(2)
              c(1) = c(5)
              c(2) = c(6)
          END IF

          tol = 0.5D0*tolx*max(ab,rel*abs(c(3)))
          tol1 = 2.0D0*x02ajf(0.0D0)*max(ab,rel*abs(c(3)))
          diff2 = 0.5D0* (c(1)-c(3))
          c(12) = diff2
          diff2 = diff2 + c(3)
          IF (c(12).NE.0.0D0) THEN
              IF (abs(c(12)).LE.tol) THEN
                  IF (c(16).GE.0.0D0) THEN
                      y = c(1)
                      i = 0

                  ELSE
                      i = 4
                  END IF

                  ind = 0
                  GO TO 90

              ELSE IF (abs(c(12)).GT.tol1) THEN
                  IF (c(17).LT.2.5D0) THEN
                      tol = tol*sign(1.0D0,c(12))
                      diff1 = (c(3)-c(5))*c(4)
                      IF (c(17).LE.1.5D0) THEN
                          diff = c(6) - c(4)

                      ELSE IF (c(7).NE.c(3) .AND. c(7).NE.c(5)) THEN
                          c(9) = (c(8)-c(4))/ (c(7)-c(3))
                          c(10) = (c(8)-c(6))/ (c(7)-c(5))
                          diff1 = c(10)*diff1
                          diff = c(9)*c(6) - c(10)*c(4)

                      ELSE
                          GO TO 70

                      END IF

                      IF (diff1.LT.0.0D0) THEN
                          diff1 = -diff1
                          diff = -diff
                      END IF

                      IF (abs(diff1).LE.c(14) .OR.
     &                    diff1.LE.diff*tol) THEN
                          c(11) = tol

                      ELSE IF (diff1.GE.c(12)*diff) THEN
                          c(11) = c(12)

                      ELSE
                          c(11) = diff1/diff
                      END IF

                  ELSE
                      c(11) = c(12)
                  END IF

                  c(7) = c(5)
                  c(8) = c(6)
                  c(5) = c(3)
                  c(6) = c(4)
                  c(3) = c(3) + c(11)
                  x = c(3)
                  y = c(1)
                  ind = 4
                  RETURN

              END IF

          END IF

70        ind = 0
          i = 5
          GO TO 90

80        y = x
          ind = 0
          i = 0

      ELSE
          i = 3
          ind = 0
      END IF

90    ifail = p01abf(ifail,i,srname,0,p01rec)
      END SUBROUTINE c05azf

C --------------------------------------------------------------------
      SUBROUTINE m01caf(rv,m1,m2,order,ifail)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01CAF SORTS A VECTOR OF REAL NUMBERS INTO ASCENDING
C     OR DESCENDING ORDER.
C
C     M01CAF IS BASED ON SINGLETON'S IMPLEMENTATION OF THE
C     'MEDIAN-OF-THREE' QUICKSORT ALGORITHM, BUT WITH TWO
C     ADDITIONAL MODIFICATIONS. FIRST, SMALL SUBFILES ARE
C     SORTED BY AN INSERTION SORT ON A SEPARATE FINAL PASS.
C     SECOND, IF A SUBFILE IS PARTITIONED INTO TWO VERY
C     UNBALANCED SUBFILES, THE LARGER OF THEM IS FLAGGED FOR
C     SPECIAL TREATMENT: BEFORE IT IS PARTITIONED, ITS END-
C     POINTS ARE SWAPPED WITH TWO RANDOM POINTS WITHIN IT;
C     THIS MAKES THE WORST CASE BEHAVIOUR EXTREMELY UNLIKELY.
C
C     THE MAXIMUM LENGTH OF A SMALL SUBFILE IS DEFINED BY THE
C     VARIABLE MINQIK, SET TO 15.
C
C     THE ROUTINE ASSUMES THAT THE NUMBER OF ELEMENTS TO BE
C     SORTED DOES NOT EXCEED MINQIK*2**MAXSTK.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      INTEGER maxstk
      PARAMETER (maxstk=40)
      CHARACTER*6 srname
      PARAMETER (srname='M01CAF')
      INTEGER minqik
      PARAMETER (minqik=15)
C     ..
C     .. Scalar Arguments ..
      INTEGER ifail,m1,m2
      CHARACTER order*1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION rv(m2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,rand,x
      INTEGER i,i1,i2,ierr,ir1,ir2,ir3,istk,j,j1,j2,k,leng,nrec
C     ..
C     .. Local Arrays ..
      INTEGER ihigh(maxstk),ilow(maxstk)
      CHARACTER p01rec(2)*80
C     ..
C     .. External Functions ..
cc      INTEGER p01abf
cc      EXTERNAL p01abf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,mod
C     ..
C     .. Save statement ..
      SAVE ir1,ir2,ir3
C     ..
C     .. Data statements ..
      DATA ir1,ir2,ir3/15223,17795,28707/
C     ..
C
C       CHECK THE PARAMETERS AND DECIDE IF QUICKSORT IS NEEDED.
C
      IF (m2.LT.1 .OR. m1.LT.1 .OR. m1.GT.m2) THEN
          ierr = 1
          WRITE (p01rec,FMT=9000) m1,m2
          nrec = 2

      ELSE IF (order.NE.'A' .AND. order.NE.'a' .AND. order.NE.'D' .AND.
     &         order.NE.'d') THEN
          ierr = 2
          WRITE (p01rec,FMT=9010) order
          nrec = 1

      ELSE IF (m1.EQ.m2) THEN
          ierr = 0

      ELSE
          ierr = 0
          leng = m2 - m1 + 1
          IF (leng.GT.minqik) THEN
C
C           INITIALISE AND START QUICKSORT ON THE WHOLE VECTOR.
C
              istk = 0
              i = m1
              j = m2
10            CONTINUE
C
C           IF THE PREVIOUS PASS WAS BAD, CHANGE THE END VALUES AT
C           RANDOM.
C
              IF (i.LT.0) THEN
                  i = -i
                  ir1 = 171*mod(ir1,177) - 2* (ir1/177)
                  ir2 = 172*mod(ir2,176) - 35* (ir2/176)
                  ir3 = 170*mod(ir3,178) - 63* (ir3/178)
                  IF (ir1.LT.0) ir1 = ir1 + 30269
                  IF (ir2.LT.0) ir2 = ir2 + 30307
                  IF (ir3.LT.0) ir3 = ir3 + 30323
                  rand = mod(dble(ir1)/30269.0D0+dble(ir2)/30307.0D0+
     &                   dble(ir3)/30323.0D0,1.0D0)
                  k = i + rand* (j-i)
                  x = rv(i)
                  rv(i) = rv(k)
                  rv(k) = x
                  k = i + j - k
                  x = rv(k)
                  rv(k) = rv(j)
                  rv(j) = x
              END IF
C
C           CALCULATE A MEDIAN BY SINGLETONS METHOD.
C
              k = (i+j)/2
              IF (rv(i).GT.rv(j)) THEN
                  x = rv(i)
                  rv(i) = rv(j)
                  rv(j) = x
              END IF

              a = rv(k)
              IF (a.LT.rv(i)) THEN
                  rv(k) = rv(i)
                  rv(i) = a
                  a = rv(k)

              ELSE IF (a.GT.rv(j)) THEN
                  rv(k) = rv(j)
                  rv(j) = a
                  a = rv(k)
              END IF
C
C           SPLIT THE VECTOR INTO TWO ASCENDING PARTS.  THIS IS WHERE
C           THE TIME IS SPENT.
C
              i1 = i
              j1 = j
20            CONTINUE
              i1 = i1 + 1
              IF (rv(i1).LT.a) THEN
                  GO TO 20

              ELSE
30                CONTINUE
                  j1 = j1 - 1
                  IF (rv(j1).GT.a) GO TO 30
                  IF (i1.LT.j1) THEN
                      x = rv(i1)
                      rv(i1) = rv(j1)
                      rv(j1) = x
                      GO TO 20

                  END IF

              END IF
C
C           STACK ONE SUBFILE, IF APPROPRIATE, AND CARRY ON.
C
              i2 = i1 - i
              j2 = j - j1
              IF (j2.LE.i2) THEN
                  IF (i2.GT.minqik) THEN
C
C                 TEST FOR VERY UNBALANCED SUBFILES
C                 ( THE DETAILS OF THE TEST ARE FAIRLY ARBITRARY.)
C
                      IF (5* (j2+5).LT.i2) i = -i
                      IF (j2.LE.minqik) THEN
                          j = i1 - 1

                      ELSE
                          istk = istk + 1
                          ilow(istk) = i
                          ihigh(istk) = i1 - 1
                          i = j1 + 1
                      END IF

                      GO TO 10

                  ELSE IF (istk.GT.0) THEN
                      i = ilow(istk)
                      j = ihigh(istk)
                      istk = istk - 1
                      GO TO 10

                  END IF
C
C              DEAL WITH THE CASE WHEN THE SECOND PART IS LARGER.
C
              ELSE IF (j2.GT.minqik) THEN
C
C                 TEST FOR VERY UNBALANCED SUBFILES
C                 ( THE DETAILS OF THE TEST ARE FAIRLY ARBITRARY.)
C
                  IF (5* (i2+5).LT.j2) j1 = - (j1+2)
                  IF (i2.LE.minqik) THEN
                      i = j1 + 1

                  ELSE
                      istk = istk + 1
                      ilow(istk) = j1 + 1
                      ihigh(istk) = j
                      j = i1 - 1
                  END IF

                  GO TO 10

              ELSE IF (istk.GT.0) THEN
                  i = ilow(istk)
                  j = ihigh(istk)
                  istk = istk - 1
                  GO TO 10

              END IF

          END IF
C
C           TIDY UP AND DO AN ASCENDING INSERTION SORT.
C
          DO 50 i = m1 + 1,m2
              a = rv(i)
              j = i - 1
              IF (a.LT.rv(j)) THEN
40                CONTINUE
                  rv(j+1) = rv(j)
                  j = j - 1
                  IF (j.GE.m1) THEN
                      IF (a.LT.rv(j)) GO TO 40
                  END IF

                  rv(j+1) = a
              END IF

50        CONTINUE
C
C           REVERSE THE ORDER IF NECESSARY AND RETURN.
C
          IF ((order.EQ.'D') .OR. (order.EQ.'d')) THEN
              DO 60 i = m1, (m1+m2-1)/2
                  i1 = m1 + m2 - i
                  x = rv(i)
                  rv(i) = rv(i1)
                  rv(i1) = x
60            CONTINUE
          END IF
C
      END IF

      IF (ierr.NE.0) THEN
          ifail = p01abf(ifail,ierr,srname,nrec,p01rec)

      ELSE
          ifail = 0
      END IF
C
9000  FORMAT (' ** On entry, one or more of the following parameter va',
     &       'lues is illegal',/,'    M1 =',I16,'  M2 =',I16)
9010  FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)

      END SUBROUTINE m01caf

C --------------------------------------------------------------------

      INTEGER FUNCTION p01abf(ifail,ierror,srname,nrec,rec)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER ierror,ifail,nrec
      CHARACTER srname* (*)
C     ..
C     .. Array Arguments ..
      CHARACTER rec(*)* (*)
C     ..
C     .. Local Scalars ..
      INTEGER i,nerr
      CHARACTER mess*72
C     ..
C     .. External Subroutines ..
cc      EXTERNAL p01abz,x04aaf,x04baf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,mod
C     ..
      IF (ierror.NE.0) THEN
C        Abnormal exit from calling routine
          IF (ifail.EQ.-1 .OR. ifail.EQ.0 .OR. ifail.EQ.-13 .OR.
     &        (ifail.GT.0.AND.mod(ifail/10,10).NE.0)) THEN
C           Noisy exit
              CALL x04aaf(0,nerr)
              DO 10 i = 1,nrec
                  CALL x04baf(nerr,rec(i))
10            CONTINUE
              IF (ifail.NE.-13) THEN
                  WRITE (mess,FMT=9000) srname,ierror
                  CALL x04baf(nerr,mess)
                  IF (abs(mod(ifail,10)).NE.1) THEN
C                 Hard failure
                      CALL x04baf(nerr,
     &                     ' ** NAG hard failure - execution terminated'
     &                            )
                      CALL p01abz

                  ELSE
C                 Soft failure
                      CALL x04baf(nerr,
     &                         ' ** NAG soft failure - control returned'
     &                            )
                  END IF

              END IF

          END IF

      END IF

      p01abf = ierror
C
9000  FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     &       ' =',I6)

      END FUNCTION p01abf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x01aaf(x)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x01aaf = 3.14159265358979323846264338328D0
      END FUNCTION x01aaf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02agf(x)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE SMALLEST POSITIVE FLOATING-POINT NUMBER  R
C     EXACTLY REPRESENTABLE ON THE COMPUTER SUCH THAT  -R, 1.0/R,
C     AND -1.0/R CAN ALL BE COMPUTED WITHOUT OVERFLOW OR UNDERFLOW.
C     ON MANY MACHINES THE CORRECT VALUE CAN BE DERIVED FROM THOSE
C     OF X02AAF, X02ABF AND X02ACF AS FOLLOWS
C
C     IF (X02ABF(X)*X02ACF(X).GE.1.0) X02AGF = X02ABF(X)
C     IF (X02ABF(X)*X02ACF(X).LT.1.0)
C    *                            X02AGF = (1.0+X02AAF(X))/X02ACF(X)
C
C     THE CORRECT VALUE SHOULD BE DEFINED AS A CONSTANT,
C     POSSIBLY IN SOME BINARY, OCTAL OR HEXADECIMAL REPRESENTATION,
C     AND INSERTED INTO THE ASSIGNMENT STATEMENT BELOW.
C
C     X IS A DUMMY ARGUMENT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. External Functions ..
cc      DOUBLE PRECISION x02amf
cc      EXTERNAL x02amf
C     ..
      x02agf = x02amf(x)
      END FUNCTION x02agf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02ahf(x)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     * MAXIMUM ARGUMENT FOR SIN AND COS *
C     RETURNS THE LARGEST POSITIVE REAL NUMBER MAXSC SUCH THAT
C     SIN(MAXSC) AND COS(MAXSC) CAN BE SUCCESSFULLY COMPUTED
C     BY THE COMPILER SUPPLIED SIN AND COS ROUTINES.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x02ahf = 2.147483648D9
      END FUNCTION x02ahf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02ajf(x)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x02ajf = 1.110223024625157D-16
      END FUNCTION x02ajf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02akf(x)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x02akf = 1d-300
c      x02akf = 2.9387358770557188D-39
      END FUNCTION x02akf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02alf(x)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1 - B**(-P)) * B**EMAX  (THE LARGEST POSITIVE MODEL
C     NUMBER)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x02alf = 1d300
c      x02alf = 1.7014118346046923D38
      END FUNCTION x02alf

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION x02amf(x)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
      x02amf = 1d-300
c      x02amf = 5.8774717541114401D-39
      END FUNCTION x02amf

C --------------------------------------------------------------------

      SUBROUTINE x04aaf(i,nerr)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER i,nerr
C     ..
C     .. Local Scalars ..
      INTEGER nerr1
C     ..
C     .. Save statement ..
      SAVE nerr1
C     ..
C     .. Data statements ..
      DATA nerr1/6/
C     ..
      IF (i.EQ.0) nerr = nerr1
      IF (i.EQ.1) nerr1 = nerr
      END SUBROUTINE x04aaf

C --------------------------------------------------------------------

      SUBROUTINE x04baf(nout,rec)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER nout
      CHARACTER rec* (*)
C     ..
C     .. Local Scalars ..
      INTEGER i
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC len
C     ..
      IF (nout.GE.0) THEN
C        Remove trailing blanks
          DO 10 i = len(rec),2,-1
              IF (rec(i:i).NE.' ') GO TO 20
10        CONTINUE
C        Write record to external file
20        WRITE (nout,FMT=9000) rec(1:i)
      END IF
C
9000  FORMAT (A)

      END SUBROUTINE x04baf

C --------------------------------------------------------------------

      SUBROUTINE p01abz
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
      STOP

      END SUBROUTINE p01abz
C     MARK 11.5

      end module MARCOMOD

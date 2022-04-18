C =====================================================================
C ======================= SUBROUTINES FOR EIGENVALUES =================
C =====================================================================

C First revision: 27th April 1995

      SUBROUTINE SLEUTH(A,B,ELAM,ERR,K,SYM,SL4COF,SL4BCS,TOL,NMESH,
     &                  IMATCH,NXTRAP,WORK,IWORK,WORKS,IFAIL)
C     .. Parameters ..
      INTEGER NCOARS,IELAM3
      PARAMETER (NCOARS=40,IELAM3=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,ELAM,ERR,TOL
      INTEGER IFAIL,IMATCH,IWORK,K,NMESH,NXTRAP
      LOGICAL SYM
C     ..
C     .. Array Arguments ..
C REMARK: WORKS must be of length at least 27*NCOARS+56
      DOUBLE PRECISION WORK(0:IWORK,1:6),WORKS(27*NCOARS+56)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS,SL4COF
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALOC,BLOC,CORR1,CORR2,EL,ELAM0,ELAM1,ELAM2,EPS,
     &                 EU,MACHPT,MCHEPS,OSC,OSCL,TEORAT,TOL1,TOLI,XSHIFT
      INTEGER COUNT,I,IFO,II,IM1ST,IMESH,IREFIN,IY1,IY2,J,JMAX,KNTMAX,
     &        KR,KXTRAP,MAXTRP,MULTI,NM1ST,NMAX,NREC
      LOGICAL RXTRAP
      CHARACTER*6 SRNAME
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ELAM3(IELAM3,IELAM3)
      CHARACTER*80 REC(2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X02AJF
      INTEGER P01ABF
      EXTERNAL X02AJF,P01ABF
C     ..
C     .. External Subroutines ..
      EXTERNAL COARS4,GETVEC,MESH4,SOL4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
      IFO = IFAIL

      MCHEPS = X02AJF(1.D0)
      KXTRAP = IELAM3
      MAXTRP = IELAM3

C Error Trapping
      IF (TOL.LE.MCHEPS**0.6D0 .OR. A.GE.B .OR. ABS(IFAIL).GT.1 .OR.
     &    ABS(ERR).LE.0.D0 .OR. K.LT.0 .OR. IWORK.LT.NCOARS) THEN
C We have an input error.
          IF (ABS(IFAIL).GT.1) THEN
              WRITE (REC,FMT=9000) IFAIL

 9000         FORMAT (' ** Illegal initial IFAIL value of',I8,/,
     &               ' ** IFAIL reset to 1')

              IFAIL = 1
              IFO = 1
              NREC = 2
              GO TO 180

          END IF

          IFAIL = 1
          IF (TOL.LE.MCHEPS**0.6D0) THEN
              WRITE (REC,FMT=9010) MCHEPS**0.6D0,TOL

 9010         FORMAT (' ** On entry TOL was',G18.8,/,
     &               ' ** TOL must be greater than',G18.8)

              NREC = 2
              GO TO 180

          END IF

          IF (A.GE.B) THEN
              WRITE (REC,FMT=9020) A,B

 9020         FORMAT (' ** A must be less than B; you had',/,' ** A=',
     &               G18.8,' and B=',G18.8)

              NREC = 2
              GO TO 180

          END IF

          IF (ERR.EQ.0.D0) THEN
              WRITE (REC,FMT=9030)

 9030         FORMAT (' ** ERR must be non-zero on entry')

              NREC = 1
              GO TO 180

          END IF

          IF (K.LT.0) THEN
              WRITE (REC,FMT=9040) K

 9040         FORMAT (' ** K must be non-negative; you had',/,' ** K=',
     &               I8)

              NREC = 2
              GO TO 180

          END IF

          IF (IWORK.LT.20) THEN
              WRITE (REC,FMT=9050) IWORK

 9050         FORMAT (' ** IWORK must be at least 20;',/,
     &               ' ** you had IWORK=',I8)

              NREC = 2
              GO TO 180

          END IF

      END IF

      MCHEPS = X02AJF(0.D0)
C First we get a very crude initial approximation to the eigenvalue

      IF (SYM) THEN
          XSHIFT = 0.5D0* (A+B)
      ELSE
          XSHIFT = 0.D0
      END IF

      NMESH = NCOARS
      ALOC = (A-XSHIFT)/ (1.D0+ABS(A-XSHIFT))
      BLOC = (B-XSHIFT)/ (1.D0+ABS(B-XSHIFT))
      CALL COARS4(NMESH,ALOC,BLOC,WORK(0,2),3)
C Now unmap the mesh back into the original interval
      DO 10 I = 0,NMESH
          WORK(I,2) = WORK(I,2)/ (1.D0-ABS(WORK(I,2))) + XSHIFT
   10 CONTINUE

      DO 20 I = 1,NMESH
          CALL SL4COF((WORK(I-1,2)+WORK(I,2))/2.D0,WORK(I,3),WORK(I,4),
     &                WORK(I,5),WORK(I,6))
   20 CONTINUE

      IMATCH = NMESH/2
      ELAM1 = ELAM
      EPS = ABS(ERR)
      TOL1 = SQRT(MCHEPS)
      KNTMAX = 100
      IFAIL = 0

      CALL SOL4(WORK(0,2),NMESH,IMATCH,WORK(1,3),WORK(1,4),WORK(1,5),
     &          WORK(1,6),ELAM1,EPS,EL,EU,K,TOL1,KNTMAX,SL4BCS,IFAIL)
      IF (IFAIL.NE.0) THEN
          IF (IFAIL.EQ.2) THEN
              WRITE (REC,FMT=9060)

 9060         FORMAT (' ** On first call to SOL4, cannot do zero-count',
     &               /,' ** correction at centre: unsuitable U and V.')

          ELSE IF (IFAIL.EQ.3) THEN
              WRITE (REC,FMT=9070)

 9070         FORMAT (' ** On first call to SOL4, cannot do zero-count',
     &               /,
     &        ' ** correction at centre: overflow computing detU, detV.'
     &               )

          ELSE IF (IFAIL.EQ.4) THEN

              WRITE (REC,FMT=9080) IFAIL

 9080         FORMAT (
     &              ' ** On first call to SOL4, an error occurred while'
     &               ,/,' ** processing a Theta matrix; IFAIL=',I8)

          ELSE IF (IFAIL.EQ.6) THEN
              WRITE (REC,FMT=9090)

 9090         FORMAT (' ** On first call to SOL4, cannot do zero-count',
     &               /,
     &          ' ** correction on a mesh interval: unsuitable U and V.'
     &               )

          ELSE IF (IFAIL.EQ.7) THEN
              WRITE (REC,FMT=9100)

 9100         FORMAT (' ** On first call to SOL4, cannot do zero-count',
     &               /,' ** correction on a mesh interval: overflow',
     &               ' computing detU, detV.')

          ELSE IF (IFAIL.EQ.9) THEN
              WRITE (REC,FMT=9110)

 9110         FORMAT (' ** On first call to SOL4, no eigenvalue can be',
     &               /,' ** found; ERR is too small.')

          ELSE IF (IFAIL.EQ.10) THEN
              WRITE (REC,FMT=9120)

 9120         FORMAT (
     &               ' ** On first call to SOL4, an eigenvalue has been'
     &               ,/,
     &               ' ** bracketed but cannot be accurately located.')

          END IF

          NREC = 2
          ELAM = ELAM1
          WORK(0,3) = ELAM1
          GO TO 180

      END IF

      ELAM0 = ELAM1

      IY1 = 19* (NMESH+2) + 11
      IY2 = IY1 + 4* (NMESH+1)
C WORKS must be of length at least 27*NMESH+65

      CALL GETVEC(WORK(0,2),NMESH,IMATCH,WORK(1,3),WORK(1,4),WORK(1,5),
     &            WORK(1,6),ELAM1,MULTI,WORKS,WORKS(IY1),WORKS(IY2),
     &            SL4BCS,IFAIL)

C Update the matchpoint
      IF (.NOT.SYM) THEN
          I = IMATCH
C          OSC = (WORK(I,4) + SQRT(MAX(0.D0,WORK(I,4)**2+4.D0*WORK(I,
C     &          3)* (ELAM1*WORK(I,6)-WORK(I,5)))))/WORK(I,3)
C          DO 30 I = 1,NMESH - 1
C            OSCL = (WORK(I,4) + SQRT(MAX(0.D0,WORK(I,4)**2+4.D0*WORK(I,
C     &               3)* (ELAM1*WORK(I,6)-WORK(I,5)))))/WORK(I,3)

          OSC = 4.D0*WORK(I,3)* (ELAM1*WORK(I,6)-WORK(I,5))/WORK(I,3)
          DO 30 I = 1,NMESH - 1
              OSCL = 4.D0*WORK(I,3)* (ELAM1*WORK(I,6)-WORK(I,5))/
     &               WORK(I,3)
              IF (OSCL.GT.OSC) THEN
                  OSC = OSCL
                  IMATCH = I
              END IF

   30     CONTINUE
      END IF

      MACHPT = WORK(IMATCH,2)

C Now use this eigenvalue approximation to assist in meshing.

C Store the matchpoint and the number of mesh-intervals of the coarse
C initial mesh:
      IM1ST = IMATCH
      NM1ST = NMESH

C Set meshing tolerance for automatic meshing:
C      TOLI = 1.D-2
C      TOLI = 0.5D0*SQRT(TOL)
C      TOLI = TOL
      TOLI = 1.D0
      COUNT = 0
      NMAX = IWORK
   40 IFAIL = 0

      CALL MESH4(MACHPT,A,ELAM1,SL4COF,WORK(0,1),NMAX,0,IMESH,TOLI,
     &           WORKS(IY1),WORKS(IY2),MULTI,IM1ST,NM1ST,WORK(0,2),
     &           IWORK,IFAIL)
      IF (IFAIL.NE.0) THEN
          IF (COUNT.LT.3) THEN
              TOLI = TOLI*5.D0
              COUNT = COUNT + 1
              GO TO 40
          END IF
          GO TO 60
      END IF

      IF (SYM) THEN
C Generate the second half of the mesh by reflection
          NMESH = 2*IMESH
          DO 50 I = 1,IMESH
              WORK(IMESH+I,1) = 2.D0*MACHPT - WORK(IMESH-I,1)
   50     CONTINUE
          GO TO 80
      END IF

      CALL MESH4(MACHPT,B,ELAM1,SL4COF,WORK(0,1),NMAX,IMESH,NMESH,TOLI,
     &           WORKS(IY1),WORKS(IY2),MULTI,IM1ST,NM1ST,WORK(0,2),
     &           IWORK,IFAIL)
      IF (IFAIL.NE.0) THEN
          IF (COUNT.LT.3) THEN
              TOLI = TOLI*5.D0
              COUNT = COUNT + 1
              GO TO 40

          END IF

      END IF

   60 IF (IFAIL.EQ.0) THEN
          NMESH = NMESH + IMESH
          IMATCH = IMESH
      ELSE
C The meshing has failed. Adopt contingency measures:
          NMESH = NCOARS
          IMATCH = NMESH/2
          IF (SYM) THEN
              XSHIFT = 0.5D0* (A+B)

          ELSE
              XSHIFT = 0.D0
          END IF

          ALOC = (A-XSHIFT)/ (1.D0+ABS(A-XSHIFT))
          BLOC = (B-XSHIFT)/ (1.D0+ABS(B-XSHIFT))
          CALL COARS4(NMESH,ALOC,BLOC,WORK(0,1),3)
C Now unmap the mesh back into the original interval
          DO 70 I = 0,NMESH
              WORK(I,1) = WORK(I,1)/ (1.D0-ABS(WORK(I,1))) + XSHIFT
   70     CONTINUE
      END IF

C Evaluate the coefficients
   80 DO 90 I = 1,NMESH
          CALL SL4COF((WORK(I-1,1)+WORK(I,1))/2.D0,WORK(I,2),WORK(I,3),
     &                WORK(I,4),WORK(I,5))
   90 CONTINUE

      IREFIN = 0
      RXTRAP = .FALSE.
      EPS = ABS(ERR)

C Update the matchpoint again.
  100 IF (SYM) THEN
          IMATCH = NMESH/2
      ELSE
          I = IMATCH
C          OSC = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &          2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

          OSC = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/WORK(I,2)

          DO 110 I = 1,NMESH - 1
C            OSCL = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &               2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

              OSCL = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/
     &               WORK(I,2)

              IF (OSCL.GT.OSC) THEN
                  OSC = OSCL
                  IMATCH = I
              END IF

  110     CONTINUE
      END IF
C We are now ready to commence computation of decent eigenvalue
C approximations.

      IFAIL = 0
      TOL1 = MCHEPS**0.75D0

      CALL SOL4(WORK(0,1),NMESH,IMATCH,WORK(1,2),WORK(1,3),WORK(1,4),
     &          WORK(1,5),ELAM1,EPS,EL,EU,K,TOL1,KNTMAX,SL4BCS,IFAIL)

      WORK(0,3) = ELAM1
      WRITE (6,FMT=*) 'NMESH,ELAM1:',NMESH,ELAM1

      IF (IFAIL.NE.0) THEN
          IF (IFAIL.EQ.2) THEN
              WRITE (REC,FMT=9130)

 9130         FORMAT (' ** On 2nd call to SOL4, cannot do zero-count',/,
     &               ' ** correction at centre: unsuitable U and V.')

          ELSE IF (IFAIL.EQ.3) THEN
              WRITE (REC,FMT=9140)

 9140         FORMAT (' ** On 2nd call to SOL4, cannot do zero-count',/,
     &        ' ** correction at centre: overflow computing detU, detV.'
     &               )

          ELSE IF (IFAIL.EQ.4) THEN

              WRITE (REC,FMT=9150) IFAIL

 9150         FORMAT (' ** On 2nd call to SOL4, an error occurred while'
     &               ,/,' ** processing a Theta matrix; IFAIL=',I8)

          ELSE IF (IFAIL.EQ.6) THEN
              WRITE (REC,FMT=9160)

 9160         FORMAT (' ** On 2nd call to SOL4, cannot do zero-count',/,
     &          ' ** correction on a mesh interval: unsuitable U and V.'
     &               )

          ELSE IF (IFAIL.EQ.7) THEN
              WRITE (REC,FMT=9170)

 9170         FORMAT (' ** On 2nd call to SOL4, cannot do zero-count',/,
     &               ' ** correction on a mesh interval: overflow',
     &               ' computing detU, detV.')

          ELSE IF (IFAIL.EQ.9) THEN
              WRITE (REC,FMT=9180)

 9180         FORMAT (' ** On 2nd call to SOL4, no eigenvalue can be',/,
     &               ' ** found; ERR is too small.')

          ELSE IF (IFAIL.EQ.10) THEN
              WRITE (REC,FMT=9190)

 9190         FORMAT (' ** On 2nd call to SOL4, an eigenvalue has been',
     &               /,' ** bracketed but cannot be accurately located.'
     &               )

          END IF

          NREC = 2
          ELAM = ELAM1
          GO TO 180

      END IF

C Update the matchpoint:
      IF (.NOT.SYM) THEN
          I = IMATCH
C          OSC = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &          2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

          OSC = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/WORK(I,2)

          DO 120 I = 1,NMESH - 1
C            OSCL = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &               2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

              OSCL = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/
     &               WORK(I,2)

              IF (OSCL.GT.OSC) THEN
                  OSC = OSCL
                  IMATCH = I
              END IF

  120     CONTINUE
      END IF

      EPS = ABS(ELAM1-ELAM0)
      ELAM0 = ELAM1
      IF (EPS.LT.0.1D0*MAX(1.D0,ABS(ELAM1))) THEN
C The error seems to be small enough to justify commencing
C Richardson extrapolation
          ELAM3(1,1) = ELAM1
          KR = 2
          RXTRAP = .TRUE.
C         WRITE (6,FMT=*) 'Started Richardson Extrapolation'
      END IF
      EPS = MAX(0.1D0,EPS)
  130 DO 140 II = NMESH,1,-1
          WORK(2*II,1) = WORK(II,1)
          WORK(2*II-1,1) = 0.5D0* (WORK(II,1)+WORK(II-1,1))
  140 CONTINUE
      NMESH = 2*NMESH

      DO 150 II = 1,NMESH
          CALL SL4COF((WORK(II-1,1)+WORK(II,1))/2.D0,WORK(II,2),
     &                WORK(II,3),WORK(II,4),WORK(II,5))
  150 CONTINUE
C Now update the matchpoint properly:

      IF (.NOT.SYM) THEN
          I = IMATCH*2
C          OSC = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &          2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

          OSC = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/WORK(I,2)

          DO 160 I = 1,NMESH - 1
C            OSCL = (WORK(I,3) + SQRT(MAX(0.D0,WORK(I,3)**2+4.D0*WORK(I,
C     &               2)* (ELAM1*WORK(I,5)-WORK(I,4)))))/WORK(I,2)

              OSCL = 4.D0*WORK(I,2)* (ELAM1*WORK(I,5)-WORK(I,4))/
     &               WORK(I,2)

              IF (OSCL.GT.OSC) THEN
                  OSC = OSCL
                  IMATCH = I
              END IF
  160     CONTINUE
      ELSE
          IMATCH = IMATCH*2
      END IF

      IF (.NOT.RXTRAP) THEN
          IREFIN = IREFIN + 1
          IF (IREFIN.LT.10 .AND. 2*NMESH.LT.IWORK) GO TO 100
          IFAIL = 8
          ELAM = ELAM1
          WRITE (REC,FMT=9200)
 9200     FORMAT (' ** The eigenvalue approximations are varying',/,
     &           ' ** wildly as the mesh is refined')
          NREC = 2
          GO TO 180
      END IF

      TOL1 = MCHEPS**0.75D0
      JMAX = MIN(KXTRAP,KR-1)
      ELAM2 = ELAM3(KR-1,JMAX)
      EPS = ABS(ERR)
      IF (KR.GT.2) EPS = ABS(ELAM3(KR-1,1)-ELAM3(KR-2,1))
      IF (KR.GT.2) TOL1 = MAX(TOL1,0.01D0*
     &                    (EPS/ (1.D0+ABS(ELAM3(KR-2,1))))**2)

      IFAIL = 0

      CALL SOL4(WORK(0,1),NMESH,IMATCH,WORK(1,2),WORK(1,3),WORK(1,4),
     &          WORK(1,5),ELAM2,EPS,EL,EU,K,TOL1,KNTMAX,SL4BCS,IFAIL)

      WORK(0,3) = ELAM2

      IF (IFAIL.NE.0) THEN
          IF (IFAIL.EQ.2) THEN
              WRITE (REC,FMT=9210)

 9210         FORMAT (' ** On later call to SOL4, cannot do zero-count',
     &               /,' ** correction at centre: unsuitable U and V.')

          ELSE IF (IFAIL.EQ.3) THEN
              WRITE (REC,FMT=9220)

 9220         FORMAT (' ** On later call to SOL4, cannot do zero-count',
     &               /,
     &        ' ** correction at centre: overflow computing detU, detV.'
     &               )

          ELSE IF (IFAIL.EQ.4) THEN

              WRITE (REC,FMT=9230) IFAIL

 9230         FORMAT (
     &              ' ** On later call to SOL4, an error occurred while'
     &               ,/,' ** processing a Theta matrix; IFAIL=',I8)

          ELSE IF (IFAIL.EQ.6) THEN
              WRITE (REC,FMT=9240)

 9240         FORMAT (' ** On later call to SOL4, cannot do zero-count',
     &               /,
     &          ' ** correction on a mesh interval: unsuitable U and V.'
     &               )

          ELSE IF (IFAIL.EQ.7) THEN
              WRITE (REC,FMT=9250)

 9250         FORMAT (' ** On later call to SOL4, cannot do zero-count',
     &               /,' ** correction on a mesh interval: overflow',
     &               'computing detU, detV.')

          ELSE IF (IFAIL.EQ.9) THEN
              WRITE (REC,FMT=9260)

 9260         FORMAT (' ** On later call to SOL4, no eigenvalue can be',
     &               /,' ** found; ERR is too small.')

          ELSE IF (IFAIL.EQ.10) THEN
              WRITE (REC,FMT=9270)

 9270         FORMAT (
     &               ' ** On later call to SOL4, an eigenvalue has been'
     &               ,/,
     &               ' ** bracketed but cannot be accurately located.')

          END IF

          NREC = 2
          GO TO 180

      END IF

      ELAM3(KR,1) = ELAM2

      JMAX = MIN(KXTRAP,KR)

C This loop may do one more Richardson extrapolation than is strictly
C justified, but what the hell -- if EPS passes the error test it
C will probably be OK.

      DO 170 J = 2,JMAX
          EPS = (ELAM3(KR,J-1)-ELAM3(KR-1,J-1))/
     &          (2.D0** (2* (J-1))-1.D0)
          ELAM3(KR,J) = ELAM3(KR,J-1) + EPS
  170 CONTINUE

C Check the error:

C      ELAM2 = ELAM3(KR,KR)
      ELAM2 = ELAM3(KR,JMAX)
      EPS = ABS(EPS)

      ELAM = ELAM2
      ERR = EPS
      NXTRAP = KR
      IF (EPS.LT.TOL*MAX(1.D0,ABS(ELAM))) THEN
C We have got our approximation
          IFAIL = 0
          RETURN

      END IF

C Now look through the Richardson extrapolation table and see if there
C is evidence of growing Richardson corrections. If the Richardson
C corrections are growing then put an upper limit on the number of
C extrapolations to be attempted.

      IF (KR.GT.2) THEN
C First check that the mesh-size is small enough to justify
C extrapolation
          CORR2 = ELAM3(KR,1) - ELAM3(KR-1,1)
          IF (ABS(CORR2).GT.0.1D0*MAX(1.D0,ABS(ELAM3(KR,1)))) THEN
C              WRITE (6,FMT=*)
C     &        'The mesh appears too coarse: extrapolation restarted'
C The mesh is too coarse for extrapolation to be justified yet.
              IF (KR.GT.MAXTRP) THEN
                  IFAIL = 8
                  WRITE (REC,FMT=9280)
 9280             FORMAT (
     &               ' ** More than 10 extrapolations would be required'
     &                   ,/,' ** to obtain the required accuracy')
                  NREC = 2
                  GO TO 180
              END IF
C              ELAM3(1,1) = ELAM3(KR,1)
C              MAXTRP = MAXTRP - KR
C              KR = 2
C              EPS = 2.D0*ABS(CORR2)
              RXTRAP = .FALSE.
              ELAM1 = ELAM3(1,1)
              MAXTRP = MAXTRP - KR
              GO TO 130
          END IF
      END IF

      IF (JMAX.GT.2 .AND. KR.GT.2) THEN
          CORR2 = (ELAM3(KR,JMAX-2)-ELAM3(KR-1,JMAX-2))/
     &            MAX(1.D0,ABS(ELAM3(KR,1)))
          CORR1 = (ELAM3(KR-1,JMAX-2)-ELAM3(KR-2,JMAX-2))/
     &            MAX(1.D0,ABS(ELAM3(KR,1)))
          TEORAT = 0.5D0** (2* (JMAX-1))
          IF (ABS(CORR1).LT.X02AJF(0.D0) .AND. ABS(CORR2).GT.TOL) THEN
C Since CORR2 did not satisfy the error test we are in a situation
C where CORR2 is `large' and CORR1 is `small'. This indicates
C that Richardson extrapolation is probably not justified.
              KXTRAP = JMAX
C              WRITE (6,FMT=*) 'Richardson table curtailed (1)'
          END IF
          IF (CORR2/CORR1.GT.2.D0*TEORAT) THEN
C The last step of Richardson extrapolation was probably
C not justified. Fix KXTRAP accordingly.
              KXTRAP = JMAX
C              WRITE (6,FMT=*) 'Richardson table curtailed (2)'
          END IF
      END IF

C      WRITE (6,FMT=*) 'KXTRAP=',KXTRAP

      IF (KR.GE.MAXTRP) THEN
          IFAIL = 8

          WRITE (REC,FMT=9290)

 9290     FORMAT (' ** More than 10 extrapolations would be required',/,
     &           ' ** to obtain the required accuracy')

          NREC = 2
          GO TO 180

      END IF

      IF (2*NMESH.GT.IWORK) THEN
          IFAIL = 11
          WRITE (REC,FMT=9300)

 9300     FORMAT (' ** IWORK needs to be at least twice as large',/,
     &           ' ** to obtain the required accuracy')

          NREC = 2
          GO TO 180

      END IF

      KR = KR + 1
      GO TO 130

C Error handling:
  180 IFAIL = P01ABF(IFO,IFAIL,SRNAME,NREC,REC)
      RETURN

      END



C -------------------------------------------------------------------


      SUBROUTINE SL4ARR(A,B,NMESH,XMESH,PP,SP,QP,WP,SL4COF)
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B
      INTEGER NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PP(NMESH),QP(NMESH),SP(NMESH),WP(NMESH),
     &                 XMESH(0:NMESH)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4COF
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
      XMESH(0) = A
      DO 10 I = 1,NMESH
          XMESH(I) = A + (B-A)*DBLE(I)/DBLE(NMESH)
          CALL SL4COF((XMESH(I)+XMESH(I-1))*0.5D0,PP(I),SP(I),QP(I),
     &                WP(I))
   10 CONTINUE

      RETURN

      END


C -------------------------------------------------------------------


      SUBROUTINE SOL4(XMESH,NMESH,IMATCH,P,S,Q,W,ELAM,EPSO,EL,EU,K,TOL,
     &                KNTMAX,SL4BCS,IFAIL)

C This routine solves a 4th order SL problem with piecewise constant
C coefficients to within a tolerance TOL

C     .. Scalar Arguments ..
      DOUBLE PRECISION EL,ELAM,EPSO,EU,TOL
      INTEGER IFAIL,IMATCH,K,KNTMAX,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION P(NMESH),Q(NMESH),S(NMESH),W(NMESH),
     &                 XMESH(0:NMESH)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DL,DU,E1,E2,EPS,SHIFT,TOLOC
      INTEGER IFLAG,IND,KNT
      LOGICAL DUAS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(17)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION SL4MIS
      EXTERNAL SL4MIS
C     ..
C     .. External Subroutines ..
      EXTERNAL C05AZF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,MAX,MIN
C     ..
      IFAIL = 0
      EPS = EPSO
      SHIFT = DBLE(K) + 5.0D-1
      IF (EPS.EQ.0.D0) THEN
          EPS = 0.1D0*MAX(1.D0,ABS(ELAM))
      END IF

      EL = ELAM - EPS
      EU = ELAM + EPS
      KNT = 0
      DUAS = .FALSE.

   10 DL = SL4MIS(XMESH,NMESH,IMATCH,P,S,Q,W,EL,SHIFT,SL4BCS,IFAIL)
      IF (IFAIL.NE.0) RETURN
C      WRITE (6,FMT=*) 'EL, DL:',EL,DL

      IF (DL.GT.0.D0) THEN
          KNT = KNT + 1
          EU = EL
          DU = DL
          EPS = 2.D0*EPS
          EL = EL - EPS
          DUAS = .TRUE.
          GO TO 10

      END IF

      IF (DUAS) GO TO 30
      DU = SL4MIS(XMESH,NMESH,IMATCH,P,S,Q,W,EU,SHIFT,SL4BCS,IFAIL)
      IF (IFAIL.NE.0) RETURN
C      WRITE (6,FMT=*) 'EU, DU:',EU,DU

   20 IF (DL*DU.GT.0.D0) THEN
C We have not yet bracketed an eigenvalue
          KNT = KNT + 1
          IF (KNT.GT.KNTMAX) THEN
              IFAIL = 9
              RETURN

          END IF

          EPS = 2.D0*EPS
          IF (DU.LT.0.D0) THEN
              EL = EU
              DL = DU
              EU = EU + EPS
              DU = SL4MIS(XMESH,NMESH,IMATCH,P,S,Q,W,EU,SHIFT,SL4BCS,
     &             IFAIL)
              IF (IFAIL.NE.0) RETURN
C              WRITE (6,FMT=*) 'EU, DU (2):',EU,DU

          ELSE
              EU = EL
              DU = DL
              EL = EL - EPS
              DL = SL4MIS(XMESH,NMESH,IMATCH,P,S,Q,W,EL,SHIFT,SL4BCS,
     &             IFAIL)
              IF (IFAIL.NE.0) RETURN
C              WRITE (6,FMT=*) 'EL, DL (2):',EL,DL
          END IF

          GO TO 20

      END IF

C If we are here it means that we have bracketed an eigenvalue and we
C can now think about locating it accurately with C05AZF

   30 KNT = 0
      IND = -1
      TOLOC = TOL
      C(1) = DU
      IFLAG = 1
      E1 = EL
      E2 = EU

   40 CALL C05AZF(E1,E2,DL,TOLOC,1,C,IND,IFLAG)

      IF (IND.EQ.0) GO TO 50

      KNT = KNT + 1
      IF (KNT.GT.KNTMAX) THEN
          IFAIL = 10
          RETURN

      END IF

      DL = SL4MIS(XMESH,NMESH,IMATCH,P,S,Q,W,E1,SHIFT,SL4BCS,IFAIL)
      IF (IFAIL.NE.0.D0) RETURN

      IF (DL.GT.0.D0) EU = MIN(EU,E1)
      IF (DL.LT.0.D0) EL = MAX(EL,E1)

C      WRITE (6,FMT=*) 'ELAM D(ELAM),NMESH:',E1,DL,NMESH

      IF (ABS(EU-EL).LT.TOLOC*MAX(1.D0,ABS(EL))) GO TO 50

      GO TO 40

   50 ELAM = E1
      RETURN

      END


C -------------------------------------------------------------------


      DOUBLE PRECISION FUNCTION SL4MIS(XMESH,NMESH,IMATCH,PP,SP,QP,WP,
     &                 ELAM,SHIFT,SL4BCS,IFAIL)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,SHIFT
      INTEGER IFAIL,IMATCH,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PP(1:NMESH),QP(1:NMESH),SP(1:NMESH),WP(1:NMESH),
     &                 XMESH(0:NMESH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ARGC,OMEGA1,OMEGA2,XEND,XO
      INTEGER I,IEND,NL,NR,SIGMA
      LOGICAL ISING,PRN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION UC(2,2),UL(2,2),UR(2,2),VC(2,2),VL(2,2),VR(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X01AAF,X02AJF
      EXTERNAL X01AAF,X02AJF
C     ..
C     .. External Subroutines ..

      EXTERNAL CORRCT,F06YAF,SPTH4,ZCOUNT
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,MAX,MIN
C     ..
      IFAIL = 0
      PRN = .FALSE.
      SL4MIS = 1.D0/X02AJF(0.D0)
C Get the boundary condition at the left-hand endpoint:
      CALL SL4BCS(0,ISING,XMESH(0),UL,VL,ELAM)
      NL = 0

C Now step across the mesh from XMESH(0) to XMESH(IMATCH),
C carrying the matrices UL, VL such that Theta = (VL+iUL).inv(VL-iUL)
C and carrying ARGL which is Argdet(Theta).
      DO 10 I = 1,IMATCH
          XO = XMESH(I-1)
          XEND = XMESH(I)
          IEND = 2
          IF (I.EQ.1) IEND = 0
C The routine ZCOUNT does the stepping across the interval
C [XMESH(I-1),XMESH(I)] on which the coefficients P, S, Q and
C W have the constant values PP(I), SP(I), QP(I) and WP(I).
C The eigenparameter has the value ELAM. PRN is a logical
C variable which controls printing of intermediate results
C and IFAIL is an integer which is used to flag errors.

          CALL ZCOUNT(UL,VL,PP(I),SP(I),QP(I),WP(I),ELAM,XO,XEND,NL,
     &                IEND,IFAIL)
          IF (IFAIL.NE.0) RETURN
   10 CONTINUE


C Get the boundary condition at the right-hand endpoint:
      CALL SL4BCS(1,ISING,XMESH(NMESH),UR,VR,ELAM)
      NR = 0

      DO 20 I = NMESH,IMATCH + 1,-1
          IEND = 2
          IF (I.EQ.NMESH) IEND = 1
          CALL ZCOUNT(UR,VR,PP(I),SP(I),QP(I),WP(I),ELAM,XMESH(I),
     &                XMESH(I-1),NR,IEND,IFAIL)
          IF (IFAIL.NE.0) RETURN
   20 CONTINUE

C Compute the central correction sigma:

      IEND = 2
      IF (IMATCH.EQ.0) IEND = 0
      IF (IMATCH.EQ.NMESH) IEND = 1
      CALL CORRCT(UL,VL,UR,VR,SIGMA,IEND,IFAIL)
      IF (IFAIL.NE.0) THEN
          IFAIL = IFAIL + 1
          RETURN
      END IF

      SL4MIS = DBLE(NL+NR+SIGMA) - SHIFT

C Compute the phase angles of Theta_{R}^{*}Theta_{L}
      CALL F06YAF('T','N',2,2,2,1.d0,VR,2,UL,2,0.d0,UC,2)
      CALL F06YAF('T','N',2,2,2,-1.D0,UR,2,VL,2,1.D0,UC,2)

      CALL F06YAF('T','N',2,2,2,1.D0,VR,2,VL,2,0.D0,VC,2)
      CALL F06YAF('T','N',2,2,2,1.D0,UR,2,UL,2,1.D0,VC,2)

      CALL SPTH4(UC,VC,OMEGA1,OMEGA2,ARGC,'L',IFAIL)
      IF (IFAIL.NE.0) THEN
          IFAIL = 4
          RETURN
      END IF

      IF (SL4MIS.GE.0.D0 .AND. SL4MIS.LT.1.D0) THEN
          SL4MIS = MIN(OMEGA1,OMEGA2)
      ELSE IF (SL4MIS.LT.0.D0 .AND. SL4MIS.GT.-1.D0) THEN
          SL4MIS = MAX(OMEGA1,OMEGA2)*0.5D0/X01AAF(0.D0) - 1.D0
      END IF

      RETURN

      END


C -------------------------------------------------------------------


      SUBROUTINE ZCOUNT(U,V,P,S,Q,W,ELAM,XO,XEND,N,IEND,IFAIL)

C This subroutine integrates the matrices U,V of a 4th order
C Sturm-Liouville problem across an interval on which the
C coefficients are constant, and calculates the nullity count N.

C The formula for the contribution of the current interval to N is
C
C             N = NO + SIGMA
C
C where NO is the nullity count across (XO, XEND) of the solution
C UO, satisfying Dirichlet boundary conditions at XEND, and SIGMA
C is the correction parameter CORR at XO, calculated by the
C subroutine CORRCT.
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,W,XEND,XO
      INTEGER IEND,IFAIL,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,B,C,DISC1,DISC2,DUMMY,M,SCAL
      INTEGER I,IFLAG,J,NO,SIGMA
      LOGICAL LINTSP,LNEG,LPOS,SPCIAL
C     ..
C     .. Local Arrays ..

      DOUBLE PRECISION CU(2,2),CV(2,2),RNORM(2,2),SU(2,2),SV(2,2),
     &                 UO(2,2),UTEMP(2,2),VO(2,2),VTEMP(2,2),Y1(4),
     &                 Y2(4),Y3(4),Y4(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL CORRCT,F06YAF,FUNSOL,INTNLM,LARGLM,NEGLM,SL4CNT
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC ABS,ATAN2,COS,MAX,SIN,SQRT
C     ..

C Copy the matrices into temporary storage for later
C use.
      DO 20 J = 1,2
          DO 10 I = 1,2
              UTEMP(I,J) = U(I,J)
              VTEMP(I,J) = V(I,J)
   10     CONTINUE
   20 CONTINUE
C
C      CALL FUNSOL(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4,SCAL)

      DISC1 = 4.D0*P* (ELAM*W-Q)
      DISC2 = S*S + DISC1

      LPOS = DISC1 .GE. 0.D0
      LINTSP = DISC2 .GE. 0.D0 .AND. S - SQRT(ABS(DISC2)) .GE. 0.D0
      LNEG = DISC2 .LT. 0.D0

      IF (LPOS .OR. LINTSP) THEN
          SCAL = ABS(XEND-XO)*SQRT(0.5D0* (S+SQRT(ABS(DISC2)))/P)

      ELSE IF (LNEG) THEN
          SCAL = ABS(XEND-XO)*0.5D0*SQRT((S+SQRT(ABS(DISC1)))/P)

      ELSE
          SCAL = 0.D0
      END IF

      SPCIAL = SCAL .GT. 1.D-1

      IF (.NOT.SPCIAL) THEN
C In this case we must accept the output from FUNSOL as correct
          CALL FUNSOL(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4,SCAL)

          CU(1,1) = Y1(1)
          CU(1,2) = Y2(1)
          CU(2,1) = Y1(2)
          CU(2,2) = Y2(2)

          SU(1,1) = Y3(1)
          SU(1,2) = Y4(1)
          SU(2,1) = Y3(2)
          SU(2,2) = Y4(2)

          SV(1,1) = Y1(3)
          SV(1,2) = Y2(3)
          SV(2,1) = Y1(4)
          SV(2,2) = Y2(4)

          CV(1,1) = Y3(3)
          CV(1,2) = Y4(3)
          CV(2,1) = Y3(4)
          CV(2,2) = Y4(4)
C
C Now do the update of U,V from xo to xend, using the formulae
C               u(xend) = cu.u(xo)+su.v(xo),
C               v(xend) = sv.u(xo)+cv.v(xo).
C
          CALL F06YAF('N','N',2,2,2,1.D0,CU,2,UTEMP,2,0.D0,UO,2)
          CALL F06YAF('N','N',2,2,2,1.D0,SU,2,VTEMP,2,1.D0,UO,2)

          CALL F06YAF('N','N',2,2,2,1.D0,CV,2,VTEMP,2,0.D0,VO,2)
          CALL F06YAF('N','N',2,2,2,1.D0,SV,2,UTEMP,2,1.D0,VO,2)
C
C Now do some postmultiplication to avoid ill-conditioning
C of these matrices.
C
          B = -UO(1,1)*VO(2,2) + UO(1,2)*VO(2,1) - UO(2,2)*VO(1,1) +
     &        UO(2,1)*VO(1,2)
          AC = UO(1,1)*UO(2,2) - UO(2,1)*UO(1,2) - VO(1,1)*VO(2,2) +
     &         VO(2,1)*VO(1,2)
          C = UO(1,1)*UO(2,2) - UO(2,1)*UO(1,2) + VO(1,1)*VO(2,2) -
     &        VO(2,1)*VO(1,2)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

          IF (C.GE.0.D0) THEN
              M = ATAN2(B,AC)

          ELSE
              M = ATAN2(-B,-AC)
          END IF

          C = COS(M/2.D0)
          AC = SIN(M/2.D0)

C Form the normalising matrix:

          DO 40 J = 1,2
              DO 30 I = 1,2
                  RNORM(I,J) = UO(I,J)*C - VO(I,J)*AC
   30         CONTINUE
   40     CONTINUE

C Form the inverse of the normalising matrix:

          C = RNORM(1,1)*RNORM(2,2) - RNORM(1,2)*RNORM(2,1)
          AC = RNORM(1,1)/C
          RNORM(1,1) = RNORM(2,2)/C
          RNORM(2,2) = AC
          RNORM(1,2) = -RNORM(1,2)/C
          RNORM(2,1) = -RNORM(2,1)/C

C Peform the postmultiplication: U = UO.RNORM,  V = VO.RNORM

          CALL F06YAF('N','N',2,2,2,1.D0,UO,2,RNORM,2,0.D0,U,2)
          CALL F06YAF('N','N',2,2,2,1.D0,VO,2,RNORM,2,0.D0,V,2)

C
C Now set the U and V (evaluated at X = 0) for Dirichlet conditions
C at X =  XEND.  (Thus U(XEND) = 0,  V(XEND) = I.)
          UO(1,1) = -SU(1,1)
          UO(1,2) = SU(1,2)
          UO(2,1) = SU(2,1)
          UO(2,2) = -SU(2,2)

          VO(1,1) = CV(1,1)
          VO(1,2) = -CV(1,2)
          VO(2,1) = -CV(2,1)
          VO(2,2) = CV(2,2)


      ELSE
          IF (LPOS) THEN
              CALL LARGLM(ELAM,P,S,Q,W,U,V,UO,VO,XO,XEND)
          ELSE IF (LINTSP) THEN
              CALL INTNLM(ELAM,P,S,Q,W,U,V,UO,VO,XO,XEND)
          ELSE IF (LNEG) THEN
              CALL NEGLM(ELAM,P,S,Q,W,U,V,UO,VO,XO,XEND)
          END IF

      END IF
C Normalise U and V so that the greatest element of the two is
C equal to 1.
      DUMMY = 0.D0
      DO 60 J = 1,2
          DO 50 I = 1,2
              DUMMY = MAX(DUMMY,ABS(U(I,J)),ABS(V(I,J)))
   50     CONTINUE
   60 CONTINUE

      DO 80 J = 1,2
          DO 70 I = 1,2
              U(I,J) = U(I,J)/DUMMY
              V(I,J) = V(I,J)/DUMMY
   70     CONTINUE
   80 CONTINUE

C
C We have now updated U, V from XO to XEND.
C Now we will calculate the correction parameter SIGMA
C
      IF (XO.LE.XEND) THEN
          CALL CORRCT(UTEMP,VTEMP,UO,VO,SIGMA,IEND,IFLAG)
      ELSE
          CALL CORRCT(UO,VO,UTEMP,VTEMP,SIGMA,IEND,IFLAG)
      END IF
      IF (IFLAG.NE.0) THEN
          IF (IFLAG.EQ.1) IFAIL = 6
          IF (IFLAG.EQ.2) IFAIL = 7
          RETURN
      END IF
C
C Next find the nullity count for UO.
C
      CALL SL4CNT(XEND-XO,P,S,Q,W,ELAM,NO)
C
C Now we can update the total nullity count N:
C
      N = N + NO + SIGMA

      END

C --------------------------------------------------------------------

      SUBROUTINE CORRCT(UL,VL,UR,VR,CORR,IEND,IFLAG)
C
C  This subroutine calculates the correction term "corr"
C  in the formula  N(lambda) = No(lambda) + corr.
C  The correction term is calculated at an endpoint xo of
C  a mesh interval.  At this point we have input data (UL,VL)
C  and (UR,VR).  There are 4 cases to consider, depending on
C  whether det(UL) and det(UR) are or are not zero.  It is convenient
C  to consider a 5th case, when UL or UR is the zero matrix.  Thus,
C  we have the following cases:
C    (1)  det(UL) and det(UR) are not zero.
C    (2)  One of the matrices UL or UR is zero.
C    (3)  det(UL) = 0;  UL and det(UR) are not zero.
C    (4)  det(UL) and UR are not zero;   det(UR) = 0.
C    (5)  det(UL) = det(UR) =0;  UL and UR are not zero.
C  The integer IFLAG is normally 0, but equals 1 if
C  det(DIFF) > 0, while trace(DIFF) = 0,
C  where DIFF = WL - WR, and equals 2 if the matrices U
C  and V do not fit any of the cases.
C
C  IEND tells CORRCT if it is being called at the end of an
C  interval.
C  IEND = 0 means left endpoint
C  IEND = 1 means right endpoint
C  IEND = 2 means not an endpoint.
C
C--------------------------------------------------------------------
C
C
C     .. Scalar Arguments ..
      INTEGER CORR,IEND,IFLAG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION UL(2,2),UR(2,2),VL(2,2),VR(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,ADL,ADR,AL,AR,DET,DUL,DUR,DVL,DVR,F,G,KL,KR,ML,
     &                 MR,NL,NR,SGN,SGNL,SGNR,SL,SR,TRC
      LOGICAL KN,MN,ZL,ZR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DIFF(2,2),VUAL(2,2),VUAR(2,2)
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
C     ..
      IFLAG = 0

      DUL = UL(1,1)*UL(2,2) - UL(1,2)*UL(2,1)
      DUR = UR(1,1)*UR(2,2) - UR(1,2)*UR(2,1)
C
      DVL = VL(1,1)*VL(2,2) - VL(1,2)*VL(2,1)
      DVR = VR(1,1)*VR(2,2) - VR(1,2)*VR(2,1)
C

      VUAL(1,1) = UL(2,2)*VL(1,1) - UL(2,1)*VL(1,2)
      VUAL(1,2) = UL(2,2)*VL(2,1) - UL(2,1)*VL(2,2)
      VUAL(2,1) = VUAL(1,2)
      VUAL(2,2) = UL(1,1)*VL(2,2) - UL(1,2)*VL(2,1)
C
      VUAR(1,1) = UR(2,2)*VR(1,1) - UR(2,1)*VR(1,2)
      VUAR(1,2) = UR(2,2)*VR(2,1) - UR(2,1)*VR(2,2)
      VUAR(2,1) = VUAR(1,2)
      VUAR(2,2) = UR(1,1)*VR(2,2) - UR(1,2)*VR(2,1)
C
      SGNL = SIGN(1.D0,DUL)
      SGNR = SIGN(1.D0,DUR)
C
C  Case (1):  det(UL) and det(UR) are not zero.

      IF (DUL.NE.0.D0 .AND. DUR.NE.0.D0) THEN

          DIFF(1,1) = VUAL(1,1)*ABS(DUR)*SGNL - VUAR(1,1)*ABS(DUL)*SGNR
          DIFF(1,2) = VUAL(1,2)*ABS(DUR)*SGNL - VUAR(1,2)*ABS(DUL)*SGNR
          DIFF(2,1) = DIFF(1,2)
          DIFF(2,2) = VUAL(2,2)*ABS(DUR)*SGNL - VUAR(2,2)*ABS(DUL)*SGNR
C
          DET = DIFF(1,1)*DIFF(2,2) - DIFF(1,2)*DIFF(2,1)
          TRC = DIFF(1,1) + DIFF(2,2)


          IF (DET.LT.0.D0) THEN
              CORR = 1
              RETURN
C
          ELSE IF (DET.GT.0.D0) THEN

              IF (TRC.LT.0.D0) THEN
                  CORR = 2
                  RETURN

              ELSE IF (TRC.GT.0.D0) THEN
                  CORR = 0
                  RETURN

              ELSE IF (TRC.EQ.0.D0) THEN
                  IFLAG = 1
                  RETURN

              END IF

          ELSE IF (DET.EQ.0.D0) THEN

              IF (TRC.LT.0.D0) THEN
                  CORR = 1
                  RETURN

              ELSE
                  CORR = 0
                  RETURN

              END IF
          END IF
      END IF
C
C
C  Case (2):   UL = 0 or UR = 0.
C
C

      ZL = UL(1,1) .EQ. 0.D0 .AND. UL(1,2) .EQ. 0.D0 .AND.
     &     UL(2,1) .EQ. 0.D0 .AND. UL(2,2) .EQ. 0.D0
      ZR = UR(1,1) .EQ. 0.D0 .AND. UR(1,2) .EQ. 0.D0 .AND.
     &     UR(2,1) .EQ. 0.D0 .AND. UR(2,2) .EQ. 0.D0
C
C
      IF (ZL .OR. ZR) THEN

          IF (IEND.EQ.0) THEN

              IF (ZL) THEN
                  CORR = 0
                  RETURN

              ELSE IF (DUL.EQ.0.D0) THEN
                  CORR = 1
                  RETURN

              ELSE
                  CORR = 2
                  RETURN

              END IF

          ELSE IF (IEND.EQ.1) THEN

              IF (ZR) THEN
                  CORR = 0
                  RETURN

              ELSE IF (DUR.EQ.0.D0) THEN
                  CORR = 1
                  RETURN

              ELSE
                  CORR = 2
                  RETURN

              END IF

          ELSE IF (IEND.EQ.2) THEN
              CORR = 2
              RETURN

          END IF
      END IF
C
C
      KL = VUAL(1,1)
      ML = VUAL(1,2)
      NL = VUAL(2,2)
C
      KR = VUAR(1,1)
      MR = VUAR(1,2)
      NR = VUAR(2,2)
C
C
C  Case (3):  det(UL) = 0;  UL and det(UR) are not zero.
C
C
      IF (DUL.EQ.0 .AND. DUR.NE.0) THEN
          SGN = SIGN(1.D0,KL+NL)
          A = ABS(KL+NL)
          ADR = ABS(DUR)

          IF (KL.NE.0.D0) THEN
              F = (KR*ML*ML-2*MR*KL*ML+NR*KL*KL)*A*SGNR
              G = DVL* (KL*KL+ML*ML)*ADR*SGN

              IF (F.LE.G) THEN

                  IF (IEND.EQ.0) THEN
                      CORR = 0
                      RETURN

                  ELSE
                      CORR = 1
                      RETURN

                  END IF

              ELSE

                  IF (IEND.EQ.0) THEN
                      CORR = 1
                      RETURN

                  ELSE
                      CORR = 2
                      RETURN

                  END IF

              END IF

          ELSE
              F = KR*A*SGNR
              G = DVL*ADR*SGN

              IF (F.LE.G) THEN

                  IF (IEND.EQ.0) THEN
                      CORR = 0
                      RETURN

                  ELSE
                      CORR = 1
                      RETURN

                  END IF

              ELSE

                  IF (IEND.EQ.0) THEN
                      CORR = 1
                      RETURN

                  ELSE
                      CORR = 2
                      RETURN

                  END IF

              END IF
          END IF

      END IF
C
C
C  Case (4):  det(UL) and UR are not zero;  det (UR) = 0.
C
C
      IF (DUL.NE.0 .AND. DUR.EQ.0) THEN
          SGN = SIGN(1.D0,KR+NR)
          A = ABS(KR+NR)
          ADL = ABS(DUL)

          IF (KR.NE.0.D0) THEN
              F = (KL*MR*MR-2*ML*KR*MR+NL*KR*KR)*A*SGNL
              G = DVR* (KR*KR+MR*MR)*ADL*SGN

              IF (F.GE.G) THEN

                  IF (IEND.EQ.1) THEN
                      CORR = 0
                      RETURN

                  ELSE
                      CORR = 1
                      RETURN

                  END IF

              ELSE

                  IF (IEND.EQ.1) THEN
                      CORR = 1
                      RETURN

                  ELSE
                      CORR = 2
                      RETURN

                  END IF

              END IF

          ELSE
              F = KL*A*SGNL
              G = DVR*ADL*SGN

              IF (F.GE.G) THEN

                  IF (IEND.EQ.1) THEN
                      CORR = 0
                      RETURN

                  ELSE
                      CORR = 1
                      RETURN

                  END IF

              ELSE

                  IF (IEND.EQ.1) THEN
                      CORR = 1
                      RETURN

                  ELSE
                      CORR = 2
                      RETURN

                  END IF

              END IF
          END IF

      END IF
C
C
C  Case (5):  det(UL) = det(UR) = 0;  UL and UR are not zero.
C
C
      IF (DUL.EQ.0.D0 .AND. DUR.EQ.0.D0) THEN
          KN = KL*NR .EQ. KR*NL
          MN = ML*NR .EQ. MR*NL

          AL = ABS(KL+NL)
          SL = SIGN(1.D0,KL+NL)
          AR = ABS(KR+NR)
          SR = SIGN(1.D0,KR+NR)

          F = DVL*AR*SL
          G = DVR*AL*SR

          IF (KN .AND. MN .AND. F.GE.G) THEN

              IF (IEND.EQ.2) THEN
                  CORR = 1
                  RETURN

              ELSE
                  CORR = 0
                  RETURN

              END IF

          ELSE

              IF (IEND.EQ.2) THEN
                  CORR = 2
                  RETURN

              ELSE
                  CORR = 1
                  RETURN

              END IF
          END IF
      END IF

      IFLAG = 2

      RETURN

      END


C --------------------------------------------------------------------

      SUBROUTINE LARGLM(ELAM,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,AOMX,B,BETA,C,COSHCS,COSHNX,COSHSN,COSOMX,
     &                 DETCU,DETCV,DETSU,DETSV,DETU,DETV,DISC,FAC,H,H2,
     &                 H2SNSN,M,NO,NO2,NU,NU2,NU3,NU4,NU6,NUX,OMEGA,
     &                 OMEGA2,OMEGA3,OMEGA4,OMEGA6,OMX,PNO,PNO2,PNORAD,
     &                 PRAD,PRAD2,RAD,RAD2,SINHCS,SINHNX,SINHSN,SINOMX,
     &                 TRUVA
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(2,2),CUSVA(2,2),SUCVA(2,2),T(2,2),UUA(2,2),
     &                 UVA(2,2),VUA(2,2),VVA(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION PHI
      EXTERNAL PHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,SIN,SQRT
C     ..
C     .. External Subroutines ..

C     ..
      DISC = SQRT(S*S+4.D0*P* (ELAM*W-Q))
      NU = SQRT(0.5d0* (S+DISC)/P)
      OMEGA = SQRT(0.5d0* (-S+DISC)/P)

      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)
      H = XEND - XO
      H2 = H*H
      NUX = ABS(NU*H)
      FAC = EXP(-NUX)
      IF (NUX.GE.0.1D0) THEN
          SINHNX = (1.D0-FAC*FAC)/ (2.D0*NUX)
      ELSE
          SINHNX = PHI(NUX*NUX)*FAC
      END IF

      COSHNX = (1.D0+FAC*FAC)/2.D0
      OMX = OMEGA*H
      AOMX = ABS(OMX)
      IF (AOMX.GE.0.1D0) THEN
          SINOMX = SIN(OMX)/OMX
      ELSE
          SINOMX = PHI(-OMX*OMX)
      END IF

      COSOMX = COS(OMX)

      COSHSN = H*COSHNX*SINOMX
      SINHCS = H*SINHNX*COSOMX
      COSHCS = COSHNX*COSOMX
      SINHSN = SINHNX*SINOMX
      H2SNSN = H2*SINHSN
      OMEGA2 = OMEGA*OMEGA
      NU2 = NU*NU
      OMEGA3 = OMEGA*OMEGA2
      NU3 = NU*NU2
      OMEGA4 = OMEGA2*OMEGA2
      NU4 = NU2*NU2
      OMEGA6 = OMEGA2*OMEGA4
      NU6 = NU2*NU4
      RAD = OMEGA2 + NU2
      RAD2 = RAD*RAD
      PRAD = P*RAD
      PRAD2 = PRAD*PRAD
      NO = NU*OMEGA
      PNO = P*NO
      PNORAD = PNO*RAD
      NO2 = NO*NO
      PNO2 = P*NO2

C Compute SU.Adjoint(CV).exp(-NUX):

      SUCVA(1,1) = (COSHSN-SINHCS)/PRAD
      SUCVA(1,2) = ((OMEGA2-NU2)* (FAC-COSHCS)+2.D0*NO2*H2SNSN)/
     &             (PRAD*RAD)
      SUCVA(2,1) = SUCVA(1,2)
      SUCVA(2,2) = (NU2*SINHCS+OMEGA2*COSHSN)/PRAD

C Compute CU.Adjoint(SV).exp(-NUX):

      CUSVA(1,1) = P* (NU4*SINHCS-OMEGA4*COSHSN)/RAD
      CUSVA(1,2) = -PNO2* ((NU4+OMEGA4)*H2SNSN+
     &             (NU2-OMEGA2)* (FAC-COSHCS))/RAD2
      CUSVA(2,1) = CUSVA(1,2)
      CUSVA(2,2) = -PNO2* (NU2*COSHSN+OMEGA2*SINHCS)/RAD

C Compute A = U.Adjoint(V):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = V(1,1)*U(2,2) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.A.Adjoint(SV))*exp(-NUX):

      T(1,1) = A(1,1)*COSHCS + 2.D0*A(1,2)* (OMEGA2*COSHSN+NU2*SINHCS)/
     &         RAD + A(2,2)*H2SNSN

      T(1,2) = (A(1,1)*NO2* (SINHCS-COSHSN)+
     &         A(1,2)* ((OMEGA2-NU2)* (FAC* (OMEGA2-NU2)+
     &         2.D0*NO2*H2SNSN)+4.D0*NO2*COSHCS)/RAD+
     &         A(2,2)* (NU2*SINHCS+OMEGA2*COSHSN))/RAD

      T(2,1) = T(1,2)


      T(2,2) = -A(1,1)*NO2*H2SNSN + 2.D0*A(1,2)*NO2* (SINHCS-COSHSN)/
     &         RAD + A(2,2)*COSHCS

C Compute U.Adjoint(V) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      DETCU = ((NU4+OMEGA4)*COSHCS+NO2* ((NU2-OMEGA2)*H2SNSN+2.D0*FAC))/
     &        RAD2
      DETSU = ((NU2-OMEGA2)*H2SNSN+2.D0* (FAC-COSHCS))/PRAD2
      DETCV = DETCU
      DETSV = P*PNO* ((OMEGA6-NU6)*NO*H2SNSN+
     &        2.D0*NU3*OMEGA3* (FAC-COSHCS))/RAD2

C Compute trace(CU.A.Adjoint(SU)).exp(-NUX):

      BETA = A(1,1)* (OMEGA2*COSHSN+NU2*SINHCS)/PRAD -
     &       2.D0*A(1,2)* ((NU2-OMEGA2)* (FAC-COSHCS)-2.D0*NO2*H2SNSN)/
     &       (P*RAD2) + A(2,2)* (COSHSN-SINHCS)/PRAD

      UUA(1,1) = DETU*DETCU + DETV*DETSU + BETA
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)

C Compute trace(CV.A.Adjoint(SV)).exp(-NUX):

      BETA = (-A(1,1)*PNO2* (NU2*COSHSN+OMEGA2*SINHCS)-
     &       2.D0*PNO2*A(1,2)* ((NU2-OMEGA2)* (FAC-COSHCS)+ (NU4+
     &       OMEGA4)*H2SNSN)/RAD-A(2,2)*P* (OMEGA4*COSHSN-NU4*SINHCS))/
     &       RAD

      VVA(1,1) = DETU*DETSV + DETV*DETCV + BETA
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

      DO 20 J = 1,2
          DO 10 I = 1,2
              U(I,J) = UUA(I,J)*COS(M/2.D0) - UVA(I,J)*SIN(M/2.D0)
              V(I,J) = VUA(I,J)*COS(M/2.D0) - VVA(I,J)*SIN(M/2.D0)
   10     CONTINUE
   20 CONTINUE


C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

C NOTE: in the following formulae we want CV and SU for BACKWARDS
C integration. We exploit the fact that the diagonal terms of
C SUCVA are odd functions of the direction of integration while
C the off-diagonal terms are even functions of the direction
C of integration, while the determinants which appear in the
C formulae are even functions of the direction of integration.

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END

C -------------------------------------------------------------------

      SUBROUTINE INTNLM(ELAM,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,AEHH,AETAH,ARHH,B,BETA,C,CMSN,CNCM,CNSM,
     &                 COSHMX,COSHNX,DETCU,DETCV,DETSU,DETSV,DETU,DETV,
     &                 DISC,ETA,ETA2,ETAH,ETASQ,FAC,FACL,H,H2,M,MU,MU2,
     &                 MU3,MU4,MU6,MUX,NM,NU,NU2,NU3,NU4,NU6,NUX,PNM,
     &                 PNMR,PRAD,PRAD2,RAD,RAD2,RADH2,RH,RH2,RHH,RHSQ,
     &                 SINHMX,SINHNX,SN2ERT,SN2HRT,SNERAT,SNHR1,SNHR2,
     &                 SNHRAT,SNSM,TRUVA
      INTEGER I,J
      LOGICAL SMALL
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(2,2),CUSVA(2,2),SUCVA(2,2),T(2,2),UUA(2,2),
     &                 UVA(2,2),VUA(2,2),VVA(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION PHI,SNHDIF
      EXTERNAL PHI,SNHDIF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,MIN,SIGN,SIN,SQRT
C     ..
      DISC = SQRT(S*S+4.D0*P* (ELAM*W-Q))
      NU2 = 0.5D0* (S+DISC)/P
      MU2 = NU2 - DISC/P
      NU = SQRT(NU2)
      MU = SQRT(MU2)
      NU3 = NU*NU2
      MU3 = MU*MU2
      NU4 = NU2*NU2
      MU4 = MU2*MU2
      NU6 = NU3*NU3
      MU6 = MU3*MU3

      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)

      H = XEND - XO
      H2 = H*H
      NUX = ABS(NU*H)
      SINHNX = (1.D0-EXP(-2.D0*NUX))/2.D0
      COSHNX = 1.D0 - SINHNX
      SINHNX = SINHNX*SIGN(1.D0,H)

      MUX = ABS(MU*H)
      SINHMX = (1.D0-EXP(-2.D0*MUX))/2.D0
      COSHMX = 1.D0 - SINHMX
      SINHMX = SINHMX*SIGN(1.D0,H)

      FAC = EXP(-NUX-MUX)

      CMSN = COSHMX*SINHNX
      CNSM = COSHNX*SINHMX
      CNCM = COSHNX*COSHMX
      SNSM = SINHNX*SINHMX

C If any of the following quantities is small then we are
C in trouble; Taylor series expansions have to be used.

C      RAD = NU**2-MU**2
      RAD = DISC/P
      RAD2 = RAD*RAD
      PRAD = DISC
      PRAD2 = DISC*DISC
      NM = NU*MU
      PNM = P*NM
      PNMR = PNM*RAD

      SMALL = H2*MIN(ABS(RAD),ABS(PRAD),ABS(PNMR),PNM) .LT. 1.D-6

C Compute SU.Adjoint(CV).exp(-(NUX+MUX)):

      IF (SMALL) THEN
C Compute certain quantities which will be jolly useful
C when calculating the matrix elements
          ETA2 = NU + MU
          ETA = 0.5D0*ETA2
          ETASQ = ETA*ETA
          ETAH = ETA*H
          AETAH = ABS(ETAH)
          RH2 = NU - MU
          RH = 0.5D0*RH2
          RHSQ = RH*RH
          RHH = RH*H
          ARHH = ABS(RHH)

          IF (ABS(RHH).LT.0.1D0) THEN
              SNHRAT = PHI(RHH*RHH)
              SNHRAT = SNHRAT*SNHRAT*FAC
              SN2HRT = PHI(4.D0*RHH*RHH)*FAC

          ELSE
              FACL = EXP(-MUX)
              SNHRAT = 0.5D0* (1.D0-EXP(-2.D0*ARHH))/ARHH
              SNHRAT = SNHRAT*FACL
              SNHRAT = SNHRAT*SNHRAT
              SN2HRT = (1.D0-EXP(-4.D0*ARHH))/ (4.D0*ARHH)
              SN2HRT = SN2HRT*FACL*FACL
          END IF

          IF (ABS(ETAH).LT.0.1D0) THEN
              SNERAT = PHI(ETAH*ETAH)
              SNERAT = SNERAT*SNERAT*FAC
              SN2ERT = PHI(4.D0*ETAH*ETAH)*FAC

          ELSE
              SNERAT = 0.5D0* (1.D0-EXP(-2.D0*AETAH))/AETAH
              SNERAT = SNERAT*SNERAT
              SN2ERT = (1.D0-EXP(-4.D0*AETAH))/ (4.D0*AETAH)
          END IF

          IF (MUX.LT.0.1D0) THEN
              SNHR1 = PHI(MUX*MUX)*EXP(-MUX)

          ELSE
              SNHR1 = SIGN(1.D0,H)*SINHMX/MUX
          END IF

          IF (NUX.LT.0.1D0) THEN
              SNHR2 = PHI(NUX*NUX)*EXP(-NUX)

          ELSE
              SNHR2 = SIGN(1.D0,H)*SINHNX/NUX
          END IF

      END IF

      IF (.NOT.SMALL) THEN
          SUCVA(1,1) = (NU*CNSM-MU*CMSN)/PNMR
          SUCVA(1,2) = - (S* (FAC-CNCM)+2.D0*PNM*SNSM)/PRAD2
          SUCVA(2,1) = SUCVA(1,2)
          SUCVA(2,2) = (NU*CMSN-MU*CNSM)/PRAD

      ELSE
          RADH2 = ABS(RAD)*H2
          AEHH = ABS(ETASQ-RHSQ)*H2
          IF (RADH2.GE.0.01D0 .AND. AEHH.LT.0.01D0) THEN

              SUCVA(1,1) = H* (COSHNX*SNHR1-COSHMX*SNHR2)/PRAD

          ELSE IF (RADH2.LT.0.01D0 .AND. AEHH.GE.0.01D0) THEN

              SUCVA(1,1) = H* (SN2ERT-SN2HRT)/ (2.D0*PNM)

          ELSE

              SUCVA(1,1) = 2.D0* (H**3)*FAC*
     &                     SNHDIF(2.D0*AETAH,2.D0*ARHH)/P
          END IF

          SUCVA(1,2) = 0.25D0*H*H* (SNERAT+SNHRAT)/P
          SUCVA(2,1) = SUCVA(1,2)
          SUCVA(2,2) = 0.5D0*H* (SN2ERT+SN2HRT)/P

      END IF

C Compute CU.Adjoint(SV).exp(-NUX-MUX):

      IF (.NOT.SMALL) THEN
          CUSVA(1,1) = P* (NU3*CMSN-MU3*CNSM)/RAD
          CUSVA(1,2) = PNM* (NM* (NU2+MU2)* (FAC-CNCM)+ (NU4+MU4)*SNSM)/
     &                 RAD2
          CUSVA(2,1) = CUSVA(1,2)
          CUSVA(2,2) = PNM* (NU3*CNSM-MU3*CMSN)/RAD

      ELSE
          FACL = 0.25D0*PNM*H

          CUSVA(1,1) = P*H*0.5D0* ((ETASQ+3.D0*RHSQ)*SN2HRT+
     &                 (RHSQ+3.D0*ETASQ)*SN2ERT)


          CUSVA(1,2) = 0.25D0*PNM* (H2* (RHSQ*SNERAT-ETASQ*SNHRAT)+
     &                 3.D0*SNSM)
          CUSVA(2,1) = CUSVA(1,2)

          CUSVA(2,2) = 2.D0*FACL* ((RHSQ+3.D0*ETASQ)*SN2ERT-
     &                 (ETASQ+3.D0*RHSQ)*SN2HRT)
      END IF

C Compute A = U.Adjoint(V):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = V(1,1)*U(2,2) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.A.Adjoint(SV))*exp(-NUX-MUX):

      IF (.NOT.SMALL) THEN

          T(1,1) = A(1,1)*CNCM + 2.D0*A(1,2)* (NU*CMSN-MU*CNSM)/RAD +
     &             A(2,2)*SNSM/NM

          T(1,2) = (A(1,1)*NM* (NU*CNSM-MU*CMSN)+
     &             A(1,2)* ((FAC* (NU2+MU2)+2.D0*NM*SNSM)* (NU2+MU2)-
     &             4.D0*NM*NM*CNCM)/RAD+A(2,2)* (NU*CMSN-MU*CNSM))/RAD
          T(2,1) = T(1,2)

          T(2,2) = NM* (A(1,1)*SNSM+2.D0*A(1,2)* (NU*CNSM-MU*CMSN)/
     &             RAD) + A(2,2)*CNCM

      ELSE

          T(1,1) = A(1,1)*CNCM + 2.D0*P*A(1,2)*SUCVA(2,2) +
     &             A(2,2)*H2*SNHR1*SNHR2

          T(1,2) = A(1,1)*PNM*NM*SUCVA(1,1) +
     &             0.5D0*A(1,2)* (FAC+CNCM-H2*
     &             (ETASQ*SNHRAT+RHSQ*SNERAT)) + A(2,2)*P*SUCVA(2,2)
          T(2,1) = T(1,2)

          T(2,2) = NM* (A(1,1)*SNSM+2.D0*A(1,2)*PNM*SUCVA(1,1)) +
     &             A(2,2)*CNCM

      END IF

C Compute U.Adjoint(V) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      IF (.NOT.SMALL) THEN

          DETCU = ((NU4+MU4)*CNCM-NM* ((NU2+MU2)*SNSM+2.D0*NM*FAC))/RAD2
          DETCV = DETCU
          DETSU = (2.D0*NM* (FAC-CNCM)+ (NU2+MU2)*SNSM)/ (PRAD2*NM)
          DETSV = P*P*NM* (2.D0*NU3*MU3* (FAC-CNCM)+ (NU6+MU6)*SNSM)/
     &            RAD2

      ELSE

          DETCU = 0.25D0* (H2* (ETASQ*SNHRAT+RHSQ*SNERAT)+FAC+3.D0*CNCM)
          DETCV = DETCU
          RADH2 = ABS(RAD)*H2
          AEHH = ABS(ETASQ-RHSQ)*H2
          IF (RADH2.GE.0.01D0 .AND. AEHH.LT.0.01D0) THEN
              FACL = H2*SNHR1*SNHR2

              DETSU = (2.D0* (FAC-CNCM)+ (NU2+MU2)*FACL)/PRAD2

          ELSE IF (RADH2.LT.0.01D0 .AND. AEHH.GE.0.01D0) THEN
              FACL = H*0.5D0/P

              DETSU = (SNERAT-SNHRAT)*FACL*FACL/ (ETASQ-RHSQ)

          ELSE
              FACL = 0.5D0*H2/P

              DETSU = (SQRT(SNERAT)-SQRT(SNHRAT))*SNHDIF(AETAH,ARHH)*
     &                FACL*FACL
          END IF

          DETSV = (P*PNM/8.D0)* (2.D0*H*H*
     &            (RHSQ*RHSQ*SNERAT-ETASQ*ETASQ*SNHRAT)+
     &            3.D0*NM* (CNCM-FAC)+15.D0* (ETASQ+RHSQ)*SNSM)
      END IF


C Compute trace(CU.A.Adjoint(SU)).exp(-NUX-MUX):

C      IF (.NOT.SMALL) THEN
C
C          BETA = A(1,1)* (NU*CMSN-MU*CNSM)/PRAD +
C     &           2.D0*A(1,2)* ((NU2+MU2)* (CNCM-FAC)-2.D0*NM*SNSM)/
C     &           (P*RAD2) + A(2,2)* (NU*CNSM-MU*CMSN)/PNMR
C
C      ELSE

      BETA = A(1,1)*SUCVA(2,2) + 2.D0*A(1,2)*SUCVA(1,2) +
     &       A(2,2)*SUCVA(1,1)
C      END IF

      UUA(1,1) = DETU*DETCU + DETV*DETSU + BETA
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)


C Compute trace(CV.A.Adjoint(SV)).exp(-NUX-MUX):

C      IF (.NOT.SMALL) THEN
C
C          BETA = (A(1,1)*PNM* (NU3*CNSM-MU3*CMSN)+
C     &           2.D0*A(1,2)*PNM* ((NU4+MU4)*SNSM+NM* (NU2+MU2)* (FAC-
C     &           CNCM))/RAD+A(2,2)*P* (NU3*CMSN-MU3*CNSM))/RAD
C
C      ELSE

      BETA = A(1,1)*CUSVA(2,2) + 2.D0*A(1,2)*CUSVA(1,2) +
     &       A(2,2)*CUSVA(1,1)

C      END IF

      VVA(1,1) = DETU*DETSV + DETV*DETCV + BETA
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

      DO 20 J = 1,2
          DO 10 I = 1,2
              U(I,J) = UUA(I,J)*COS(M/2.D0) - UVA(I,J)*SIN(M/2.D0)
              V(I,J) = VUA(I,J)*COS(M/2.D0) - VVA(I,J)*SIN(M/2.D0)
   10     CONTINUE
   20 CONTINUE


C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF


C NOTE: in the following formulae we want CV and SU for BACKWARDS
C integration. We exploit the fact that the diagonal terms of
C SUCVA are odd functions of the direction of integration while
C the off-diagonal terms are even functions of the direction
C of integration, while the determinants which appear in the
C formulae are even functions of the direction of integration.

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END
C -------------------------------------------------------------------

      SUBROUTINE NEGLM(ELAM,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,ALFA,ALFA2,AX,B,BETA,BETA2,BX,C,COSBX,COSCOS,
     &                 COSHAX,CSHCSH,D,DETCU,DETCV,DETSU,DETSV,DETU,
     &                 DETV,FAC,H,M,PRAD,RAD,RADA,RADB,RADM,SINBX,
     &                 SINCOS,SINHAX,SINSIN,SNHCSH,SNHSNH,TR,TRUVA
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(2,2),CUSVA(2,2),SUCVA(2,2),T(2,2),UUA(2,2),
     &                 UVA(2,2),VUA(2,2),VVA(2,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,SIGN,SIN,SQRT
C     ..
C     .. External Functions ..

      DOUBLE PRECISION PHI
      EXTERNAL PHI
C     ..
      D = 2.D0*SQRT(P* (Q-ELAM*W))
      ALFA = 0.5D0*SQRT((S+D)/P)
      BETA = 0.5D0*SQRT((-S+D)/P)

      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)
      H = XEND - XO
      AX = ABS(ALFA*H)
      FAC = EXP(-AX)
      IF (AX.GE.0.1D0) THEN
          SINHAX = (1.D0-FAC*FAC)/2.D0
          COSHAX = 1.D0 - SINHAX
          SINHAX = SINHAX*SIGN(1.D0,H)/ALFA

      ELSE
          SINHAX = H*PHI(AX*AX)*FAC
          COSHAX = (1.D0+FAC*FAC)/2.D0
      END IF

      BX = BETA*H
      IF (ABS(BX).GE.0.1D0) THEN
          SINBX = SIN(BX)*FAC/BETA

      ELSE
          SINBX = H*PHI(-BX*BX)*FAC
      END IF

      COSBX = COS(BX)*FAC

      SNHCSH = SINHAX*COSHAX
      SNHSNH = SINHAX*SINHAX
      CSHCSH = COSHAX*COSHAX

      SINCOS = SINBX*COSBX
      SINSIN = SINBX*SINBX
      COSCOS = COSBX*COSBX

      ALFA2 = ALFA*ALFA
      BETA2 = BETA*BETA
      RAD = ALFA2 + BETA2
      PRAD = P*RAD
      RADM = ALFA2 - BETA2
      RADA = 3.D0*ALFA2 - BETA2
      RADB = 3.D0*BETA2 - ALFA2

C Compute SU.Adjoint(CV).exp(-AX):

      SUCVA(1,1) = 0.5D0* (SNHCSH-SINCOS)/PRAD
      SUCVA(1,2) = 0.25D0* (SNHSNH+SINSIN)/P
      SUCVA(2,1) = SUCVA(1,2)
      SUCVA(2,2) = 0.5D0* (SNHCSH+SINCOS)/P

C Compute CU.Adjoint(SV).exp(-AX):

      CUSVA(1,1) = 0.5D0* (RADA*SNHCSH-RADB*SINCOS)/P
      CUSVA(1,2) = -0.25D0*RAD* (SNHSNH*BETA2-3.D0*CSHCSH+SINSIN*ALFA2+
     &             3.D0*COSCOS)/P
      CUSVA(2,1) = CUSVA(1,2)
      CUSVA(2,2) = 0.5D0*PRAD* (RADA*SNHCSH+RADB*SINCOS)

C Compute A = U.Adjoint(V)(xo):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = U(2,2)*V(1,1) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.Adjoint(A).Adjoint(SV))*exp(-AX):

      T(1,1) = A(1,1)* (CSHCSH-SINSIN*BETA2) + A(1,2)* (SNHCSH+SINCOS) +
     &         A(2,2)* (CSHCSH-COSCOS)/RAD

      T(1,2) = 0.5D0*RAD*A(1,1)* (SNHCSH-SINCOS) +
     &         0.5D0*A(1,2)* (CSHCSH+COSCOS+BETA2*SNHSNH-ALFA2*SINSIN) +
     &         0.5D0*A(2,2)* (SNHCSH+SINCOS)

      T(2,1) = T(1,2)

      T(2,2) = A(1,1)*RAD* (CSHCSH-COSCOS) +
     &         A(1,2)*RAD* (SNHCSH-SINCOS) +
     &         A(2,2)* (CSHCSH-SINSIN*BETA2)

C Compute U.Adjoint(V)(xend) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      DETCU = 0.25D0*RADM* (SNHSNH+SINSIN) + 0.5D0* (CSHCSH+COSCOS)
      DETSU = 0.25* (SNHSNH-SINSIN)/ (P*PRAD)
      DETCV = DETCU
      DETSV = 0.25D0*P*PRAD* (RADA*RADA*SNHSNH-RADB*RADB*SINSIN)

C Compute trace(CU.A.Adjoint(SU)).exp(-AX):

      TR = 0.5D0* (A(1,1)* (SNHCSH+SINCOS)+A(1,2)* (SNHSNH+SINSIN))/P +
     &     0.5D0*A(2,2)* (SNHCSH-SINCOS)/PRAD

      UUA(1,1) = DETU*DETCU + DETV*DETSU + TR
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)

C Compute trace(CV.Adjoint(A).Adjoint(SU)).exp(-AX):

      TR = 0.5D0*P*RADA* (RAD* (A(1,1)*SNHCSH+A(1,2)*SNHSNH)+
     &     A(2,2)*SNHCSH) + 0.5D0*P*RADB*
     &     (RAD* (A(1,1)*SINCOS+A(1,2)*SINSIN)-A(2,2)*SINCOS)

      VVA(1,1) = DETU*DETSV + DETV*DETCV + TR
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

      U(1,1) = UUA(1,1)*COS(M/2.D0) - UVA(1,1)*SIN(M/2.D0)
      U(1,2) = -UVA(1,2)*SIN(M/2.D0)
      U(2,1) = -UVA(2,1)*SIN(M/2.D0)
      U(2,2) = UUA(2,2)*COS(M/2.D0) - UVA(2,2)*SIN(M/2.D0)

      V(1,1) = VUA(1,1)*COS(M/2.D0) - VVA(1,1)*SIN(M/2.D0)
      V(1,2) = VUA(1,2)*COS(M/2.D0)
      V(2,1) = VUA(2,1)*COS(M/2.D0)
      V(2,2) = VUA(2,2)*COS(M/2.D0) - VVA(2,2)*SIN(M/2.D0)


C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END


C -------------------------------------------------------------------


      SUBROUTINE FUNSOL(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4,SCAL)
C     .. Local Scalars ..
      DOUBLE PRECISION ALFA,ALFA1,ALFAH,ALFASQ,BETA,BETA1,BETAH,BETASQ,
     &                 COSBX,COSHAX,COSHHX,COSHNX,COSHX,COSNX,COSOMX,D,
     &                 DISC,ETA,ETAH,ETASQ,FAC,FACETA,FACMU,FACNU,FACRH,
     &                 H,MU,MUH,NU,NUH,NUSQ,OMEGA,OMEGA1,OMEGA2,OMH1,
     &                 OMH2,OMSQ,RAD,RH,RHH,RHSQ,SINBX,SINHAX,SINHHX,
     &                 SINHMX,SINHNX,SINHX,SINNX,SINOM1,SINOM2,SINOMX,
     &                 SIZE,YTEMP2,YTEMP3
      INTEGER I
C     ..
C     .. External Functions ..
C This subroutine finds the fundamental matrix for a fourth order
C Sturm-Liouville equation with constant coefficients
C
C py'''' - sy'' + qy = elam wy,
C
C on an interval [xo,xend]

C First measure the size of the coefficients.
      DOUBLE PRECISION PHI
      EXTERNAL PHI
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,SCAL,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION Y1(4),Y2(4),Y3(4),Y4(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL FUNMAT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,EXP,SIN,SQRT
C     ..
      SCAL = 0.D0
      H = XEND - XO
      SIZE = ABS(H)* (ABS((ELAM*W-Q)/P)+ABS(S/P))
      IF (ABS(SIZE).LT.0.02D0) THEN
C We can easily use the matrix exponential to get our fundamental
C solutions.
          CALL FUNMAT(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4)
          SCAL = 0.D0

      ELSE
C We have to do the solutions the hard way. There are four cases
C to consider, each of which will have to be subdivided yet
C further.

          ALFA1 = 0.5D0*S/P
          DISC = S*S + 4.D0*P* (ELAM*W-Q)
          BETA1 = 0.5D0*SQRT(ABS(DISC))/P

          IF (DISC.LT.0.D0) THEN

C We have complex roots. We must compute their real and imaginary
C parts. This is tha case of large negative lambda.
              D = 2.D0*SQRT(P* (Q-ELAM*W))
              ALFA = 0.5D0*SQRT((S+D)/P)
              BETA = 0.5D0*SQRT((D-S)/P)
              ALFAH = ABS(ALFA*H)
              SCAL = ALFAH
              BETAH = BETA*H
              FAC = EXP(-ALFAH)
              COSHAX = (1.D0+FAC*FAC)/2.D0
              IF (ALFAH.GT.0.1D0) THEN
                  SINHAX = (1.D0-FAC*FAC)/ (ALFAH+ALFAH)

              ELSE
                  SINHAX = PHI(ALFAH*ALFAH)*FAC
              END IF

              SINHAX = H*SINHAX
              COSBX = COS(BETAH)
              IF (ABS(BETAH).GT.0.1D0) THEN
                  SINBX = SIN(BETAH)/BETAH

              ELSE
                  SINBX = PHI(-BETAH*BETAH)
              END IF

              SINBX = H*SINBX

              BETASQ = BETA*BETA
              ALFASQ = ALFA*ALFA

              Y1(1) = COSHAX*COSBX + (BETASQ-ALFASQ)*SINHAX*SINBX/2.D0
              Y2(1) = (COSBX*SINHAX+COSHAX*SINBX)/2.D0
              Y3(1) = (COSBX*SINHAX-COSHAX*SINBX)/
     &                (2.D0*P* (ALFASQ+BETASQ))
              Y4(1) = SINHAX*SINBX/ (2.D0*P)

          ELSE IF (DISC.GE.0.D0 .AND. ALFA1.GT.BETA1) THEN

C We have real roots only. This is the case of intermediate
C lambda and positive s.
              NU = SQRT(ALFA1+BETA1)
              MU = SQRT(ALFA1-BETA1)

              ETA = 0.5D0* (NU+MU)
              RH = 0.5D0* (NU-MU)
              ETASQ = ETA*ETA
              RHSQ = RH*RH
              ETAH = ABS(ETA*H)
              RHH = ABS(RH*H)
              SCAL = ETAH + RHH
C This means that SCAL = ABS(NU*H)
              FACETA = EXP(-ETAH)
              FACRH = EXP(-RHH)
              COSHNX = (1.D0+FACETA*FACETA)/2.D0
              COSHHX = (1.D0+FACRH*FACRH)/2.D0
              IF (ETAH.GT.0.1D0) THEN
                  SINHNX = (1.D0-FACETA*FACETA)/ (ETAH+ETAH)

              ELSE
                  SINHNX = PHI(ETAH*ETAH)*FACETA
              END IF

              SINHNX = H*SINHNX
              IF (RHH.GT.0.1D0) THEN
                  SINHHX = (1.D0-FACRH*FACRH)/ (RHH+RHH)

              ELSE
                  SINHHX = PHI(RHH*RHH)*FACRH
              END IF

              SINHHX = H*SINHHX
              Y1(1) = COSHNX*COSHHX - 0.5D0* (ETASQ+RHSQ)*SINHNX*SINHHX
              Y2(1) = 0.5D0* (COSHNX*SINHHX+COSHHX*SINHNX)
              Y4(1) = 0.5D0*SINHNX*SINHHX/P
              IF (ETA*RH.GT.0.1D0) THEN
C This means that NU**2-MU**2 is not small

                  NUH = ABS(NU*H)
                  MUH = ABS(MU*H)
                  FACNU = EXP(-NUH)
                  FACMU = EXP(-MUH)
                  IF (NUH.GT.0.1D0) THEN
                      SINHNX = (1.D0-FACNU*FACNU)/ (NUH+NUH)

                  ELSE
                      SINHNX = PHI(NUH*NUH)*FACNU
                  END IF

                  SINHNX = H*SINHNX

                  IF (MUH.GT.0.1D0) THEN
                      SINHMX = (1.D0-FACMU*FACMU)/ (MUH+MUH)

                  ELSE
                      SINHMX = PHI(MUH*MUH)*FACMU
                  END IF

                  SINHMX = H*SINHMX*EXP(- (NUH-MUH))

                  Y3(1) = (SINHMX-SINHNX)/ (P* (NU**2-MU**2))

              ELSE

                  Y3(1) = (COSHHX*SINHNX-COSHNX*SINHHX)/
     &                    (2.D0*P* (ETASQ-RHSQ))
              END IF

          ELSE IF (DISC.GE.0.D0 .AND. ALFA1.LE.BETA1 .AND.
     &             ALFA1+BETA1.GE.0.D0) THEN

C We have two real roots and two complex roots. This is
C the case of large positive lambda.
              NUSQ = ALFA1 + BETA1
              OMSQ = BETA1 - ALFA1
              NU = SQRT(NUSQ)
              OMEGA = SQRT(OMSQ)
              RAD = 2.D0*BETA1
              NUH = ABS(NU*H)
              SCAL = NUH
              FAC = EXP(-NUH)
              IF (NUH.GT.0.1D0) THEN
                  SINHNX = (1.D0-FAC*FAC)/ (NUH+NUH)

              ELSE
                  SINHNX = PHI(NUH*NUH)*FAC
              END IF

              SINHNX = H*SINHNX

              OMH1 = OMEGA*H
              IF (ABS(OMH1).GT.0.1D0) THEN
                  SINOMX = SIN(OMH1)/OMH1

              ELSE
                  SINOMX = PHI(-OMH1*OMH1)
              END IF

              SINOMX = H*SINOMX*FAC
              COSOMX = COS(OMH1)*FAC
              COSHNX = (1.D0+FAC*FAC)/2.D0

              Y1(1) = (NUSQ*COSOMX+OMSQ*COSHNX)/RAD
              Y2(1) = (NUSQ*SINHNX+OMSQ*SINOMX)/RAD
              RAD = P*RAD
              Y3(1) = (SINOMX-SINHNX)/RAD
              Y4(1) = (COSHNX-COSOMX)/RAD



          ELSE IF (DISC.GE.0.D0 .AND. ALFA1.LE.BETA1 .AND.
     &             ALFA1+BETA1.LE.0.D0) THEN

C This case is the case of purely complex roots; it is the
C case of intermediate lambda and negative s.
              SCAL = 0.D0
              OMEGA1 = SQRT(ABS(ALFA1+BETA1))
              OMEGA2 = SQRT(ABS(ALFA1-BETA1))
              ETA = 0.5D0* (OMEGA1+OMEGA2)
              RH = 0.5D0* (OMEGA2-OMEGA1)
              ETASQ = ETA*ETA
              RHSQ = RH*RH
              ETAH = ETA*H
              RHH = RH*H
              COSNX = COS(ETAH)
              COSHX = COS(RHH)
              IF (ABS(ETAH).GT.0.1D0) THEN
                  SINNX = SIN(ETAH)/ETAH

              ELSE
                  SINNX = PHI(-ETAH*ETAH)
              END IF

              SINNX = H*SINNX
              IF (ABS(RHH).GT.0.1D0) THEN
                  SINHX = SIN(RHH)/RHH

              ELSE
                  SINHX = PHI(-RHH*RHH)
              END IF

              SINHX = H*SINHX

              Y1(1) = COSNX*COSHX + (ETASQ+RHSQ)*SINNX*SINHX/2.D0
              Y2(1) = (COSNX*SINHX+COSHX*SINNX)/2.D0
              Y4(1) = SINNX*SINHX/ (2.D0*P)
              IF (ETA*RH.GT.0.1D0) THEN
C OMEGA2**2-OMEGA1**2 is not small
                  OMH1 = OMEGA1*H
                  OMH2 = OMEGA2*H
                  IF (ABS(OMH1).GT.0.1D0) THEN
                      SINOM1 = SIN(OMH1)/OMH1

                  ELSE
                      SINOM1 = PHI(-OMH1*OMH1)
                  END IF

                  SINOM1 = H*SINOM1
                  IF (ABS(OMH2).GT.0.1D0) THEN
                      SINOM2 = SIN(OMH2)/OMH2

                  ELSE
                      SINOM2 = PHI(-OMH2*OMH2)
                  END IF

                  SINOM2 = H*SINOM2
                  Y3(1) = (SINOM2-SINOM1)/ (2.D0*P*BETA1)

              ELSE
                  Y3(1) = - (COSHX*SINNX-COSNX*SINHX)/
     &                    (2.D0*P* (ETASQ-RHSQ))
              END IF

          END IF

          DO 10 I = 1,3
              Y1(I+1) = (Q-ELAM*W)*Y3(I)
              Y2(I+1) = Y1(I) + S*Y4(I)
              Y3(I+1) = -Y4(I)
              Y4(I+1) = Y2(I)/P
   10     CONTINUE

          YTEMP2 = Y1(3)
          YTEMP3 = Y1(4)
          Y1(4) = P*YTEMP2
          Y1(3) = -P*YTEMP3 + S*Y1(2)

          YTEMP2 = Y2(3)
          YTEMP3 = Y2(4)
          Y2(4) = P*YTEMP2
          Y2(3) = -P*YTEMP3 + S*Y2(2)

          YTEMP2 = Y3(3)
          YTEMP3 = Y3(4)
          Y3(4) = P*YTEMP2
          Y3(3) = -P*YTEMP3 + S*Y3(2)

          YTEMP2 = Y4(3)
          YTEMP3 = Y4(4)
          Y4(4) = P*YTEMP2
          Y4(3) = -P*YTEMP3 + S*Y4(2)

      END IF

      RETURN

      END


C -------------------------------------------------------------------


C This subroutine finds the exponential of a small matrix
C arising from a 4th order SL problem

      SUBROUTINE FUNMAT(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION Y1(4),Y2(4),Y3(4),Y4(4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H,SCALAR
      INTEGER I
C     ..
C     .. Local Arrays ..

      DOUBLE PRECISION A(4,4),MNEST(4,4),STORE(4,4)
C     ..
C     .. External Subroutines ..
      EXTERNAL F06QFF,F06QHF,NEST
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC REAL
C     ..
      H = XEND - XO
C Set up the matrix A.
C First set all elements to zero:
      CALL F06QHF('G',4,4,0.D0,0.D0,A,4)
C Now fill in the non-zero elements:
      A(1,2) = H
      A(2,4) = H/P
      A(3,1) = H* (Q-ELAM*W)
      A(4,2) = H*S
      A(4,3) = -H
C Now do the nested multiplication. Start by setting
C mnest to be the identity.
      CALL F06QHF('G',4,4,0.D0,1.D0,MNEST,4)

      DO 10 I = 9,1,-1
          SCALAR = 1.D0/REAL(I)
          CALL NEST(MNEST,4,4,SCALAR,A,4,STORE,4)
   10 CONTINUE

C Now recover the matrix exponential into the vectors y1,y2,y3 and
C y4.


      CALL F06QFF('G',4,1,MNEST(1,1),4,Y1,4)
      CALL F06QFF('G',4,1,MNEST(1,2),4,Y2,4)
      CALL F06QFF('G',4,1,MNEST(1,3),4,Y3,4)
      CALL F06QFF('G',4,1,MNEST(1,4),4,Y4,4)

      RETURN

      END


C -------------------------------------------------------------------


      SUBROUTINE NEST(MNEST,IMNEST,N,SCALAR,FACTOR,IFACT,STORE,ISTOR)
C This routine performs the operation
C
C mnest = Identity + scalar.factor.mnest,
C
C where scalar is a real number, factor and mnest are nxn matrices,
C imnest is the first dimension of mnest as declared in the
C calling program and ifact is the first dimension of factor
C as declared in the calling program.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION SCALAR
      INTEGER IFACT,IMNEST,ISTOR,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FACTOR(IFACT,N),MNEST(IMNEST,N),STORE(ISTOR,N)
C     ..
C     .. External Subroutines ..
      EXTERNAL F06QFF,F06QHF,F06YAF
C     ..
      CALL F06QFF('G',N,N,MNEST,IMNEST,STORE,ISTOR)
      CALL F06QHF('G',N,N,0.D0,1.D0,MNEST,IMNEST)
      CALL F06YAF('N','N',N,N,N,SCALAR,FACTOR,IFACT,STORE,ISTOR,1.D0,
     &            MNEST,IMNEST)
      RETURN

      END


C ---------------------------------------------------------------------


      SUBROUTINE SPTH4(U,V,OMEGA1,OMEGA2,ARG,LR,IFAIL)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ARG,OMEGA1,OMEGA2
      INTEGER IFAIL
      CHARACTER LR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,BLOC,C,DEFALT,DISC,PI,R,SGNA,SGNB,SGNC,TWOPI
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X01AAF
      EXTERNAL X01AAF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,MAX,MIN,SIGN,SQRT
C     ..
C     .. External Subroutines ..

C     ..
      IFAIL = 0
      PI = X01AAF(0.D0)
      TWOPI = 2.D0*PI

      IF (LR.NE.'R') THEN
          DEFALT = 0.D0

      ELSE
          DEFALT = TWOPI
      END IF

      A = V(1,1)*V(2,2) - V(1,2)*V(2,1)
      B = U(1,2)*V(2,1) + U(2,1)*V(1,2) - U(1,1)*V(2,2) - U(2,2)*V(1,1)
      C = U(1,1)*U(2,2) - U(1,2)*U(2,1)

      R = MAX(ABS(A),ABS(B),ABS(C))
      A = A/R
      B = -B/R
      BLOC = B
      C = C/R

      IF (B.EQ.0.D0) THEN
          IF (C.EQ.0.D0) THEN
              OMEGA1 = DEFALT
              OMEGA2 = DEFALT

          ELSE IF (A.EQ.0.D0) THEN
              OMEGA1 = PI
              OMEGA2 = OMEGA1

          ELSE
              OMEGA2 = ATAN2(-C,A)
              IF (OMEGA1.LT.0.D0) OMEGA1 = OMEGA1 + TWOPI
              IF (OMEGA2.LT.0.D0) OMEGA2 = OMEGA2 + TWOPI
          END IF

      ELSE
          IF (A.EQ.0.D0) THEN
              OMEGA1 = PI
              IF (C.EQ.0.D0) THEN
                  OMEGA2 = DEFALT

              ELSE
                  OMEGA2 = 2.D0*ATAN2(C,B)
                  IF (OMEGA2.LT.0.D0) OMEGA2 = OMEGA2 + TWOPI
              END IF

          ELSE
              IF (C.EQ.0.D0) THEN
                  OMEGA1 = DEFALT
                  OMEGA2 = 2.D0*ATAN2(B,A)
                  IF (OMEGA2.LT.0.D0) OMEGA2 = OMEGA2 + TWOPI

              ELSE
C We are now in the situation where A, B and C are all non-zero.
                  SGNA = SIGN(1.D0,A)
                  SGNB = SIGN(1.D0,B)
                  SGNC = SIGN(1.D0,C)
                  B = ABS(B)
                  DISC = SQRT(ABS(B*B-4.D0*A*C))
                  IF (SGNA*SGNC.GT.0.D0) THEN
                      OMEGA1 = 2.D0*ATAN2(SGNB*ABS(B-DISC),2.D0*A)

                  ELSE
                      OMEGA1 = 2.D0*ATAN2(-SGNB*ABS(B-DISC),2.D0*A)
                  END IF

                  OMEGA2 = 2.D0*ATAN2(SGNB*ABS(B+DISC),2.D0*A)
                  IF (OMEGA1.LT.0.D0) OMEGA1 = OMEGA1 + TWOPI
                  IF (OMEGA2.LT.0.D0) OMEGA2 = OMEGA2 + TWOPI
              END IF

          END IF

      END IF

      ARG = OMEGA1
      OMEGA1 = MIN(OMEGA1,OMEGA2)
      OMEGA2 = MAX(ARG,OMEGA2)
      ARG = OMEGA1 + OMEGA2


      B = BLOC

      RETURN

      END


C -----------------------------------------------------------------


      SUBROUTINE SL4CNT(XEND,P,S,Q,W,ELAM,NZ)
C This routine counts the number of zeros of det(U(x)), where
C U is the 2x2 U-matrix of the Hamiltonian system associated
C with the 4th order SLDE (py'')''-(sy')'+qy = elam w y satisfying
C the initial condition U(0) = 0. The zeros counted are those
C which live in the range (0,xend] (or [xend,0) if xend<0).
C
C IMPORTANT NOTE: det(U(x)) is an EVEN function of x so we can
C always assume that xend>0.
C
C     .. Scalar Arguments ..
C The values of xend, p, s, q, w, and elam are all supplied as
C input. The value of nz -- the number of zeros -- is output.

      DOUBLE PRECISION ELAM,P,Q,S,W,XEND
      INTEGER NZ
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CS,CSH,DET,DIFF,DISC,NU,OMEGA,OMEGA1,OMEGA2,PI,
     &                 SDISC,SN,SNH,SUB,T,X
      INTEGER ICASE,IZ
C     ..
C     .. External Functions ..
      DOUBLE PRECISION PHI,X01AAF
      EXTERNAL PHI,X01AAF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,EXP,INT,MOD,SIN,SQRT
C     ..
      PI = X01AAF(0.D0)

C The first thing to do is to classify the case that we are in.

      DISC = S*S + 4.D0*P* (ELAM*W-Q)
      SDISC = SQRT(ABS(DISC))


      IF (DISC.LT.0.D0) THEN
C We are in the case of `large negative lambda'.
          ICASE = 1

      ELSE IF (DISC.GE.0.D0 .AND. S-SDISC.GE.0.D0) THEN
C This is the case of `Intermediate lambda and positive s'.
          ICASE = 2

      ELSE IF (DISC.GE.0.D0 .AND. S+SDISC.LE.0.D0) THEN
C This is the case of `Intermediate lambda and negative s'.
          OMEGA1 = SQRT(0.5D0* (SDISC-S)/P)
          OMEGA2 = SQRT(0.5D0* (-SDISC-S)/P)
          ICASE = 3

      ELSE IF (DISC.GE.0.D0 .AND. SDISC.GE.ABS(S)) THEN
C This is the case of `large positive lambda'.
          NU = SQRT(0.5D0* (S+SDISC)/P)
          OMEGA = SQRT(0.5D0* (SDISC-S)/P)
          ICASE = 4
      END IF

C Now we treat each case according to the theory in our paper.
C Cases 1 and 2 are a real scoosh, while 3 and 4 are much
C harder. Case 3 is the worst of the lot by a long way.

C Exploit the fact that det(U(x)) is an even function of x:
      X = ABS(XEND)

      GO TO (10,20,30,40) ICASE

C Case 1:
   10 NZ = 0
      RETURN

C Case 2:
   20 NZ = 0
      RETURN

C Case 3:
   30 IF (OMEGA1.EQ.OMEGA2) THEN
          NZ = 0

      ELSE
          IF (X*ABS(OMEGA1-OMEGA2).LE.2.D0*PI) THEN
              NZ = 0
C A correction is necessary under certain circumstances
              IF (X*ABS(OMEGA1+OMEGA2).GE.0.1D0) THEN
                  DIFF = ((OMEGA1-OMEGA2)**2)*
     &                   (1.D0-COS((OMEGA1+OMEGA2)*X)) -
     &                   ((OMEGA1+OMEGA2)**2)* (1.D0-
     &                   COS((OMEGA1-OMEGA2)*X))

              ELSE
                  DIFF = - (OMEGA1*OMEGA2* ((OMEGA1+OMEGA2)*X)**2*
     &                   ((OMEGA1-OMEGA2)*X)**2)/6.D0
              END IF
C              IF (X*ABS(OMEGA1-OMEGA2).EQ.2.D0*PI .OR.
C     &            DIFF.LE.0.D0) NZ = 1
              IF (X*ABS(OMEGA1-OMEGA2).EQ.2.D0*PI .OR.
     &            DIFF.GT.0.D0) NZ = 1

          ELSE
              IZ = INT(X*ABS(OMEGA1-OMEGA2)*0.5D0/PI)
              NZ = 2*IZ - 1
C A correction is necessary under certain circumstances
              IF (((OMEGA1-OMEGA2)**2)* (1.D0-COS((OMEGA1+
     &            OMEGA2)*X)).LE. ((OMEGA1+OMEGA2)**2)*
     &            (1.D0-COS((OMEGA1-OMEGA2)*X))) THEN
                  NZ = NZ + 1

              ELSE IF (((OMEGA1-OMEGA2)**2)*
     &                 (1.D0-COS((OMEGA1+OMEGA2)*X)).GE.
     &                 ((OMEGA1+OMEGA2)**2)* (1.D0-COS((OMEGA1-
     &                 OMEGA2)*X)) .AND. (OMEGA1+OMEGA2)*
     &                 SIN((OMEGA1+OMEGA2)*X)* (1.D0-COS((OMEGA1-
     &                 OMEGA2)*X)).GE. (OMEGA1-OMEGA2)*
     &                 SIN((OMEGA1-OMEGA2)*X)* (1.D0-COS((OMEGA1+
     &                 OMEGA2)*X))) THEN
                  NZ = NZ + 2
              END IF

          END IF

      END IF

      RETURN

C Case 4:
   40 IF (OMEGA.EQ.0.D0 .OR. X*OMEGA.LE.PI) THEN
C There are no zeros in this case
          NZ = 0
C ???
          T = NU*X
          IF (T.GT.0.1D0) THEN
              SNH = 0.5D0* (1.D0-EXP(-2.D0*T))/T
              CSH = 0.5D0* (1.D0+EXP(-2.D0*T))
              SUB = EXP(-T)

          ELSE
              SNH = PHI(T*T)
              CSH = 0.5D0* (EXP(T)+EXP(-T))
              SUB = 1.D0
          END IF

          T = OMEGA*X
          CS = COS(T)
          IF (T.GT.0.1D0) THEN
              SN = SIN(T)/T

          ELSE
              SN = PHI(-T*T)
          END IF

          DET = CS*CSH - SUB - 0.5D0* (NU**2-OMEGA**2)*SN*SNH

C ???
      ELSE
          NZ = INT(X*OMEGA/PI) - 1

C We may have to increment nz by 1 to get the zero-count right
          T = NU*X
          IF (T.GT.0.1D0) THEN
              SNH = 0.5D0* (1.D0-EXP(-2.D0*T))/T
              CSH = 0.5D0* (1.D0+EXP(-2.D0*T))
              SUB = EXP(-T)

          ELSE
              SNH = PHI(T*T)
              CSH = 0.5D0* (EXP(T)+EXP(-T))
              SUB = 1.D0
          END IF

          T = OMEGA*X
          CS = COS(T)
          IF (T.GT.0.1D0) THEN
              SN = SIN(T)/T

          ELSE
              SN = PHI(-T*T)
          END IF

          DET = CS*CSH - SUB - 0.5D0* (NU**2-OMEGA**2)*SN*SNH

C If the determinant has not changed sign and the number of zeros in
C previous intervals is odd then there must be a further zero to make
C the total number of zeros even.
          IF (DET.LE.0.D0 .AND. MOD(NZ,2).NE.0) THEN
              NZ = NZ + 1

          ELSE IF (DET.GE.0.D0 .AND. MOD(NZ,2).EQ.0) THEN
              NZ = NZ + 1
          END IF
C If the determinant has changed sign and the number of zeros in
C previous intervals is even then there must be a further zero to make
C the total number of zeros odd.



      END IF


      RETURN

      END


C =====================================================================
C ============== SUBROUTINES FOR EIGENFUNCTION APPROXIMATION ==========
C =====================================================================


      SUBROUTINE EFN(XMESH,P,S,Q,W,ELAM,X,YSTOR1,IMATCH,NMESH,Y,HOT)
C
C This subroutines takes as input the arrays containing the nodal values
C of the eigenfunction, and uses Hermite interpolation to determine the
C values of the eigenfunction and of its quasiderivatives at internodal
C points.
C
C The parameter HOT should be set to .FALSE. on the first call. It will
C be reset to .TRUE. by EFN; on subsequent calls for computation of the
C same eigenfunction, HOT should not be changed by the user. This allows
C the routine to reduce the amount of searching required to find a mesh
C interval containing the point at which the eigenfunction is to be
C evaluated.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,X
      INTEGER IMATCH,NMESH
      LOGICAL HOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION P(1:NMESH),Q(1:NMESH),S(1:NMESH),W(1:NMESH),
     &                 XMESH(0:NMESH),Y(1:4),YSTOR1(4,0:NMESH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPS,LSTVLU,X1,X2,XLOC
      INTEGER I,J,K,L,M
C     ..
C     .. External Subroutines ..
      EXTERNAL AEFP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X02AJF
      EXTERNAL X02AJF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Save statement ..
      SAVE EPS,LSTVLU,M
C     ..

      IF (.NOT.HOT) THEN
          J = 1
          K = NMESH
          L = 1
          LSTVLU = XMESH(0)
          EPS = 100.D0*X02AJF(0.D0)
          HOT = .TRUE.
          GO TO 10
      END IF

      IF (X.GT.LSTVLU) THEN
          J = M
          K = NMESH
          L = 1

      ELSE
          J = M
          K = 1
          L = -1
      END IF

   10 XLOC = MIN(X,XMESH(NMESH))
      XLOC = MAX(XLOC,XMESH(0))

      DO 20 I = J,K,L
          X1 = XMESH(I-1)
          X2 = XMESH(I)
          IF (XLOC.LT.X2 .AND. XLOC.GE.X1) THEN
              IF (ABS(XLOC-X2).LT.EPS) THEN
C Just assign the value at X2
                  Y(4) = YSTOR1(4,I)
                  Y(3) = YSTOR1(3,I)
                  Y(2) = YSTOR1(2,I)
                  Y(1) = YSTOR1(1,I)
              ELSE IF (ABS(XLOC-X1).LT.EPS) THEN
C Just assign the value at X1
                  Y(4) = YSTOR1(4,I-1)
                  Y(3) = YSTOR1(3,I-1)
                  Y(2) = YSTOR1(2,I-1)
                  Y(1) = YSTOR1(1,I-1)
              ELSE
C Do Hermite interpolation
                  CALL AEFP(X1,X2,XLOC,YSTOR1(4,I-1),
     &                      (S(I)*YSTOR1(2,I-1)-YSTOR1(3,I-1)),
     &                      YSTOR1(4,I), (S(I)*YSTOR1(2,I)-YSTOR1(3,I)),
     &                      Y(4))
                  CALL AEFP(X1,X2,XLOC,YSTOR1(1,I-1),YSTOR1(2,I-1),
     &                      YSTOR1(1,I),YSTOR1(2,I),Y(1))
                  CALL AEFP(X1,X2,XLOC,YSTOR1(3,I-1),
     &                      (Q(I)-ELAM*W(I))*YSTOR1(1,I-1),YSTOR1(3,I),
     &                      (Q(I)-ELAM*W(I))*YSTOR1(1,I-1),Y(3))
                  CALL AEFP(X1,X2,XLOC,YSTOR1(2,I-1),
     &                      (YSTOR1(4,I-1)/P(I)),YSTOR1(2,I),
     &                      (YSTOR1(4,I)/P(I)),Y(2))
              END IF
              M = I
              LSTVLU = XLOC
              GO TO 30

          END IF

   20 CONTINUE

   30 CONTINUE

      RETURN

      END


C ---------------------------------------------------------------------


      SUBROUTINE AEFP(XO,XEND,INTP,FXO,FDXO,FXEND,FDXEND,Y)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FDXEND,FDXO,FXEND,FXO,INTP,XEND,XO,Y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H,H2,H3,XA,XA2,XB,XB2
C     ..

      H = XEND - XO
      H2 = H*H
      H3 = H*H2
      XA = INTP - XEND
      XB = INTP - XO
      XA2 = XA*XA
      XB2 = XB*XB

      Y = FXO* ((XA2/H2)+ ((2.D0*XB*XA2)/H3)) +
     &    FXEND* ((XB2/H2)- ((2.D0*XA*XB2)/H3)) + FDXO* (XB*XA2/H2) +
     &    FDXEND* (XA*XB2/H2)


      RETURN

      END


C ----------------------------------------------------------------------


      SUBROUTINE SL4EFN(X,Y1,Y2,HOT,WORK,IWORK,NMESH,IMATCH,ELAM,MULTI,
     &                  WORKS,YSTOR1,YSTOR2,SL4BCS,IFAIL)
C This subroutine computes eigenfunctions of fourth-order
C Sturm-Liouville problems.

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,X
      INTEGER IFAIL,IMATCH,IWORK,MULTI,NMESH
      LOGICAL HOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(0:IWORK,1:6),WORKS(19*NMESH+48),Y1(1:4),
     &                 Y2(1:4),YSTOR1(1:4,0:NMESH),YSTOR2(1:4,0:NMESH)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ARGL,EPS,FAC,NUX,NUX2
      INTEGER FMTCH1,FNMSH1,I,IAL,IAR,IDL1,IDL2,IDR1,IDR2,IFO,IMTCH1,
     &        INDEX,INDEX2,INL,INR,IRL,IRR,IUL,IUR,IV,IVL,IVR,IW,IX,IZL,
     &        IZR,IZVL,IZVR,J,NMSH1,NREC,TMTCH1,TNMSH1
      LOGICAL PRN
      CHARACTER*6 SRNAME
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RINV(2,2),U(2,2),ULOC(2,2),V(2,2),VLOC(2,2),
     &                 YLOC(4)
      CHARACTER*80 REC(2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X02AJF
      INTEGER P01ABF
      EXTERNAL X02AJF,P01ABF
C     ..
C     .. External Subroutines ..
      EXTERNAL F06PAF,GETVEC,THTNT1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,MAX,MIN
C     ..
C     .. Save statement ..
      SAVE IMTCH1,TMTCH1,FMTCH1,NMSH1,TNMSH1,FNMSH1,INL,INR,IRL,IRR,IUL,
     &     IUR,IVL,IVR,IV,IW,IZL,IZR,IZVL,IZVR,IAL,IAR,IDL1,IDL2,IDR1,
     &     IDR2
C     ..
      SRNAME = 'SL4EFN'

      IFO = IFAIL
      IF (ABS(IFO).GT.1) IFO = 0

      IFAIL = 0

      EPS = 1000.D0*X02AJF(0.D0)

      IF (.NOT.HOT) THEN

          IMTCH1 = IMATCH + 1
          TMTCH1 = 2*IMTCH1
          FMTCH1 = 4*IMTCH1
          NMSH1 = NMESH - IMATCH + 1
          TNMSH1 = 2*NMSH1
          FNMSH1 = 4*NMSH1

          INL = 1
          INR = INL + IMTCH1
          IRL = INR + NMSH1
          IRR = IRL + FMTCH1
          IUL = IRR + FNMSH1
          IUR = IUL + FMTCH1
          IVL = IUR + FNMSH1
          IVR = IVL + FMTCH1
          IV = IVR + FNMSH1
          IW = IV + TNMSH1
          IZL = IW + TMTCH1
          IZR = IZL + TMTCH1
          IZVL = IZR + TNMSH1
          IZVR = IZVL + TMTCH1
          IAL = IZVR + TNMSH1
          IAR = IAL + 1
          IDL1 = IAR + 1
          IDL2 = IDL1 + 2
          IDR1 = IDL2 + 2
          IDR2 = IDR1 + 2

C We must start by evaluating the eigenfunction at all the mesh-points.
C NOTE: the unextrapolated value of the eigenvalue is hidden in
C WORK(0,3).
C

          CALL GETVEC(WORK(0,1),NMESH,IMATCH,WORK(1,2),WORK(1,3),
     &                WORK(1,4),WORK(1,5),WORK(0,3),MULTI,WORKS,YSTOR1,
     &                YSTOR2,SL4BCS,IFAIL)
          IF (IFAIL.NE.0) GO TO 250

C We have VSTORE and WSTORE in WORKS(IV..) and WORKS(IW..). If the
C eigenfunction in question has multiplicity greater than 1, then we
C must recover the second copy of VSTOR and WSTOR.

          IF (MULTI.EQ.2) THEN
              WORKS(IZVL+2*IMATCH) = WORKS(IDL1)
              WORKS(IZVL+2*IMATCH+1) = WORKS(IDL1+1)
              DO 10 I = IMATCH,1,-1
                  CALL F06PAF('N',2,2,1.D0,WORKS(IRL+4*I),2,
     &                        WORKS(IZVL+2*I),1,0.D0,
     &                        WORKS(IZVL+2* (I-1)),1)
   10         CONTINUE
              WORKS(IZVR) = WORKS(IDR1)
              WORKS(IZVR+1) = WORKS(IDR1+1)
              DO 20 I = IMATCH,NMESH
                  CALL F06PAF('N',2,2,1.D0,WORKS(IRR+4* (I-IMATCH)),2,
     &                        WORKS(IZVR+2* (I-IMATCH)),1,0.D0,
     &                        WORKS(IZVR+2* (I-IMATCH+1)),1)
   20         CONTINUE
          END IF

      END IF

C We are now ready to get the eigenfunction value at the point X
C specified by the user.


      IF (X.LT.WORK(0,1)-EPS .OR. X.GT.WORK(NMESH,1)+EPS) THEN
          IFAIL = 1
          GO TO 250
      END IF
      X = MAX(X,WORK(0,1))
      X = MIN(X,WORK(NMESH,1))


      IX = 0
   30 IF (WORK(IX+1,1).GE.X .OR. IX+1.GE.NMESH) GO TO 40
      IX = IX + 1
      GO TO 30

C Now we have found an interval [WORK(IX,1),WORK(IX+1,1)] containing
C the point X. We integrate over the interval to find the value of the
C eigenfunction at X from the endpoint values.

   40 IF (X.LT.WORK(IX,1)+1.D2*X02AJF(0.D0)) THEN
          DO 50 I = 1,4
              Y1(I) = YSTOR1(I,IX)
   50     CONTINUE
          IF (MULTI.EQ.2) THEN
              DO 60 I = 1,4
                  Y2(I) = YSTOR2(I,IX)
   60         CONTINUE
          END IF
      ELSE IF (X.GT.WORK(IX+1,1)-1.D2*X02AJF(0.D0)) THEN
          DO 70 I = 1,4
              Y1(I) = YSTOR1(I,IX+1)
   70     CONTINUE
          IF (MULTI.EQ.2) THEN
              DO 80 I = 1,4
                  Y2(I) = YSTOR2(I,IX+1)
   80         CONTINUE
          END IF
      ELSE
C In this case we are not close to one of the endpoints so we must do
C some work. First, we need the integers which indicate the position
C in the storage array WORKS of various essential animals:

          PRN = .FALSE.
          ARGL = 0.D0
          NUX = 0.D0
          NUX2 = 0.D0
          IFAIL = 0

          IF (IX.LT.IMATCH) THEN
C We get our eigenfunction value with a forward integration.
              INDEX = IUL + 4*IX
              INDEX2 = IVL + 4*IX
              DO 100 J = 1,2
                  DO 90 I = 1,2
                      U(I,J) = WORKS(INDEX)
                      V(I,J) = WORKS(INDEX2)
                      INDEX = INDEX + 1
                      INDEX2 = INDEX2 + 1
   90             CONTINUE
  100         CONTINUE

              CALL THTNT1(U,V,RINV,ARGL,NUX,WORK(IX+1,2),WORK(IX+1,3),
     &                    WORK(IX+1,4),WORK(IX+1,5),WORK(0,3),
     &                    WORK(IX,1),X,PRN,IFAIL)
              IF (IFAIL.NE.0) GO TO 250

              DO 120 J = 1,2
                  DO 110 I = 1,2
                      ULOC(I,J) = U(I,J)
                      VLOC(I,J) = V(I,J)
  110             CONTINUE
  120         CONTINUE

              CALL THTNT1(ULOC,VLOC,RINV,ARGL,NUX2,WORK(IX+1,2),
     &                    WORK(IX+1,3),WORK(IX+1,4),WORK(IX+1,5),
     &                    WORK(0,3),X,WORK(IX+1,1),PRN,IFAIL)
              IF (IFAIL.NE.0) GO TO 250

C
C We now have all the information that we need in order to recover
C the eigenfunction. We use the formula
C
C / y(x)  \  = U.RINV.WSTOR_{IX+1}.
C \ y'(x) /             ALFAL.EXP(-NUX2-SUM_{I=IX+2}^{IMATCH}NLSTOR(I))
C
C / -(py''(x))'+sy'(x)  \  = V.RINV.WSTOR_{IX+1}.ALFAL.
C \  py''(x)            /      EXP(-NUX2-SUM_{I=IX+2}^{IMATCH}NLSTOR(I))

C
C valid for x < xmesh(IMATCH). WSTOR is in WORKS(IW+2*(IX+1)). ALFAL is
C in WORKS(IAL). NLSTOR(I) is in WORKS(INL+I).
C
              IF (MULTI.NE.2) THEN

                  CALL F06PAF('N',2,2,1.D0,RINV,2,WORKS(IW+2* (IX+1)),1,
     &                        0.D0,YLOC(1),1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y1(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y1(3),1)
C Now do the scaling.
                  NUX = NUX2
                  IF (IX+1.LT.IMATCH) THEN
                      DO 130 I = IX + 2,IMATCH
                          NUX = NUX + WORKS(INL+I)
  130                 CONTINUE
                  END IF
                  FAC = EXP(-NUX)*WORKS(IAL)
                  DO 140 I = 1,4
                      Y1(I) = Y1(I)*FAC
  140             CONTINUE

              ELSE
C We must also get the second eigenfunction.

                  CALL F06PAF('N',2,2,1.D0,RINV,2,WORKS(IW+2* (IX+1)),1,
     &                        0.D0,YLOC(1),1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y2(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y2(3),1)

                  CALL F06PAF('N',2,2,1.D0,RINV,2,WORKS(IZVL+2* (IX+1)),
     &                        1,0.D0,YLOC(1),1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y1(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y1(3),1)

C Now do the scaling.
                  NUX = NUX2
                  IF (IX+1.LT.IMATCH) THEN
                      DO 150 I = IX + 2,IMATCH
                          NUX = NUX + WORKS(INL+I)
  150                 CONTINUE
                  END IF
                  FAC = EXP(-NUX)*WORKS(IAL)
                  DO 160 I = 1,4
                      Y1(I) = Y1(I)*FAC
                      Y2(I) = Y2(I)*FAC
  160             CONTINUE

              END IF

          ELSE
C We get our eigenfunction value with a backward integration.
              INDEX = IUR + 4* (IX+1-IMATCH)
              INDEX2 = IVR + 4* (IX+1-IMATCH)
              DO 180 J = 1,2
                  DO 170 I = 1,2
                      U(I,J) = WORKS(INDEX)
                      V(I,J) = WORKS(INDEX2)
                      INDEX = INDEX + 1
                      INDEX2 = INDEX2 + 1
  170             CONTINUE
  180         CONTINUE

              CALL THTNT1(U,V,RINV,ARGL,NUX,WORK(IX+1,2),WORK(IX+1,3),
     &                    WORK(IX+1,4),WORK(IX+1,5),WORK(0,3),
     &                    WORK(IX+1,1),X,PRN,IFAIL)
              IF (IFAIL.NE.0) GO TO 250

              DO 200 J = 1,2
                  DO 190 I = 1,2
                      ULOC(I,J) = U(I,J)
                      VLOC(I,J) = V(I,J)
  190             CONTINUE
  200         CONTINUE

              CALL THTNT1(ULOC,VLOC,RINV,ARGL,NUX2,WORK(IX+1,2),
     &                    WORK(IX+1,3),WORK(IX+1,4),WORK(IX+1,5),
     &                    WORK(0,3),X,WORK(IX,1),PRN,IFAIL)
              IF (IFAIL.NE.0) GO TO 250

C
C We now have all the information that we need in order to recover
C the eigenfunction. We use the formula
C
C / y(x)  \  = U.RINV.VSTOR_{IX}.ALFAR.
C \ y'(x) /           EXP(-NUX2-SUM_{I=IMATCH}^{IX-1}NRSTOR(I))
C
C / -(py''(x))'+sy'(x)  \  = V.RINV.VSTOR_{IX}.ALFAR.
C \  py''(x)            /   EXP(-NUX2-SUM_{I=IMATCH}^{IX-1}NRSTOR(I))
C
C valid for x > xmesh(IMATCH). VSTOR is in WORKS(IV+2*(IX-IMATCH)).
C ALFAR is in WORKS(IAR). NRSTOR(I) is in WORKS(INR+I-1).
C
              IF (MULTI.NE.2) THEN

                  CALL F06PAF('N',2,2,1.D0,RINV,2,
     &                        WORKS(IV+2* (IX-IMATCH)),1,0.D0,YLOC(1),1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y1(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y1(3),1)
C Now do the scaling.
                  NUX = NUX2
                  IF (IX.GT.IMATCH) THEN
                      DO 210 I = 0,IX - IMATCH - 1
                          NUX = NUX + WORKS(INR+I)
  210                 CONTINUE
                  END IF
                  FAC = EXP(-NUX)*WORKS(IAR)
                  DO 220 I = 1,4
                      Y1(I) = Y1(I)*FAC
  220             CONTINUE

              ELSE
C We must also get the second eigenfunction.

                  CALL F06PAF('N',2,2,1.D0,RINV,2,
     &                        WORKS(IV+2* (IX-IMATCH)),1,0.D0,YLOC(1),1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y2(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y2(3),1)

                  CALL F06PAF('N',2,2,1.D0,RINV,2,
     &                        WORKS(IZVR+2* (IX-IMATCH)),1,0.D0,YLOC(1),
     &                        1)
                  CALL F06PAF('N',2,2,1.D0,U,2,YLOC(1),1,0.D0,Y1(1),1)
                  CALL F06PAF('N',2,2,1.D0,V,2,YLOC(1),1,0.D0,Y1(3),1)

C Now do the scaling.
                  NUX = NUX2
                  IF (IX.GT.IMATCH) THEN
                      DO 230 I = 0,IX - IMATCH - 1
                          NUX = NUX + WORKS(INR+I)
  230                 CONTINUE
                  END IF
                  FAC = EXP(-NUX)*WORKS(IAR)
                  DO 240 I = 1,4
                      Y1(I) = Y1(I)*FAC
                      Y2(I) = Y2(I)*FAC
  240             CONTINUE

              END IF
          END IF

      END IF

      HOT = .TRUE.

  250 IF (IFAIL.EQ.1) THEN
          WRITE (REC,FMT=9000)
 9000     FORMAT (' ** The X which you specified was out of range')
          NREC = 1
      ELSE IF (ABS(IFAIL).GT.1) THEN
          IFAIL = 2
          WRITE (REC,FMT=9010)
 9010     FORMAT (' ** Some Theta matrices have arisen whose',/,
     &           ' ** eigenvalues cannot be accurately determined')
          NREC = 2
      END IF

      IFAIL = P01ABF(IFO,IFAIL,SRNAME,NREC,REC)

      RETURN
      END

C ----------------------------------------------------------------------

      SUBROUTINE GETVCN(WORK,IWORK,NMESH,IMATCH,ELAM,MULTI,WORKS,YSTOR1,
     &                  YSTOR2,SL4BCS,IFAIL)
C INPUT VARIABLES
C ---------------
C WORK is input from a previous call to SLEUTH, and contains the
C information about the mesh and the coefficient functions PP, SP, QP
C and WP. It must be of length (0:IWORK,1:6), where IWORK is the same
C as on the call to SLEUTH.
C
C NMESH is the value returned by SLEUTH, as is IMATCH.
C ELAM is the eigenvalue approximation returned by SLEUTH
C
C OUTPUT VARIABLES
C ----------------
C MULTI is the multiplicity of ELAM as an eigenvalue of the problem;
C MULTI = 1 or 2.
C
C YSTOR1 is aan array of length (1:4,0:NMESH) such that
C YSTOR1(1,i) contains y(x(i)) for i = 0,..,NMESH
C YSTOR1(2,i) contains y'(x(i)) for i = 0,..,NMESH
C YSTOR1(3,i) contains -(py'')'(x(i))+sy'(x(i)) for i = 0,..,NMESH
C YSTOR1(4,i) contains py''(x(i)) for i = 0,..,NMESH
C
C YSTOR2 contains the same information for the second eigenfunction
C if MULTI = 2.
C
C WORKSPACE
C ---------
C WORKS is an array of length (1:19*NMESH+48). It is used as workspace.
C
C ROUTINES SUPPLIED BY THE USER
C -----------------------------
C SL4BCS is the same routine as was used by SLEUTH.
C
C

C Remark: WORK(0,3) contains, hidden, the last un-extrapolated
C eigenvalue approximation to be calculated.

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM
      INTEGER IFAIL,IMATCH,IWORK,MULTI,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(0:IWORK,1:6),WORKS(19*NMESH+48),
     &                 YSTOR1(1:4,0:NMESH),YSTOR2(1:4,0:NMESH)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C     .. External Subroutines ..
      EXTERNAL GETVEC
C     ..

      CALL GETVEC(WORK(0,1),NMESH,IMATCH,WORK(1,2),WORK(1,3),WORK(1,4),
     &            WORK(1,5),WORK(0,3),MULTI,WORKS,YSTOR1,YSTOR2,SL4BCS,
     &            IFAIL)

      RETURN
      END

C ---------------------------------------------------------------------

      SUBROUTINE GETVEC(XMESH,NMESH,IMATCH,PP,SP,QP,WP,ELAM,MULTI,WORKS,
     &                  YSTOR1,YSTOR2,SL4BCS,IFAIL)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM
      INTEGER IFAIL,IMATCH,MULTI,NMESH
C     ..
C     .. Array Arguments ..
C WORKS must be of length at least 19*NMESH+48
      DOUBLE PRECISION PP(1:NMESH),QP(1:NMESH),SP(1:NMESH),
     &                 WORKS(19*NMESH+48),WP(1:NMESH),XMESH(0:NMESH),
     &                 YSTOR1(4,0:NMESH),YSTOR2(4,0:NMESH)
C     ..
C     .. External Subroutines ..
      EXTERNAL GETVC1
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C IMPORTANT NOTE: The array WORKS must be of length at least
C 19*(NMESH+2)+10

C     .. Local Scalars ..
      INTEGER FMTCH1,FNMSH1,IAL,IAR,IDL1,IDL2,IDR1,IDR2,IMTCH1,INL,INR,
     &        IRL,IRR,IUL,IUR,IV,IVL,IVR,IW,IZL,IZR,IZVL,IZVR,NMSH1,
     &        TMTCH1,TNMSH1
C     ..
      IFAIL = 0

      IMTCH1 = IMATCH + 1
      TMTCH1 = 2*IMTCH1
      FMTCH1 = 4*IMTCH1
      NMSH1 = NMESH - IMATCH + 1
      TNMSH1 = 2*NMSH1
      FNMSH1 = 4*NMSH1

      INL = 1
      INR = INL + IMTCH1
      IRL = INR + NMSH1
      IRR = IRL + FMTCH1
      IUL = IRR + FNMSH1
      IUR = IUL + FMTCH1
      IVL = IUR + FNMSH1
      IVR = IVL + FMTCH1
      IV = IVR + FNMSH1
      IW = IV + TNMSH1
      IZL = IW + TMTCH1
      IZR = IZL + TMTCH1
      IZVL = IZR + TNMSH1
      IZVR = IZVL + TMTCH1
      IAL = IZVR + TNMSH1
      IAR = IAL + 1
      IDL1 = IAR + 1
      IDL2 = IDL1 + 2
      IDR1 = IDL2 + 2
      IDR2 = IDR1 + 2

      CALL GETVC1(XMESH,NMESH,IMATCH,PP,SP,QP,WP,ELAM,MULTI,WORKS(INL),
     &            WORKS(INR),WORKS(IRL),WORKS(IRR),WORKS(IUL),
     &            WORKS(IUR),WORKS(IVL),WORKS(IVR),WORKS(IV),WORKS(IW),
     &            WORKS(IZL),WORKS(IZR),WORKS(IZVL),WORKS(IZVR),YSTOR1,
     &            YSTOR2,WORKS(IAL),WORKS(IAR),WORKS(IDL1),WORKS(IDL2),
     &            WORKS(IDR1),WORKS(IDR2),SL4BCS,IFAIL)

      RETURN

      END

C ---------------------------------------------------------------------


      SUBROUTINE GETVC1(XMESH,NMESH,IMATCH,PP,SP,QP,WP,ELAM,MULTI,
     &                  NLSTOR,NRSTOR,RLSTOR,RRSTOR,ULSTOR,URSTOR,
     &                  VLSTOR,VRSTOR,VSTORE,WSTORE,ZLSTOR,ZRSTOR,
     &                  ZVLSTR,ZVRSTR,YSTOR1,YSTOR2,ALFAL,ALFAR,DL1,DL2,
     &                  DR1,DR2,SL4BCS,IFAIL)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ALFAL,ALFAR,ELAM
      INTEGER IFAIL,IMATCH,MULTI,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DL1(2),DL2(2),DR1(2),DR2(2),NLSTOR(0:IMATCH),
     &                 NRSTOR(IMATCH:NMESH),PP(1:NMESH),QP(1:NMESH),
     &                 RLSTOR(2,2,0:IMATCH),RRSTOR(2,2,IMATCH:NMESH),
     &                 SP(1:NMESH),ULSTOR(2,2,0:IMATCH),
     &                 URSTOR(2,2,IMATCH:NMESH),VLSTOR(2,2,0:IMATCH),
     &                 VRSTOR(2,2,IMATCH:NMESH),VSTORE(2,IMATCH:NMESH),
     &                 WP(1:NMESH),WSTORE(2,0:IMATCH),XMESH(0:NMESH),
     &                 YSTOR1(4,0:NMESH),YSTOR2(4,0:NMESH),
     &                 ZLSTOR(2,0:IMATCH),ZRSTOR(2,IMATCH:NMESH),
     &                 ZVLSTR(2,0:IMATCH),ZVRSTR(2,IMATCH:NMESH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ARGL,ARGR,BETAL,BETAR,DUMMY2,EXCONT,MAXUC,NUX,
     &                 RATIO,XEND,XO
      INTEGER I,J,JJ,K
      LOGICAL ISING,LAST,PRN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(4),RINV(2,2),UC(2,2),UL(2,2),UR(2,2),
     &                 VL(2,2),VR(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X02AJF
      EXTERNAL X02AJF
C     ..
C     .. External Subroutines ..

      EXTERNAL F06PAF,F06QFF,F06QHF,F06YAF,THTNT1
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4BCS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,LOG,MAX,SQRT
C     ..

      IFAIL = 0
      PRN = .FALSE.
C Get the boundary condition at the left-hand endpoint:
      CALL SL4BCS(0,ISING,XMESH(0),UL,VL,ELAM)

      CALL F06QFF('G',2,2,UL,2,ULSTOR(1,1,0),2)
      CALL F06QFF('G',2,2,VL,2,VLSTOR(1,1,0),2)
      CALL F06QHF('G',2,2,0.D0,1.D0,RLSTOR(1,1,0),2)
      NLSTOR(0) = 0.D0

C Now step across the mesh from XMESH(0) to XMESH(IMATCH),
C carrying the matrices UL, VL such that Theta = (VL+iUL).inv(VL-iUL)
C ARGL is not important so we set it (arbitarily) to zero.

      ARGL = 0.D0
      DO 30 I = 1,IMATCH

          XO = XMESH(I-1)
          XEND = XMESH(I)

C The routine THTNT1 does the stepping across the interval
C [XMESH(I-1),XMESH(I)] on which the coefficients P, S, Q and
C W have the constant values PP(I), SP(I), QP(I) and WP(I).
C The eigenparameter has the value ELAM. PRN is a logical
C variable which controls printing of intermediate results
C and IFAIL is an integer which is used to flag errors.


          CALL THTNT1(UL,VL,RINV,ARGL,NUX,PP(I),SP(I),QP(I),WP(I),ELAM,
     &                XO,XEND,PRN,IFAIL)
          IF (IFAIL.NE.0) RETURN


          DO 20 K = 1,2
              DO 10 J = 1,2
                  ULSTOR(J,K,I) = UL(J,K)
                  VLSTOR(J,K,I) = VL(J,K)
                  RLSTOR(J,K,I) = RINV(J,K)
   10         CONTINUE
   20     CONTINUE
          NLSTOR(I) = NUX

   30 CONTINUE

C Get the boundary condition at the right-hand endpoint:

      CALL SL4BCS(1,ISING,XMESH(NMESH),UR,VR,ELAM)
      CALL F06QFF('G',2,2,UR,2,URSTOR(1,1,NMESH),2)
      CALL F06QFF('G',2,2,VR,2,VRSTOR(1,1,NMESH),2)
      CALL F06QHF('G',2,2,0.D0,1.D0,RRSTOR(1,1,NMESH),2)
      NRSTOR(NMESH) = 0.D0
      ARGR = 0.D0

      DO 60 I = NMESH,IMATCH + 1,-1
          XO = XMESH(I)
          XEND = XMESH(I-1)
          CALL THTNT1(UR,VR,RINV,ARGR,NUX,PP(I),SP(I),QP(I),WP(I),ELAM,
     &                XO,XEND,PRN,IFAIL)
          IF (IFAIL.NE.0) RETURN
          DO 50 K = 1,2
              DO 40 J = 1,2
                  URSTOR(J,K,I-1) = UR(J,K)
                  VRSTOR(J,K,I-1) = VR(J,K)
                  RRSTOR(J,K,I-1) = RINV(J,K)
   40         CONTINUE
   50     CONTINUE
          NRSTOR(I-1) = NUX
   60 CONTINUE

C Now match things up in the centre

      CALL F06YAF('T','N',2,2,2,1.d0,VR,2,UL,2,0.D0,UC,2)
      CALL F06YAF('T','N',2,2,2,-1.D0,UR,2,VL,2,1.D0,UC,2)

C Now determine the rank deficiency of UC and its left and
C right null-vectors.

      MAXUC = 0.0D0
      DO 80 I = 1,2
          DO 70 J = 1,2
              MAXUC = MAX(MAXUC,ABS(UC(I,J)))
   70     CONTINUE
   80 CONTINUE

      IF (MAXUC.LT. (10.D0*X02AJF(0.D0))) THEN
          MULTI = 2
          DL1(1) = 1.D0
          DL1(2) = 0.D0
          DL2(1) = 0.D0
          DL2(2) = 1.D0
          DR1(1) = 1.D0
          DR1(2) = 0.D0
          DR2(1) = 0.D0
          DR2(1) = 1.D0

      ELSE IF ((UC(1,2)**2+UC(1,1)**2).GT. (UC(2,2)**2+UC(2,1)**2)) THEN

          MULTI = 1
          DL1(1) = -UC(1,2)
          DL1(2) = UC(1,1)

      ELSE

          MULTI = 1
          DL1(1) = -UC(2,2)
          DL1(2) = UC(2,1)

      END IF

      IF (MAXUC.GT. (10.D0*X02AJF(0.D0))) THEN

          IF ((UC(2,1)**2+UC(1,1)**2).GT. (UC(2,2)**2+UC(1,2)**2)) THEN

              DR1(1) = -UC(2,1)
              DR1(2) = UC(1,1)

          ELSE

              DR1(1) = -UC(2,2)
              DR1(2) = UC(1,2)

          END IF

      END IF


C Now we have DL and DR we can calc Wi and Vi

      LAST = .TRUE.
      IF (MULTI.EQ.2) LAST = .FALSE.

   90 IF (MULTI.EQ.2 .AND. LAST) THEN
C Set up conditions for second eigenfunction
          WSTORE(1,IMATCH) = DL2(1)
          WSTORE(2,IMATCH) = DL2(2)

      ELSE
C Set up conditions for first eigenfunction
          WSTORE(1,IMATCH) = DL1(1)
          WSTORE(2,IMATCH) = DL1(2)
      END IF

      DO 100 I = IMATCH,1,-1
          CALL F06PAF('N',2,2,1.D0,RLSTOR(1,1,I),2,WSTORE(1,I),1,0.D0,
     &                WSTORE(1,I-1),1)
C Scale to avoid possible accumulation to overflow:
          DUMMY2 = SQRT(WSTORE(1,I-1)**2+WSTORE(2,I-1)**2)
          IF (DUMMY2.GT.1.D0/SQRT(X02AJF(0.D0))) THEN
C             DO 101 J = I-1,IMATCH
              J = I - 1
              WSTORE(1,J) = WSTORE(1,J)/DUMMY2
              WSTORE(2,J) = WSTORE(2,J)/DUMMY2
              NLSTOR(I) = NLSTOR(I) - LOG(DUMMY2)
C  101        CONTINUE
          END IF
  100 CONTINUE

      IF (MULTI.EQ.2 .AND. LAST) THEN
C Set up conditions for second eigenfunction
          VSTORE(1,IMATCH) = DR2(1)
          VSTORE(2,IMATCH) = DR2(2)
      ELSE
C Set up conditions for first eigenfunction
          VSTORE(1,IMATCH) = DR1(1)
          VSTORE(2,IMATCH) = DR1(2)
      END IF

      DO 110 I = IMATCH,NMESH - 1
          CALL F06PAF('N',2,2,1.D0,RRSTOR(1,1,I),2,VSTORE(1,I),1,0.D0,
     &                VSTORE(1,I+1),1)
C Scale to avoid possible accumulation to overflow:
          DUMMY2 = SQRT(VSTORE(1,I+1)**2+VSTORE(2,I+1)**2)
          IF (DUMMY2.GT.1.D0/SQRT(X02AJF(0.D0))) THEN
C             DO 111 J = IMATCH,I+1
              J = I + 1
              VSTORE(1,J) = VSTORE(1,J)/DUMMY2
              VSTORE(2,J) = VSTORE(2,J)/DUMMY2
              NRSTOR(J) = NRSTOR(J) - LOG(DUMMY2)
C  111        CONTINUE
          END IF
  110 CONTINUE


C Calculate Zl and Zr

      CALL F06PAF('N',2,2,1.D0,ULSTOR(1,1,IMATCH),2,WSTORE(1,IMATCH),1,
     &            0.D0,ZLSTOR(1,IMATCH),1)
      CALL F06PAF('N',2,2,1.D0,VLSTOR(1,1,IMATCH),2,WSTORE(1,IMATCH),1,
     &            0.D0,ZVLSTR(1,IMATCH),1)

      EXCONT = 0.D0
      DO 120 I = IMATCH - 1,0,-1
          EXCONT = EXCONT + NLSTOR(I+1)
          CALL F06PAF('N',2,2,EXP(-EXCONT),ULSTOR(1,1,I),2,WSTORE(1,I),
     &                1,0.D0,ZLSTOR(1,I),1)
          CALL F06PAF('N',2,2,EXP(-EXCONT),VLSTOR(1,1,I),2,WSTORE(1,I),
     &                1,0.D0,ZVLSTR(1,I),1)
  120 CONTINUE

      CALL F06PAF('N',2,2,1.D0,URSTOR(1,1,IMATCH),2,VSTORE(1,IMATCH),1,
     &            0.D0,ZRSTOR(1,IMATCH),1)
      CALL F06PAF('N',2,2,1.D0,VRSTOR(1,1,IMATCH),2,VSTORE(1,IMATCH),1,
     &            0.D0,ZVRSTR(1,IMATCH),1)

      EXCONT = 0.D0
      DO 130 I = IMATCH + 1,NMESH
          EXCONT = EXCONT + NRSTOR(I-1)
          CALL F06PAF('N',2,2,EXP(-EXCONT),URSTOR(1,1,I),2,VSTORE(1,I),
     &                1,0.D0,ZRSTOR(1,I),1)
          CALL F06PAF('N',2,2,EXP(-EXCONT),VRSTOR(1,1,I),2,VSTORE(1,I),
     &                1,0.D0,ZVRSTR(1,I),1)
C      WRITE (14,FMT=*)
C      WRITE (14,FMT=*) 'URSTOR:'
C      WRITE (14,FMT=*) URSTOR(1,1,I),URSTOR(1,2,I)
C      WRITE (14,FMT=*) URSTOR(2,1,I),URSTOR(2,2,I)
C      WRITE (14,FMT=*) 'VSTORE:'
C      WRITE (14,FMT=*) VSTORE(1,I),VSTORE(2,I)
C      WRITE (14,FMT=*)
C      WRITE (14,FMT=*) 'I,EXCONT,NRSTOR(I-1),ZRSTOR(1,I)',I,EXCONT,
C     & NRSTOR(I-1),ZRSTOR(1,I)
  130 CONTINUE


C Now compute BL,BR

      BETAL = 0.D0
      DO 140 I = 1,IMATCH
          BETAL = BETAL + WP(I)* (((ZLSTOR(1,I)+ZLSTOR(1,I-1))*0.5D0)**
     &            2)* (XMESH(I)-XMESH(I-1))
  140 CONTINUE

      BETAR = 0.D0
      DO 150 I = IMATCH + 1,NMESH
          BETAR = BETAR + WP(I)* (((ZRSTOR(1,I)+ZRSTOR(1,I-1))*0.5D0)**
     &            2)* (XMESH(I)-XMESH(I-1))
C      WRITE (14,FMT=*) 'I,BETAR-so-far:',I,BETAR
C      WRITE (14,FMT=*) 'WP(I),ZRSTOR(1,I),ZRSTOR(1,I-1),H:',
C     & WP(I),ZRSTOR(1,I),ZRSTOR(1,I-1),XMESH(I)-XMESH(I-1)
  150 CONTINUE

C Now calculat alfal and alfar

      DUMMY(1) = ZLSTOR(1,IMATCH)**2 + ZRSTOR(1,IMATCH)**2
      DUMMY(2) = ZLSTOR(2,IMATCH)**2 + ZRSTOR(2,IMATCH)**2
      DUMMY(3) = ZVLSTR(1,IMATCH)**2 + ZVRSTR(1,IMATCH)**2
      DUMMY(4) = ZVLSTR(2,IMATCH)**2 + ZVRSTR(2,IMATCH)**2

      J = 1
      DUMMY2 = DUMMY(J)
      DO 160 JJ = 2,4
          IF (DUMMY(JJ).GT.DUMMY2) THEN
              J = JJ
              DUMMY2 = DUMMY(JJ)
          END IF
  160 CONTINUE

      IF (J.LE.2) THEN
          DUMMY(1) = ZLSTOR(J,IMATCH)
          DUMMY(2) = ZRSTOR(J,IMATCH)
      ELSE
          DUMMY(1) = ZVLSTR(J-2,IMATCH)
          DUMMY(2) = ZVRSTR(J-2,IMATCH)
      END IF

      IF (DUMMY(1).GE.DUMMY(2)) THEN
          RATIO = DUMMY(2)/DUMMY(1)
          ALFAR = 1.D0/SQRT(BETAR+BETAL*RATIO**2)
          ALFAL = ALFAR*RATIO
      ELSE
          RATIO = DUMMY(1)/DUMMY(2)
          ALFAL = 1.D0/SQRT(BETAL+BETAR*RATIO**2)
          ALFAR = ALFAL*RATIO
      END IF

C      WRITE (14,FMT=*) 'ALFAL,ALFAR:',ALFAL,ALFAR
C      WRITE (14,FMT=*) 'DUMMY(1),DUMMY(2),DUMMY(3),DUMMY(4):',DUMMY(1),
C     & DUMMY(2),DUMMY(3),DUMMY(4)
C      WRITE (14,FMT=*) 'J:',J
C      IF (J.LE.2) THEN
C      WRITE (14,FMT=*) 'ZLSTOR(J,IMATCH),ZRSTOR(J,IMATCH):',
C     & ZLSTOR(J,IMATCH),ZRSTOR(J,IMATCH)
C      ELSE
C      WRITE (14,FMT=*) 'ZVLSTR(J-2,IMATCH),ZVRSTR(J-2,IMATCH):',
C     & ZVLSTR(J-2,IMATCH),ZVRSTR(J-2,IMATCH)
C      END IF
C      WRITE (14,FMT=*) 'BETAL,BETAR:',BETAL,BETAR


C Now compute yl(i) and yr(i)

      IF (MULTI.EQ.2 .AND. LAST) THEN
          DO 170 I = 0,IMATCH - 1
              YSTOR2(1,I) = ALFAL*ZLSTOR(1,I)
              YSTOR2(2,I) = ALFAL*ZLSTOR(2,I)
              YSTOR2(3,I) = ALFAL*ZVLSTR(1,I)
              YSTOR2(4,I) = ALFAL*ZVLSTR(2,I)
  170     CONTINUE

          DO 180 I = IMATCH,NMESH
              YSTOR2(1,I) = ALFAR*ZRSTOR(1,I)
              YSTOR2(2,I) = ALFAR*ZRSTOR(2,I)
              YSTOR2(3,I) = ALFAR*ZVRSTR(1,I)
              YSTOR2(4,I) = ALFAR*ZVRSTR(2,I)
  180     CONTINUE

      ELSE
          DO 190 I = 0,IMATCH - 1
              YSTOR1(1,I) = ALFAL*ZLSTOR(1,I)
              YSTOR1(2,I) = ALFAL*ZLSTOR(2,I)
              YSTOR1(3,I) = ALFAL*ZVLSTR(1,I)
              YSTOR1(4,I) = ALFAL*ZVLSTR(2,I)
C              WRITE (14,FMT=*) 'YSTOR1,ZLSTOR:',YSTOR1(1,I),ZLSTOR(1,I)
  190     CONTINUE

          DO 200 I = IMATCH,NMESH
              YSTOR1(1,I) = ALFAR*ZRSTOR(1,I)
              YSTOR1(2,I) = ALFAR*ZRSTOR(2,I)
              YSTOR1(3,I) = ALFAR*ZVRSTR(1,I)
              YSTOR1(4,I) = ALFAR*ZVRSTR(2,I)
C              WRITE (14,FMT=*) 'YSTOR1,ZRSTOR:',YSTOR1(1,I),ZRSTOR(1,I)
  200     CONTINUE

          IF (MULTI.EQ.2) THEN
C LAST must be .FALSE., i.e. we have done only one
C eigenfunction.
              LAST = .TRUE.
              GO TO 90

          END IF

      END IF

      RETURN

      END


C ----------------------------------------------------------------------

      SUBROUTINE THTNT1(U,V,RINV,ARGDET,NUX,P,S,Q,W,ELAM,XO,XEND,PRN,
     &                  IFAIL)

C This subroutine integrates the theta matrix of a 4th order
C Sturm-Liouville problem across an interval on which the
C coefficients are constant, updating argdet(theta) in the
C process. The theta matrix itself is represented by two
C 2x2 matrices u and v such that theta = (v+iu)(v-iu)^{-1}.

C First we update the theta matrix (i.e. the matrices u and
C v) so that we will have all the necessary information for
C updating argdet(theta).

C The formulae which we use for doing the update are as follows.
C If xo < xend, then
C
C Argdet(Theta(xend)) = Argdet(Theta(xo))
C          + (Sum of Phase Angles of Theta(xend) in [0,2pi))
C          - (Sum of Phase Angles of Theta_{0}^{*}Theta(xo) in [0,2pi))
C          - (Sum of Phase Angles of Theta_{0}(xo) in [0,2pi))
C          + 2(n+2)pi,
C
C where n is the number of zeros of det(Uo) in [xo,xend) (remember
C that Theta_{0}(xend) = I).
C If xo > xend, then
C
C Argdet(Theta(xend)) = Argdet(Theta(xo))
C          - (Sum of Phase Angles of Theta^{*}(xend) in [0,2pi))
C          + (Sum of Phase Angles of Theta^{*}Theta_{0}(xo) in [0,2pi))
C          - (Sum of Phase Angles of Theta_{0}(xo) in (0,2pi])
C          - 2.n.pi
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ARGDET,ELAM,NUX,P,Q,S,W,XEND,XO
      INTEGER IFAIL
      LOGICAL PRN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RINV(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,ARG0,ARGEND,ARGTHO,B,C,D,DISC1,DISC2,DUMMY,EC,
     &                 M,OMEGA1,OMEGA2,SCAL,TWOPI
      INTEGER I,J,N
      LOGICAL LINTSP,LNEG,LPOS,SPCIAL
C     ..
C     .. Local Arrays ..

      DOUBLE PRECISION CU(2,2),CV(2,2),RNORM(2,2),SU(2,2),SV(2,2),
     &                 ULOC(2,2),UO(2,2),UTEMP(2,2),VLOC(2,2),VO(2,2),
     &                 VTEMP(2,2),Y1(4),Y2(4),Y3(4),Y4(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION X01AAF
      EXTERNAL X01AAF
C     ..
C     .. External Subroutines ..
      EXTERNAL F06YAF,FUNSOL,ITNLM1,LRGLM1,NEGLM1,SL4CNT,SPTH4
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC ABS,ATAN2,COS,DBLE,LOG,MAX,SIN,SQRT
C     ..
      TWOPI = 2.D0*X01AAF(0.D0)



C Copy the matrices into temporary storage for later
C use.
      DO 20 J = 1,2
          DO 10 I = 1,2
              UTEMP(I,J) = U(I,J)
              VTEMP(I,J) = V(I,J)
   10     CONTINUE
   20 CONTINUE
C

      DISC1 = 4.D0*P* (ELAM*W-Q)
      DISC2 = S*S + DISC1

      LPOS = DISC1 .GE. 0.D0
      LINTSP = DISC2 .GE. 0.D0 .AND. S - SQRT(ABS(DISC2)) .GE. 0.D0
      LNEG = DISC2 .LT. 0.D0

      IF (LPOS .OR. LINTSP) THEN
          SCAL = ABS(XEND-XO)*SQRT(0.5D0* (S+SQRT(ABS(DISC2)))/P)

      ELSE IF (LNEG) THEN
          SCAL = ABS(XEND-XO)*0.5D0*SQRT((S+SQRT(ABS(DISC1)))/P)

      ELSE
          SCAL = 0.D0
      END IF

      SPCIAL = SCAL .GT. 1.D-1

      IF (.NOT.SPCIAL) THEN
C In this case we must accept the output from FUNSOL as correct

C          WRITE (8,FMT=*) 'Doing simple case'
          CALL FUNSOL(XO,XEND,P,S,Q,W,ELAM,Y1,Y2,Y3,Y4,SCAL)

          CU(1,1) = Y1(1)
          CU(1,2) = Y2(1)
          CU(2,1) = Y1(2)
          CU(2,2) = Y2(2)

          SU(1,1) = Y3(1)
          SU(1,2) = Y4(1)
          SU(2,1) = Y3(2)
          SU(2,2) = Y4(2)

          SV(1,1) = Y1(3)
          SV(1,2) = Y2(3)
          SV(2,1) = Y1(4)
          SV(2,2) = Y2(4)

          CV(1,1) = Y3(3)
          CV(1,2) = Y4(3)
          CV(2,1) = Y3(4)
          CV(2,2) = Y4(4)
C
C Now do the update of Theta from xo to xend, using the formulae
C               u(xend) = cu.u(xo)+su.v(xo),
C               v(xend) = sv.u(xo)+cv.v(xo).
C
          CALL F06YAF('N','N',2,2,2,1.D0,CU,2,UTEMP,2,0.D0,UO,2)
          CALL F06YAF('N','N',2,2,2,1.D0,SU,2,VTEMP,2,1.D0,UO,2)

          CALL F06YAF('N','N',2,2,2,1.D0,CV,2,VTEMP,2,0.D0,VO,2)
          CALL F06YAF('N','N',2,2,2,1.D0,SV,2,UTEMP,2,1.D0,VO,2)
C
C Now do some postmultiplication to avoid ill-conditioning
C of these matrices.
C
          B = -UO(1,1)*VO(2,2) + UO(1,2)*VO(2,1) - UO(2,2)*VO(1,1) +
     &        UO(2,1)*VO(1,2)
          AC = UO(1,1)*UO(2,2) - UO(2,1)*UO(1,2) - VO(1,1)*VO(2,2) +
     &         VO(2,1)*VO(1,2)
          C = UO(1,1)*UO(2,2) - UO(2,1)*UO(1,2) + VO(1,1)*VO(2,2) -
     &        VO(2,1)*VO(1,2)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

          IF (C.GE.0.D0) THEN
              M = ATAN2(B,AC)

          ELSE
              M = ATAN2(-B,-AC)
          END IF

          C = COS(M/2.D0)
          AC = SIN(M/2.D0)

C Form the normalising matrix which is equal to the return velaue of R:

          DO 40 J = 1,2
              DO 30 I = 1,2
                  RNORM(I,J) = UO(I,J)*C - VO(I,J)*AC
   30         CONTINUE
   40     CONTINUE

C In this case a uniform scaling is used by the routine FUNSOL, so

          NUX = SCAL

C Form the inverse of the normalising matrix:

          D = RNORM(1,1)*RNORM(2,2) - RNORM(1,2)*RNORM(2,1)
          EC = RNORM(1,1)/D
          RINV(1,1) = RNORM(2,2)/D
          RINV(2,2) = EC
          RINV(1,2) = -RNORM(1,2)/D
          RINV(2,1) = -RNORM(2,1)/D

C
C Peform the postmultiplication:

          CALL F06YAF('N','N',2,2,2,1.D0,UO,2,RINV,2,0.D0,U,2)
          CALL F06YAF('N','N',2,2,2,1.D0,VO,2,RINV,2,0.D0,V,2)

C U and V have now been postmultiplied and V+iU should be
C reasonably well-conditioned.

C Now set the U and V for Theta_{0}
          UO(1,1) = -SU(1,1)
          UO(1,2) = SU(1,2)
          UO(2,1) = SU(2,1)
          UO(2,2) = -SU(2,2)

          VO(1,1) = CV(1,1)
          VO(1,2) = -CV(1,2)
          VO(2,1) = -CV(2,1)
          VO(2,2) = CV(2,2)
C          WRITE (8,FMT=*) XO,XEND,'Used FUNSOL'
      ELSE
          IF (LPOS) THEN
              CALL LRGLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,UO,VO,XO,XEND)
C              WRITE (8,FMT=*) XO,XEND,'Used LRGLM1'
          ELSE IF (LINTSP) THEN
              CALL ITNLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,UO,VO,XO,XEND)
C              WRITE (8,FMT=*) XO,XEND,'Used ITNLM1'
          ELSE IF (LNEG) THEN
              CALL NEGLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,UO,VO,XO,XEND)
C              WRITE (8,FMT=*) XO,XEND,'Used NEGLM1'
          END IF

      END IF
C Normalise U and V so that the greatest element of the two is
C equal to 1.
      DUMMY = 0.D0
      DO 60 J = 1,2
          DO 50 I = 1,2
              DUMMY = MAX(DUMMY,ABS(U(I,J)),ABS(V(I,J)))
   50     CONTINUE
   60 CONTINUE

      DO 80 J = 1,2
          DO 70 I = 1,2
              U(I,J) = U(I,J)/DUMMY
              V(I,J) = V(I,J)/DUMMY
   70     CONTINUE
   80 CONTINUE

C Include this into the value of NUX
      NUX = NUX + LOG(DUMMY)

C Now all we have to do is get the eigenvalues of lots of
C different matrices.

C Get the eigenvalues of Theta_{0}^{*}Theta(xo).

      CALL F06YAF('T','N',2,2,2,1.D0,VO,2,UTEMP,2,0.D0,ULOC,2)
      CALL F06YAF('T','N',2,2,2,-1.D0,UO,2,VTEMP,2,1.D0,ULOC,2)

      CALL F06YAF('T','N',2,2,2,1.D0,VO,2,VTEMP,2,0.D0,VLOC,2)
      CALL F06YAF('T','N',2,2,2,1.D0,UO,2,UTEMP,2,1.D0,VLOC,2)

      CALL SPTH4(ULOC,VLOC,OMEGA1,OMEGA2,ARG0,'L',IFAIL)
      IF (IFAIL.NE.0) THEN
          IFAIL = 5
          RETURN

      END IF

      IF (XO.GT.XEND) THEN
C We actually want the phase-angles of the hermitian conjugate
C of Theta_{0}^{*}Theta(xo)
          IF (OMEGA1.NE.0.D0) OMEGA1 = TWOPI - OMEGA1
          IF (OMEGA2.NE.0.D0) OMEGA2 = TWOPI - OMEGA2
          ARG0 = -OMEGA1 - OMEGA2
      END IF

C Get the phase angles of Theta(xend)

      CALL SPTH4(U,V,OMEGA1,OMEGA2,ARGEND,'L',IFAIL)
      IF (IFAIL.NE.0) THEN
          IFAIL = 6
          RETURN

      END IF

      IF (XO.GT.XEND) THEN
C We actually want the phase-angles of the hermitian conjugate
C of Theta(xend)
          IF (OMEGA1.NE.0.D0) OMEGA1 = TWOPI - OMEGA1
          IF (OMEGA2.NE.0.D0) OMEGA2 = TWOPI - OMEGA2
          ARGEND = -OMEGA1 - OMEGA2
      END IF

C Get phase-angles of Theta_{0}(xo), normalised to lie in (0,2pi]

      CALL SPTH4(UO,VO,OMEGA1,OMEGA2,ARGTHO,'R',IFAIL)
      IF (IFAIL.NE.0) THEN
          IFAIL = 7
          RETURN

      END IF

C Do zero count:
      CALL SL4CNT(XEND-XO,P,S,Q,W,ELAM,N)

C Do the update of ARGDET:

      IF (XO.LE.XEND) THEN
          ARGDET = ARGDET + ARGEND - ARG0 + DBLE(N+2)*TWOPI - ARGTHO
      ELSE
          ARGDET = ARGDET + ARGEND - ARG0 - DBLE(N)*TWOPI - ARGTHO
      END IF

      END

C --------------------------------------------------------------------

      SUBROUTINE LRGLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,NUX,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),RINV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,AOMX,B,BETA,C,COSHCS,COSHNX,COSHSN,COSOMX,
     &                 DETCU,DETCV,DETR,DETSU,DETSV,DETU,DETV,DISC,FAC,
     &                 H,H2,H2SNSN,M,NO,NO2,NU,NU2,NU3,NU4,NU6,OMEGA,
     &                 OMEGA2,OMEGA3,OMEGA4,OMEGA6,OMX,PNO,PNO2,PNORAD,
     &                 PRAD,PRAD2,RAD,RAD2,SINHCS,SINHNX,SINHSN,SINOMX,
     &                 TRUVA,ZCOSMT,ZSINMT
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(2,2),CU(2,2),CUSVA(2,2),R(2,2),SUCVA(2,2),
     &                 SV(2,2),T(2,2),TEMP(2,2),UT(2,2),UUA(2,2),
     &                 UVA(2,2),VT(2,2),VUA(2,2),VVA(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION PHI
      EXTERNAL PHI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,SIN,SQRT
C     ..
C     .. External Subroutines ..

      EXTERNAL F06YAF
C     ..
      DISC = SQRT(S*S+4.D0*P* (ELAM*W-Q))
      NU = SQRT(0.5d0* (S+DISC)/P)
      OMEGA = SQRT(0.5d0* (-S+DISC)/P)



      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)
      H = XEND - XO
      H2 = H*H
      NUX = ABS(NU*H)
      FAC = EXP(-NUX)
      IF (NUX.GE.0.1D0) THEN
          SINHNX = (1.D0-FAC*FAC)/ (2.D0*NUX)

      ELSE
          SINHNX = PHI(NUX*NUX)*FAC
      END IF

      COSHNX = (1.D0+FAC*FAC)/2.D0
      OMX = OMEGA*H
      AOMX = ABS(OMX)
      IF (AOMX.GE.0.1D0) THEN
          SINOMX = SIN(OMX)/OMX

      ELSE
          SINOMX = PHI(-OMX*OMX)
      END IF

      COSOMX = COS(OMX)

      COSHSN = H*COSHNX*SINOMX
      SINHCS = H*SINHNX*COSOMX
      COSHCS = COSHNX*COSOMX
      SINHSN = SINHNX*SINOMX
      H2SNSN = H2*SINHSN
      OMEGA2 = OMEGA*OMEGA
      NU2 = NU*NU
      OMEGA3 = OMEGA*OMEGA2
      NU3 = NU*NU2
      OMEGA4 = OMEGA2*OMEGA2
      NU4 = NU2*NU2
      OMEGA6 = OMEGA2*OMEGA4
      NU6 = NU2*NU4
      RAD = OMEGA2 + NU2
      RAD2 = RAD*RAD
      PRAD = P*RAD
      PRAD2 = PRAD*PRAD
      NO = NU*OMEGA
      PNO = P*NO
      PNORAD = PNO*RAD
      NO2 = NO*NO
      PNO2 = P*NO2


C Compute the CU,CV,SU,SV matrices here

      CU(1,1) = (NU2*COSOMX*FAC+OMEGA2*COSHNX)/RAD
      CU(1,2) = H* (NU2*SINHNX+OMEGA2*SINOMX*FAC)/RAD
      CU(2,1) = H* (-NU2*OMEGA2*SINOMX*FAC+OMEGA2*NU2*SINHNX)/RAD
      CU(2,2) = (NU2*COSHNX+OMEGA2*COSOMX*FAC)/RAD

      SU(1,1) = H* (SINOMX*FAC-SINHNX)/PRAD
      SU(1,2) = (COSHNX-COSOMX*FAC)/PRAD
      SU(2,1) = -SU(1,2)
      SU(2,2) = H* (NU2*SINHNX+OMEGA2*SINOMX*FAC)/PRAD

      CV(1,1) = (NU2*COSOMX*FAC+OMEGA2*COSHNX)/RAD
      CV(1,2) = H* (NO* (NO*SINOMX*FAC-NO*SINHNX))/RAD
      CV(2,1) = H* (-OMEGA2*SINOMX*FAC-NU2*SINHNX)/RAD
      CV(2,2) = (NU2*COSHNX+OMEGA2*COSOMX*FAC)/RAD

      SV(1,1) = H* (-PNO* (NU2*NO*SINOMX*FAC+OMEGA2*NO*SINHNX))/RAD
      SV(1,2) = (-P*NU2*OMEGA2* (COSHNX-COSOMX*FAC))/RAD
      SV(2,1) = -SV(1,2)
      SV(2,2) = H* (P* (NU4*SINHNX-OMEGA4*SINOMX*FAC))/RAD

C Now we can compute UT and VT


      CALL F06YAF('N','N',2,2,2,1.D0,CU,2,U,2,0.D0,UT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,SU,2,V,2,1.D0,UT,2)

      CALL F06YAF('N','N',2,2,2,1.D0,SV,2,U,2,0.D0,VT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,CV,2,V,2,1.D0,VT,2)



C Compute SU.Adjoint(CV).exp(-NUX):

      SUCVA(1,1) = (COSHSN-SINHCS)/PRAD
      SUCVA(1,2) = ((OMEGA2-NU2)* (FAC-COSHCS)+2.D0*NO2*H2SNSN)/
     &             (PRAD*RAD)
      SUCVA(2,1) = SUCVA(1,2)
      SUCVA(2,2) = (NU2*SINHCS+OMEGA2*COSHSN)/PRAD

C Compute CU.Adjoint(SV).exp(-NUX):

      CUSVA(1,1) = P* (NU4*SINHCS-OMEGA4*COSHSN)/RAD
      CUSVA(1,2) = -PNO2* ((NU4+OMEGA4)*H2SNSN+
     &             (NU2-OMEGA2)* (FAC-COSHCS))/RAD2
      CUSVA(2,1) = CUSVA(1,2)
      CUSVA(2,2) = -PNO2* (NU2*COSHSN+OMEGA2*SINHCS)/RAD

C Compute A = U.Adjoint(V):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = V(1,1)*U(2,2) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.A.Adjoint(SV))*exp(-NUX):

      T(1,1) = A(1,1)*COSHCS + 2.D0*A(1,2)* (OMEGA2*COSHSN+NU2*SINHCS)/
     &         RAD + A(2,2)*H2SNSN

      T(1,2) = (A(1,1)*NO2* (SINHCS-COSHSN)+
     &         A(1,2)* ((OMEGA2-NU2)* (FAC* (OMEGA2-NU2)+
     &         2.D0*NO2*H2SNSN)+4.D0*NO2*COSHCS)/RAD+
     &         A(2,2)* (NU2*SINHCS+OMEGA2*COSHSN))/RAD

      T(2,1) = T(1,2)


      T(2,2) = -A(1,1)*NO2*H2SNSN + 2.D0*A(1,2)*NO2* (SINHCS-COSHSN)/
     &         RAD + A(2,2)*COSHCS

C Compute U.Adjoint(V) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      DETCU = ((NU4+OMEGA4)*COSHCS+NO2* ((NU2-OMEGA2)*H2SNSN+2.D0*FAC))/
     &        RAD2
      DETSU = ((NU2-OMEGA2)*H2SNSN+2.D0* (FAC-COSHCS))/PRAD2
      DETCV = DETCU
      DETSV = P*PNO* ((OMEGA6-NU6)*NO*H2SNSN+
     &        2.D0*NU3*OMEGA3* (FAC-COSHCS))/RAD2

C Compute trace(CU.A.Adjoint(SU)).exp(-NUX):

      BETA = A(1,1)* (OMEGA2*COSHSN+NU2*SINHCS)/PRAD -
     &       2.D0*A(1,2)* ((NU2-OMEGA2)* (FAC-COSHCS)-2.D0*NO2*H2SNSN)/
     &       (P*RAD2) + A(2,2)* (COSHSN-SINHCS)/PRAD

      UUA(1,1) = DETU*DETCU + DETV*DETSU + BETA
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)

C Compute trace(CV.A.Adjoint(SV)).exp(-NUX):

      BETA = (-A(1,1)*PNO2* (NU2*COSHSN+OMEGA2*SINHCS)-
     &       2.D0*PNO2*A(1,2)* ((NU2-OMEGA2)* (FAC-COSHCS)+ (NU4+
     &       OMEGA4)*H2SNSN)/RAD-A(2,2)*P* (OMEGA4*COSHSN-NU4*SINHCS))/
     &       RAD

      VVA(1,1) = DETU*DETSV + DETV*DETCV + BETA
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF


C Now we have got the value of M we can compute R.

      ZCOSMT = COS(M/2.D0)
      ZSINMT = SIN(M/2.D0)


      DO 20 J = 1,2
          DO 10 I = 1,2
              R(I,J) = UT(I,J)*ZCOSMT - VT(I,J)*ZSINMT
   10     CONTINUE
   20 CONTINUE

C Compute the DET R and RINV

      DETR = ZCOSMT**2*UUA(1,1) + ZSINMT**2*VVA(1,1) -
     &       ZSINMT*ZCOSMT* (UVA(1,1)+UVA(2,2))


      RINV(1,1) = R(2,2)/DETR
      RINV(1,2) = -R(1,2)/DETR
      RINV(2,1) = -R(2,1)/DETR
      RINV(2,2) = R(1,1)/DETR

      CALL F06YAF('N','N',2,2,2,1.D0,UT,2,RINV,2,0.D0,TEMP,2)

      DO 40 J = 1,2
          DO 30 I = 1,2
              U(I,J) = (UUA(I,J)*ZCOSMT-UVA(I,J)*ZSINMT)/DETR
              V(I,J) = (VUA(I,J)*ZCOSMT-VVA(I,J)*ZSINMT)/DETR
   30     CONTINUE
   40 CONTINUE

      NUX = 0.D0

C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

C NOTE: in the following formulae we want CV and SU for BACKWARDS
C integration. We exploit the fact that the diagonal terms of
C SUCVA are odd functions of the direction of integration while
C the off-diagonal terms are even functions of the direction
C of integration, while the determinants which appear in the
C formulae are even functions of the direction of integration.

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END

C -------------------------------------------------------------------

      SUBROUTINE ITNLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,NUX,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),RINV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..


      DOUBLE PRECISION AC,ADJUST,AEHH,AETAH,ARHH,B,BETA,C,CEHCHH,CEHSHH,
     &                 CMSN,CNCM,CNSM,COSHEH,COSHHH,COSHMX,COSHNX,DETCU,
     &                 DETCV,DETR,DETSU,DETSV,DETU,DETV,DISC,ETA,ETA2,
     &                 ETAH,ETASQ,ETSMRS,ETSRHS,FAC,FACL,H,H2,M,MU,MU2,
     &                 MU3,MU4,MU6,MUX,NM,NU,NU2,NU3,NU4,NU6,PNM,PNMR,
     &                 PRAD,PRAD2,RAD,RAD2,RADH2,RH,RH2,RHH,RHSQ,SEHCHH,
     &                 SEHSHH,SINHMX,SINHNX,SN2ERT,SN2HRT,SNERAT,SNERNS,
     &                 SNHR1,SNHR2,SNHRAT,SNHRNS,SNSM,TRUVA,ZCOSMT,
     &                 ZSINMT
      INTEGER I,J
      LOGICAL SMALL
C     ..
C     .. Local Arrays ..

C      DOUBLE PRECISION CUT(2,2),CVT(2,2),SUT(2,2),SVT(2,2)
      DOUBLE PRECISION A(2,2),CU(2,2),CUSVA(2,2),R(2,2),SUCVA(2,2),
     &                 SV(2,2),T(2,2),UT(2,2),UUA(2,2),UVA(2,2),VT(2,2),
     &                 VUA(2,2),VVA(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION PHI,SNHDIF
      EXTERNAL PHI,SNHDIF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,MIN,SIGN,SIN,SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL F06YAF
C     ..

C      WRITE (8,FMT=*) 'Debugging output from ITNLM1'

      DISC = SQRT(S*S+4.D0*P* (ELAM*W-Q))
      NU2 = 0.5D0* (S+DISC)/P
      MU2 = NU2 - DISC/P
      NU = SQRT(NU2)
      MU = SQRT(MU2)
      NU3 = NU*NU2
      MU3 = MU*MU2
      NU4 = NU2*NU2
      MU4 = MU2*MU2
      NU6 = NU3*NU3
      MU6 = MU3*MU3

      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)

      H = XEND - XO
      H2 = H*H
      NUX = ABS(NU*H)
      SINHNX = (1.D0-EXP(-2.D0*NUX))/2.D0
      COSHNX = 1.D0 - SINHNX
      SINHNX = SINHNX*SIGN(1.D0,H)

      MUX = ABS(MU*H)
      SINHMX = (1.D0-EXP(-2.D0*MUX))/2.D0
      COSHMX = 1.D0 - SINHMX
      SINHMX = SINHMX*SIGN(1.D0,H)

      FAC = EXP(-NUX-MUX)

      CMSN = COSHMX*SINHNX
      CNSM = COSHNX*SINHMX
      CNCM = COSHNX*COSHMX
      SNSM = SINHNX*SINHMX

C If any of the following quantities is small then we are
C in trouble; Taylor series expansions have to be used.

C      RAD = NU**2-MU**2
      RAD = DISC/P
      RAD2 = RAD*RAD
      PRAD = DISC
      PRAD2 = DISC*DISC
      NM = NU*MU
      PNM = P*NM
      PNMR = PNM*RAD

      SMALL = H2*MIN(ABS(RAD),ABS(PRAD),ABS(PNMR),PNM) .LT. 1.D-6

C      WRITE (8,FMT=*) 'SMALL is ',SMALL


C Compute SU.Adjoint(CV).exp(-(NUX+MUX)):

      IF (SMALL) THEN
C Compute certain quantities which will be jolly useful
C when calculating the matrix elements
          ETA2 = NU + MU
          ETA = 0.5D0*ETA2
          ETASQ = ETA*ETA
          ETAH = ETA*H
          AETAH = ABS(ETAH)
          RH2 = NU - MU
          RH = 0.5D0*RH2
          RHSQ = RH*RH
          RHH = RH*H
          ARHH = ABS(RHH)

          COSHEH = 0.5D0* (1.D0+EXP(-2.D0*AETAH))
          COSHHH = 0.5D0* (1.D0+EXP(-2.D0*ARHH))
          CEHCHH = COSHEH*COSHHH

          IF (ARHH.LT.0.1D0) THEN
              SNHRAT = PHI(RHH*RHH)
              SNHRNS = SNHRAT*EXP(-ARHH)
              SNHRAT = SNHRAT*SNHRAT*FAC
              SN2HRT = PHI(4.D0*RHH*RHH)*FAC
          ELSE
              FACL = EXP(-MUX)
              SNHRAT = 0.5D0* (1.D0-EXP(-2.D0*ARHH))/ARHH
              SNHRNS = SNHRAT
              SNHRAT = SNHRAT*FACL
              SNHRAT = SNHRAT*SNHRAT
              SN2HRT = (1.D0-EXP(-4.D0*ARHH))/ (4.D0*ARHH)
              SN2HRT = SN2HRT*FACL*FACL
          END IF

          IF (AETAH.LT.0.1D0) THEN
              SNERAT = PHI(ETAH*ETAH)
              SNERNS = SNERAT*EXP(-AETAH)
              SNERAT = SNERAT*SNERAT*FAC
              SN2ERT = PHI(4.D0*ETAH*ETAH)*FAC

          ELSE
              SNERAT = 0.5D0* (1.D0-EXP(-2.D0*AETAH))/AETAH
              SNERNS = SNERAT
              SNERAT = SNERAT*SNERAT
              SN2ERT = (1.D0-EXP(-4.D0*AETAH))/ (4.D0*AETAH)
          END IF
C          WRITE (8,FMT=*) 'Check COSHEH:',COSH(ETAH)*EXP(-AETAH)-COSHEH
C          WRITE (8,FMT=*) 'Check COSHRH:',COSH(RHH)*EXP(-ARHH)-COSHHH
C          WRITE (8,FMT=*) 'Check SNHRNS:',SINH(RHH)*EXP(-ARHH)-
C     &                     SNHRNS*RHH
C          WRITE (8,FMT=*) 'Check SNERNS:',SINH(ETAH)*EXP(-ARHH)-
C     &                     SNERNS*ETAH

          CEHSHH = COSHEH*SNHRNS*H
          SEHSHH = SNERNS*SNHRNS*H2
          SEHCHH = COSHHH*SNERNS*H

          IF (MUX.LT.0.1D0) THEN
              SNHR1 = PHI(MUX*MUX)*EXP(-MUX)
          ELSE
              SNHR1 = SIGN(1.D0,H)*SINHMX/MUX
          END IF

          IF (NUX.LT.0.1D0) THEN
              SNHR2 = PHI(NUX*NUX)*EXP(-NUX)
          ELSE
              SNHR2 = SIGN(1.D0,H)*SINHNX/NUX
          END IF

      END IF

C Here we can compute the SU,SV,CU,CV matrices

      ADJUST = EXP(-NUX+MUX)

      IF (.NOT.SMALL) THEN
          CU(1,1) = (NU2*COSHMX*ADJUST-MU2*COSHNX)/RAD
          CU(1,2) = (NU*SINHNX-MU*SINHMX*ADJUST)/RAD
          CU(2,1) = NM* (NU*SINHMX*ADJUST-MU*SINHNX)/RAD
          CU(2,2) = (NU2*COSHNX-MU2*COSHMX*ADJUST)/RAD

          SU(1,1) = (((SINHMX)/MU)*ADJUST- ((SINHNX)/NU))/PRAD
          SU(1,2) = (COSHNX-COSHMX*ADJUST)/PRAD
          SU(2,1) = -SU(1,2)
          SU(2,2) = (NU*SINHNX-MU*SINHMX*ADJUST)/PRAD

          SV(1,1) = PNM* (NU3*SINHMX*ADJUST-MU3*SINHNX)/RAD
          SV(1,2) = PNM*NM* (COSHNX-COSHMX*ADJUST)/RAD
          SV(2,1) = -SV(1,2)
          SV(2,2) = (P* (NU3*SINHNX-MU3*SINHMX*ADJUST))/RAD

          CV(1,1) = (NU2*COSHMX*ADJUST-MU2*COSHNX)/RAD
          CV(1,2) = -NM* (NU*SINHMX*ADJUST-MU*SINHNX)/RAD
          CV(2,1) = - (NU*SINHNX-MU*SINHMX*ADJUST)/RAD
          CV(2,2) = (NU2*COSHNX-MU2*COSHMX*ADJUST)/RAD

      ELSE

          ETSRHS = ETASQ + RHSQ
          ETSMRS = ETASQ - RHSQ

          CU(1,1) = CEHCHH - 0.5D0*ETSRHS*SEHSHH
          CU(1,2) = 0.5D0* (CEHSHH+SEHCHH)
          CU(2,1) = 0.5D0*ETSMRS* (SEHCHH-CEHSHH)
          CU(2,2) = CEHCHH + 0.5D0*ETSRHS*SEHSHH

          SU(1,2) = (0.5D0*SEHSHH)/P
          SU(2,1) = -SU(1,2)
          SU(2,2) = (0.5D0* (CEHSHH+SEHCHH))/P

          CV(1,1) = CU(1,1)
          CV(1,2) = -CU(2,1)
          CV(2,1) = -CU(1,2)
          CV(2,2) = CU(2,2)

          SV(1,1) = 0.5D0*PNM* ((3.D0*ETASQ+RHSQ)*SEHCHH-
     &              (ETASQ+3.D0*RHSQ)*CEHSHH)
          SV(1,2) = 0.5D0*PNM*NM*SEHSHH
          SV(2,1) = -SV(1,2)
          SV(2,2) = 0.5D0*P* ((ETASQ+3.0D0*RHSQ)*CEHSHH+
     &              (RHSQ+3.0D0*ETASQ)*SEHCHH)

          RADH2 = ABS(RAD)*H2
          AEHH = ABS(ETASQ-RHSQ)*H2

          IF (RADH2.GE.0.01D0 .AND. AEHH.LT.0.01D0) THEN
C              WRITE (8,FMT=*) 'Case 1 computation of SU(1,1)'
              SU(1,1) = (1/PRAD)* (SNHR1*ADJUST-SNHR2)*H

          ELSE IF (RADH2.LT.0.01D0 .AND. AEHH.GE.0.01D0) THEN
C              WRITE (8,FMT=*) 'Case 2 computation of SU(1,1)'

              SU(1,1) = H* (SEHCHH-CEHSHH)/ (2.D0*PNM)

          ELSE
C              WRITE (8,FMT=*) 'Case 3 computation of SU(1,1)'

              SU(1,1) = ((H**3)/P)*SNHDIF(NUX,MUX)*EXP(-NUX)

          END IF


      END IF

C Now we can compute UT and VT


      CALL F06YAF('N','N',2,2,2,1.D0,CU,2,U,2,0.D0,UT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,SU,2,V,2,1.D0,UT,2)

      CALL F06YAF('N','N',2,2,2,1.D0,SV,2,U,2,0.D0,VT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,CV,2,V,2,1.D0,VT,2)


      IF (.NOT.SMALL) THEN
          SUCVA(1,1) = (NU*CNSM-MU*CMSN)/PNMR
          SUCVA(1,2) = - (S* (FAC-CNCM)+2.D0*PNM*SNSM)/PRAD2
          SUCVA(2,1) = SUCVA(1,2)
          SUCVA(2,2) = (NU*CMSN-MU*CNSM)/PRAD

      ELSE
          RADH2 = ABS(RAD)*H2
          AEHH = ABS(ETASQ-RHSQ)*H2
          IF (RADH2.GE.0.01D0 .AND. AEHH.LT.0.01D0) THEN

              SUCVA(1,1) = H* (COSHNX*SNHR1-COSHMX*SNHR2)/PRAD

          ELSE IF (RADH2.LT.0.01D0 .AND. AEHH.GE.0.01D0) THEN

              SUCVA(1,1) = H* (SN2ERT-SN2HRT)/ (2.D0*PNM)

          ELSE

              SUCVA(1,1) = 2.D0* (H**3)*FAC*
     &                     SNHDIF(2.D0*AETAH,2.D0*ARHH)/P
          END IF

          SUCVA(1,2) = 0.25D0*H*H* (SNERAT+SNHRAT)/P
          SUCVA(2,1) = SUCVA(1,2)
          SUCVA(2,2) = 0.5D0*H* (SN2ERT+SN2HRT)/P

      END IF

C Compute CU.Adjoint(SV).exp(-NUX-MUX):

      IF (.NOT.SMALL) THEN
          CUSVA(1,1) = P* (NU3*CMSN-MU3*CNSM)/RAD
          CUSVA(1,2) = PNM* (NM* (NU2+MU2)* (FAC-CNCM)+ (NU4+MU4)*SNSM)/
     &                 RAD2
          CUSVA(2,1) = CUSVA(1,2)
          CUSVA(2,2) = PNM* (NU3*CNSM-MU3*CMSN)/RAD

      ELSE
          FACL = 0.25D0*PNM*H

          CUSVA(1,1) = P*H*0.5D0* ((ETASQ+3.D0*RHSQ)*SN2HRT+
     &                 (RHSQ+3.D0*ETASQ)*SN2ERT)


          CUSVA(1,2) = 0.25D0*PNM* (H2* (RHSQ*SNERAT-ETASQ*SNHRAT)+
     &                 3.D0*SNSM)
          CUSVA(2,1) = CUSVA(1,2)

          CUSVA(2,2) = 2.D0*FACL* ((RHSQ+3.D0*ETASQ)*SN2ERT-
     &                 (ETASQ+3.D0*RHSQ)*SN2HRT)
      END IF

C Compute A = U.Adjoint(V):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = V(1,1)*U(2,2) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.A.Adjoint(SV))*exp(-NUX-MUX):

      IF (.NOT.SMALL) THEN

          T(1,1) = A(1,1)*CNCM + 2.D0*A(1,2)* (NU*CMSN-MU*CNSM)/RAD +
     &             A(2,2)*SNSM/NM

          T(1,2) = (A(1,1)*NM* (NU*CNSM-MU*CMSN)+
     &             A(1,2)* ((FAC* (NU2+MU2)+2.D0*NM*SNSM)* (NU2+MU2)-
     &             4.D0*NM*NM*CNCM)/RAD+A(2,2)* (NU*CMSN-MU*CNSM))/RAD
          T(2,1) = T(1,2)

          T(2,2) = NM* (A(1,1)*SNSM+2.D0*A(1,2)* (NU*CNSM-MU*CMSN)/
     &             RAD) + A(2,2)*CNCM

      ELSE

          T(1,1) = A(1,1)*CNCM + 2.D0*P*A(1,2)*SUCVA(2,2) +
     &             A(2,2)*H2*SNHR1*SNHR2

          T(1,2) = A(1,1)*PNM*NM*SUCVA(1,1) +
     &             0.5D0*A(1,2)* (FAC+CNCM-H2*
     &             (ETASQ*SNHRAT+RHSQ*SNERAT)) + A(2,2)*P*SUCVA(2,2)
          T(2,1) = T(1,2)

          T(2,2) = NM* (A(1,1)*SNSM+2.D0*A(1,2)*PNM*SUCVA(1,1)) +
     &             A(2,2)*CNCM

      END IF

C Compute U.Adjoint(V) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      IF (.NOT.SMALL) THEN

          DETCU = ((NU4+MU4)*CNCM-NM* ((NU2+MU2)*SNSM+2.D0*NM*FAC))/RAD2
          DETCV = DETCU
          DETSU = (2.D0*NM* (FAC-CNCM)+ (NU2+MU2)*SNSM)/ (PRAD2*NM)
          DETSV = P*P*NM* (2.D0*NU3*MU3* (FAC-CNCM)+ (NU6+MU6)*SNSM)/
     &            RAD2

      ELSE

          DETCU = 0.25D0* (H2* (ETASQ*SNHRAT+RHSQ*SNERAT)+FAC+3.D0*CNCM)
          DETCV = DETCU
          RADH2 = ABS(RAD)*H2
          AEHH = ABS(ETASQ-RHSQ)*H2
          IF (RADH2.GE.0.01D0 .AND. AEHH.LT.0.01D0) THEN
              FACL = H2*SNHR1*SNHR2

              DETSU = (2.D0* (FAC-CNCM)+ (NU2+MU2)*FACL)/PRAD2

          ELSE IF (RADH2.LT.0.01D0 .AND. AEHH.GE.0.01D0) THEN
              FACL = H*0.5D0/P

              DETSU = (SNERAT-SNHRAT)*FACL*FACL/ (ETASQ-RHSQ)

          ELSE
              FACL = 0.5D0*H2/P

              DETSU = FACL*FACL*SNHDIF(AETAH,ARHH)*FAC

          END IF

          DETSV = (P*PNM/8.D0)* (2.D0*H*H*
     &            (RHSQ*RHSQ*SNERAT-ETASQ*ETASQ*SNHRAT)+
     &            3.D0*NM* (CNCM-FAC)+15.D0* (ETASQ+RHSQ)*SNSM)
      END IF


C Compute trace(CU.A.Adjoint(SU)).exp(-NUX-MUX):

      BETA = A(1,1)*SUCVA(2,2) + 2.D0*A(1,2)*SUCVA(1,2) +
     &       A(2,2)*SUCVA(1,1)
C      END IF

      UUA(1,1) = DETU*DETCU + DETV*DETSU + BETA
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)


C Compute trace(CV.A.Adjoint(SV)).exp(-NUX-MUX):

      BETA = A(1,1)*CUSVA(2,2) + 2.D0*A(1,2)*CUSVA(1,2) +
     &       A(2,2)*CUSVA(1,1)

C      END IF

      VVA(1,1) = DETU*DETSV + DETV*DETCV + BETA
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)
      ELSE
          M = ATAN2(-B,-AC)
      END IF

C Now we have got the value of M we can compute R.

      ZCOSMT = COS(M/2.D0)
      ZSINMT = SIN(M/2.D0)


      DO 20 J = 1,2
          DO 10 I = 1,2
              R(I,J) = UT(I,J)*ZCOSMT - VT(I,J)*ZSINMT
   10     CONTINUE
   20 CONTINUE

C Compute the DET R and RINV
      DETR = ZCOSMT**2*UUA(1,1) + ZSINMT**2*VVA(1,1) -
     &       ZSINMT*ZCOSMT* (UVA(1,1)+UVA(2,2))

C      WRITE (8,FMT=*) 'DETR:',DETR
C
C      WRITE (8,FMT=*) 'DETR*ADJUST-DET(R),NUX,MUX:',DETR*ADJUST-
C     &                 R(1,1)*R(2,2)+R(1,2)*R(2,1),
C     &                 NUX,MUX
C      NUX = NUX-MUX
C      NUX = 0.D0
      NUX = MUX
C      NUX = NUX + 2.D0*MUX
      RINV(1,1) = R(2,2)/DETR
      RINV(1,2) = -R(1,2)/DETR
      RINV(2,1) = -R(2,1)/DETR
      RINV(2,2) = R(1,1)/DETR

      DO 40 J = 1,2
          DO 30 I = 1,2
              U(I,J) = (UUA(I,J)*ZCOSMT-UVA(I,J)*ZSINMT)/DETR
              V(I,J) = (VUA(I,J)*ZCOSMT-VVA(I,J)*ZSINMT)/DETR
   30     CONTINUE
   40 CONTINUE

C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF


C NOTE: in the following formulae we want CV and SU for BACKWARDS
C integration. We exploit the fact that the diagonal terms of
C SUCVA are odd functions of the direction of integration while
C the off-diagonal terms are even functions of the direction
C of integration, while the determinants which appear in the
C formulae are even functions of the direction of integration.

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END


C -------------------------------------------------------------------


      SUBROUTINE NEGLM1(ELAM,NUX,RINV,P,S,Q,W,U,V,SU,CV,XO,XEND)

C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,NUX,P,Q,S,W,XEND,XO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CV(2,2),RINV(2,2),SU(2,2),U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,ALFA,ALFA2,ALFA4,AX,B,BETA,BETA2,BETA4,BX,C,
     &                 CAXCBX,CAXSBX,COSBX,COSBX1,COSCOS,COSHAX,CSHCSH,
     &                 D,DETCU,DETCV,DETR,DETSU,DETSV,DETU,DETV,FAC,H,M,
     &                 PRAD,RAD,RAD2,RADA,RADB,RADM,SAXCBX,SAXSBX,SINBX,
     &                 SINBX1,SINCOS,SINHAX,SINSIN,SNHCSH,SNHSNH,TR,
     &                 TRUVA,ZCOSMT,ZSINMT
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(2,2),CU(2,2),CUSVA(2,2),R(2,2),SUCVA(2,2),
     &                 SV(2,2),T(2,2),UT(2,2),UUA(2,2),UVA(2,2),VT(2,2),
     &                 VUA(2,2),VVA(2,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN2,COS,EXP,SIGN,SIN,SQRT
C     ..
C     .. External Functions ..

      DOUBLE PRECISION PHI
      EXTERNAL PHI
C     ..
C     .. External Subroutines ..
      EXTERNAL F06YAF
C     ..
      D = 2.D0*SQRT(P* (Q-ELAM*W))
      ALFA = 0.5D0*SQRT((S+D)/P)
      BETA = 0.5D0*SQRT((-S+D)/P)

      DETU = U(1,1)*U(2,2) - U(2,1)*U(1,2)
      DETV = V(1,1)*V(2,2) - V(2,1)*V(1,2)


      H = XEND - XO
      AX = ABS(ALFA*H)
      FAC = EXP(-AX)
      NUX = AX
      IF (AX.GE.0.1D0) THEN
          SINHAX = (1.D0-FAC*FAC)/2.D0
          COSHAX = 1.D0 - SINHAX
          SINHAX = SINHAX*SIGN(1.D0,H)/ALFA

      ELSE
          SINHAX = H*PHI(AX*AX)*FAC
          COSHAX = (1.D0+FAC*FAC)/2.D0
      END IF

      BX = BETA*H
      IF (ABS(BX).GE.0.1D0) THEN
          SINBX1 = SIN(BX)/BETA

      ELSE
          SINBX1 = H*PHI(-BX*BX)
      END IF

      SINBX = SINBX1*FAC

      COSBX1 = COS(BX)
      COSBX = COSBX1*FAC

      CAXCBX = COSHAX*COSBX1
      CAXSBX = COSHAX*SINBX1
      SAXSBX = SINHAX*SINBX1
      SAXCBX = SINHAX*COSBX1

      SNHCSH = (SINHAX*COSHAX)
      SNHSNH = (SINHAX*SINHAX)
      CSHCSH = (COSHAX*COSHAX)

      SINCOS = (SINBX*COSBX)
      SINSIN = (SINBX*SINBX)
      COSCOS = (COSBX*COSBX)

      ALFA2 = ALFA*ALFA
      ALFA4 = ALFA2*ALFA2
      BETA2 = BETA*BETA
      BETA4 = BETA2*BETA2
      RAD = ALFA2 + BETA2
      RAD2 = RAD*RAD
      PRAD = P*RAD
      RADM = ALFA2 - BETA2
      RADA = 3.D0*ALFA2 - BETA2
      RADB = 3.D0*BETA2 - ALFA2
C Compute the CU,CV,SU,SV matrices here.

      CU(1,1) = -RADM*SAXSBX*0.5D0 + CAXCBX
      CU(1,2) = 0.5D0*SAXCBX + 0.5D0*CAXSBX
      CU(2,1) = -0.5D0*RAD* (CAXSBX-SAXCBX)
      CU(2,2) = CAXCBX + 0.5D0*RADM*SAXSBX

      SU(1,1) = (SAXCBX-CAXSBX)/ (2.0D0*PRAD)
      SU(1,2) = SAXSBX/ (2.0D0*P)
      SU(2,1) = -SU(1,2)
      SU(2,2) = (CAXSBX+SAXCBX)/ (2.0D0*P)

      SV(1,1) = (((-2.D0*ALFA4+2.D0*BETA4+RAD2)*CAXSBX)+
     &          ((-2.D0*BETA4+2.D0*ALFA4+RAD2)*SAXCBX))* (0.5D0*P)
      SV(1,2) = RAD2*SAXSBX*0.5D0*P
      SV(2,1) = -SV(1,2)
      SV(2,2) = (RADA*SAXCBX-RADB*CAXSBX)* (0.5D0*P)

      CV(1,1) = -RADM*SAXSBX*0.5D0 + CAXCBX
      CV(1,2) = (RAD*CAXSBX-RAD*SAXCBX)/2.0D0
      CV(2,1) = -0.5D0* (CAXSBX+SAXCBX)
      CV(2,2) = (RADM*SAXSBX*0.5D0) + CAXCBX

C Now we can compute UT and VT


      CALL F06YAF('N','N',2,2,2,1.D0,CU,2,U,2,0.D0,UT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,SU,2,V,2,1.D0,UT,2)

      CALL F06YAF('N','N',2,2,2,1.D0,SV,2,U,2,0.D0,VT,2)
      CALL F06YAF('N','N',2,2,2,1.D0,CV,2,V,2,1.D0,VT,2)


C Compute SU.Adjoint(CV).exp(-AX):

      SUCVA(1,1) = 0.5D0* (SNHCSH-SINCOS)/PRAD
      SUCVA(1,2) = 0.25D0* (SNHSNH+SINSIN)/P
      SUCVA(2,1) = SUCVA(1,2)
      SUCVA(2,2) = 0.5D0* (SNHCSH+SINCOS)/P

C Compute CU.Adjoint(SV).exp(-AX):

      CUSVA(1,1) = 0.5D0* (RADA*SNHCSH-RADB*SINCOS)/P
      CUSVA(1,2) = -0.25D0*RAD* (SNHSNH*BETA2-3.D0*CSHCSH+SINSIN*ALFA2+
     &             3.D0*COSCOS)/P
      CUSVA(2,1) = CUSVA(1,2)
      CUSVA(2,2) = 0.5D0*PRAD* (RADA*SNHCSH+RADB*SINCOS)

C Compute A = U.Adjoint(V)(xo):

      A(1,1) = U(1,1)*V(2,2) - U(1,2)*V(2,1)
      A(1,2) = U(1,2)*V(1,1) - U(1,1)*V(1,2)
      A(2,1) = A(1,2)
      A(2,2) = U(2,2)*V(1,1) - U(2,1)*V(1,2)

C Compute T = (CU.A.Adjoint(CV) + SU.Adjoint(A).Adjoint(SV))*exp(-AX):

      T(1,1) = A(1,1)* (CSHCSH-SINSIN*BETA2) + A(1,2)* (SNHCSH+SINCOS) +
     &         A(2,2)* (CSHCSH-COSCOS)/RAD

      T(1,2) = 0.5D0*RAD*A(1,1)* (SNHCSH-SINCOS) +
     &         0.5D0*A(1,2)* (CSHCSH+COSCOS+BETA2*SNHSNH-ALFA2*SINSIN) +
     &         0.5D0*A(2,2)* (SNHCSH+SINCOS)

      T(2,1) = T(1,2)

      T(2,2) = A(1,1)*RAD* (CSHCSH-COSCOS) +
     &         A(1,2)*RAD* (SNHCSH-SINCOS) +
     &         A(2,2)* (CSHCSH-SINSIN*BETA2)

C Compute U.Adjoint(V)(xend) = (CU.Uo+SU.Vo).Adjoint(SV.Uo+CV.Vo):

      UVA(1,1) = DETU*CUSVA(1,1) + DETV*SUCVA(1,1) + T(1,1)
      UVA(1,2) = DETU*CUSVA(1,2) + DETV*SUCVA(1,2) + T(1,2)
      UVA(2,1) = DETU*CUSVA(2,1) + DETV*SUCVA(2,1) + T(2,1)
      UVA(2,2) = DETU*CUSVA(2,2) + DETV*SUCVA(2,2) + T(2,2)

      VUA(1,1) = UVA(2,2)
      VUA(2,2) = UVA(1,1)
      VUA(1,2) = -UVA(1,2)
      VUA(2,1) = -UVA(2,1)

      DETCU = 0.25D0*RADM* (SNHSNH+SINSIN) + 0.5D0* (CSHCSH+COSCOS)
      DETSU = 0.25* (SNHSNH-SINSIN)/ (P*PRAD)
      DETCV = DETCU
      DETSV = 0.25D0*P*PRAD* (RADA*RADA*SNHSNH-RADB*RADB*SINSIN)

C Compute trace(CU.A.Adjoint(SU)).exp(-AX):

      TR = 0.5D0* (A(1,1)* (SNHCSH+SINCOS)+A(1,2)* (SNHSNH+SINSIN))/P +
     &     0.5D0*A(2,2)* (SNHCSH-SINCOS)/PRAD

      UUA(1,1) = DETU*DETCU + DETV*DETSU + TR
      UUA(1,2) = 0.D0
      UUA(2,1) = 0.D0
      UUA(2,2) = UUA(1,1)

C Compute trace(CV.Adjoint(A).Adjoint(SU)).exp(-AX):

      TR = 0.5D0*P*RADA* (RAD* (A(1,1)*SNHCSH+A(1,2)*SNHSNH)+
     &     A(2,2)*SNHCSH) + 0.5D0*P*RADB*
     &     (RAD* (A(1,1)*SINCOS+A(1,2)*SINSIN)-A(2,2)*SINCOS)

      VVA(1,1) = DETU*DETSV + DETV*DETCV + TR
      VVA(1,2) = 0.D0
      VVA(2,1) = 0.D0
      VVA(2,2) = VVA(1,1)

      TRUVA = UVA(1,1) + UVA(2,2)

      AC = UUA(1,1) - VVA(1,1)
      B = -TRUVA
      C = UUA(1,1) + VVA(1,1)

C We now choose M to maximise the function
C
C ABS(AC.cos(M) + BC.sin(M) + C):

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

C Now we have got the value of M we can compute R.


      ZCOSMT = COS(M/2.D0)
      ZSINMT = SIN(M/2.D0)

      DO 20 J = 1,2
          DO 10 I = 1,2
              R(I,J) = (UT(I,J)*ZCOSMT-VT(I,J)*ZSINMT)
   10     CONTINUE
   20 CONTINUE

C Compute the DET R and RINV

      DETR = ZCOSMT**2*UUA(1,1) + ZSINMT**2*VVA(1,1) -
     &       ZSINMT*ZCOSMT* (UVA(1,1)+UVA(2,2))



      RINV(1,1) = R(2,2)/DETR
      RINV(1,2) = -R(1,2)/DETR
      RINV(2,1) = -R(2,1)/DETR
      RINV(2,2) = R(1,1)/DETR

      U(1,1) = (UUA(1,1)*ZCOSMT-UVA(1,1)*ZSINMT)/DETR
      U(1,2) = (-UVA(1,2)*ZSINMT)/DETR
      U(2,1) = (-UVA(2,1)*ZSINMT)/DETR
      U(2,2) = (UUA(2,2)*ZCOSMT-UVA(2,2)*ZSINMT)/DETR

      V(1,1) = (VUA(1,1)*ZCOSMT-VVA(1,1)*ZSINMT)/DETR
      V(1,2) = (VUA(1,2)*ZCOSMT)/DETR
      V(2,1) = (VUA(2,1)*ZCOSMT)/DETR
      V(2,2) = (VUA(2,2)*ZCOSMT-VVA(2,2)*ZSINMT)/DETR

C That's done the very tricky renormalisation. Now do the SU and CV
C matrices which are the U and V matrices for Theta_{0}.

      TRUVA = SUCVA(1,1) + SUCVA(2,2)

      AC = (DETSU-DETCV)
      B = -TRUVA
      C = (DETSU+DETCV)

      IF (C.GE.0.D0) THEN
          M = ATAN2(B,AC)

      ELSE
          M = ATAN2(-B,-AC)
      END IF

      SU(1,1) = DETSU*COS(M/2.D0) + SUCVA(1,1)*SIN(M/2.D0)
      SU(1,2) = -SUCVA(1,2)*SIN(M/2.D0)
      SU(2,1) = -SUCVA(2,1)*SIN(M/2.D0)
      SU(2,2) = DETSU*COS(M/2.D0) + SUCVA(2,2)*SIN(M/2.D0)

      CV(1,1) = -SUCVA(2,2)*COS(M/2.D0) - DETCV*SIN(M/2.D0)
      CV(1,2) = -SUCVA(1,2)*COS(M/2.D0)
      CV(2,1) = -SUCVA(2,1)*COS(M/2.D0)
      CV(2,2) = -SUCVA(1,1)*COS(M/2.D0) - DETCV*SIN(M/2.D0)

      RETURN

      END


C =====================================================================
C ==================== SUBROUTINES FOR AUTOMATIC MESHING ==============
C =====================================================================


      SUBROUTINE MESH4(XO,XEND,ELAM,SL4COF,XMESH,NMESH,IC,IMESH,TOL,
     &                 YSTOR1,YSTOR2,MULTI,IM1ST,CNMESH,WORK,IWORK,
     &                 IFAIL)
C
C The purpose of this subroutine is to compute a mesh which can be
C used as an initial mesh for solving a 4th-order Sturm-Liouville
C eigenvalue problems.
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,TOL,XEND,XO
      INTEGER CNMESH,IC,IFAIL,IM1ST,IMESH,IWORK,MULTI,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(0:IWORK,1:5),XMESH(0:NMESH),
     &                 YSTOR1(4,0:CNMESH),YSTOR2(4,0:CNMESH)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4COF
C     ..
C     .. Parameters ..
      DOUBLE PRECISION HLF,ONE,SAFE,TRD,FIVE,SAFE3,P008
      PARAMETER (HLF=5.d-1,ONE=1.d0,SAFE=0.8D0,TRD=ONE/3.d0,FIVE=5.d0,
     &          SAFE3=7.29d-1,P008=8.d-3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DI,EPS,ERREST,H,HMAX,RATIO,TEND,TO,TOLOC,TOLOW,
     &                 X1,X3
      INTEGER I,ICNT,IDI,IRF,J,KNTR
      LOGICAL PHASE1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION GETERR,X02AJF
      EXTERNAL GETERR,X02AJF
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,LOG,MAX,MIN,SIGN
C     ..
      TO = XO/ (1.D0+ABS(XO))
      TEND = XEND/ (1.D0+ABS(XEND))
      IFAIL = 0
      EPS = X02AJF(0.D0)*1.D2
      IF (TO.EQ.TEND) THEN
          IMESH = 0
          XMESH(IC) = XO
          RETURN
      ELSE
          DI = SIGN(1.D0,TEND-TO)
          IDI = DI
      END IF

      IF (IDI.GT.0) THEN
          I = IC

      ELSE
          I = NMESH
      END IF

      XMESH(I) = TO

      ICNT = 0
      KNTR = 0
      IRF = INT(LOG(X02AJF(0.D0))*HLF/LOG(SAFE))
      TOLOC = TOL
      TOLOW = TOLOC*SAFE
      PHASE1 = .TRUE.
C Set initial stepsize:
      H = MIN(ABS(TEND-TO),1.D0)/20.D0
C Set maximum stepsize:
      HMAX = ABS(TEND-TO)/20.D0

   10 X3 = XMESH(I) + H*DI
      X1 = XMESH(I)
      IF (IDI.GT.0) X3 = MIN(X3,TEND)
      IF (IDI.LT.0) X3 = MAX(X3,TEND)

C Get an error estimate.

      ERREST = GETERR(X1,X3,SL4COF,ELAM,YSTOR1,YSTOR2,MULTI,IM1ST,
     &         CNMESH,WORK,IWORK)

      IF (ERREST.GT.TOLOC) THEN
          RATIO = SAFE
          IF (TOLOC.LT.ERREST*SAFE3) RATIO = (TOLOC/ERREST)**TRD
          H = H*RATIO
          ICNT = ICNT + 1
          IF (ICNT.LT.IRF) GO TO 10
          IFAIL = 10
          RETURN

      END IF

      IF (ERREST.LT.TOLOW .AND. ICNT.EQ.0) THEN
          RATIO = FIVE
          IF (ERREST.GT.TOLOW*P008) RATIO = (TOLOW/ERREST)**TRD
          H = MIN(HMAX,H*RATIO)
          IF (PHASE1 .AND. H.LT.HMAX) THEN
              KNTR = KNTR + 1
              IF (KNTR.LT.5) GO TO 10
          END IF

      END IF

      PHASE1 = .FALSE.
      ICNT = 0
      I = I + IDI
      XMESH(I) = X3
      IF ((IDI.GT.0.AND.X3.LT.TEND-EPS) .OR.
     &    (IDI.LT.0.AND.X3.GT.TEND+EPS)) THEN
C We have not yet reached the end of the range. Meshing must continue
C so we must check that there is space to hold at least another one
C mesh-point.
          IF (IDI.GT.0 .AND. I.GE.NMESH) THEN
              IFAIL = 11
              RETURN

          END IF

          IF (IDI.LT.0 .AND. I.LE.IC) THEN
              IFAIL = 11
              RETURN

          END IF

          GO TO 10

      END IF

C If we are here, it means that the meshing has finished. We must
C move the mesh into the place where the calling routine expects
C to find it.

      IF (IDI.LT.0 .AND. I.GT.IC) THEN
          DO 20 J = 0,NMESH - I
              XMESH(IC+J) = XMESH(I+J)
   20     CONTINUE
      END IF

      IF (IDI.GT.0) THEN
          IMESH = I - IC
C Transform mesh back to original coordinates
          DO 30 J = IC,I
              XMESH(J) = XMESH(J)/ (1.D0-ABS(XMESH(J)))
   30     CONTINUE
      ELSE
          IMESH = NMESH - I
C Transform mesh back to original coordinates
          DO 40 J = IC,IMESH + IC
              XMESH(J) = XMESH(J)/ (1.D0-ABS(XMESH(J)))
   40     CONTINUE
      END IF


      RETURN

      END


C ---------------------------------------------------------------------


      DOUBLE PRECISION FUNCTION GETERR(X1,X3,SL4COF,ELAM,YSTOR1,YSTOR2,
     &                 MULTI,IMATCH,NMESH,WORK,IWORK)
C     .. External Subroutines ..
      EXTERNAL EFN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION Y(1:4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,X1,X3
      INTEGER IMATCH,IWORK,MULTI,NMESH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(0:IWORK,5),YSTOR1(4,0:NMESH),
     &                 YSTOR2(4,0:NMESH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DDPY21,DDPY23,DY21,DY23,ERRST,FAC1,FAC3,P1,P2,P3,
     &                 Q1,Q2,Q3,S1,S2,S3,W1,W2,W3,X2,Y21,Y23
      LOGICAL HOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SL4COF
C     ..
      X2 = (X3+X1)*0.5D0

C At the moment we only form the error estimate using one of the
C eigenfunctions, in the case of an eigenvalue of multiplicity 2.
      FAC1 = 1.D0/ (1.D0-ABS(X1))**2
      FAC3 = 1.D0/ (1.D0-ABS(X3))**2

      CALL SL4COF(X1/ (1.D0-ABS(X1)),P1,S1,Q1,W1)
      CALL SL4COF(X2/ (1.D0-ABS(X2)),P2,S2,Q2,W2)
      CALL SL4COF(X3/ (1.D0-ABS(X3)),P3,S3,Q3,W3)

      HOT = .FALSE.

      CALL EFN(WORK(0,1),WORK(1,2),WORK(1,3),WORK(1,4),WORK(1,5),ELAM,
     &         X1/ (1.D0-ABS(X1)),YSTOR1,IMATCH,NMESH,Y,HOT)

      Y21 = Y(1)**2
      DY21 = Y(2)**2
      DDPY21 = Y(4)**2

      Y21 = MAX(Y21,ABS(Y(1)))
      DY21 = MAX(DY21,ABS(Y(2)))
      DDPY21 = MAX(DDPY21,ABS(Y(4)))

      CALL EFN(WORK(0,1),WORK(1,2),WORK(1,3),WORK(1,4),WORK(1,5),ELAM,
     &         X3/ (1.D0-ABS(X3)),YSTOR1,IMATCH,NMESH,Y,HOT)

      Y23 = Y(1)**2
      DY23 = Y(2)**2
      DDPY23 = Y(4)**2

      Y23 = MAX(Y23,ABS(Y(1)))
      DY23 = MAX(DY23,ABS(Y(2)))
      DDPY23 = MAX(DDPY23,ABS(Y(4)))

      ERRST = ((Q1-Q2)-ELAM* (W1-W2))*Y21 + (S1-S2)*DY21 +
     &        ((1.D0/P2)- (1.D0/P1))*DDPY21

      ERRST = FAC3* (((Q3-Q2)-ELAM* (W3-W2))*Y23+ (S3-S2)*DY23+
     &        ((1.D0/P2)- (1.D0/P3))*DDPY23) + FAC1*ERRST

      ERRST = ERRST* ((X3-X1)/6.D0)

      GETERR = ABS(ERRST)

      RETURN

      END


C ---------------------------------------------------------------


      SUBROUTINE COARS4(N,A,B,XMESH,ISING)

C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B
      INTEGER ISING,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XMESH(0:N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
      XMESH(0) = A
      XMESH(N) = B

      IF (ISING.EQ.0) THEN
C     THE PROBLEM IS REGULAR, USE A UNIFORM MESH

          H = (B-A)/DBLE(N)
          DO 10 I = 1,N - 1
              XMESH(I) = XMESH(I-1) + H
   10     CONTINUE

      ELSE IF (ISING.EQ.1) THEN
C     X = A IS SINGULAR BUT X = B IS NOT

          H = (B-A)/DBLE(N-4)
          XMESH(1) = XMESH(0) + H/16.D0
          XMESH(2) = XMESH(1) + H/16.D0
          XMESH(3) = XMESH(2) + H/8.D0
          XMESH(4) = XMESH(3) + H/4.D0
          XMESH(5) = XMESH(4) + H/2.D0
          DO 20 I = 6,N - 1
              XMESH(I) = XMESH(I-1) + H
   20     CONTINUE

      ELSE IF (ISING.EQ.2) THEN
C     X=A IS REGULAR BUT X = B IS SINGULAR

          H = (B-A)/DBLE(N-4)

          XMESH(N-1) = XMESH(N) - H/16.D0
          XMESH(N-2) = XMESH(N-1) - H/16.D0
          XMESH(N-3) = XMESH(N-2) - H/8.D0
          XMESH(N-4) = XMESH(N-3) - H/4.D0
          XMESH(N-5) = XMESH(N-4) - H/2.D0
          DO 30 I = N - 6,1,-1
              XMESH(I) = XMESH(I+1) - H
   30     CONTINUE

      ELSE IF (ISING.EQ.3) THEN
C     BOTH ENDS OF THE INTERVAL ARE SINGULAR POINTS.

          H = (B-A)/DBLE(N-8)

          XMESH(1) = XMESH(0) + H/16.D0
          XMESH(2) = XMESH(1) + H/16.D0
          XMESH(3) = XMESH(2) + H/8.D0
          XMESH(4) = XMESH(3) + H/4.D0
          XMESH(5) = XMESH(4) + H/2.D0

          XMESH(N-1) = XMESH(N) - H/16.D0
          XMESH(N-2) = XMESH(N-1) - H/16.D0
          XMESH(N-3) = XMESH(N-2) - H/8.D0
          XMESH(N-4) = XMESH(N-3) - H/4.D0
          XMESH(N-5) = XMESH(N-4) - H/2.D0

          DO 40 I = 6,N - 6
              XMESH(I) = XMESH(I-1) + H
   40     CONTINUE

      ELSE
C     THERE IS AN ERROR IN THE VALUE OF ISING
      END IF

      RETURN

      END


C =====================================================================
C ========================== AUXILIARY SUBROUTINES ====================
C =====================================================================


      DOUBLE PRECISION FUNCTION PHI(V)
C     .. Parameters ..
      DOUBLE PRECISION A0,A1,A2,A3,A4,A5
      PARAMETER (A0=1.D0,A1=A0/6.D0,A2=A1/20.D0,A3=A2/42.D0,A4=A3/72.D0,
     &          A5=A4/110.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION V
C     ..
      PHI = A0 + V* (A1+V* (A2+V* (A3+V* (A4+V*A5))))
      RETURN

      END


C ---------------------------------------------------------------------


      DOUBLE PRECISION FUNCTION SNHDIF(ETA,H)
C This function computes the quantity
C      (sinh(eta)/eta - sinh(h)/h)/(eta**2-h**2)
C in the case where eta and h are both small.
C
C     .. Parameters ..
      DOUBLE PRECISION FAC1,FAC2,FAC3,FAC4
      PARAMETER (FAC1=1.D0/6.D0,FAC2=FAC1/20.D0,FAC3=FAC2/42.D0,
     &          FAC4=FAC3/72.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ETA,H
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ETA2,ETA4,ETA6,H2,H4,H6
C     ..
      H2 = H*H
      H4 = H2*H2
      H6 = H4*H2
      ETA2 = ETA*ETA
      ETA4 = ETA2*ETA2
      ETA6 = ETA4*ETA2

      SNHDIF = FAC1 + FAC2* (ETA2+H2) + FAC3* (ETA4+H4+ETA2*H2) +
     &         FAC4* (ETA6+H6+ETA4*H2+ETA2*H4)

      RETURN

      END

C ---------------------------------------------------------------------

      SUBROUTINE MONIT(A,IA,N,IPRN)
C     .. Scalar Arguments ..
      INTEGER IA,IPRN,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(IA,N)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
      DO 10 J = 1,N
          WRITE (IPRN,FMT=*) (REAL(A(J,I)),I=1,N)
   10 CONTINUE
      RETURN

      END

c
c
      SUBROUTINE DABMOD(N,NCAPMD,MCAP,DEPS,IROUT,DA,DB,DXIR,DXII,IS,
     +                  DALPHA,DBETA,NCAPD,KOUNTD,IERRD,DBE,DX,DW,DE,
     +                  DP0,DP1,DP2)
c      double precision theta
c      common/par/theta
c     .. Scalar Arguments ..
      DOUBLE PRECISION DEPS
      INTEGER IERRD,IROUT,KOUNTD,MCAP,N,NCAPD,NCAPMD
c     ..
c     .. Array Arguments ..  The upper bounds DXII, DXIR and IS are all
c                            MCAP, but this doesn't work in Fortran 77
c                            if MCAP = 0.
      DOUBLE PRECISION DA(NCAPMD),DALPHA(N),DB(NCAPMD),DBE(N),DBETA(N),
     +                 DE(NCAPMD),DP0(NCAPMD),DP1(NCAPMD),DP2(NCAPMD),
     +                 DW(NCAPMD),DX(NCAPMD),DXII(*),DXIR(*)
      INTEGER IS(*)
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION DEPSM,DP
      INTEGER IED,IERRGD,IMU,INCAP,K,MU
c     ..
c     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
c     ..
c     .. External Subroutines ..
      EXTERNAL DGAUSS,DLANCZ,DSTI
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC ABS
c     ..
      integer i1mach, nout
      nout = i1mach(2)
      IERRD = 0
      DEPSM = D1MACH(3)
      INCAP = 1
      DO 10 K = 1,N
          DALPHA(K) = DA(K)
          DBETA(K) = DB(K)
   10 CONTINUE
      IF (MCAP.GT.0) THEN
          KOUNTD = -1
          DO 20 K = 1,N
              DBETA(K) = 0.D0
   20     CONTINUE
          NCAPD = (2*N-1)/2
   30     DO 40 K = 1,N
              DBE(K) = DBETA(K)
   40     CONTINUE
          KOUNTD = KOUNTD + 1
          IF (KOUNTD.GT.1) INCAP = 2** (KOUNTD/5)*N
          NCAPD = NCAPD + INCAP
          IF (NCAPD.GT.NCAPMD) THEN
              IERRD = NCAPMD
              RETURN

          END IF

          CALL DGAUSS(NCAPD,DA,DB,DEPSM,DX,DW,IERRGD,DE)
          IF (IERRGD.NE.0) THEN
              WRITE (NOUT,FMT=9000) IERRGD
              IERRD = 1
              RETURN

          END IF

          DO 60 K = 1,NCAPD
              DP = 1.D0
              IMU = 0
              DO 50 MU = 1,MCAP
                  IF (IMU.EQ.0) THEN
                      IF (DXII(MU).EQ.0.D0) THEN
                          DP = ((1.D0+DXIR(MU)*DX(K))**IS(MU))*DP

                      ELSE
                          DP = (((1.D0+DXIR(MU)*DX(K))**2+
     +                         (DXII(MU)*DX(K))**2)**IS(MU))*DP
                          IMU = 1
                      END IF

                  ELSE
                      IMU = 0
                  END IF

   50         CONTINUE
              DW(K) = DW(K)/DP
c          dw(k)=dsqrt(1.d0+.5d0*theta*dx(k))*dw(k)/dp
   60     CONTINUE
          IF (IROUT.EQ.1) THEN
              CALL DSTI(N,NCAPD,DX,DW,DALPHA,DBETA,IED,DP0,DP1,DP2)
              IF (IED.NE.0) THEN
                  IERRD = 2
                  RETURN

              END IF

          ELSE
              CALL DLANCZ(N,NCAPD,DX,DW,DALPHA,DBETA,IED,DP0,DP1)
              IF (IED.NE.0) THEN
                  IERRD = 2
                  RETURN

              END IF

          END IF

          DO 70 K = 1,N
              IF (ABS(DBETA(K)-DBE(K)).GT.DEPS*DBETA(K)) GO TO 30
   70     CONTINUE
      END IF

      RETURN

 9000 FORMAT (1X,'ierrgd in dgauss=',I5)
      END
c
c
      SUBROUTINE DGAUSS(N,DALPHA,DBETA,DEPS,DZERO,DWEIGH,IERR,DE)
c
c This is a double-precision version of the routine  gauss.
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION DEPS
      INTEGER IERR,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DALPHA(N),DBETA(N),DE(N),DWEIGH(N),DZERO(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DB,DC,DF,DG,DP,DR,DS
      INTEGER I,II,J,K,L,M,MML
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
C     ..
      IF (N.LT.1) THEN
          IERR = -1
          RETURN

      END IF

      IERR = 0
      DZERO(1) = DALPHA(1)
      IF (DBETA(1).LT.0.D0) THEN
          IERR = -2
          RETURN

      END IF

      DWEIGH(1) = DBETA(1)
      IF (N.EQ.1) RETURN
      DWEIGH(1) = 1.D0
      DE(N) = 0.D0
      DO 10 K = 2,N
          DZERO(K) = DALPHA(K)
          IF (DBETA(K).LT.0.D0) THEN
              IERR = -2
              RETURN

          END IF

          DE(K-1) = SQRT(DBETA(K))
          DWEIGH(K) = 0.D0
   10 CONTINUE
      DO 80 L = 1,N
          J = 0
   20     DO 30 M = L,N
              IF (M.EQ.N) GO TO 40
              IF (ABS(DE(M)).LE.DEPS* (ABS(DZERO(M))+ABS(DZERO(M+
     +            1)))) GO TO 40
   30     CONTINUE
   40     DP = DZERO(L)
          IF (M.EQ.L) GO TO 80
          IF (J.EQ.30) GO TO 120
          J = J + 1
          DG = (DZERO(L+1)-DP)/ (2.D0*DE(L))
          DR = SQRT(DG*DG+1.D0)
          DG = DZERO(M) - DP + DE(L)/ (DG+SIGN(DR,DG))
          DS = 1.D0
          DC = 1.D0
          DP = 0.D0
          MML = M - L
          DO 70 II = 1,MML
              I = M - II
              DF = DS*DE(I)
              DB = DC*DE(I)
              IF (ABS(DF).LT.ABS(DG)) GO TO 50
              DC = DG/DF
              DR = SQRT(DC*DC+1.D0)
              DE(I+1) = DF*DR
              DS = 1.D0/DR
              DC = DC*DS
              GO TO 60

   50         DS = DF/DG
              DR = SQRT(DS*DS+1.D0)
              DE(I+1) = DG*DR
              DC = 1.D0/DR
              DS = DS*DC
   60         DG = DZERO(I+1) - DP
              DR = (DZERO(I)-DG)*DS + 2.D0*DC*DB
              DP = DS*DR
              DZERO(I+1) = DG + DP
              DG = DC*DR - DB
              DF = DWEIGH(I+1)
              DWEIGH(I+1) = DS*DWEIGH(I) + DC*DF
              DWEIGH(I) = DC*DWEIGH(I) - DS*DF
   70     CONTINUE
          DZERO(L) = DZERO(L) - DP
          DE(L) = DG
          DE(M) = 0.D0
          GO TO 20

   80 CONTINUE
      DO 100 II = 2,N
          I = II - 1
          K = I
          DP = DZERO(I)
          DO 90 J = II,N
              IF (DZERO(J).GE.DP) GO TO 90
              K = J
              DP = DZERO(J)
   90     CONTINUE
          IF (K.EQ.I) GO TO 100
          DZERO(K) = DZERO(I)
          DZERO(I) = DP
          DP = DWEIGH(I)
          DWEIGH(I) = DWEIGH(K)
          DWEIGH(K) = DP
  100 CONTINUE
      DO 110 K = 1,N
          DWEIGH(K) = DBETA(1)*DWEIGH(K)*DWEIGH(K)
  110 CONTINUE
      RETURN

  120 IERR = L
      RETURN

      END
c
c
      SUBROUTINE DGCHRS(N,IOPT,DA,DB,DX,DHP,DHN,DALPHA,DBETA,IERRD)
C     .. Scalar Arguments ..
      DOUBLE PRECISION DHN,DHP,DX
      INTEGER IERRD,IOPT,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DA(N),DALPHA(N),DB(N),DBETA(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DE,DEH,DQ,DQH
      INTEGER K
C     ..
      IERRD = 0
      IF (IOPT.EQ.1) THEN
          DALPHA(1) = DX - DB(1)/DHP
          DBETA(1) = -DHP
          DQ = -DB(1)/DHP
          DO 10 K = 2,N
              DE = DA(K-1) - DX - DQ
              DBETA(K) = DQ*DE
              DQ = DB(K)/DE
              DALPHA(K) = DQ + DE + DX
   10     CONTINUE
          RETURN

      ELSE IF (IOPT.EQ.2) THEN
          DALPHA(1) = DX* (DHP+DHN)/ (DHP-DHN)
          DBETA(1) = - (DHP-DHN)/ (2.D0*DX)
          DQ = -DB(1)/DHP
          DQH = -DHP/DBETA(1)
          DE = 0.D0
          DO 20 K = 2,N
              DEH = DQ + DE + 2.D0*DX - DQH
              DBETA(K) = DQH*DEH
              DE = DA(K-1) - DX - DQ
              DQH = DQ*DE/DEH
              DALPHA(K) = DQH + DEH - DX
              DQ = DB(K)/DE
   20     CONTINUE
          RETURN

      ELSE
          IERRD = 1
          RETURN

      END IF

      END
c
c
      SUBROUTINE DGQRAT(N,MCAP,DALPHA,DBETA,DXIR,DXII,IS,DZG,DWG,DCONST,
     +                  IERRD,DE)
c     .. Scalar Arguments ..
      DOUBLE PRECISION DCONST
      INTEGER IERRD,MCAP,N
c     ..
c     .. Array Arguments ..  The upper bounds DXII, DXIR and IS are all
c                            MCAP, but this doesn't work in Fortran 77
c                            if MCAP = 0.
      DOUBLE PRECISION DALPHA(N+1),DBETA(N+1),DE(N),DWG(N),DXII(*),
     +                 DXIR(*),DZG(N)
      INTEGER IS(*)
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION DEPSM,DP
      INTEGER IERRGD,IMU,K,MU
c     ..
c     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
c     ..
c     .. External Subroutines ..
      EXTERNAL DGAUSS
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC DBLE
c     ..
      integer i1mach, nout
      nout = i1mach(2)
      IERRD = 0
      DEPSM = D1MACH(3)
      CALL DGAUSS(N,DALPHA,DBETA,DEPSM,DZG,DWG,IERRGD,DE)
      IF (IERRGD.NE.0) THEN
          WRITE (nout,FMT=9000) IERRGD
          IERRD = 1
          RETURN

      END IF

      DCONST = DBETA(1)
      DO 20 K = 1,N
          DCONST = DBETA(K+1)*DCONST/DBLE(2*K* (2*K-1))
          IF (MCAP.GT.0) THEN
              DP = 1.D0
              IMU = 0
              DO 10 MU = 1,MCAP
                  IF (IMU.EQ.0) THEN
                      IF (DXII(MU).EQ.0.D0) THEN
                          DP = ((1.D0+DXIR(MU)*DZG(K))**IS(MU))*DP

                      ELSE
                          DP = (((1.D0+DXIR(MU)*DZG(K))**2+
     +                         (DXII(MU)*DZG(K))**2)**IS(MU))*DP
                          IMU = 1
                      END IF

                  ELSE
                      IMU = 0
                  END IF

   10         CONTINUE
              DWG(K) = DP*DWG(K)
          END IF

   20 CONTINUE
      RETURN

 9000 FORMAT (1X,'ierrgd in dgauss=',I5)
      END
c
c
      SUBROUTINE DKNUM(N,NU0,NUMAX,DX,DY,DEPS,DA,DB,DRHOR,DRHOI,NU,IERR,
     +                 DROLDR,DROLDI)
c
c This is a double-precision version of the routine  knum.
c
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION DEPS,DX,DY
      INTEGER IERR,N,NU,NU0,NUMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DA(NUMAX),DB(NUMAX),DRHOI(*),DRHOR(*),DROLDI(*),
     +                 DROLDR(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DDEN,DRI,DRR,DT
      INTEGER J,J1,K,NP1
C     ..
      IERR = 0
      NP1 = N + 1
      IF (NU0.GT.NUMAX) THEN
          IERR = NU0
          RETURN

      END IF

      IF (NU0.LT.NP1) NU0 = NP1
      NU = NU0 - 5
      DO 10 K = 1,NP1
          DRHOR(K) = 0.D0
          DRHOI(K) = 0.D0
   10 CONTINUE
   20 NU = NU + 5
      IF (NU.GT.NUMAX) THEN
          IERR = NUMAX
          GO TO 60

      END IF

      DO 30 K = 1,NP1
          DROLDR(K) = DRHOR(K)
          DROLDI(K) = DRHOI(K)
   30 CONTINUE
      DRR = 0.D0
      DRI = 0.D0
      DO 40 J = 1,NU
          J1 = NU - J + 1
          DDEN = (DX-DA(J1)-DRR)**2 + (DY-DRI)**2
          DRR = DB(J1)* (DX-DA(J1)-DRR)/DDEN
          DRI = -DB(J1)* (DY-DRI)/DDEN
          IF (J1.LE.NP1) THEN
              DRHOR(J1) = DRR
              DRHOI(J1) = DRI
          END IF

   40 CONTINUE
      DO 50 K = 1,NP1
          IF ((DRHOR(K)-DROLDR(K))**2+ (DRHOI(K)-DROLDI(K))**2.GT.
     +        DEPS* (DRHOR(K)**2+DRHOI(K)**2)) GO TO 20
   50 CONTINUE
   60 IF (N.EQ.0) RETURN
      DO 70 K = 2,NP1
          DT = DRHOR(K)*DRHOR(K-1) - DRHOI(K)*DRHOI(K-1)
          DRHOI(K) = DRHOR(K)*DRHOI(K-1) + DRHOI(K)*DRHOR(K-1)
          DRHOR(K) = DT
   70 CONTINUE
      RETURN

      END
c
c
      SUBROUTINE DLANCZ(N,NCAP,DX,DW,DALPHA,DBETA,IERR,DP0,DP1)
c
c This is a double-precision version of the routine  lancz.
c
C     .. Scalar Arguments ..
      INTEGER IERR,N,NCAP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DALPHA(N),DBETA(N),DP0(NCAP),DP1(NCAP),DW(NCAP),
     +                 DX(NCAP)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DGAM,DPI,DRHO,DSIG,DT,DTK,DTMP,DTSIG,DXLAM
      INTEGER I,K
C     ..
      IF (N.LE.0 .OR. N.GT.NCAP) THEN
          IERR = 1
          RETURN

      ELSE
          IERR = 0
      END IF

      DO 10 I = 1,NCAP
          DP0(I) = DX(I)
          DP1(I) = 0.D0
   10 CONTINUE
      DP1(1) = DW(1)
      DO 30 I = 1,NCAP - 1
          DPI = DW(I+1)
          DGAM = 1.D0
          DSIG = 0.D0
          DT = 0.D0
          DXLAM = DX(I+1)
          DO 20 K = 1,I + 1
              DRHO = DP1(K) + DPI
              DTMP = DGAM*DRHO
              DTSIG = DSIG
              IF (DRHO.LE.0.D0) THEN
                  DGAM = 1.D0
                  DSIG = 0.D0

              ELSE
                  DGAM = DP1(K)/DRHO
                  DSIG = DPI/DRHO
              END IF

              DTK = DSIG* (DP0(K)-DXLAM) - DGAM*DT
              DP0(K) = DP0(K) - (DTK-DT)
              DT = DTK
              IF (DSIG.LE.0.D0) THEN
                  DPI = DTSIG*DP1(K)

              ELSE
                  DPI = (DT**2)/DSIG
              END IF

              DTSIG = DSIG
              DP1(K) = DTMP
   20     CONTINUE
   30 CONTINUE
      DO 40 K = 1,N
          DALPHA(K) = DP0(K)
          DBETA(K) = DP1(K)
   40 CONTINUE
      RETURN

      END
c
c
      SUBROUTINE DRECUR(N,IPOLY,DAL,DBE,DA,DB,IDERR)
c
c This is a double-precision version of the routine  recur.
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION DAL,DBE
      INTEGER IDERR,IPOLY,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DA(N),DB(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DAL2,DALPBE,DBE2,DKM1,DLMACH,DT
      INTEGER K
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,DGAMMA,DLGA
      EXTERNAL D1MACH,DGAMMA,DLGA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,EXP,LOG,SQRT
C     ..
      IF (N.LT.1) THEN
          IDERR = 3
          RETURN

      END IF

      DLMACH = LOG(D1MACH(2))
      IDERR = 0
      DO 10 K = 1,N
          DA(K) = 0.D0
   10 CONTINUE
      IF (IPOLY.EQ.1) THEN
          DB(1) = 2.D0
          IF (N.EQ.1) RETURN
          DO 20 K = 2,N
              DKM1 = DBLE(K-1)
              DB(K) = 1.D0/ (4.D0-1.D0/ (DKM1*DKM1))
   20     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.2) THEN
          DA(1) = .5D0
          DB(1) = 1.D0
          IF (N.EQ.1) RETURN
          DO 30 K = 2,N
              DA(K) = .5D0
              DKM1 = DBLE(K-1)
              DB(K) = .25D0/ (4.D0-1.D0/ (DKM1*DKM1))
   30     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.3) THEN
          DB(1) = 4.D0*ATAN(1.D0)
          IF (N.EQ.1) RETURN
          DB(2) = .5D0
          IF (N.EQ.2) RETURN
          DO 40 K = 3,N
              DB(K) = .25D0
   40     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.4) THEN
          DB(1) = 2.D0*ATAN(1.D0)
          IF (N.EQ.1) RETURN
          DO 50 K = 2,N
              DB(K) = .25D0
   50     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.5) THEN
          DB(1) = 4.D0*ATAN(1.D0)
          DA(1) = .5D0
          IF (N.EQ.1) RETURN
          DO 60 K = 2,N
              DB(K) = .25D0
   60     CONTINUE
          RETURN

      ELSE IF (IPOLY.EQ.6) THEN
          IF (DAL.LE.-1.D0 .OR. DBE.LE.-1.D0) THEN
              IDERR = 1
              RETURN

          ELSE
              DALPBE = DAL + DBE
              DA(1) = (DBE-DAL)/ (DALPBE+2.D0)
              DT = (DALPBE+1.D0)*LOG(2.D0) + DLGA(DAL+1.D0) +
     +             DLGA(DBE+1.D0) - DLGA(DALPBE+2.D0)
              IF (DT.GT.DLMACH) THEN
                  IDERR = 2
                  DB(1) = D1MACH(2)

              ELSE
                  DB(1) = EXP(DT)
              END IF

              IF (N.EQ.1) RETURN
              DAL2 = DAL*DAL
              DBE2 = DBE*DBE
              DA(2) = (DBE2-DAL2)/ ((DALPBE+2.D0)* (DALPBE+4.D0))
              DB(2) = 4.D0* (DAL+1.D0)* (DBE+1.D0)/
     +                ((DALPBE+3.D0)* (DALPBE+2.D0)**2)
              IF (N.EQ.2) RETURN
              DO 70 K = 3,N
                  DKM1 = DBLE(K-1)
                  DA(K) = .25D0* (DBE2-DAL2)/
     +                    (DKM1*DKM1* (1.D0+.5D0*DALPBE/DKM1)*
     +                    (1.D0+.5D0* (DALPBE+2.D0)/DKM1))
                  DB(K) = .25D0* (1.D0+DAL/DKM1)* (1.D0+DBE/DKM1)*
     +                    (1.D0+DALPBE/DKM1)/ ((1.D0+.5D0* (DALPBE+
     +                    1.D0)/DKM1)* (1.D0+.5D0* (DALPBE-1.D0)/DKM1)*
     +                    (1.D0+.5D0*DALPBE/DKM1)**2)
   70         CONTINUE
              RETURN

          END IF

      ELSE IF (IPOLY.EQ.7) THEN
          IF (DAL.LE.-1.D0) THEN
              IDERR = 1
              RETURN

          ELSE
              DA(1) = DAL + 1.D0
              DB(1) = DGAMMA(DAL+1.D0,IDERR)
              IF (IDERR.EQ.2) DB(1) = D1MACH(2)
              IF (N.EQ.1) RETURN
              DO 80 K = 2,N
                  DKM1 = DBLE(K-1)
                  DA(K) = 2.D0*DKM1 + DAL + 1.D0
                  DB(K) = DKM1* (DKM1+DAL)
   80         CONTINUE
              RETURN

          END IF

      ELSE IF (IPOLY.EQ.8) THEN
          DB(1) = SQRT(4.D0*ATAN(1.D0))
          IF (N.EQ.1) RETURN
          DO 90 K = 2,N
              DB(K) = .5D0*DBLE(K-1)
   90     CONTINUE
          RETURN

      ELSE
          IDERR = 4
      END IF

      END

      DOUBLE PRECISION FUNCTION DLGA(DX)
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DC,DP,DS,DT,DY
      REAL DPREC,Y,Y0
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DBDEN(8),DBNUM(8)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,EXP,LOG,LOG10,REAL
C     ..
C     .. Data statements ..
c
c This routine evaluates the logarithm of the gamma function by a
c combination of recurrence and asymptotic approximation.
c
c The entries in the next data statement are the numerators and
c denominators, respectively, of the quantities B[16]/(16*15),
c B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
c numbers.
c
      DATA DBNUM/-3.617D3,1.D0,-6.91D2,1.D0,-1.D0,1.D0,-1.D0,1.D0/,
     +     DBDEN/1.224D5,1.56D2,3.6036D5,1.188D3,1.68D3,1.26D3,3.6D2,
     +     1.2D1/
C     ..
c
c The quantity  dprec  in the next statement is the number of decimal
c digits carried in double-precision floating-point arithmetic.
c
      DPREC = -LOG10(REAL(D1MACH(3)))
      DC = .5D0*LOG(8.D0*ATAN(1.D0))
      DP = 1.D0
      DY = DX
      Y = REAL(DY)
c
c The quantity  y0  below is the threshold value beyond which asymptotic
c evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
c and I.A. Stegun,Handbook of Mathematical Functions''. The constants
c are .12118868... = ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
c
      Y0 = EXP(.121189*DPREC+.053905)
   10 IF (Y.GT.Y0) GO TO 20
      DP = DY*DP
      DY = DY + 1.D0
      Y = REAL(DY)
      GO TO 10

   20 DT = 1.D0/ (DY*DY)
c
c The right-hand side of the next assignment statement is B[18]/(18*17).
c
      DS = 4.3867D4/2.44188D5
      DO 30 I = 1,8
          DS = DT*DS + DBNUM(I)/DBDEN(I)
   30 CONTINUE
      DLGA = (DY-.5D0)*LOG(DY) - DY + DC + DS/DY - LOG(DP)
      RETURN

      END

      DOUBLE PRECISION FUNCTION DGAMMA(DX,IDERR)
c
c This evaluates the gamma function for real positive  dx, using the
c function subroutine  dlga.
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
      INTEGER IDERR
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DLMACH,DT
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,DLGA
      EXTERNAL D1MACH,DLGA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,LOG
C     ..
      DLMACH = LOG(D1MACH(2))
      IDERR = 0
      DT = DLGA(DX)
      IF (DT.GE.DLMACH) THEN
          IDERR = 2
          DGAMMA = D1MACH(2)
          RETURN

      ELSE
          DGAMMA = EXP(DT)
          RETURN

      END IF

      END
c
c
      SUBROUTINE DSTI(N,NCAP,DX,DW,DALPHA,DBETA,IERR,DP0,DP1,DP2)
c
c This is a double-precision version of the routine  sti.
c
C     .. Scalar Arguments ..
      INTEGER IERR,N,NCAP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DALPHA(N),DBETA(N),DP0(NCAP),DP1(NCAP),DP2(NCAP),
     +                 DW(NCAP),DX(NCAP)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DHUGE,DSUM0,DSUM1,DSUM2,DT,DTINY
      INTEGER K,M,NM1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      DTINY = 10.D0*D1MACH(1)
      DHUGE = .1D0*D1MACH(2)
      IERR = 0
      IF (N.LE.0 .OR. N.GT.NCAP) THEN
          IERR = 1
          RETURN

      END IF

      NM1 = N - 1
      DSUM0 = 0.D0
      DSUM1 = 0.D0
      DO 10 M = 1,NCAP
          DSUM0 = DSUM0 + DW(M)
          DSUM1 = DSUM1 + DW(M)*DX(M)
   10 CONTINUE
      DALPHA(1) = DSUM1/DSUM0
      DBETA(1) = DSUM0
      IF (N.EQ.1) RETURN
      DO 20 M = 1,NCAP
          DP1(M) = 0.D0
          DP2(M) = 1.D0
   20 CONTINUE
      DO 40 K = 1,NM1
          DSUM1 = 0.D0
          DSUM2 = 0.D0
          DO 30 M = 1,NCAP
              IF (DW(M).EQ.0.D0) GO TO 30
              DP0(M) = DP1(M)
              DP1(M) = DP2(M)
              DP2(M) = (DX(M)-DALPHA(K))*DP1(M) - DBETA(K)*DP0(M)
              IF (ABS(DP2(M)).GT.DHUGE .OR. ABS(DSUM2).GT.DHUGE) THEN
                  IERR = K
                  RETURN

              END IF

              DT = DW(M)*DP2(M)*DP2(M)
              DSUM1 = DSUM1 + DT
              DSUM2 = DSUM2 + DT*DX(M)
   30     CONTINUE
          IF (ABS(DSUM1).LT.DTINY) THEN
              IERR = -K
              RETURN

          END IF

          DALPHA(K+1) = DSUM2/DSUM1
          DBETA(K+1) = DSUM1/DSUM0
          DSUM0 = DSUM1
   40 CONTINUE
      RETURN

      END
c
c
      INTEGER FUNCTION NU0HER(N,Z,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Hermite measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL EPS
      INTEGER N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,LOG,REAL,SQRT
C     ..
      NU0HER = 2.* (SQRT(.5*REAL(N+1))+.25*LOG(1./EPS)/ABS(AIMAG(Z)))**2
      RETURN

      END
c
c
      INTEGER FUNCTION NU0JAC(N,Z,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Jacobi measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL EPS
      INTEGER N
C     ..
C     .. Local Scalars ..
      REAL ANGLE,PI,R,X,X2,Y,Y2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,ATAN,COS,LOG,REAL,SIN,SQRT
C     ..
      PI = 4.*ATAN(1.)
      X = REAL(Z)
      Y = ABS(AIMAG(Z))
      IF (X.LT.1.) THEN
          IF (X.LT.-1.) ANGLE = .5* (2.*PI+ATAN(Y/ (X-1.))+
     +                          ATAN(Y/ (X+1.)))
          IF (X.EQ.-1.) ANGLE = .5* (1.5*PI-ATAN(.5*Y))
          IF (X.GT.-1.) ANGLE = .5* (PI+ATAN(Y/ (X-1.))+ATAN(Y/ (X+1.)))

      ELSE
          IF (X.EQ.1.) ANGLE = .5* (.5*PI+ATAN(.5*Y))
          IF (X.GT.1.) ANGLE = .5* (ATAN(Y/ (X-1.))+ATAN(Y/ (X+1.)))
      END IF

      X2 = X*X
      Y2 = Y*Y
      R = ((X2-Y2-1.)**2+4.*X2*Y2)**.25
      R = SQRT((X+R*COS(ANGLE))**2+ (Y+R*SIN(ANGLE))**2)
      NU0JAC = REAL(N+1) + .5*LOG(1./EPS)/LOG(R)
      RETURN

      END
c
c
      INTEGER FUNCTION NU0LAG(N,Z,AL,EPS)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Laguerre measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
C     .. Scalar Arguments ..
      COMPLEX Z
      REAL AL,EPS
      INTEGER N
C     ..
C     .. Local Scalars ..
      REAL PHI,PI,X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AIMAG,ATAN,COS,LOG,REAL,SQRT
C     ..
      PI = 4.*ATAN(1.)
      X = REAL(Z)
      Y = AIMAG(Z)
      PHI = .5*PI
      IF (Y.LT.0.) PHI = 1.5*PI
      IF (X.EQ.0.) GO TO 10
      PHI = ATAN(Y/X)
      IF (Y.GT.0. .AND. X.GT.0.) GO TO 10
      PHI = PHI + PI
      IF (X.LT.0.) GO TO 10
      PHI = PHI + PI
   10 NU0LAG = (SQRT(REAL(N+1)+.5* (AL+1.))+
     +         LOG(1./EPS)/ (4.* (X*X+Y*Y)**.25*COS(.5* (PHI-PI))))**2 -
     +         .5* (AL+1.)
      RETURN

      END

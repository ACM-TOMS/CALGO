C  PROGRAM COUPDR

C     **********
C     MARCH 1, 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
C     VERSION 1.2
C     **********

      PROGRAM COUPDR

C     This program is for the purpose of indicating how SLCOUP
C     can be called for the purpose of obtaining eigenvalues
C     of Sturm-Liouville problems with regular and singular
C     coupled boundary conditions.
C
C     The call is of the form:
C
C     CALL SLCOUP(A, B, INTAB, P0ATA, QFATA, P0ATB, QFATB,
C    1    A1, A2, B1, B2, NUMEIG, EIG, TOL, IFLAG,
C    2    CPFUN, NCA, NCB, ALFA, K11, K12, K21, K22)
C
C     DECLARE THE NEEDED VARIABLES:

C
C     DECLARE THE NEEDED EXTERNALS:

C
C     SET THE PARAMETERS DEFINING THE INTERVAL OF
C       DEFINITION OF THE DIFFERENTIAL EQUATION:
C       (BESSEL'S EQUATION WITH NU = 0.75)
C
C     .. Scalars in Common ..
      REAL A1S,A2S,AA,ASAV,B1S,B2S,BB,BSAV,DTHDAA,DTHDBB,
     +                 EIGSAV,EPSMIN,FA,FB,GQA,GQB,
     +                 GWA,GWB,HPI,LPQA,
     +                 LPQB,LPWA,LPWB,P0ATAS,P0ATBS,PI,QFATAS,QFATBS,
     +                 TMID,TSAVEL,TSAVER,TWOPI,Z
      INTEGER IND,INTSAV,ISAVE,MDTHZ,MFS,MLS,MMWD,T21
      LOGICAL ADDD,PR
C     ..
C     .. Arrays in Common ..
      REAL TEE(100),TT(7,2),YS(200),YY(7,3,2),ZEE(100)
      INTEGER JAY(100),MMW(100),NT(2)
C     ..
C     .. Local Scalars ..
      REAL A,A1,A2,ALFA,B,B1,B2,EIG,K11,K12,K21,K22,P0ATA,
     +                 P0ATB,QFATA,QFATB,TOL
      INTEGER I,ICPFUN,IFLAG,INTAB,NCA,NCB,NP,NUMEIG
C     ..
C     .. Local Arrays ..
      REAL CPFUN(100),X(100)
C     ..
C     .. External Subroutines ..
      EXTERNAL SLCOUP
C     ..
C     .. Common blocks ..
C     COMMON /ALBE/LPWA,LPQA,LPWB,LPQB
      COMMON /ALBE/LPWA,LPQA,FA,GWA,GQA,LPWB,LPQB,FB,GWB,GQB
      COMMON /BCDATA/A1S,A2S,P0ATAS,QFATAS,B1S,B2S,P0ATBS,QFATBS
      COMMON /DATADT/ASAV,BSAV,INTSAV
      COMMON /DATAF/EIGSAV,IND
      COMMON /LP/MFS,MLS
      COMMON /PASS/YS,MMW,MMWD
      COMMON /PIE/PI,TWOPI,HPI
      COMMON /PRIN/PR,T21
      COMMON /RNDOFF/EPSMIN
      COMMON /TDATA/AA,TMID,BB,DTHDAA,DTHDBB,MDTHZ,ADDD
      COMMON /TEEZ/TEE
      COMMON /TEMP/TT,YY,NT
      COMMON /TSAVE/TSAVEL,TSAVER,ISAVE
      COMMON /Z1/Z
      COMMON /ZEEZ/JAY,ZEE
C     ..
      T21 = 21
      OPEN (T21,FILE='test.out')
      PR = .true.
      A = 0.0D0
      B = 1.0D0

C  ----------------------------------------------------------------C
C     SET THE PARAMETERS DEFINING THE COUPLED BOUNDARY CONDITIONS:
C  ----------------------------------------------------------------C

      A1 = 1.0D0
      A2 = 0.0D0
      B1 = 1.0D0
      B2 = 0.0D0

C     SET THE REMAINING PARAMETERS NEEDED BY SLCOUP:

      P0ATA = -1.0D0
      QFATA = -1.0D0
      P0ATB = -1.0D0
      QFATB = 1.0D0
      INTAB = 1
      NUMEIG = 3
      EIG = 0.0D0
      TOL = 1.D-5
      NCA = 3
      NCB = 1
      ICPFUN = 0

C  ----------------------------------------------------------------C
C     SET THE COUPLED BOUNDARY CONDITION PARAMETERS:
C     (THESE ARE THE "PERIODIC" CONDITIONS.)
C  ----------------------------------------------------------------C
      ALFA = 0.0D0
      K11 = 1.0D0
      K12 = 0.0D0
      K21 = 0.0D0
      K22 = 1.0D0

C ---------------------------------------------------------C
C  If eigenfunction values are wanted, then:

      NP = 42
      DO 10 I = 2,NP - 1
          X(I-1) = A + (I-1)* (B-A)/ (NP-1)
   10 CONTINUE
      ICPFUN = NP - 2

      DO 15 I = 1,ICPFUN
          CPFUN(I) = X(I)
   15 CONTINUE
C  ----------------------------------------------------------------C
C     CALL SLCOUP:
C  ----------------------------------------------------------------C

      CALL SLCOUP(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +            EIG,TOL,IFLAG,ICPFUN,CPFUN,NCA,NCB,ALFA,K11,K12,K21,
     +            K22)

C     WRITE OUT THE RETURNED EIGENVALUE, IT'S ESTIMATED ACCURACY,
C     AND THE IFLAG STATUS:

      WRITE (*,FMT=*) ' NUMEIG, EIG, TOL, IFLAG = ',NUMEIG,EIG,TOL,IFLAG

C
      IF (ICPFUN.GT.0) THEN
          WRITE (*,FMT=*) ' EIGENFUNCTION '
          DO 30 I = 1,ICPFUN
              WRITE (*,FMT=*) X(I),CPFUN(I)
   30     CONTINUE
      END IF
C
      CLOSE (T21)
      STOP
      END
C
      REAL FUNCTION P(X)

C     .. Scalar Arguments ..
      REAL X
C     ..
      P = 1.0D0
      RETURN
      END
C
      REAL FUNCTION Q(X)

C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Local Scalars ..
      REAL NU
C     ..
      NU = 0.75D0
      Q = (NU*NU-0.25D0)/X**2
      RETURN
      END
C
      REAL FUNCTION W(X)

C     .. Scalar Arguments ..
      REAL X
C     ..
      W = 1.0D0
      RETURN
      END
C
      SUBROUTINE UV(X,U,PUP,V,PVP,HU,HV)

C     .. Scalar Arguments ..
      REAL HU,HV,PUP,PVP,U,V,X
C     ..
C     .. Local Scalars ..
      REAL NU
C     ..
      NU = 0.75D0
      U = X** (NU+0.5D0)
      PUP = (NU+0.5D0)*X** (NU-0.5D0)
      V = X** (-NU+0.5D0)
      PVP = (-NU+0.5D0)*X** (-NU-0.5D0)
      HU = 0.0D0
      HV = 0.0D0
      RETURN
      END

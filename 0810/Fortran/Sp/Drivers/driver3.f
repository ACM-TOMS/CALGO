C  PROGRAM SEPDR

C     **********
C     MARCH 1, 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
C     VERSION 1.2
C     **********

      PROGRAM SEPDR

C     This program is for the purpose of indicating how SLEIGN2
C     can be called for the purpose of obtaining eigenvalues
C     of Sturm-Liouville problems with regular and singular
C     separated boundary conditions.
C
C     The call is of the form:
C
C     CALL SLEIGN(A, B, INTAB, P0ATA, QFATA, P0ATB, QFATB,
C    1    A1, A2, B1, B2, NUMEIG, EIG, TOL, IFLAG,
C    2    ISLFUN, SLFUN, NCA, NCB)
C
C     SET THE PARAMETERS DEFINING THE INTERVAL OF
C       DEFINITION OF THE DIFFERENTIAL EQUATION:
C       (BESSEL'S EQUATION WITH NU = 0.75)
C
C     .. Scalars in Common ..
      INTEGER T21
      LOGICAL PR
C     ..
C     .. Local Scalars ..
      REAL A,A1,A2,B,B1,B2,EIG,P0ATA,P0ATB,QFATA,QFATB,TOL
      INTEGER I,IFLAG,INTAB,ISLFUN,NCA,NCB,NP,NUMEIG
C     ..
C     .. Local Arrays ..
      REAL SLFUN(100),X(100)
C     ..
C     .. External Subroutines ..
      EXTERNAL SLEIGN
C     ..
C     .. Common blocks ..
      COMMON /PRIN/PR,T21
C     ..
      T21 = 21
      OPEN (T21,FILE='test.out')
      PR = .TRUE.
C
      A = 0.0D0
      B = 1.0D0

C ------------------------------------------------------------------
C     SET THE PARAMETERS DEFINING THE SEPARATED BOUNDARY CONDITIONS:
C ------------------------------------------------------------------

      A1 = 1.0D0
      A2 = 0.0D0
      B1 = 1.0D0
      B2 = 0.0D0

C     SET THE REMAINING PARAMETERS NEEDED BY SLEIGN2:

      P0ATA = -1.0D0
      QFATA = -1.0D0
      P0ATB = -1.0D0
      QFATB = 1.0D0
      INTAB = 1
      NUMEIG = 2
      EIG = 0.0D0
      TOL = 1.D-5
      ISLFUN = 0
      NCA = 3
      NCB = 1

C ---------------------------------------------------------
C  If eigenfunction values are wanted, then:

      NP = 22
      DO 10 I = 2,NP - 1
          X(I-1) = A + (I-1)* (B-A)/ (NP-1)
   10 CONTINUE
      ISLFUN = NP - 2

      DO 15 I = 1,ISLFUN
          SLFUN(9+I) = X(I)
   15 CONTINUE
C ---------------------------------------------------------

C     NOW CALL SLEIGN2:

      CALL SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,NUMEIG,
     +            EIG,TOL,IFLAG,ISLFUN,SLFUN,NCA,NCB)

C     WRITE OUT THE RETURNED EIGENVALUE, IT'S ESTIMATED ACCURACY,
C     AND THE IFLAG STATUS:

      WRITE (*,FMT=*) ' NUMEIG, EIG, TOL, IFLAG = ',NUMEIG,EIG,TOL,IFLAG

C ---------------------------------------------------------
      IF (ISLFUN.GT.0) THEN
C     IN THIS CASE, THE EIGENFUNCTION WAS COMPUTED :
C
          WRITE (*,FMT=*) ' EIGENFUNCTION = '
          DO 25 I = 1,ISLFUN
              WRITE (*,FMT=*) X(I),SLFUN(9+I)
   25     CONTINUE
      END IF

C ---------------------------------------------------------
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

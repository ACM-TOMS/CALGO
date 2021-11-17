C
C STENMIN MINIMIZES AN UNCONSTRAINED NONLINEAR FUNCTION IN N
C UNKNOWNS WHERE THE HESSIAN IS LARGE AND SPARSE, USING TENSOR
C METHODS.
C
C EXAMPLE OF USE FOR STENMIN.  THE TEST PROBLEM IS
C THE BROYDEN TRIDIAGONAL (SOURCE:  BUCKLEY#78 (P. 42)).
C
C  ALI BOUARICHA, OCTOBER 1994.
C  MCS DIVISION, ARGONNE NATIONAL LAB.
C
      INTEGER          NMAX, N, NZ, LIRN, LICN, ILIM, IPR, METHOD
      INTEGER          GRDFLG, HSNFLG, NDIGIT, MSG, LWRK, LIWRK
      INTEGER          TERMCD, INFORM, I
      REAL             FSCALE, GRADTL, STEPTL, FPLS, STEPMX, ONE
      PARAMETER      ( NMAX = 10000, LIRN = 50000, LICN = 500000 )
      PARAMETER      ( LIWRK = 2 * LIRN + 12 * NMAX + 2 )
      PARAMETER      ( LWRK = 7 * NMAX )
      INTEGER          IRN ( LIRN  ), ICN ( LICN )
      INTEGER          IWRK( LIWRK )
      REAL             X   ( NMAX ),  TYPX( NMAX ), XPLS( NMAX )
      REAL             GPLS( NMAX ),  HESS( LICN ), WRK ( LWRK )
      REAL             HTV ( NMAX )
      EXTERNAL         FCN, UGRAD, UHESS
      DATA ONE / 1.0E0 /

C RUN BROYDEN PROBLEM WITH N = 10000.

      N = 10000

C  COMPUTE THE STANDARD STARTING POINT.

      DO 10 I = 1, N
         X(I) = -ONE
 10   CONTINUE

C SET THE DEFAULT VALUES.

      CALL STDFLT(N, TYPX, FSCALE, GRADTL, STEPTL, ILIM, STEPMX,
     *            IPR, METHOD, GRDFLG, HSNFLG, NDIGIT, INFORM, MSG)

      GRADTL = 1.0E-3
      GRDFLG = 2
      HSNFLG = 2

C  CALL THE SPARSE OPTIMIZER.

      CALL STUMCS(N,X,NZ,IRN,LIRN,ICN,LICN,FCN,UGRAD,
     *  UHESS,TYPX,FSCALE,GRADTL,STEPTL,ILIM,STEPMX,IPR,
     *  METHOD,GRDFLG,HSNFLG,NDIGIT,MSG,XPLS,FPLS,GPLS,
     *  HESS,WRK,LWRK,IWRK,LIWRK,TERMCD,HTV,INFORM)

      STOP
      END

C
C THE FOLLOWING IS A SUBROUTINE FOR THE BROYDEN TRIDIAGONAL
C PROBLEM (SOURCE:  BUCKLEY#78 (P. 42))
C
        SUBROUTINE FCN(N, X, F)
        INTEGER N
        REAL X(N), F
C
C  LOCAL VARIABLES
C
        INTEGER I
        REAL ONE, TWO, THREE
        DATA ONE, TWO, THREE / 1.0E0, 2.0E0, 3.0E0 /
C
        F = ((THREE - TWO * X(1)) * X(1) - TWO * X(2) + ONE) ** 2
        DO 10 I = 2, N-1
           F = F + ((THREE - TWO * X(I)) * X(I) - X(I-1) -
     *         TWO * X(I+1) + ONE) ** 2
 10     CONTINUE
        F = F + ((THREE - TWO * X(N)) * X(N) - X(N-1) + ONE) ** 2
        RETURN
        END
C
C THE FOLLOWING IS A SUBROUTINE FOR THE GRADIENT OF THE BROYDEN
C TRIDIAGONAL PROBLEM
C
        SUBROUTINE UGRAD(N, X, G)
        INTEGER N
        REAL X(N), G(N)
C
C LOCAL VARIABLES
C
        INTEGER I
        REAL RL, RM, RR, ONE, TWO, THREE, FOUR
        DATA ONE, TWO, THREE, FOUR/ 1.0E0, 2.0E0, 3.0E0, 4.0E0/
C
        RL = (THREE - TWO * X(1)) * X(1) - TWO * X(2) + ONE
        RR = (THREE - TWO * X(2)) * X(2) - X(1) - TWO * X(3) + ONE
        G(1) = TWO * (RL * (THREE - FOUR * X(1)) - RR)
        DO 10 I = 2, N-1
           IF(I .NE. 2) THEN
              RL = (THREE - TWO * X(I-1)) * X(I-1) - X(I-2) -
     *              TWO * X(I) + ONE
           ENDIF
           RM = (THREE - TWO * X(I)) * X(I) - X(I-1) -
     *            TWO * X(I+1) + ONE
           IF(I .EQ. N-1) THEN
              RR = (THREE - TWO * X(N)) * X(N) - X(N-1) + ONE
           ELSE
              RR = (THREE - TWO * X(I+1)) * X(I+1) - X(I) -
     *              TWO * X(I+2) + ONE
           ENDIF
        G(I) = -TWO * (TWO * RL - RM * (THREE - FOUR * X(I)) + RR)
 10     CONTINUE
        G(N) = -TWO * (TWO * RM - RR * (THREE - FOUR * X(N)))
        RETURN
        END
C
C THE FOLLOWING IS A SUBROUTINE FOR THE HESSIAN OF THE BROYDEN
C TRIDIAGONAL PROBLEM
C
        SUBROUTINE UHESS(N,X,NZ,LICN,HESS,IRN,ICN)
        INTEGER N, NZ, LICN
        INTEGER IRN(*), ICN(LICN)
        REAL X(N), HESS(LICN)
C
C LOCAL VARIABLES
C
        INTEGER I
        REAL RL,RM,RR
        REAL ONE,TWO,THREE,FOUR
        DATA ONE, TWO, THREE, FOUR/1.0E0, 2.0E0, 3.0E0, 4.0E0/
C
        NZ = 1
        RL = (THREE - TWO * X(1)) * X(1) - TWO * X(2) + ONE
        HESS(NZ) = TWO * ((THREE - FOUR * X(1))**2 -
     *            FOUR * RL + ONE)
        IRN(NZ) = 1
        ICN(NZ) = 1
        DO 10 I = 2, N-1
           IF(I .NE. 2) THEN
              NZ = NZ + 1
              HESS(NZ) = FOUR
              IRN(NZ) = I
              ICN(NZ) = I-2
           ENDIF
           NZ = NZ + 1
           HESS(NZ) = -TWO * (TWO * (THREE - FOUR * X(I-1)) +
     *                ONE * (THREE - FOUR * X(I)))
           IRN(NZ) = I
           ICN(NZ) = I-1
           RM = (THREE - TWO * X(I)) * X(I) - X(I-1) -
     *          TWO * X(I+1) + ONE
           NZ = NZ + 1
           HESS(NZ) = -TWO * (-FOUR - (THREE - FOUR * X(I))**2 +
     *                FOUR * RM  - ONE)
           IRN(NZ) = I
           ICN(NZ) = I
 10     CONTINUE
        RR = (THREE - TWO * X(N)) * X(N) - X(N-1) + ONE
        NZ = NZ + 1
        HESS(NZ) = FOUR
        IRN(NZ) = N
        ICN(NZ) = N-2
        NZ = NZ + 1
        HESS(NZ) = -TWO * (TWO * (THREE - FOUR * X(N-1)) +
     *              THREE - FOUR * X(N))
        IRN(NZ) = N
        ICN(NZ) = N-1
        NZ = NZ + 1
        HESS(NZ) = TWO * (FOUR + (THREE - FOUR * X(N))**2 - FOUR * RR)
        IRN(NZ) = N
        ICN(NZ) = N
        RETURN
        END



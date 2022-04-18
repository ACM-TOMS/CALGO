C-----------------------------------------------------------------------
C CONSTANTS:
        INTEGER                 KCOL
C                               KCOL IS THE NUMBER OF COLLOCATION POINTS
C                               TO BE USED IN EACH SUBINTERVAL, WHICH IS
C                               EQUAL TO THE DEGREE OF THE PIECEWISE
C                               POLYNOMIALS MINUS ONE.
C                               1 < KCOL < 11.
        PARAMETER              (KCOL = 2)
        INTEGER                 NPDE
C                               NUMBER OF PDES
        PARAMETER              (NPDE = 1)
        INTEGER                 NINTMX
C                               MAXIMAL NUMBER OF INTEVALS ALLOWED
        PARAMETER              (NINTMX = 2000)
        INTEGER                 MAXVEC
C                               THE DIMENSION OF THE VECTOR OF
C                               BSPLINE COEFFICIENTS
        PARAMETER              (MAXVEC = NPDE*(NINTMX*KCOL+2))
        INTEGER                 LRP
C                               SEE THE COMMENT FOR RPAR
        PARAMETER              (LRP=74+24*NPDE*NPDE*NINTMX*KCOL+8*
     +                              NPDE*NPDE*NINTMX*KCOL*KCOL+29*
     +                              NPDE*NINTMX*KCOL+61*NPDE+14*KCOL
     +                              +35*NINTMX*KCOL+35*NINTMX+21*NINTMX
     +                              *NPDE+8*NINTMX*KCOL*KCOL+37*NPDE
     +                              *NPDE+12*NPDE*NPDE*NINTMX)
C
        INTEGER                 LIP
C                               SEE THE COMMENT FOR IPAR
        PARAMETER              (LIP =100+3*NPDE*(NINTMX*(2*KCOL+1)+4))
        INTEGER                 LCP
C                               SEE THE COMMENT FOR CPAR
        PARAMETER              (LCP =NPDE*(4+(2*KCOL+1)*NINTMX)+NPDE
     +                               *NPDE*(8+NINTMX*(2*KCOL*(KCOL+3)
     +                               +3)))
        INTEGER                 LENWRK
C                               THE DIMENSION OF ARRAY WORK WHEN WE
C                               CALL VALUES
        PARAMETER              (LENWRK =(KCOL+2)+KCOL*(NINTMX+1)+4)
        INTEGER                 NUOUT
C                               THE DIMENSION OF UOUT
        PARAMETER              (NUOUT = NPDE*101) 
        DOUBLE PRECISION        XA
C                               THE LEFT BOUNDARY POINT
        PARAMETER              (XA = 0.0D0)
        DOUBLE PRECISION        XB
C                               THE RIGHT BOUNDARY POINT
        PARAMETER              (XB = 1.0D0)
C-----------------------------------------------------------------------

        DOUBLE PRECISION        T0
C                               T0 < TOUT IS THE INITIAL TIME.
C
        DOUBLE PRECISION        TOUT
C                               TOUT IS THE DESIRED FINAL OUTPUT TIME.
C
        DOUBLE PRECISION        ATOL(NPDE)
C                               ATOL IS THE ABSOLUTE ERROR TOLERANCE
C                               REQUEST AND IS A SCALAR QUANTITY IF
C                               MFLAG(2) = 0.
C
        DOUBLE PRECISION        RTOL(NPDE)
C                               RTOL IS THE RELATIVE ERROR TOLERANCE
C                               REQUEST AND IS A SCALAR QUANTITY IF
C                               MFLAG(2) = 0.
C
        INTEGER                 NINT
C                               NINT IS THE NUMBER OF SUBINTERVALS
C                               DEFINED BY THE SPATIAL MESH X. 
C
        DOUBLE PRECISION        X(NINTMX+1)
C                               X IS THE SPATIAL MESH WHICH DIVIDES THE
C                               INTERVAL [X_A,X_B] AS: X_A = X(1) < 
C                               X(2) < X(3) < ... < X(NINT+1) = X_B.
C
        INTEGER                 MFLAG(6)
C                               THIS VECTOR OF USER INPUT DETERMINES
C                               THE INTERACTION OF BACOL WITH DASSL.
C
C       WORK STORAGE:
        DOUBLE PRECISION        RPAR(LRP)
C                               RPAR IS A FLOATING POINT WORK ARRAY
C                               OF SIZE LRP.
C
        INTEGER                 IPAR(LIP)
C                               IPAR IS AN INTEGER WORK ARRAY
C                               OF SIZE LIP.
C
        DOUBLE COMPLEX          CPAR(LCP)
C                               CPAR IS A COMPLEX WORK ARRAY
C                               OF SIZE LCP.
C
C-----------------------------------------------------------------------
        DOUBLE PRECISION        Y(MAXVEC)
C                               ON SUCCESSFUL RETURN FROM BACOL, Y IS 
C                               THE VECTOR OF BSPLINE
C                               COEFFICIENTS AT THE CURRENT TIME T0.
C
        INTEGER                 IDID
C                               IDID IS THE BACOL EXIT STATUS FLAG
C                               WHICH IS BASED ON THE EXIT STATUS FROM
C                               DASSL ON ERROR CHECKING PERFORMED BY
C                               BACOL ON INITIALIZATION.
C-----------------------------------------------------------------------
        DOUBLE PRECISION        EXACTU(NPDE)
C                               EXACT SOLUTION AT CERTAIN POINT
        DOUBLE PRECISION        UOUT(NUOUT)
C                               THE APPROXIMATION SOLUTIONS AT A SET
C                               OF POINTS
        DOUBLE PRECISION        VALWRK(LENWRK)
C                               VALWRK IS A WORK ARRAY IN VALUES
        DOUBLE PRECISION        XOUT(101)
C                               XOUT IS A SET OF SPATIAL POINTS FOR
C                               OUTPUT
        DOUBLE PRECISION        COEFF
C                               COEFF IS THE COEFFCIENT OF UXX IN THE
C                               BURGERS' EQUATION
        COMMON /BURGER/         COEFF
        INTEGER                 I,J
C-----------------------------------------------------------------------
C SUBROUTINES CALLED:
C                               BACOLR
C                               VALUES
C                               TRUU
C-----------------------------------------------------------------------

C     SET THE REMAINING INPUT PARAMETERS.
      T0 = 0.0D0
      TOUT = 0.2D0
      ATOL(1) = 1.D-6
      RTOL(1) = ATOL(1)
      NINT = 10
      COEFF = 1.D-3

C     DEFINE THE MESH BASED ON A UNIFORM STEP SIZE.
      X(1) = XA
      DO 10 I = 2, NINT
         X(I) = XA + ((I-1) * (XB - XA)) / NINT
   10 CONTINUE
      X(NINT+1) = XB

C     INITIALIZE THE MFLAG VECTOR.
      DO 20 I = 1, 6
         MFLAG(I) = 0
   20 CONTINUE
      MFLAG(5) = 1

      WRITE(6,'(/A)') 'THE INPUT IS  ' 
      WRITE(6,'(/A, I3, A, I4, 2(A, E8.2))') 'KCOL =', KCOL, ', NINT =',
     &  NINT, ', ATOL(1) =', ATOL(1), ', RTOL(1) =', RTOL(1)

      XOUT(1) = XA
      DO 30 I = 2, 100
         XOUT(I) = XA + DBLE(I - 1) * (XB - XA)/100.D0
   30 CONTINUE
      XOUT(101) = XB

      DO 50 J = 1, 5
         WRITE(6, '(/A, E8.2, A, E8.2)') 'EPS =',COEFF,',  TOUT = ',TOUT

         CALL BACOLR(T0, TOUT, ATOL, RTOL, NPDE, KCOL, NINTMX, NINT, X,
     &              MFLAG, RPAR, LRP, IPAR, LIP, CPAR, LCP, Y, IDID)

C        CHECK FOR AN ERROR FROM BACOLR.
         WRITE(6,'(/A, I5)') 'IDID =', IDID
         IF (IDID .LT. 1) GOTO 100

         CALL VALUES(KCOL, XOUT, NINT, X, NPDE, 101, 0, UOUT, Y, VALWRK)

         WRITE(6,'(/A)') 'THE OUTPUT IS  ' 
         WRITE(6,'(/A, I3, A, I4)') 'KCOL =', KCOL, ', NINT =', NINT       
         WRITE(6,'(/A, 2A)') '         XOUT       ', 
     &                       '       UOUT       ', '      EXACTU'
         DO 40 I = 1, 101
            CALL TRUU(TOUT, XOUT(I), EXACTU, NPDE) 
            WRITE(6, '(/3E18.4)') XOUT(I), UOUT(I), EXACTU(1)
   40    CONTINUE

         IDID = 1
         MFLAG(1) = 1
         TOUT = TOUT + 0.2D0
            
   50 CONTINUE
         
      GOTO 9999

  100 CONTINUE
      WRITE(6,'(A)') 'CANNOT PROCEED DUE TO ERROR FROM BACOLR.'

 9999 STOP
      END
C-----------------------------------------------------------------------
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO DEFINE THE INFORMATION ABOUT THE 
C       PDE REQUIRED TO FORM THE ANALYTIC JACOBIAN MATRIX FOR THE DAE
C       OR ODE SYSTEM. ASSUMING THE PDE IS OF THE FORM
C                        UT = F(T, X, U, UX, UXX)
C       THIS ROUTINE RETURNS THE JACOBIANS D(F)/D(U), D(F)/D(UX), AND
C       D(F)/D(UXX).
C
C       COMPATIBLE WITH MSCPDE AND PDECOL/EPDCOL.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DFDU(NPDE,NPDE)
C                               DFDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUX(NPDE,NPDE)
C                               DFDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUXX(NPDE,NPDE)
C                               DFDUXX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SECOND SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DFDU(1,1) = -UX(1)
      DFDUX(1,1) = -U(1)
      DFDUXX(1,1) = COEFF
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C
C       COMPATIBLE ONLY WITH MSCPDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C 
C COMMON VARIABLES:
        DOUBLE PRECISION COEFF
        COMMON /BURGER/ COEFF
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      A1 = (0.5D0 - 4.95D0 * T) * 0.5D-1 / COEFF
      A2 = (0.5D0 - 0.75D0 * T) * 0.25D0 / COEFF
      A3 = 0.1875D0 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C     
      DBDU(1,1) = 1.0D0
      DBDUX(1,1) = 0.0D0
      DBDT(1) = -(0.24D-1 *EXPA1*EXPA2+0.22275D0 *EXPA1*EXPA3
     *           +0.9375D-1 * EXPA2*EXPA3)/
     *           COEFF/(EXPA1+EXPA2+EXPA3)**2
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DIFBXB(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE RIGHT SPATIAL END POINT 1 = XB. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C
C       COMPATIBLE ONLY WITH MSCPDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C
C COMMON VARIABLES:
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      A1 = (-0.5D0 - 4.95D0 * T) * 0.5D-1 / COEFF
      A2 = (-0.5D0 - 0.75D0 * T) * 0.25D0 / COEFF
      A3 = - 0.3125D0 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DBDU(1,1) = 1.0D0
      DBDUX(1,1) = 0.0D0
      DBDT(1) = -(0.24D-1 *EXPA1*EXPA2+0.22275D0 *EXPA1*EXPA3
     *           +0.9375D-1 * EXPA2*EXPA3)/
     *           COEFF/(EXPA1+EXPA2+EXPA3)**2
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       LEFT SPATIAL END POINT X = XA.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XA).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XA).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE LEFT BOUNDARY POINT.
C-----------------------------------------------------------------------
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------
C
      A1 = (0.5D0 - 4.95D0 * T) * 0.5D-1 / COEFF
      A2 = (0.5D0 - 0.75D0 * T) * 0.25D0 / COEFF
      A3 = 0.1875D0 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)

      BVAL(1) = U(1)-(0.1D0*EXPA1+0.5D0*EXPA2+EXPA3)/(EXPA1+EXPA2+EXPA3)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       RIGHT SPATIAL END POINT X = XB.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XB).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XB).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE RIGHT BOUNDARY POINT.
C-----------------------------------------------------------------------
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------
      A1 = (-0.5D0 - 4.95D0 * T) * 0.5D-1 / COEFF
      A2 = (-0.5D0 - 0.75D0 * T) * 0.25D0 / COEFF
      A3 = - 0.3125D0 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)

      BVAL(1) = U(1)-(0.1D0*EXPA1+0.5D0*EXPA2+EXPA3)/(EXPA1+EXPA2+EXPA3)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C
C       COMPATIBLE WITH MSCPDE AND PDECOL/EPDCOL.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
C-----------------------------------------------------------------------
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      FVAL(1) = COEFF*UXX(1) - U(1)*UX(1)
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TRUU(T, X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION PROVIDES THE EXACT SOLUTION OF THE PDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
C OUTPUT:
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE EXACT SOLUTION AT THE 
C                               POINT (T,X).
C-----------------------------------------------------------------------
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------
      A1 = (-X + 0.5D0 - 4.95D0 * T) * 0.5D-1 / COEFF
      A2 = (-X + 0.5D0 - 0.75D0 * T) * 0.25D0 / COEFF
      A3 = (-X + 0.375D0) * 0.5 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)

      U(1) = (0.1D0*EXPA1+0.5D0*EXPA2+EXPA3)/(EXPA1+EXPA2+EXPA3)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UINIT(X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL 
C       CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0 
C       AT THE SPATIAL COORDINATE X.
C
C       COMPATIBLE WITH MSCPDE AND PDECOL/EPDCOL.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        DOUBLE PRECISION        X
C                               THE SPATIAL COORDINATE.
C
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
C OUTPUT:
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS VECTOR OF INITIAL VALUES OF
C                               THE UNKNOWN FUNCTION AT T = T0 AND THE
C                               GIVEN VALUE OF X.
C-----------------------------------------------------------------------
      DOUBLE PRECISION COEFF
      COMMON /BURGER/ COEFF
C-----------------------------------------------------------------------
C
C LOCAL VARIABLES
        DOUBLE PRECISION A1, A2, A3
        DOUBLE PRECISION EXPA1, EXPA2, EXPA3
        DOUBLE PRECISION TEMP
C-----------------------------------------------------------------------
C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
C
C-----------------------------------------------------------------------
      A1 = (-X + 0.5D0) * 0.5D-1 / COEFF
      A2 = (-X + 0.5D0) * 0.25D0 / COEFF
      A3 = (-X + 0.375D0) * 0.5 / COEFF
      EXPA1 = 0.D0
      EXPA2 = 0.D0
      EXPA3 = 0.D0
      TEMP = MAX(A1, A2, A3)
      IF ((A1-TEMP) .GE. -35.D0) EXPA1 = EXP(A1-TEMP)
      IF ((A2-TEMP) .GE. -35.D0) EXPA2 = EXP(A2-TEMP)
      IF ((A3-TEMP) .GE. -35.D0) EXPA3 = EXP(A3-TEMP)

      U(1) = (0.1D0*EXPA1+0.5D0*EXPA2+EXPA3)/(EXPA1+EXPA2+EXPA3)
C
      RETURN
      END

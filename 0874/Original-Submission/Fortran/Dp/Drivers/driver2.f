C-----------------------------------------------------------------------
C CONSTANTS:
        INTEGER                 KCOL
C                               KCOL IS THE NUMBER OF COLLOCATION POINTS
C                               TO BE USED IN EACH SUBINTERVAL, WHICH IS
C                               EQUAL TO THE DEGREE OF THE PIECEWISE
C                               POLYNOMIALS MINUS ONE.
C                               1 < KCOL < 11.
        PARAMETER              (KCOL = 5)
        INTEGER                 NPDE
C                               NUMBER OF PDES
        PARAMETER              (NPDE = 4)
        INTEGER                 NINTMX
C                               MAXIMAL NUMBER OF INTEVALS ALLOWED
        PARAMETER              (NINTMX = 3000)
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
        PARAMETER              (XA = -30.0D0)
        DOUBLE PRECISION        XB
C                               THE RIGHT BOUNDARY POINT
        PARAMETER              (XB = 90.0D0)
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
        DOUBLE PRECISION        UOUT(NUOUT)
C                               THE APPROXIMATION SOLUTIONS AT A SET
C                               OF POINTS
        DOUBLE PRECISION        VALWRK(LENWRK)
C                               VALWRK IS A WORK ARRAY IN VALUES
        DOUBLE PRECISION        XOUT(101)
C                               XOUT IS A SET OF SPATIAL POINTS FOR
C                               OUTPUT
        DOUBLE PRECISION        EXACTU(NPDE)
        INTEGER                 I,II
C-----------------------------------------------------------------------
C SUBROUTINES CALLED:
C                               BACOLR
C                               VALUES
C-----------------------------------------------------------------------

C     SET THE REMAINING INPUT PARAMETERS.
      T0 = 0.0D0
      TOUT = 1.0D+0
      ATOL(1) = 1.D-6
      RTOL(1) = ATOL(1)
      NINT = 10

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

      WRITE(6,'(/A)') 'THE INPUT IS  '
      WRITE(6,'(/A, I3, A, I4, 2(A, E8.2))') 'KCOL =', KCOL, ', NINT =',
     &  NINT, ', ATOL(1) =', ATOL(1), ', RTOL(1) =', RTOL(1)
      WRITE(6,'(/A, E8.2)') 'TOUT = ', TOUT

      CALL BACOLR(T0, TOUT, ATOL, RTOL, NPDE, KCOL, NINTMX, NINT, X,
     &           MFLAG, RPAR, LRP, IPAR, LIP, CPAR, LCP, Y, IDID)

C     CHECK FOR AN ERROR FROM BACOL.
      WRITE(6,'(/A, I5)') 'IDID =', IDID
      IF (IDID .LT. 1) GOTO 100

      XOUT(1) = XA
      DO 30 I = 2, 100
         XOUT(I) = XA + DBLE(I - 1) * (XB - XA)/100.D0
   30 CONTINUE
      XOUT(101) = XB

      CALL VALUES(KCOL, XOUT, NINT, X, NPDE, 101, 0, UOUT, Y, VALWRK)
      
      WRITE(6,'(/A)') 'THE OUTPUT IS  '
      WRITE(6,'(/A, I3, A, I4)') 'KCOL =', KCOL, ', NINT =', NINT
      WRITE(6,'(/A, 4A)') '     XOUT    ', '     UOUT(1)    ',
     &                    '   EXACTU(1)  ', '   UOUT(2)    ',
     &                    '   EXACTU(2)'
      DO 40 I = 1, 101
         II = (I - 1) * 4
         CALL TRUU(T0, XOUT(I), EXACTU, NPDE)
         WRITE(6, '(/5E14.4)') XOUT(I), UOUT(II+1), EXACTU(1),
     &                        UOUT(II+2), EXACTU(2)
   40 CONTINUE
      WRITE(6,'(/A, 4A)') '     XOUT    ', '     UOUT(3)    ',
     &                    '   EXACTU(3)  ', '   UOUT(4)    ',
     &                    '   EXACTU(4)'
      DO 50 I = 1, 101
         II = (I - 1) * 4
         CALL TRUU(T0, XOUT(I), EXACTU, NPDE)
         WRITE(6, '(/5E14.4)') XOUT(I), UOUT(II+3), EXACTU(3),
     &                        UOUT(II+4), EXACTU(4)
   50 CONTINUE

      GOTO 9999

  100 CONTINUE
      WRITE(6,'(A)') 'CANNOT PROCEED DUE TO ERROR FROM BACOL.'

 9999 STOP
      END
C-----------------------------------------------------------------------
c     this is a problem in Ismail and Taha (2001)
c     u1 --- real part of the first PDE component
c     u2 --- imaginary part of the first PDE component
c     u3 --- real part of the second PDE component
c     u4 --- imaginary part of the second PDE component 
c
c-----------------------------------------------------------------------
      subroutine derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx, npde)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is used to define the information about the 
c       PDE required to form the analytic Jacobian matrix for the DAE
c       or ODE system. Assuming the PDE is of the form
c                        ut = f(t, x, u, ux, uxx)
c       this routine returns the Jacobians d(f)/d(u), d(f)/d(ux), and
c       d(f)/d(uxx).
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c Input:
        integer                 npde
c                               The number of PDEs in the system.
c
        double precision        t
c                               The current time coordinate.
c
        double precision        x
c                               The current spatial coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
        double precision        uxx(npde)
c                               uxx(1:npde) is the approximation of the
c                               second spatial derivative of the 
c                               solution at the point (t,x).
c
c Output:
        double precision        dfdu(npde,npde)
c                               dfdu(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dfdux(npde,npde)
c                               dfdux(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the spatial derivative of the 
c                               unknown function u.
c
        double precision        dfduxx(npde,npde)
c                               dfduxx(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the second spatial derivative of the 
c                               unknown function u.
c-----------------------------------------------------------------------

        double precision        epsi, epsi2
c-----------------------------------------------------------------------
c
c     Assign dfdu(1:npde,1:npde), dfdux(1:npde,1:npde), and
c     dfduxx(1:npde,1:npde) according to the right hand side of the PDE 
c     in terms of u(1:npde), ux(1:npde), uxx(1:npde).
c
      epsi = 2.d0/3.d0
      epsi2 = 2.d0 * epsi

      dfdu(1,1) = - 2.d0 * u(1) * u(2)
      dfdu(1,2) = - u(1) * u(1) - 3.d0 * u(2) * u(2) - epsi 
     &            * (u(3) * u(3) + u(4) * u(4))
      dfdu(1,3) = - epsi2 * u(3) * u(2)
      dfdu(1,4) = - epsi2 * u(4) * u(2)
      dfdu(2,1) = 3.d0 * u(1) * u(1) + u(2) * u(2) + epsi
     &            * (u(3) * u(3) + u(4) * u(4))
      dfdu(2,2) = 2.d0 * u(2) * u(1)
      dfdu(2,3) = epsi2 * u(3) * u(1)
      dfdu(2,4) = epsi2 * u(4) * u(1)
      dfdu(3,1) = - epsi2 * u(1) * u(4) 
      dfdu(3,2) = - epsi2 * u(2) * u(4)
      dfdu(3,3) = - 2.d0 * u(3) * u(4) 
      dfdu(3,4) = - 3.d0 * u(4) * u(4) - u(3) * u(3) - epsi
     &            * (u(1) * u(1) + u(2) * u(2))
      dfdu(4,1) = epsi2 * u(1) * u(3) 
      dfdu(4,2) = epsi2 * u(1) * u(3) 
      dfdu(4,3) = 3.d0 * u(3) * u(3) + u(4) * u(4) + epsi
     &            * (u(1) * u(1) + u(2) * u(2))
      dfdu(4,4) = 2.d0 * u(4) * u(3)
c
      dfdux(1,1) = - 0.5d0
      dfdux(1,2) = 0.d0
      dfdux(1,3) = 0.d0
      dfdux(1,4) = 0.d0
      dfdux(2,1) = 0.d0
      dfdux(2,2) = - 0.5d0
      dfdux(2,3) = 0.d0
      dfdux(2,4) = 0.d0
      dfdux(3,1) = 0.d0
      dfdux(3,2) = 0.d0
      dfdux(3,3) = 0.5d0
      dfdux(3,4) = 0.d0
      dfdux(4,1) = 0.d0
      dfdux(4,2) = 0.d0
      dfdux(4,3) = 0.d0
      dfdux(4,4) = 0.5d0
c
      dfduxx(1,1) = 0.d0
      dfduxx(1,2) = - 0.5d0
      dfduxx(1,3) = 0.d0
      dfduxx(1,4) = 0.d0
      dfduxx(2,1) = 0.5d0
      dfduxx(2,2) = 0.d0
      dfduxx(2,3) = 0.d0
      dfduxx(2,4) = 0.d0
      dfduxx(3,1) = 0.d0
      dfduxx(3,2) = 0.d0
      dfduxx(3,3) = 0.d0
      dfduxx(3,4) = - 0.5d0
      dfduxx(4,1) = 0.d0
      dfduxx(4,2) = 0.d0
      dfduxx(4,3) = 0.5d0
      dfduxx(4,4) = 0.d0
c
      return
      end
c-----------------------------------------------------------------------
      subroutine difbxa(t, u, ux, dbdu, dbdux, dbdt, npde)
c-----------------------------------------------------------------------
c Purpose:
c       The subroutine is used to define the differentiated boundary 
c       conditions at the left spatial end point x = xa. For the 
c       boundary condition equation
c                              b(t, u, ux) = 0
c       the partial derivatives db/du, db/dux, and db/dt are supplied 
c       by this routine.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c Input:
        integer                 npde
c                               The number of PDEs in the system.
c
        double precision        t
c                               The current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
c Output:
        double precision        dbdu(npde,npde)
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dbdux(npde,npde)
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the 
c                               unknown function u.
c
        double precision        dbdt(npde)
c                               dbdt(i) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to time t.
c-----------------------------------------------------------------------
c
c     Assign dbdu(1:npde,1:npde), dbdu(1:npde,1:npde), and dbdt(1:npde)
c     according to the right boundary condition equation in terms of 
c     u(1:npde), ux(1:npde), uxx(1:npde).
c
      dbdu(1,1) = 0.d0
      dbdu(1,2) = 0.d0
      dbdu(1,3) = 0.d0
      dbdu(1,4) = 0.d0
      dbdu(2,1) = 0.d0
      dbdu(2,2) = 0.d0
      dbdu(2,3) = 0.d0
      dbdu(2,4) = 0.d0
      dbdu(3,1) = 0.d0
      dbdu(3,2) = 0.d0
      dbdu(3,3) = 0.d0
      dbdu(3,4) = 0.d0
      dbdu(4,1) = 0.d0
      dbdu(4,2) = 0.d0
      dbdu(4,3) = 0.d0
      dbdu(4,4) = 0.d0
c
      dbdux(1,1) = 1.d0
      dbdux(1,2) = 0.d0
      dbdux(1,3) = 0.d0
      dbdux(1,4) = 0.d0
      dbdux(2,1) = 0.d0
      dbdux(2,2) = 1.d0
      dbdux(2,3) = 0.d0
      dbdux(2,4) = 0.d0
      dbdux(3,1) = 0.d0
      dbdux(3,2) = 0.d0
      dbdux(3,3) = 1.d0
      dbdux(3,4) = 0.d0
      dbdux(4,1) = 0.d0
      dbdux(4,2) = 0.d0
      dbdux(4,3) = 0.d0
      dbdux(4,4) = 1.d0
c
      dbdt(1) = 0.d0
      dbdt(2) = 0.d0
      dbdt(3) = 0.d0
      dbdt(4) = 0.d0
c
      return
      end
c-----------------------------------------------------------------------
      subroutine difbxb(t, u, ux, dbdu, dbdux, dbdt, npde)
c-----------------------------------------------------------------------
c Purpose:
c       The subroutine is used to define the differentiated boundary 
c       conditions at the right spatial end point 1 = xb. For the 
c       boundary condition equation
c                              b(t, u, ux) = 0
c       the partial derivatives db/du, db/dux, and db/dt are supplied 
c       by this routine.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c Input:
        integer                 npde
c                               The number of PDEs in the system.
c
        double precision        t
c                               The current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
c Output:
        double precision        dbdu(npde,npde)
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dbdux(npde,npde)
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the 
c                               unknown function u.
c
        double precision        dbdt(npde)
c                               dbdt(i) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to time t.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Assign dbdu(1:npde,1:npde), dbdu(1:npde,1:npde), and dbdt(1:npde)
c     according to the right boundary condition equation in terms of 
c     u(1:npde), ux(1:npde), uxx(1:npde).
c
      dbdu(1,1) = 0.d0
      dbdu(1,2) = 0.d0
      dbdu(1,3) = 0.d0
      dbdu(1,4) = 0.d0
      dbdu(2,1) = 0.d0
      dbdu(2,2) = 0.d0
      dbdu(2,3) = 0.d0
      dbdu(2,4) = 0.d0
      dbdu(3,1) = 0.d0
      dbdu(3,2) = 0.d0
      dbdu(3,3) = 0.d0
      dbdu(3,4) = 0.d0
      dbdu(4,1) = 0.d0
      dbdu(4,2) = 0.d0
      dbdu(4,3) = 0.d0
      dbdu(4,4) = 0.d0
c
      dbdux(1,1) = 1.d0
      dbdux(1,2) = 0.d0
      dbdux(1,3) = 0.d0
      dbdux(1,4) = 0.d0
      dbdux(2,1) = 0.d0
      dbdux(2,2) = 1.d0
      dbdux(2,3) = 0.d0
      dbdux(2,4) = 0.d0
      dbdux(3,1) = 0.d0
      dbdux(3,2) = 0.d0
      dbdux(3,3) = 1.d0
      dbdux(3,4) = 0.d0
      dbdux(4,1) = 0.d0
      dbdux(4,2) = 0.d0
      dbdux(4,3) = 0.d0
      dbdux(4,4) = 1.d0
c
      dbdt(1) = 0.d0
      dbdt(2) = 0.d0
      dbdt(3) = 0.d0
      dbdt(4) = 0.d0
c
      return
      end
c-----------------------------------------------------------------------
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
      BVAL(1) = UX(1) 
      BVAL(2) = UX(2) 
      BVAL(3) = UX(3)
      BVAL(4) = UX(4)
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
      BVAL(1) = UX(1)
      BVAL(2) = UX(2)
      BVAL(3) = UX(3)
      BVAL(4) = UX(4)
C
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine f(t, x, u, ux, uxx, fval, npde)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine defines the right hand side vector of the 
c       npde dimensional parabolic partial differential equation 
c                        ut = f(t, x, u, ux, uxx). 
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c Input:
        integer                 npde
c                               The number of PDEs in the system.
c
        double precision        t
c                               The current time coordinate.
c
        double precision        x
c                               The current spatial coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
        double precision        uxx(npde)
c                               uxx(1:npde) is the approximation of the
c                               second spatial derivative of the 
c                               solution at the point (t,x).
c
c Output:
        double precision        fval(npde)
c                               fval(1:npde) is the right hand side
c                               vector f(t, x, u, ux, uxx) of the PDE.
c-----------------------------------------------------------------------
c
c     Assign fval(1:npde) according to the right hand side of the PDE 
c     in terms of u(1:npde), ux(1:npde), uxx(1:npde).
c
      fval(1) = - 0.5d0 * ux(1) - 0.5d0 * uxx(2) - u(2)
     &          * ((u(1) * u(1) + u(2) * u(2)) + 2.d0/3.d0
     &          * ((u(3) * u(3) + u(4) * u(4))))
      fval(2) = - 0.5d0 * ux(2) + 0.5d0 * uxx(1) + u(1)
     &          * ((u(1) * u(1) + u(2) * u(2)) + 2.d0/3.d0
     &          * ((u(3) * u(3) + u(4) * u(4))))
      fval(3) = 0.5d0 * ux(3) - 0.5d0 * uxx(4) - u(4)
     &          * ((u(3) * u(3) + u(4) * u(4)) + 2.d0/3.d0
     &          * ((u(1) * u(1) + u(2) * u(2))))
      fval(4) = 0.5d0 * ux(4) + 0.5d0 * uxx(3) + u(3)
     &          * ((u(3) * u(3) + u(4) * u(4)) + 2.d0/3.d0
     &          * ((u(1) * u(1) + u(2) * u(2))))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine uinit(x, u, npde)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is used to return the npde-vector of initial 
c       conditions of the unknown function at the initial time t = t0 
c       at the spatial coordinate x.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c Input:
        double precision        x
c                               The spatial coordinate.
c
        integer                 npde
c                               The number of PDEs in the system.
c
c Output:
        double precision        u(npde)
c                               u(1:npde) is vector of initial values of
c                               the unknown function at t = t0 and the
c                               given value of x.
c-----------------------------------------------------------------------
        double precision        tempt1, tempt2

c
c     Assign u(1:npde) the initial values of u(t0,x).
c
      
      tempt1 = sqrt(6.d0/5.d0)
      tempt2 = sqrt(2.d0)

      u(1) = tempt1 / cosh(tempt2 * x) * cos(0.5d0 * x)
      u(2) = tempt1 / cosh(tempt2 * x) * sin(0.5d0 * x)  
      u(3) = tempt1 / cosh(tempt2 * x) * cos(1.5d0 * x)  
      u(4) = tempt1 / cosh(tempt2 * x) * sin(1.5d0 * x)  
c
      return
      end
C-----------------------------------------------------------------------
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
      double precision        tempt1, tempt2

      tempt1 = sqrt(6.d0/5.d0)
      tempt2 = sqrt(2.d0)

      u(1) = tempt1 / cosh(tempt2*(x-t)) * cos(0.5d0*x+0.625d0*t)
      u(2) = tempt1 / cosh(tempt2*(x-t)) * sin(0.5d0*x+0.625d0*t)  
      u(3) = tempt1 / cosh(tempt2*(x-t)) * cos(1.5d0*x+0.625d0*t)
      u(4) = tempt1 / cosh(tempt2*(x-t)) * sin(1.5d0*x+0.625d0*t)  
C
      RETURN
      END


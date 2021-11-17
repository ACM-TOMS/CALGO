C The Reaction-Convection-Diffusion (RCD) system.
C This problem is also referred to as the Catalytic Surface Reaction 
C Model (CSRM). It is a system consisting of four PDES (npde=4).
C See F for PDE, BNDXA, BNDXB for the boundary conditions, and
C UNIT for the initial conditions.
C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
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
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      FVAL(1) = -UX(1)+n*(D1*U(3)-A1*U(1)*(1.d0-U(3)-U(4)))+(1.d0/Pe1)
     *          *UXX(1)

      FVAL(2) = -UX(2)+n*(D2*U(4)-A2*U(2)*(1.d0-U(3)-U(4)))+(1.d0/Pe1)
     *          *UXX(2)

      FVAL(3) = A1*U(1)*(1.d0-U(3)-U(4))-D1*U(3)-R*U(3)*U(4)*(1.d0-
     *          U(3)-U(4))**2+(1.d0/Pe2)*UXX(3)

      FVAL(4) = A2*U(2)*(1.d0-U(3)-U(4))-D2*U(4)-R*U(3)*U(4)*(1.d0-
     *          U(3)-U(4))**2+(1.d0/Pe2)*UXX(4) 
C
      RETURN
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
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DFDU(1,1) = -n*A1*(1-U(3)-U(4))
      DFDU(1,2) = 0.d0
      DFDU(1,3) = n*D1+n*A1*U(1)
      DFDU(1,4) = n*A1*U(1)

      DFDU(2,1) = 0.d0
      DFDU(2,2) = -n*A2+n*A2*U(3)+n*A2*U(4)
      DFDU(2,3) = n*A2*U(2)
      DFDU(2,4) = n*D2+n*A2*U(2)

      DFDU(3,1) = A1-A1*U(3)-A1*U(4)
      DFDU(3,2) = 0.d0
      DFDU(3,3) = -A1*U(1)-D1-R*U(4)+4.d0*R*U(3)*U(4)+2.d0*R*U(4)**2-
     *           4.d0*R*U(3)*U(4)**2-3.d0*R*U(3)**2*U(4)-R*U(4)**3
      DFDU(3,4) = -A1*U(1)-R*U(3)+2.d0*R*U(3)**2+4.d0*R*U(3)*U(4)-
     *           4.d0*R*U(3)**2*U(4)-R*U(3)**3-3.d0*R*U(3)*U(4)**2

      DFDU(4,1) = 0.d0
      DFDU(4,2) = A2-A2*U(3)-A2*U(4)
      DFDU(4,3) = -A2*U(2)-R*U(4)+2.d0*R*U(4)**2+4.d0*R*U(3)*U(4)-
     *           4.d0*R*U(4)**2*U(3)-2.d0*R*U(3)*U(4)-R*U(4)**3
      DFDU(4,4) = -A2*U(2)-D2-R*U(3)+2.D0*R*U(3)**2+4.D0*R*U(3)*U(4)-
     *           4.D0*R*U(3)**2*U(4)-R*U(3)**2-3.D0*R*U(3)*U(4)**2


      DO I=1,4
        DO J=1,4
            DFDUX(I,J) = 0.d0
            DFDUXX(I,J) = 0.d0
        END DO
      END DO

      DFDUX(1,1) = -1.d0
      DFDUX(2,2) = -1.d0

      DFDUXX(1,1) = 1.d0/Pe1
      DFDUXX(2,2) = 1.d0/Pe1
      DFDUXX(3,3) = 1.d0/Pe2
      DFDUXX(4,4) = 1.d0/Pe2
      
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
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
c       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
        

      DO I=1,4
        DO J=1,4
          DBDU(I,J) = 0.d0
          DBDUX(I,J) = 0.d0
        END DO
        DBDT(I) = 0.d0
      END DO

      DBDU(1,1) = -Pe1
      DBDU(2,2) = -Pe1

      DBDUX(1,1) = 1.d0
      DBDUX(2,2) = 1.d0
      DBDUX(3,3) = 1.d0
      DBDUX(4,4) = 1.d0

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
C-----------------------------------------------------------------------
C       LOOP INDICES
        INTEGER I, J
C-----------------------------------------------------------------------
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).

      DO I=1,4
        DO J=1,4
          DBDU(I,J) = 0.d0
          DBDUX(I,J) = 0.d0
        END DO
        DBDT(I) = 0.d0
        DBDUX(I,I) = 1.d0
      END DO
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
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
      
      BVAL(1) = UX(1)+Pe1*(2.d0-C-U(1))
      BVAL(2) = UX(2)+Pe1*(C-U(2))
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
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------

      BVAL(1)= UX(1)
      BVAL(2)= UX(2)
      BVAL(3)= UX(3)
      BVAL(4)= UX(4)

C
      RETURN
      END

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
C Ture solution is unknown
      U(1) = 0.d0
      U(2) = 0.d0
      U(3) = 0.d0
      U(4) = 0.d0
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
C       PROBLEM CONSTANTS
C
        DOUBLE PRECISION        A1,A2,D1,D2,R,C,n,Pe1,Pe2
        COMMON /RCDSTUFF/       A1,A2,D1,D2,R,C,n,Pe1,Pe2

C-----------------------------------------------------------------------
     
      U(1) = 2.d0-C
      U(2) = C
      U(3) = 0.d0
      U(4) = 0.d0
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine header(nout)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine writes a header describing the npde dimensional
c       parabolic partial differential equation
c                        ut = f(t, x, u, ux, uxx).
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 nout
c                               nout is the output unit number.
c-----------------------------------------------------------------------
c constants:
        double precision        t0
        parameter              (t0 = 0.0d0)
c
        double precision        xa
        parameter              (xa = 0.0d0)
c
        double precision        xb
        parameter              (xb = 1.0d0)
c-----------------------------------------------------------------------
c
      write(nout,95) 'The RCD system'
      write(nout,95) 'domain:'
      write(nout,96) '   t0 =', t0, ' < t,'
      write(nout,96) '   xa =', xa, ' <= x <= xb =', xb, ','
c
      return
   95 format(a)
   96 format(a,e13.5,a,e13.5,a,e13.5,a,e13.5,a)
      end



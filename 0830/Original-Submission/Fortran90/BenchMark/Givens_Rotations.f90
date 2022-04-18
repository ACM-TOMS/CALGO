    MODULE givens_rotations
! .. Use Statements ..
      USE import_kinds, ONLY : sp => single_precision, dp => double_precision
! ..
! .. Generic Interface Blocks ..
      INTERFACE rot
        MODULE PROCEDURE s_rot
        MODULE PROCEDURE d_rot
      END INTERFACE
      INTERFACE rotm
        MODULE PROCEDURE s_rotm
        MODULE PROCEDURE d_rotm
      END INTERFACE
! ..
    CONTAINS
! Modified by R. Hanson 24 April 2002.  Use more stable downdating.
! Ref. is from Bjorck, Stewart and Stewart, and Chambers.
! !!!!!!!!!!!!!!!!!!!!!!!!! Documentation !!!!!!!!!!!!!!!!!!!!!
! This module contains a set of routines that efficiently implement
! the calculation and application of a Givens plane rotation using
! both the straightforward and modified methods. 

! These routines are intended to replace the blas 1 routines
! drot, drotm, drotg, drotmg, srot, srotm, srotg and srotmg.

! Six user callable routines have been provided.

! Calls to d_rot are intended to be direct replacements for a pair of
! calls to the blas level 1 routines drotg and drot; i.e., the
! calculation of the rotation parameters and the application of the
! rotation to a pair of rows has been amalgamated in a single call.
! s_rot replaces calls to srotg and srot.

! In a similar fashion d_rotm and s_rotm replace pairs of calls to drotg
! and drotg, and srotm and srotmg for the modified version of the
! algorithm.

! The routine *_drot and *_drotm mirror the original blas routines in
! their calling sequences. The data to be transformed is passed to the
! routines as a pair of one dimensional, assumed size arrays and the incx
! and incy stride parameters are also available. These routines may be
! used in Fortran 77 codes that use an array element to define the start
! of an array.

! In addition we have provided a pair of generic routines rot and rotm
! that provide the functionality of the *_rot and *_rotm routines
! but allow the same name to be used for either double or single
! precision data. These routines use assumed size arrays to store
! the row data. 

! Experiments were conducted using assumed shape arrays designed to
! shorten the calling sequences by allowing the effect of the use of
! strides to be obtained by using a Fortran 90 array slice; for example
! x(1:1+(k-1)*incx:incx). Also the use of assumed shape arrays would
! allow the length of the data matrices to be determined internally using
! the size function.  Unfortunately these codes proved to be extrememly
! inefficient at run time and were abandoned.

! The generic routines are implemented as wrappers to the 
! *_rot and *_rotm routines.

! In line with the original blas routines both single and double
! precision versions of the routines are available. We have used a
! consistent naming convention for arguments throughout the routines.

! !!!!!!!!!!!!!!!!!!!!!!!!! Routine Descriptions !!!!!!!!!!!!!!!!!!!!!!!

!  S_ROT and D_ROT
!  ---------------
!  construct a Givens transformation, i.e, they compute
!        c = cos(theta) and s = sin(theta) such that either

!         ( c  s) ( w1 )   ( r )
!         (     ) (    ) = (   )
!         (-s  c) ( z1 )   ( 0 )

!  or, when dropping data, they compute the hyperbolic
!  transformation matrix of the form

!         ( c is) 
!         (-is c)

!  where i is the imaginary unit and this unit is implicitly applied.

!  When data is dropped the condition, 
!               w1**2-z1**2 = (w1-z1)*(w1+z1) > 0
!  is required.


!  These routines then apply the computed transformation to the data

!         ( x_i )
!         (     )
!         ( y_i )

!  where the ith element of the  vector x is

!     1 + (i-1)*incx,   if incx >= 0,

!  and similarly for y and incy.

!  Note that these routines do not allow negative values of incx
!  and incy, nor do they compute the z value returned by *rotg
!  In addition these routines also allow for row removal which
!  was only available with the modified Givens routines *ROTMG
!  in the original Blas level 1 library.

!  S_ROTM and D_ROTM
!  -----------------
!  construct the modified Givens rotation matrix H which
!  zeros the second component of the 2-vector

!    (w1,z1) = (sqrt(d1)*x1, sqrt(d2)*x2)'

!  where  d1 and d2 are the reciprocals of the scaling factors
!  described in the introduction to the paper and appearing
!  in the original blas level 1 routines.

!  Row removal mathematically requires d2<0 but this is flagged
!  in our code by setting d1<0 for efficiency reasons.

!  When data is dropped the condition, 
!               w1**2-z1**2 = (w1-z1)*(w1+z1) > 0
!  is required.

!  and apply it to the 2 x N 
!         ( x' )
!  matrix (    ), where the elements of the x-vector used are  
!         ( y' )

!         1 + (i-1)*incx,   if incx >= 0,

!  and similarly for y and incy.

!  If the value of 

! NOTE that these routines require the reciprocal of the d values where
! the original blas routines used the d values directly. Also these
! routines do not allow negative values of incx and incy

! ROT and ROTM
! ------------
! are Fortran 90 generic interfaces to the routines *_ROT and *_ROTM
! respectively. These routines allow the same routine name to be used 
! for either single or double precision data.

! In addition, ROTM uses a structure to return the details of the
! transformation matrix used rather than a floating point array.

! !!!!!!!!!!!!!!!!!!!!!!!!! Parameter Descriptions !!!!!!!!!!!!!!!!!!!!!!!
!
! *_ROT, *_ROTM, ROT and ROTM
! ---------------------------
!          (In what follows REAL should be replaced by DOUBLE PRECISION
!           when calling the D_ROT and D_ROTM routines. Either REAL or DOUBLE
!           PRECISION data may be used with the ROT and ROTM routines provided 
!           that all the real parameters are of the same type in each call)

! W1,Z1  - REAL (ROT only)
!          On entry, these define the two vector whose second element is
!          to be annihilated.
!          On exit, W1 is set to the value of r (see above) and Z1 is
!          set to zero.

! K      - INTEGER
!          On entry, specifies the number of elements (columns)
!          in the two data vectors x and y.
!          Unchanged on exit.

! X(*),Y(*)    - REAL assumed size arrays
!          On entry, these contain the row data to be transformed.
!          On exit, unless an error is detected in the input data,
!          these contain the transformed data.

! INCX   - INTEGER
!          On entry, INCX specifies the increment parameter used to step
!          through the array SX. Must be positive.
!          Unchanged on exit.

! INCY   - INTEGER
!          On entry, INCY specifies the increment parameter used to step
!          through the array SY. Must be positive.
!          Unchanged on exit.

! C, S   - REAL (ROT only)
!          On exit, these are elements of the plane rotation or 
!          hyperbolic matrix.

! X1, X2 - REAL (ROTM only)
!          On entry, contains the values that, when combined with the
!          scaling factors RD1 and RD2 (below) give the two vector
!          (W1, Z1) whose second element is to be annihilated.
!          On exit, X1 contains the transformed value and X2 is
!          set to zero.

! RD1, RD2 - REAL (ROTM only)
!          On entry, contain the scaling values used for the modified
!          transformation. These values are the reciprocal of the values
!          (D1 and D2) used in the original level 1 blas routines.
!          On exit, contain the scaling values for the transformed vector.

! PARAM(5)  - REAL array. (ROTM only)
!          On exit, the array param is used to control the form of the 
!          transformation matrix returned by the modified givens
!          algorithm by setting param(1) = sflag where

!             sflag = -1.0e0      sflag = 0.0e0       sflag = 1.0e0

!               (sh11  sh12)        (1.e0 sh12)        ( sh11 1.e0)
!          H =  (          )        (         )        (          )
!               (sh21  sh22)        (sh21 1.e0)        (-1.e0 sh22)

!          In addition if param(1) = -2 then H is set to I.

!          Elements 2--5 of param are used to specify the values of
!          sh11, sh21, sh12, sh22 respectively. Values of 1.0e0, -1.0e0,
!          or 0.0e0, implied by the value of param(1), are stored
!          in param but are not used as multipliers in the transformation 
!          stage.
      SUBROUTINE s_rotm(rd1,rd2,x1,y1,k,x,incx,y,incy,param)
!  Routine to implement the modified Givens transformation
!  This routine is equivalent to the two level 1 blas routines
!  srotm and srotmg.
! .. Parameters ..
        REAL (sp), PARAMETER :: one = 1E0_sp
        REAL (sp), PARAMETER :: quarter = 25E-2_sp
        REAL (sp), PARAMETER :: zero = 0E0_sp
        INTEGER, PARAMETER :: clts = 1, error = -1, rescale = 2, slec = 0
! ..
! .. Scalar Arguments ..
        REAL (sp), INTENT (INOUT) :: rd1, rd2, x1, y1
        INTEGER, INTENT (IN) :: incx, incy, k
! ..
! .. Array Arguments ..
        REAL (sp), INTENT (OUT) :: param(5)
        REAL (sp), INTENT (INOUT) :: x(*), y(*)
! ..
! .. Local Scalars ..
        REAL (sp) :: gam, gamsq, h11, h12, h21, h22, p1, p2, rgam, sx1, sy1, u
        INTEGER :: i, kx, ky, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, huge, min, sqrt, tiny
! ..
! Check input parameters; quick return if illegal
        IF (rd1==zero .OR. rd2<=zero .OR. (k>0 .AND. (incx<=0 .OR. incy<= &
            0))) THEN
          param(1) = error
          RETURN
        END IF
        n = k
        kx = 1
        ky = 1

        sx1 = x1
        sy1 = y1
        h11 = one
        h21 = zero
        h12 = zero
        h22 = one

        p2 = rd1*sy1
        p1 = rd2*sx1

! This value is dependent on the underlying
! floating-point arithmetic and need be computed once.
! To make the code thread-safe it is done every time.
! When rescaling is needed gam and rgam are computed.
!       gamsq = sqrt(min(huge(one),one/tiny(one))*quarter)
        gamsq = 4.611686E+18
! If |c| >= |s|:
        IF (p1*sx1>=abs(p2*sy1)) THEN
          param(1) = slec
! There is no transformation necessary.
          IF (p2==zero) GO TO 10
          h21 = -sy1/sx1
          h12 = p2/p1
          u = one - h12*h21
          sx1 = sx1*u
          rd1 = rd1*u
          rd2 = rd2*u
! Check if rescaling is required on either diagonal factor.
          IF (abs(rd1)+rd2>gamsq) GO TO 20
! Transformation : (1 h12)
!                  (h21 1)
          param(2:5) = (/h11, h21, h12, h22/)
! Equal increments are likely:
          IF (incx==incy) THEN
! Unit strides are usually best:
            IF (incx==1) THEN
              DO i = 1, n
                p1 = x(i) + param(4)*y(i)
      			y(i) = y(i) + param(3)*x(i)
                x(i) = p1
              END DO

            ELSE
! Non-unit but equal strides:
              DO i = 1, n
                p1 = x(kx) + param(4)*y(kx)
                y(kx) = y(kx) + param(3)*x(kx)
                x(kx) = p1
                kx = kx + incx
              END DO
            END IF

          ELSE
! Non-equal strides:
            DO i = 1, n
              p1 = x(kx) + param(4)*y(ky)
              y(ky) = y(ky) + param(3)*x(kx)
              x(kx) = p1

              kx = kx + incx
              ky = ky + incy
            END DO
          END IF

          GO TO 10

        ELSE
          h22 = sx1/sy1
          h11 = p1/p2
          h12 = one
          h21 = -one
          u = one + h11*h22
          p1 = rd2*u
          rd2 = rd1*u
          sx1 = sy1*u
          rd1 = p1
          param(1) = clts

! Check if rescaling is required on either diagonal factor.
          IF (abs(rd1)+rd2>gamsq) GO TO 20

! Transformation : ( h11 1)
!                  (-1 h22)
          param(2:5) = (/h11, h21, h12, h22/)

! Perhaps only an interchange and sign change is needed.
          IF (h11==zero .AND. h22==zero) THEN
            DO i = 1, n
              p1 = y(ky)
              y(ky) = -x(kx)
              x(kx) = p1
              kx = kx + incx
              ky = ky + incy
            END DO
            GO TO 10
          END IF

! Equal increments are likely:
          IF (incx==incy) THEN
! Unit strides are usually best:
            IF (incx==1) THEN
              DO i = 1, n
                p1 = y(i) + param(2)*x(i)
                y(i) = -x(i) + param(5)*y(i)
                x(i) = p1
              END DO
            ELSE
! Non-unit but equal strides:
              DO i = 1, n
                p1 = y(kx) + param(2)*x(kx)
                y(kx) = -x(kx) + param(5)*y(kx)
                x(kx) = p1

                kx = kx + incx
              END DO
            END IF

          ELSE
! Non-equal strides:
            DO i = 1, n
              p1 = y(ky) + param(2)*x(kx)
              y(ky) = -x(kx) + param(5)*y(ky)
              x(kx) = p1
              kx = kx + incx
              ky = ky + incy
            END DO
          END IF

        END IF

10      CONTINUE
! The form of the matrix is defined in param(1) when here.
        param(2:5) = (/h11, h21, h12, h22/)

        y1 = zero
        x1 = sx1
        RETURN

20      CONTINUE
! Compute gam and rgam if they may be needed for rescaling.
        gam = sqrt(gamsq)
        rgam = one/gam
        IF (abs(rd1)>gamsq) THEN
          rd1 = (rd1*rgam)*rgam
          sx1 = sx1*gam
          h11 = h11*rgam
          h12 = h12*rgam
          param(1) = rescale
        END IF

        IF (rd2>gamsq) THEN
          rd2 = (rd2*rgam)*rgam
          h21 = h21*rgam
          h22 = h22*rgam
          param(1) = rescale
        END IF

! Apply scaled transformation.  This is a rare event, after start.
        DO i = 1, n
          p1 = h11*x((i-1)*incx+1) + h12*y((i-1)*incy+1)
          y((i-1)*incy+1) = h21*x((i-1)*incx+1) + h22*y((i-1)*incy+1)
          x((i-1)*incx+1) = p1
        END DO
        param(2:5) = (/h11, h21, h12, h22/)
        GO TO 10

      END SUBROUTINE s_rotm

      SUBROUTINE d_rotm(rd1,rd2,x1,y1,k,x,incx,y,incy,param)
!  Routine to implement the modified Givens transformation
!  This routine is equivalent to the two level 1 blas routines
!  drotm and drotmg.
! .. Parameters ..
        REAL (dp), PARAMETER :: one = 1E0_dp
        REAL (dp), PARAMETER :: quarter = 25E-2_dp
        REAL (dp), PARAMETER :: zero = 0E0_dp
        INTEGER, PARAMETER :: clts = 1, error = -1, rescale = 2, slec = 0
! ..
! .. Scalar Arguments ..
        REAL (dp), INTENT (INOUT) :: rd1, rd2, x1, y1
        INTEGER, INTENT (IN) :: incx, incy, k
! ..
! .. Array Arguments ..
        REAL (dp), INTENT (OUT) :: param(5)
        REAL (dp), INTENT (INOUT) :: x(*), y(*)
! ..
! .. Local Scalars ..
        REAL (dp) :: gam, gamsq, h11, h12, h21, h22, p1, p2, rgam, sx1, sy1, u
        INTEGER :: i, kx, ky, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, huge, min, sqrt, tiny
! ..
! Check input parameters; quick return if illegal
        IF (rd1==zero .OR. rd2<=zero .OR. (k>0 .AND. (incx<=0 .OR. incy<= &
            0))) THEN
          param(1) = error
          RETURN
        END IF
        n = k
        kx = 1
        ky = 1

        sx1 = x1
        sy1 = y1
        h11 = one
        h21 = zero
        h12 = zero
        h22 = one

        p2 = rd1*sy1
        p1 = rd2*sx1

! Set the value of gamsq on first call to the
! routine. This value is dependent on the underlying
! floating-point arithmetic and need be computed once.
! To make the code thread-safe it is done every time.
! When rescaling is needed gam and rgam are computed.
!       gamsq = sqrt(min(huge(one),one/tiny(one))*quarter)
        gamsq = 3.3519519824856493D+153
! If |c| >= |s|:
        IF (p1*sx1>=abs(p2*sy1)) THEN
          param(1) = slec
! There is no transformation necessary.
          IF (p2==zero) GO TO 10
          h21 = -sy1/sx1
          h12 = p2/p1
          u = one - h12*h21
          sx1 = sx1*u
          rd1 = rd1*u
          rd2 = rd2*u
! Check if rescaling is required on either diagonal factor.
          IF (abs(rd1)+rd2>gamsq) GO TO 20
! Transformation : (1 h12)
!                  (h21 1)
          param(2:5) = (/h11, h21, h12, h22/)
! Equal increments are likely:
          IF (incx==incy) THEN
! Unit strides are usually best:
            IF (incx==1) THEN
              DO i = 1, n
                p1 = x(i) + param(4)*y(i)
                y(i) = y(i) + param(3)*x(i)
                x(i) = p1
              END DO

            ELSE
! Non-unit but equal strides:
              DO i = 1, n
                p1 = x(kx) + param(4)*y(kx)
                y(kx) = y(kx) + param(3)*x(kx)
                x(kx) = p1
                kx = kx + incx
              END DO
            END IF

          ELSE
! Non-equal strides:
            DO i = 1, n
              p1 = x(kx) + param(4)*y(ky)
              y(ky) = y(ky) + param(3)*x(kx)
              x(kx) = p1

              kx = kx + incx
              ky = ky + incy
            END DO
          END IF

          GO TO 10

        ELSE
! exit with an error return if rd1 < 0 at this point
          IF (rd1<zero) THEN
            param(1) = error
            RETURN
          END IF
          h22 = sx1/sy1
          h11 = p1/p2
          h12 = one
          h21 = -one
          u = one + h11*h22
          p1 = rd2*u
          rd2 = rd1*u
          sx1 = sy1*u
          rd1 = p1
          param(1) = clts

! Check if rescaling is required on either diagonal factor.
          IF (abs(rd1)+rd2>gamsq) GO TO 20

! Transformation : ( h11 1)
!                  (-1 h22)
          param(2:5) = (/h11, h21, h12, h22/)

! Perhaps only an interchange and sign change is needed.
          IF (h11==zero .AND. h22==zero) THEN
            DO i = 1, n
              p1 = y(ky)
              y(ky) = -x(kx)
              x(kx) = p1
              kx = kx + incx
              ky = ky + incy
            END DO
            GO TO 10
          END IF

! Equal increments are likely:
          IF (incx==incy) THEN
! Unit strides are usually best:
            IF (incx==1) THEN
              DO i = 1, n
                p1 = y(i) + param(2)*x(i)
                y(i) = -x(i) + param(5)*y(i)
                x(i) = p1
              END DO
            ELSE
! Non-unit but equal strides:
              DO i = 1, n
                p1 = y(kx) + param(2)*x(kx)
                y(kx) = -x(kx) + param(5)*y(kx)
                x(kx) = p1

                kx = kx + incx
              END DO
            END IF

          ELSE
! Non-equal strides:
            DO i = 1, n
              p1 = y(ky) + param(2)*x(kx)
              y(ky) = -x(kx) + param(5)*y(ky)
              x(kx) = p1
              kx = kx + incx
              ky = ky + incy
            END DO
          END IF

        END IF

10      CONTINUE
! The form of the matrix is defined in param(1) when here.
        param(2:5) = (/h11, h21, h12, h22/)

        y1 = zero
        x1 = sx1
        RETURN

20      CONTINUE
! Compute gam and rgam if they may be needed for rescaling.
        gam = sqrt(gamsq)
        rgam = one/gam
        IF (abs(rd1)>gamsq) THEN
          rd1 = (rd1*rgam)*rgam
          sx1 = sx1*gam
          h11 = h11*rgam
          h12 = h12*rgam
          param(1) = rescale
        END IF

        IF (rd2>gamsq) THEN
          rd2 = (rd2*rgam)*rgam
          h21 = h21*rgam
          h22 = h22*rgam
          param(1) = rescale
        END IF

! Apply scaled transformation.  This is a rare event, after start.
        DO i = 1, n
          p1 = h11*x((i-1)*incx+1) + h12*y((i-1)*incy+1)
          y((i-1)*incy+1) = h21*x((i-1)*incx+1) + h22*y((i-1)*incy+1)
          x((i-1)*incx+1) = p1
        END DO
        param(2:5) = (/h11, h21, h12, h22/)
        GO TO 10

      END SUBROUTINE d_rotm

      SUBROUTINE s_rot(w1,z1,k,x,incx,y,incy,c,s)
!  Routine to implement the standard Givens transformation
!  This routine is equivalent to the two level 1 blas routines
!  srot and srotg.
! .. Parameters ..
        REAL (sp), PARAMETER :: half = 5E-1_sp
        REAL (sp), PARAMETER :: one = 1E0_sp
        REAL (sp), PARAMETER :: quarter = 25E-2_sp
        REAL (sp), PARAMETER :: zero = 0E0_sp
! ..
! .. Scalar Arguments ..
        REAL (sp), INTENT (OUT) :: c, s
        REAL (sp), INTENT (INOUT) :: w1, z1
        INTEGER, INTENT (IN) :: incx, incy, k
! ..
! .. Array Arguments ..
        REAL (sp), INTENT (INOUT) :: x(*), y(*)
! ..
! .. Local Scalars ..
        REAL (sp) :: a, b, u, v
        INTEGER :: i, kx, ky, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, sqrt
! ..
! Test incx and incy are both positive
        IF (incx<=0 .OR. incy<=0) THEN
          c = zero
          s = zero
          RETURN
        END IF
        n = k
        a = w1
        b = z1
! If transformation is hyperbolic, code will drop data.
        IF (n<0) GO TO 10
        IF (abs(a)>abs(b)) THEN
! Here ABS(w1) > ABS(z1)
          u = a + a
          v = b/u
! Note that u and r have the sign of w1

          a = sqrt(quarter+v*v)
! Note that c is positive
          c = half/a
          s = v*(c+c)
          a = a*u
        ELSE
! Here ABS(w1) <= ABS(z1)
          IF (b/=zero) THEN
            u = b + b
            v = a/u
! Note that u and r have the sign of
! z1 (r is immediately stored in a)
            a = sqrt(quarter+v*v)
! Note that s is positive
            s = half/a
            c = v*(s+s)
            a = a*u
! Here w1 = z1 = 0.
          ELSE
            c = one
            s = zero
            GO TO 20
          END IF
        END IF
! Possibly B=0, so no arithmetic is needed.
        IF (c==one .AND. s==zero) GO TO 20

        kx = 1
        ky = 1
! Equal increments are likely:
        IF (incx==incy) THEN
! Unit strides are usually best:
          IF (incx==1) THEN
            DO i = 1, n
              u = c*x(i) + s*y(i)
              y(i) = -s*x(i) + c*y(i)
              x(i) = u
            END DO
          ELSE
! Non-unit but equal strides:
            DO i = 1, n
              u = c*x(kx) + s*y(kx)
              y(kx) = -s*x(kx) + c*y(kx)
              x(kx) = u

              kx = kx + incx
            END DO
          END IF

        ELSE
! Non-equal strides:
          DO i = 1, n
            u = c*x(kx) + s*y(ky)
            y(ky) = -s*x(kx) + c*y(ky)
            x(kx) = u

            kx = kx + incx
            ky = ky + incy
          END DO
        END IF

        GO TO 20

10      CONTINUE
        n = -(n+1)

! This is an error condition.
        IF(ABS(a) <= ABS(b))  THEN
          c = zero
          s = zero
          a = zero
          RETURN
        ENDIF
! Will have abs(a) > abs(b) here:
        s=b/a
        c=(one-s)*(one+s)
! May have c <= 0, even though mathmatically it should be > 0.
! This is also an error condition.
        IF(c <= zero) THEN
          c = zero
          s = zero
          a = zero
          RETURN
        ENDIF
        c=sqrt(c)
        a=c*a
! Will need reciprocal of c to apply transformation:
        b=one/c

        kx = 1
        ky = 1
! Equal increments are likely:
       IF (incx==incy) THEN
! Unit strides are usually best:
         IF (incx==1) THEN
           DO i = 1, n
             u = b*(x(i) - s*y(i))
             y(i) = -s*u + c*y(i)
             x(i) = u
           END DO
         ELSE
! Non-unit but equal strides:
            DO i = 1, n
              u = b*(x(kx) - s*y(kx))
              y(kx) = -s*u + c*y(kx)
              x(kx) = u
              kx = kx + incx
            END DO
          END IF

        ELSE
! Non-equal strides:
          DO i = 1, n
            u = b*(x(kx) - s*y(ky))
            y(ky) = -s*u + c*y(ky)
            x(kx) = u
            kx = kx + incx
            ky = ky + incy
          END DO
       END IF

20      CONTINUE
        z1 = zero
        w1 = a

      END SUBROUTINE s_rot

      SUBROUTINE d_rot(w1,z1,k,x,incx,y,incy,c,s)
!  Routine to implement the standard Givens transformation
!  This routine is equivalent to the two level 1 blas routines
!  drot and drotg.
! .. Parameters ..
        REAL (dp), PARAMETER :: half = 5E-1_dp
        REAL (dp), PARAMETER :: one = 1E0_dp
        REAL (dp), PARAMETER :: quarter = 25E-2_dp
        REAL (dp), PARAMETER :: zero = 0E0_dp
! ..
! .. Scalar Arguments ..
        REAL (dp), INTENT (OUT) :: c, s
        REAL (dp), INTENT (INOUT) :: w1, z1
        INTEGER, INTENT (IN) :: incx, incy, k
! ..
! .. Array Arguments ..
        REAL (dp), INTENT (INOUT) :: x(*), y(*)
! ..
! .. Local Scalars ..
        REAL (dp) :: a, b, u, v
        INTEGER :: i, kx, ky, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, sqrt
! ..
! Test incx and incy are both positive
        IF (incx<=0 .OR. incy<=0) THEN
          c = zero
          s = zero
          RETURN
        END IF
        n = k
        a = w1
        b = z1
! If transformation is hyperbolic, code will drop data.
        IF (n<0) GO TO 10
        IF (abs(a)>abs(b)) THEN
! Here ABS(w1) > ABS(z1)
          u = a + a
          v = b/u
! Note that u and r have the sign of W1

          a = sqrt(quarter+v*v)
! Note that c is positive
          c = half/a
          s = v*(c+c)
          a = a*u
        ELSE
! Here ABS(w1) <= ABS(z1)
          IF (b/=zero) THEN
            u = b + b
            v = a/u
! Note that u and r have the sign of
! z1 (r is immediately stored in a)
            a = sqrt(quarter+v*v)
! Note that s is positive
            s = half/a
            c = v*(s+s)
            a = a*u
! Here W1 = Z1 = 0.
          ELSE
            c = one
            s = zero
            GO TO 20
          END IF
        END IF
! Possibly B=0, so no arithmetic is needed.
        IF (c==one .AND. s==zero) GO TO 20

        kx = 1
        ky = 1
! Equal increments are likely:
        IF (incx==incy) THEN
! Unit strides are usually best:
          IF (incx==1) THEN
            DO i = 1, n
              u = c*x(i) + s*y(i)
              y(i) = -s*x(i) + c*y(i)
              x(i) = u
            END DO
          ELSE
! Non-unit but equal strides:
            DO i = 1, n
              u = c*x(kx) + s*y(kx)
              y(kx) = -s*x(kx) + c*y(kx)
              x(kx) = u

              kx = kx + incx
            END DO
          END IF

        ELSE
! Non-equal strides:
          DO i = 1, n
            u = c*x(kx) + s*y(ky)
            y(ky) = -s*x(kx) + c*y(ky)
            x(kx) = u
            kx = kx + incx
            ky = ky + incy
          END DO
        END IF

        GO TO 20

10      CONTINUE
        n = -(n+1)

        IF(ABS(a) <= ABS(b))  THEN
          c = zero
          s = zero
          a = zero
          RETURN
        ENDIF
! Will have abs(a) > abs(b) here:
        s=b/a
        c=(one-s)*(one+s)
! May have c <= 0, even though mathmatically it should be > 0.
! This is also an error condition.
        IF(c <= zero) THEN
          c = zero
          s = zero
          a = zero
          RETURN
        ENDIF
        c=sqrt(c)
        a=c*a
! Will need reciprocal of c to apply transformation:
        b=one/c

        kx = 1
        ky = 1

! Equal increments are likely:
       IF (incx==incy) THEN
! Unit strides are usually best:
         IF (incx==1) THEN
           DO i = 1, n
             u = b*(x(i) - s*y(i))
             y(i) = -s*u + c*y(i)
             x(i) = u
           END DO
         ELSE
! Non-unit but equal strides:
            DO i = 1, n
              u = b*(x(kx) - s*y(kx))
              y(kx) = -s*u + c*y(kx)
              x(kx) = u
              kx = kx + incx
            END DO
          END IF

        ELSE
! Non-equal strides:
          DO i = 1, n
            u = b*(x(kx) - s*y(ky))
            y(ky) = -s*u + c*y(ky)
            x(kx) = u
            kx = kx + incx
            ky = ky + incy
          END DO
       END IF

20      CONTINUE
        z1 = zero
        w1 = a

      END SUBROUTINE d_rot


    END MODULE givens_rotations

      program main

c***********************************************************************
c
cc TOMS392_PRB tests TOMS392.
c
c  Modified:
c
c    11 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS392_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 392, which approximates'
      write ( *, '(a)' ) '  the time evolution of a system governed'
      write ( *, '(a)' ) '  by hyperpolic PDE''s.'

      call testch

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS392_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '

      stop
      end
      SUBROUTINE TESTCH

      IMPLICIT NONE

      REAL DATA(4,81)
      REAL FM
      INTEGER I
      INTEGER IFAIL
      INTEGER M
      INTEGER N

      N = 20
C
C  GENERATE INITIAL DATA.
C
      M = 4 * N + 1
      FM = 4.0E+00 * FLOAT ( N )
      DO 100 I = 1, M
        DATA(1,I) = FLOAT ( I - 1 ) / FM
        DATA(2,I) = 0.0E+00
        DATA(3,I) = 0.0E+00
        DATA(4,I) = 2.0E+00 * EXP ( DATA(1,I) )
100   CONTINUE

      IFAIL = 0
      WRITE ( *, 900 )
900   FORMAT ( " " )

      WRITE ( *, '(A)' ) ' '
      WRITE ( *, '(A)' ) 
     &  '       X               Y               U               V'
      WRITE ( *, '(A)' ) ' '

200   DO 250 I = 1, M
        WRITE ( *, 910 ) DATA(1,I), DATA(2,I), DATA(3,I), DATA(4,I)
250   CONTINUE

910   FORMAT ( 2X,G14.6,2X,G14.6,2X,G14.6,2X,G14.6 )

      IF ( M .LE. 1 ) GO TO 300
      IF ( IFAIL .NE. 0 ) GO TO 300

      CALL CHARAC ( DATA, M, IFAIL )
      WRITE ( *, 900 )

      GO TO 200

300   CONTINUE
      WRITE ( *, 920 ) M, IFAIL
920   FORMAT ( 1X, "M =",I2," IFAIL =",I2 )
      RETURN
      END
      SUBROUTINE CHCOEF ( COEFF, XYUV )
C
C  COMPUTES COEFFICIENTS A1, A2, A3, A4, H1, B1, B2, B3, B4, H2,
C  AND STORES THEM SEQUENTIALLY IN COEFF.
C
      IMPLICIT NONE

      REAL COEFF(10)
      REAL XYUV(4)

      COEFF(1) = 1.0E+00 - XYUV(3)**2
      COEFF(2) = -XYUV(3) * XYUV(4)
      COEFF(3) = -XYUV(3) * XYUV(4)
      COEFF(4) = 1.0E+00 - XYUV(4)**2
      COEFF(5) = -4.0E+00 * XYUV(3) * EXP ( XYUV(1) )**2
      COEFF(6) = 0.0E+00
      COEFF(7) = 1.0E+00
      COEFF(8) = -1.0E+00
      COEFF(9) = 0.0E+00
      COEFF(10) = 0.0E+00

      RETURN
      END

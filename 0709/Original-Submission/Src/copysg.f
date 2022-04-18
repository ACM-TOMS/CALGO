      SUBROUTINE  ZZCOPY (N, X, INCX, Y, INCY )

C## A R G U M E N T S:
                      INTEGER INCX, INCY, N

                      REAL             X(*), Y(*)
C!!!!                 DOUBLE PRECISION X(*), Y(*)
C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.

C>RCS $HEADER$
C>RCS $LOG$

C## D E S C R I P T I O N:
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.

C     SEE FURTHER COMMENTS IN ZZAMAX.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZCOPY.
C## S U B R O U T I N E S:   MOD...INTRINSIC.
C## P A R A M E T E R S:     NONE ARE DEFINED.
C## L O C A L   D E C L:
                             INTEGER I, IX, IY, M

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO VALUES ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      IF( N .GT. 0 ) THEN

         IF ( INCX .EQ. 1 .AND. INCY .EQ. 1 ) THEN
C
C           CODE FOR BOTH INCREMENTS EQUAL TO 1
C
            M = MOD(N,7)
            IF ( M .NE. 0 ) THEN
C                               CLEAN-UP LOOP
               DO 30 I = 1,M
                  Y(I) = X(I)
   30          CONTINUE
            ENDIF
            DO 50 I = M+1,N,7
               Y(I    ) = X(I    )
               Y(I + 1) = X(I + 1)
               Y(I + 2) = X(I + 2)
               Y(I + 3) = X(I + 3)
               Y(I + 4) = X(I + 4)
               Y(I + 5) = X(I + 5)
               Y(I + 6) = X(I + 6)
   50       CONTINUE

         ELSE

C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
            IX = 1
            IY = 1
            IF ( INCX .LT. 0 ) IX = (-N+1)*INCX + 1
            IF ( INCY .LT. 0 ) IY = (-N+1)*INCY + 1
            DO 10 I = 1,N
               Y(IY) = X(IX)
               IX = IX + INCX
               IY = IY + INCY
   10       CONTINUE
         ENDIF

      ENDIF

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C                   E N D  OF ZZCOPY.
      END

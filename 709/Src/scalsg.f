      SUBROUTINE  ZZSCAL (N,A,X,INCX)

C## A R G U M E N T S:
                      INTEGER INCX,N

                      REAL             A, X(*)
C!!!!                 DOUBLE PRECISION A, X(*)
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
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.

C     SEE FURTHER COMMENTS IN ZZAMAX.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZSCAL.
C## S U B R O U T I N E S:   MOD...INTRINSIC.
C## P A R A M E T E R S:     NONE ARE DEFINED.
C## L O C A L   D E C L:
                             INTEGER I, M, NINCX

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO VALUES ARE SET.

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      IF ( N .GT. 0 ) THEN

         IF ( INCX .NE. 1 ) THEN
C
C           CODE FOR INCREMENT NOT EQUAL TO 1
C
            NINCX = N*INCX
            DO 10 I = 1,NINCX,INCX
               X(I) = A*X(I)
   10       CONTINUE

         ELSE

C        CODE FOR INCREMENT EQUAL TO 1
C
            M = MOD(N,5)
            IF ( M .NE. 0 ) THEN
C                               CLEAN-UP LOOP
               DO 30 I = 1,M
                  X(I) = A*X(I)
   30          CONTINUE
            ENDIF
            DO 50 I = M+1,N,5
               X(I  ) = A*X(I )
               X(I+1) = A*X(I+1)
               X(I+2) = A*X(I+2)
               X(I+3) = A*X(I+3)
               X(I+4) = A*X(I+4)
   50       CONTINUE
         ENDIF

      ENDIF

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C                   E N D  OF ZZSCAL.
      END

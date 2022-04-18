      SUBROUTINE BBMJDP ( DIAG, U, W, Z, PHI, VV, DELT, GAM, Y,
     - V, HV, N, NUPS, I, M, SCDIAG, IDENTY,  INNER, IW, RW, DW )

C## A R G U M E N T S:
                       INTEGER          N, NUPS, I, M, IW(*)
                       LOGICAL          IDENTY, SCDIAG

      REAL              V(N), HV(N), DIAG(N), W(N), Z(N), U(N)
C!!!! DOUBLE PRECISION  V(N), HV(N), DIAG(N), W(N), Z(N), U(N)

      REAL              PHI(0:M-1), VV(N,0:M-1), DELT(N,0:M-1)
C!!!! DOUBLE PRECISION  PHI(0:M-1), VV(N,0:M-1), DELT(N,0:M-1)

      REAL              Y(N,0:M-1), GAM(N,0:M-1)
C!!!! DOUBLE PRECISION  Y(N,0:M-1), GAM(N,0:M-1)

      EXTERNAL                INNER
      DOUBLE PRECISION DW(*), INNER
      REAL             RW(*)

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.

C>RCS $HEADER: MJDP.F,V 1.12 91/12/16 12:00:37 BUCKLEY EXP $
C>RCS $LOG:     MJDP.F,V $
C>RCS REVISION 1.12  91/12/16  12:00:37  BUCKLEY
C>RCS MINOR FIX FOR TOMS.
C>RCS
C>RCS REVISION 1.11  91/11/22  11:30:55  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.10  90/07/31  10:49:56  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:12:48  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:39:22  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:55:32  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:54:29  BUCKLEY
C>RCS INITIAL REVISION
C>RCS

C## D E S C R I P T I O N:

C     GIVEN THE QUASI NEWTON UPDATE MATRIX  H (IN ZZ' FORM) AND
C     GIVEN THE VECTOR V, THIS ROUTINE COMPUTES
C
C               HV = H * V  .
C
C     IT ALSO RETURNS THE INTERMEDIATE VALUE U = Z' * V.
C
C     IT ASSUMES HERE THAT Z REPRESENTS THE MATRIX Z[I] AT THE
C     POINT X[I].
C
C     THE DESCRIPTION OF THE ARGUMENTS IS GIVEN IN BBPOWL.
C
C------TRACES:  NOTE THAT THE TRACE PARAMETERS WILL BE THE SAME
C     FOR EACH CALL TO BBMJDP DURING ANY PARTICULAR MINIMIZATION
C     PROBLEM, SO THEY ARE SET JUST ONCE THROUGH AN ENTRY POINT.
C
C     THESE WILL BE ON THE UNIT TRU.
C     VECTORS ARE TRACED ONLY IF TRV IS TRUE AS WELL.
C
C## E N T R Y   P O I N T S: BBMJDP  THE NATURAL ENTRY POINT.
C                            BBSMJD  FOR UNCHANGING ARGUMENTS.
C## S U B R O U T I N E S:   INNER (PASSED AS ARGUMENT)
C                            MAX, MIN, SQRT ...INTRINSIC
C## P A R A M E T E R S:

      LOGICAL     DONORM,          NONORM
      PARAMETER ( DONORM = .TRUE., NONORM = .FALSE. )
      REAL              ZERO,       ONE,       TWO,       THREE
C!!!! DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      REAL              FOUR,       FIVE,      SIX,       SEVEN
C!!!! DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      REAL              EIGHT,         NINE,          TEN
C!!!! DOUBLE PRECISION  EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

C## L O C A L   D E C L:
                         INTEGER          K,  J
                         INTEGER          LOWROW, HIROW

                         REAL             ZETA, TMP
C!!!!                    DOUBLE PRECISION ZETA, TMP

C-----                   VARIABLES FOR THE ENTRY POINT.

                         LOGICAL     TR, STR, TRV, STRV
                         INTEGER     TRU, STRACN
C## S A V E:
                         SAVE   TR, TRV, TRU

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.

C##                                               E X E C U T I O N
C##                                               E X E C U T I O N
      IF ( TR ) WRITE (TRU,*) ' ***[MJDP ENTERED]***'
      IF ( TR ) WRITE (TRU,*) ' [MJDP] SCDIAG,N,NUPS,I,M,IDENTY=',
     -                                 SCDIAG,N,NUPS,I,M,IDENTY

      IF ( NUPS .EQ. 0 ) THEN
          DO 100 K = 1,N
              IF (SCDIAG) THEN
                TMP = DIAG(K)
              ELSE
                TMP = ONE
              ENDIF
              U(K)  = TMP*V(K)
              HV(K) = TMP**2 * V(K)
  100     CONTINUE

      ELSE

          IF ( TR .AND. TRV ) THEN
              CALL ZZWMAT( 0,' [MJDP] Y',ONE,Y,ZERO,Y,N,1,N,1,I,
     -                     'F',9,4,1,80,TRU)
          ENDIF
          IF (SCDIAG) THEN
            TMP = Y(N,0)*DIAG(N)
          ELSE
            TMP = Y(N,0)
          ENDIF

          DO 200 K = 1,N
              HV(K) = ZERO
  200     CONTINUE

          IF ( TR .AND. TRV ) THEN
              CALL ZZWMAT(0,' [MJDP] DELT',ONE,DELT,ZERO,DELT,N,1,N,1,I,
     -                    'F',9,4,1,80,TRU)
              CALL ZZWMAT(0,' [MJDP] GAM',ONE,GAM,ZERO,GAM,N,1,N,1,I,
     -                    'F',9,4,1,80,TRU)
          ENDIF
          CALL ZZCOPY( N, DELT(1,0), 1,  VV(1,0), 1 )
          CALL ZZSCAL( N, -TMP*GAM(N,0), VV(1,0), 1 )

          VV(N,0) = TMP + VV(N,0)
          PHI(0) = Y(N,0)**2

          DO 4000 K = N, 3-I, -1
C            K-1 IS STARTING COLUMN ALONG BOTTOM LEVEL OF TABLE

             LOWROW = MAX(0,2-K)
             HIROW  = MIN(I-1,N-K)

             IF(TR) WRITE(TRU,*) '  [MJDP] OUTER LOOP: K=',K
             DO 3000 J = LOWROW, HIROW
C               J IS CURRENT LEVEL IN TABLE, SO THAT
C               J+K-1=K-1+J IS CURRENT COLUMN AND
C               J+K         IS NEXT COLUMN TO RIGHT
C               WE ARE COMPUTING Z FOR ROW J+1

                IF(TR) WRITE(TRU,*) '  [MJDP] INNER LOOP: J=',J

                ZETA = SQRT( PHI(J)/(PHI(J)+Y(K+J-1,J))**2 )

                IF  ( J .EQ. 0 ) THEN
                   IF (SCDIAG) THEN
                       TMP = DIAG(K-1)
                   ELSE
                       TMP = ONE
                   ENDIF

                   CALL ZZCOPY ( N, DELT(1,0), 1,    W, 1 )
                   CALL ZZSCAL ( N, -GAM(K-1,0)*TMP, W, 1 )

                   W(K-1) =  TMP + W(K-1)

                ELSE IF  ( J+K .EQ. 2 ) THEN
C                                        IN COLUMN 1 SINCE J+K-1=1

                   IF (TR) WRITE(TRU,*) ' [MJDP] CREATING Z FROM ',
     -                                  'DELT SUB ', J-1
                   TMP = INNER (N,GAM(1,J),DELT(1,J-1),NONORM,IW,RW,DW)
                   CALL ZZCOPY ( N, DELT(1,J-1), 1,     W, 1 )
                   CALL ZZAXPY ( N, -TMP, DELT(1,J), 1, W, 1 )

                ELSE
                   TMP = INNER (N, GAM(1,J), Z, NONORM,IW,RW,DW)
                   CALL ZZCOPY ( N, Z, 1,               W, 1 )
                   CALL ZZAXPY ( N, -TMP, DELT(1,J), 1, W, 1 )
                ENDIF

                IF(TR) WRITE(TRU,*) '  [MJDP] COMPUTING Z'
                CALL ZZCOPY ( N, W, 1,                          Z, 1 )
                CALL ZZAXPY ( N, Y(K+J-1,J)/PHI(J), VV(1,J), 1, Z, 1 )
                CALL ZZSCAL ( N, ZETA,                          Z, 1 )


                IF ( J+K .GT. 2 ) THEN
                   IF(TR) WRITE(TRU,*) '  [MJDP] UPDATING VV, PHI',J,
     -                                ' WITH Y; LEVEL ',K+J-1
                    CALL ZZAXPY ( N, Y(K+J-1,J), W, 1, VV(1,J), 1 )
                    PHI(J) = PHI(J) + Y(K+J-1,J)**2
                ENDIF

 3000        CONTINUE

             J = HIROW + 1
             IF ( J .LT. I ) THEN
                IF(TR) WRITE(TRU,*) '  [MJDP] INITIALIZING VV',J

                TMP = INNER(N, GAM(1,J), Z, NONORM,IW,RW,DW)
                CALL ZZCOPY ( N, Z, 1,               VV(1,J), 1 )
                CALL ZZAXPY ( N, -TMP, DELT(1,J), 1, VV(1,J), 1 )
                CALL ZZSCAL ( N, Y(N,J),             VV(1,J), 1 )
                PHI(J) = Y(N,J)**2
             ELSE
                IF(TR) WRITE(TRU,*) '  [MJDP] TOP ROW: UPDATING U AT',
     -                               K+I-1
                U(K+I-1) = INNER (N, Z, V, NONORM, IW, RW, DW )
                CALL ZZAXPY ( N, U(K+I-1), Z, 1, HV, 1 )
                IF(TR) THEN
                   IF ( TRV ) THEN
                       CALL ZZWMAT(0,' [MJDP] Z',ONE,Z,ZERO,Z,N,1,N,1,1,
     -                    'F',9,4,1,80,TRU)
                       IF ( K-1+I .LE. I .AND. K-1+I .GE. 2) THEN
                           WRITE(TRU,*) '  [MJDP] CHECKING...I,K,K-1+I,'
     -                                   ,'1-K=' ,I,K,K-1+I,1-K
                           TMP = Z(1)/DELT(1,1-K)
                           CALL ZZWMAT(0,'Z-DELT1-K',ONE,Z,-TMP,
     -                      DELT(1,1-K), N,1,N,1,1,'F',9,4,1,80,TRU)
                       ENDIF
                   ENDIF
                   TMP = INNER(N, GAM(1,J), Z, NONORM,IW,RW,DW)
                   WRITE(TRU,*) ' [MJDP] ZT*GAM[',J,']=',TMP
                ENDIF
             ENDIF

 4000     CONTINUE

          IF(TR) WRITE(TRU,*) ' [MJDP] TOP ROW: UPDATING U AT 1',
     -                         ' WITH DELT SUB ',I-1
          U(1) = INNER( N, DELT(1,I-1), V, NONORM, IW, RW, DW )
          CALL ZZAXPY ( N, U(1), DELT(1,I-1), 1, HV, 1 )

      ENDIF

      GOTO 90000

C## E N T R Y  BBSMJD:
                      ENTRY  BBSMJD ( STR, STRV, STRACN)
      TR  = STR
      TRV = STRV
      TRU = STRACN
                      RETURN
C## E X I T
90000        IF ( TR ) WRITE (TRU,*) ' ===[LEAVING MJDP].'
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF BBMJDP.
                    END


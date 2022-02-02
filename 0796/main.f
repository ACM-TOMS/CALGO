C   DRIVER PROGRAM FOR TESTING SOFTWARE INVLTF  IMPLEMENTING           
C                                                                      
C   FOURIER BASED METHOD FOR THE NUMERICAL INVERSION                   
C                                                                      
C                  OF LAPLACE TRANSFORMS                              
C                                                                    
C                                                                 
C                                                              
C                                                           
C                     JANUARY , 1998      
C                                                                      
C                                                                      
C   AUTHORS: D'AMORE LUISA, GIULIANO LACCETTI, ALMERICO MURLI            
C                                                                      
C                                                                      
C
C
C
C   REFERENCES
C   ==========
C
C
C
C
C
C
C
C
C
C
C   REFERENCES
C   ==========
C
C   D'AMORE L., LACCETTI G., MURLI A.,  -    "ALGORITHM XXX: A FORTRAN
C      SOFTWARE PACKAGE FOR THE NUMERICAL INVERSION OF THE
C      LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"
C
      PROGRAM MAIN
C     ***********************************************************
C     DRIVER  PROGRAM TO TEST  THE ROUTINE INVLTF
C     FOR THE INVERSION OF A LAPLACE TRANSFORM FUNCTION.
C     THIS VERSION USES BOTH REAL AND COMPLEX DOUBLE PRECISION
C     OPERATIONS.
C     THE MAIN PROGRAM ALLOWS THE INVERSION OF A SET OF LAPLACE
C     TRANSFORM FUNCTIONS WHICH CAN BE ADDRESSED  BY A NATURAL
C     NUMBER  BETWEEN  1 AND 34.
C     FOR THE COMPLETE LIST SEE THE TABLE IN THE COMPANION PAPER.
C
C     THIS IS A SELF-CONTAINED DRIVER FOR THE INVLTF ROUTINE
C     COMPRISING A MAIN PROGRAM AND SUBPROGRAMS QDACC, BACKCF,
C     AND THE LAPLACE TRANSFORM PAIRS FZ AND FEX
C
C     THIS DRIVER ALLOWS TO OBTAIN UP TO 50 VALUES OF THE
C     INVERSE  LAPLACE FUNCTION 
C     
C     ************************************************************
C
C
C
C     .. Parameters ..
      INTEGER NMAX,FMAX,LASTFN
      PARAMETER (NMAX=550,FMAX=50,LASTFN=34)
C     ..
C     .. Scalars in Common ..
      INTEGER NFUN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EXF,SIGMA0,TOL,SSBAR
      INTEGER I,NT
      CHARACTER*100 AA,BB
      LOGICAL STOP
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX WORK(2,0:2*NMAX)
      DOUBLE PRECISION DIFABS(FMAX),DIFREL(FMAX),FZINV(FMAX),
     +                 TARRAY(FMAX),ERROR(3,FMAX),ERR(3)
      INTEGER IFAIL(FMAX),IFZEVAL(FMAX)
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX FZ
      DOUBLE PRECISION FEX
      EXTERNAL FZ,FEX
C     ..
C     .. External Subroutines ..
      EXTERNAL INVLTF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /NF/NFUN
C     ..
C     *****************************************************************
C                   SET UP THE OUTPUT
C     *****************************************************************

C
C Loop through all the tests
C
      DO 50 NFUN = 1,LASTFN

C
C Skip the tests that require the use of NAG library routines
C
C Commented out lines need to be reintroduced in the routines
C FZ and FEX after labels 9, 36, 37, 38
C
          IF (NFUN.EQ.9 .OR. NFUN.EQ.36 .OR. NFUN.EQ.37 .OR.
     +        NFUN.EQ.38) GO TO 50
          AA =
     +'*****************************************************************
     +***********************'
          BB =
     +      '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
          WRITE (*,FMT=9000) AA,AA,AA,AA
          WRITE (*,FMT=9090)
C     READ (*,FMT=*) NFUN
          WRITE (*,FMT=9050) NFUN
          WRITE (*,FMT=9110)
C     READ (*,FMT=*) SIGMA0
	  SIGMA0 = 0.0
          WRITE (*,FMT=9120) SIGMA0
	  IF (NFUN .EQ. 18) SIGMA0 = 3.0
          IF (NFUN .EQ. 23) SIGMA0 = 0.25D0
          IF (NFUN .EQ. 29) SIGMA0 = 2.D0
          WRITE (*,FMT=9080) BB,BB,BB,BB
C     *****************************************************************
C          NOW, SET INPUT PARAMETERS FOR INVLTF
C     *****************************************************************
C
C      NT IS THE NUMBER OF T-VALUES USED HERE FOR THE TEST  FUNCTIONS.
C      THE MAXIMUM ALLOWED IS THE DIMENSION OF THE ARRAY TVAL
C     *****************************************************************
C
C     IN THIS TEST PROGRAM THE INVERSE FUNCTION IS REQUESTED IN THE
C     FOLLOWING VALUES OF T:
C     T=1,20, STEP=0.5 AND T=20,65 STEP=5
C     REMARK:
C     ======
C     FOR T SMALL, THAT IS FOR  T=1,20, STEP=0.5, ALL THE RESULTS ARE
C     QUITE ACCURATE, WHILE FOR T LARGE, THAT IS FOR T=20,65 STEP=5,
C     IN SOME CASES THE RESULTS COULDN'T BE ACCURATE. IN SUCH CASES 
C     THE AUTHORS SUGGEST TO USE AN ASYMPTOTIC INVERSION METHOD. 
C     *****************************************************************

      NT = 44
      DO 10 I = 1,39
           TARRAY(I) =1.D0 + (I-1)*0.5D0
   10 CONTINUE
      TARRAY(40) = 30.
      DO 20 I = 41,44
          TARRAY(I) = 30.D0 + 5.D0* (I-40)
   20 CONTINUE
      TOL = .1D-5
      SSBAR=0.d0
      WRITE (*,FMT=9010) AA,AA
      WRITE (*,FMT=9030) AA,AA
      WRITE (*,FMT=9020) TOL
C
C ************************************************************
C                    CALL OF THE ROUTINE INVLTF
C ************************************************************
C
C
      DO 2 I=1,NT
        CALL INVLTF(TOL,TARRAY(I),FZ,SIGMA0,SSBAR,NMAX,FZINV(I),ERR,
     +            IFZEVAL(I),WORK,IFAIL(I))
        ERROR(1,I)= ERR(1)
        ERROR(2,I)= ERR(2)
        ERROR(3,I)= ERR(3)
2     CONTINUE
C
C
C
      WRITE (*,FMT=9040) AA,AA
     
      DO 30 I = 1,NT
      STOP= .FALSE.
      IF (IFAIL(I).NE. 0 )THEN
         WRITE (*,FMT=9100) I,IFAIL(I)
         STOP=.TRUE.
      ENDIF
 

      IF (.NOT. STOP) THEN    

C
              EXF = FEX(TARRAY(I))
              DIFABS(I) = ABS(EXF-FZINV(I))
              IF (EXF.NE.0.D0) THEN
                  DIFREL(I) = DIFABS(I)/ABS(EXF)

              ELSE
                  DIFREL(I) = DIFABS(I)
              END IF

              EXF = FEX(TARRAY(I))
              IF (IFAIL(I).NE.-2) THEN
                  WRITE (*,FMT=9070) TARRAY(I),FZINV(I),EXF,ERROR(1,I),
     +         DIFREL(I),ERROR(3,I),ERROR(2,I),DIFABS(I),IFZEVAL(I),
     +         IFAIL(I)

              ELSE
                  WRITE (*,FMT=9060) TARRAY(I),IFAIL(I)
              END IF
      END IF

   30     CONTINUE
   50     CONTINUE

      STOP
C
 9000 FORMAT (1X,A45,A45,/,/,15X,
     +       '           SUBROUTINE INVLTF             ',/,15X,
     +       'NUMERICAL INVERSION OF A LAPLACE TRANSFORM:',/,15X,
     +       ' THIS VERSION USES BOTH REAL AND COMPLEX',/,15X,
     +       ' DOUBLE PRECISION OPERATIONS',/,/,1X,A45,A45,/)
 9010 FORMAT (/,A45,A45,/,/,16X,'        OUTPUT',/)
 9020 FORMAT (/,1X,'TOLL --> ',E15.7)
 9030 FORMAT (1X,'T      : POINT AT WHICH THE INVERSE TRANSFORM IS',
     +       ' COMPUTED;',/,/,1X,
     +       'FEX    : EXACT VALUE OF THE INVERSE TRANSFORM;',/,/,1X,
     +       'FCAL   : COMPUTED VALUE OF THE INVERSE TRANSFORM;',/,/,1X,
     +       'ESTREL : ESTIMATED RELATIVE ERROR;',/,/,1X,
     +       'RELERR : ACTUAL RELATIVE ERROR;',/,/,1X,
     +       'ESTABS : ESTIMATED ABSOLUTE ERROR ;',/,/,1X,
     +       'ABSERR : ACTUAL ABSOLUTE ERROR;',/,/,1X,
     +       'N      : # OF FUNCTION EVALUATIONS;',/,/,1X,
     +       'IFAIL  : = 0 NO INPUT ERRORS; SUCCESSFUL RUN',
     +       '    (ACCURACY REACHED AND IFZEVAL<NMAX),',/,10X,
     +       '= 1 TOL >= 1,',/,10X,
     +       '= 2 VALT LESS THAN ZERO,',/,10X,
     +       '= -1 ACCURACY NOT REACHED AND IFZEVAL > NMAX;,',/,10X,
     +       '= -2 THE CHOICE FOR SSBAR  MAY BE NOT OPTIMAL;',/,10X,
     +       '        IN SUCH A CASE THE USER MAY SLIGHTLY',/,10X,
     +       '        CHANGE THE DEFAULT VALUE;',/,/,A45,A45,/,/)
 9040 FORMAT (/,A86,/,'   T',9X,'FCAL',9X,'FEX',8X,' ESTREL',3X,
     +       'RELERR',2X,'TRUNERR',2X,'ESTABS',3X,'ABSERR',3X,'N',2X,
     +       'IFAIL',/,A100,/,/)
 9050 FORMAT (/,' TEST FUNCTION ----->  ',I2,/)
 9060 FORMAT (F5.1,90X,I2)
 9070 FORMAT (F5.1,1X,E14.8,1X,E14.8,1X,E8.3,1X,E8.3,1X,E8.3,1X,E8.3,1X,
     +       E8.3,1X,I3,3X,I2)
 9080 FORMAT (/,/,1X,A45,A45,/,/,2X,
     +     '   THE T-VALUES AT WHICH THE INVERSE IS REQUIRED ARE T=1,20'
     +       ,' STEP=0.5 AND T=20,100 STEP=10.',/,/,1X,A45,A45)
 9090 FORMAT (1X,'TEST  FUNCTION : ')
 9100 FORMAT (/,1X,'ERROR DETECTED , I = ', I3, ' IFAIL=',I3)
 9110 FORMAT (/,' ABSCISSA OF CONVERGENCE  ---> ',/)
 9120 FORMAT (/,' ABSCISSA OF CONVERGENCE  ---> ',F5.1,/)
      END

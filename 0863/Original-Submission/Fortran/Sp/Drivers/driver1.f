      PROGRAM DRIVER1
C
C     Single Precision
C      
C     Sample program that demonstrates subroutine L2WPMA when the
C      data is read from, and the reuslts are directed to, files.
C
C     Calls subroutine L2WPMA.
C......................................................................
C
C.... N O T A T I O N, D E C L A R A T I O N AND A R G U M E N T S
C     For details see subroutine L2WPMA.
C
C.... P A R A M E T E R S
C     IX1, NX  Lower and upper bound data indices.
C     KSECT    Upper bound to number of monotonic sections. Although
C               in the parameter setting is given that KSECT=20,
C               it may well take any value that is smaller than NX.
C     IORDER   Specifies whether first monotonic section is increasing
C               or decreasing (see comments in subrtn L2WPMA).
C              = 0, first monotonic section is increasing, default
C              =-1, first monotonic section is decreasing
C     MODEWF   Integer variable, specifying the weights WF(.):
C              = 0, User defined weights, the default value
C              = 1, Weights set to unity automatically
C              = 2, automatically, due to abscissae spacing
C              = 3, automatically, due to abscissae spacing
C              = 4, automatically, due to abscissae spacing
C              = 5, only F is provided, while X(i)=i and WF(.)=1
C                                                     (automatically)
C     IFILE    Integer parameter whose default value is 0.
C               If IFILE = 1  then some printout is displayed at the
C                              end of the smoothing process.
C
C.... I N P U T ....
C     'XFWDAT'  Datafile that holds I1 and N (first record of file),
C                and the data X(I),F(I),WF(I), for I = I1,I1 + 1,...,N
C                (one record for each value of I). 'XFWDAT' is
C                provided by the user (see driver program DRIVER2).
C
C.... O U T P U T
C     1) Results from the call of subroutine L2WPMA: A best weighted
C         least squares fit with KSECTN-1 sign changes in its first
C         differences to the data X(.), F(.).
C
C     2) File 'WPMAPP' keeps the best fit and corresponding Lagrange
C         multipliers.
C     3) File 'WIVAPP' keeps the indices of the extrema of the fit
C         contained in 'WPMAPP'.
C     4) File 'WIAKN' keeps the indices of the knots of a spline
C         representation of the fit contained in 'WPMAPP'.
C     5) File 'WACTLG' keeps the indices of the active constraints
C         of the fit and the corresponding Lagrange multipliers.
C
C     SCREEN OUTPUT is monitored by parameter IFILE.
C
C......................................................................
C
C     .. Parameters ..
C        (for details see subroutine L2WPMA):
C
      INTEGER IX1,NX,NX2,KSECT,IORDER,MODEWF,IFILE
      PARAMETER (IX1 = 1,NX = 2000,NX2 = (NX - IX1 + 1)/2 + 1,
     +           KSECT = 20,IORDER = 0,MODEWF = 0,IFILE = 1)
C     ..
C     .. Local Scalars ..
      INTEGER I,I1,J,KSECTN,MODE,N,NACT,NK
C     ..
C     .. Local Arrays ..
      REAL F(IX1:NX),FT(IX1:NX),FTNEG(IX1:NX),G(0:KSECT,IX1:NX2),
     +     RG(IX1:NX2),SS(IX1:NX),X(IX1:NX),Y(IX1:NX),
     +     Z(1:NX + IX1 - 1),WF(IX1:NX),WY(IX1:NX),WFT(IX1:NX),
     +     WY1(IX1:NX),WZ(1:NX + IX1 - 1),PAR(IX1:NX)
      INTEGER INDX(IX1:NX),ITAU(0:KSECT,IX1:NX2),ITHETA(0:KSECT),
     +        IUPPER(IX1:NX2),IW(1:NX + IX1 - 1),LOWER(IX1:NX2),
     +        IAKN(IX1-1:NX - IX1 + 1),IAKNW(1:NX - IX1 + 1),
     +        IACT(1:NX-IX1)
C     ..
C     .. External Subroutines ..
      EXTERNAL L2WPMA
C     ..
C     Read data.
C
      OPEN (2,FILE = 'XFWDAT')
      REWIND 2
      READ (2,FMT=*) I1,N
      DO 10 I = I1,N
         READ (2,FMT=*) X(I),F(I),WF(I)  
C                                Whenever MODEWF = 5, file 'XFWDAT'
C                                 should contain I1,N and only F(.).
C                                 Therefore in order to use this driver
C                                 program, replace the previous
C                                 statement with: 
C                                 READ (2,FMT=*) F(I)
   10 CONTINUE
      CLOSE (2)
      PRINT 9000
C
C     Calculate best weighted least squares approximation to
C      the data with KSECTN monotonic sections.
C
      PRINT 9010
      READ *,KSECTN
C
      CALL L2WPMA(I1,N,X,F,WF,MODEWF,KSECTN,IORDER,Y,WY,NK,IAKN,NACT,
     +            IACT,PAR,ITAU,ITHETA,MODE,SS,G,RG,LOWER,IUPPER,INDX,
     +            FT,WFT,FTNEG,WY1,Z,WZ,IW,IAKNW)
C
C
      IF (MODE .EQ. 0) GO TO 80
C
C     Send results to files.
C
      OPEN (3,FILE = 'WPMAPP')
      WRITE (3,FMT=9020) I1,N
      DO 20 I = I1,N
         IF (I .EQ. I1) THEN
             WRITE (3,FMT=9030) I,X(I),F(I),WF(I),Y(I)
         ELSE
             WRITE (3,FMT=9030) I,X(I),F(I),WF(I),Y(I),PAR(I)
         END IF
   20 CONTINUE
      CLOSE (3)
      PRINT 9040
C
      OPEN (3,FILE = 'WIVAPP')
      DO 30 I = I1 - 1,KSECTN
         WRITE (3,FMT=9050) I,ITHETA(I)
   30 CONTINUE
      CLOSE (3)
C
      OPEN (3,FILE = 'WIAKN')
      DO 40 I = I1,NK
         WRITE (3,FMT=9025) I,IAKN(I),Y(IAKN(I))
   40 CONTINUE
      CLOSE (3)
C
      IF (NACT .LT. 1) GO TO 60
          OPEN (3,FILE = 'WACTLG')
      DO 50 I = 1,NACT
         WRITE (3,FMT=9025) I,IACT(I), PAR(IACT(I))
   50 CONTINUE
      CLOSE (3)
   60 CONTINUE
C
      PRINT 9080
      READ *
C
C     PRINTOUT (optional)
C      The user may direct this printout to a file.
C
      IF (IFILE .EQ. 1) THEN
          WRITE (*,FMT=9090) (J,ITHETA(J),J=I1-1,KSECTN)
          IF (NACT .GT. 1) THEN
              WRITE (*,FMT=9100) (I,IACT(I),PAR(IACT(I)),I=1,NACT)
          END IF
          WRITE (*,FMT=9110) (I,IAKN(I),Y(IAKN(I)),I=I1,NK)
C
          WRITE (*,FMT=9120)
          DO 70 I = I1,N
             IF (I .EQ. I1) THEN
                 WRITE (*,FMT=9130) I,X(I),F(I),WF(I),Y(I)
             ELSE
                 WRITE (*,FMT=9140) I,X(I),F(I),WF(I),Y(I),PAR(I)
             END IF
   70     CONTINUE
C
          WRITE (*,FMT=9150) SS(N)
      END IF
   80 READ *
C
      STOP
C
 9000 FORMAT (/5X,'X-F-W Data supply for smoothing from datafile ',
     +       ' "XFWDAT" .')
 9010 FORMAT (/5X,'Just below, input number of monotonic sections < 21'
     +       )
 9020 FORMAT (2I5,F5.2)
 9025 FORMAT (2I5,5X,E20.10)
 9030 FORMAT (I5,5E20.10)
 9040 FORMAT (//5X,'-------------------------------------------------',
     +       /5X,'Best piecewise monotonic approximation to the data',
     +       /5X,'from file "XFWDAT" is now kept in file ..."WPMAPP".')
 9050 FORMAT (2I5)
 9070 FORMAT (3I5)
 9080 FORMAT (/5X,'Indices of extrema of best approximation are now',
     +        /5X,'kept in file ............................"WIVAPP".',
     +       //5X,'Indices of spline representation of best apprxmtn',
     +       /5X,'are now kept in file ....................."WIAKN".',
     +       //5X,'Indices of active constraints of best apprxtn and ',
     +       /5X,'Lagrange multipliers are now kept in file "WACTLG".',
     +       /5X,'--------------------------------------------------')
 9090 FORMAT (//12X,'INDICES OF EXTREMA AT OPTIMUM',/12X,'Increment',
     +       6X,'Data index',/12X,'(J)',12X,'(ITHETA)',
     +       / (12X,I5,12X,I5))
 9100 FORMAT (//12X,'INDICES OF ACTIVE CONSTRAINTS AT OPTIMUM',/12X,
     +       'Increment',6X,'Active constraint',6X,'Lagrange mult',
     +       /12X,'(I)',12X,'(IACT)',17X,'(PAR(IACT))',
     +       / (12X,I5,12X,I5,12X,F12.4))
 9110 FORMAT (//12X,'INDICES OF KNOTS AT OPTIMUM',/12X,'Increment',
     +       6X,'Knot index',13X,'Spline Coeff',
     +       /12X,'(I)',12X,'(IAKN)',17X,'(Y(IAKN))'
     +       / (12X,I5,12X,I5,12X,F12.4))
 9120 FORMAT (//7X,'Data points',2X,'Measurements',2X,'Data weights',
     +       2X,'Best appxmtn',2X,'Lagrange mult'/7X,'(X)',10X,'(F)',
     +       11X,'(WF)',10X, '(Y)',11X,'(PAR)'/)
 9130 FORMAT (I5,1X,4(F12.8,1X))
 9140 FORMAT (I5,1X,4(F12.8,1X),2X,F12.8)
 9150 FORMAT (/2X,'Value of objective function at the optimum=',E20.10)
      END

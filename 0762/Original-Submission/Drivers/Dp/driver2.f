      PROGRAM tstlld
C
C     PACKAGE WHICH PRINTS THE RESULTS OF CALL TO LLDRLF TO SCREEN
C
C     .. Parameters ..
      INTEGER bigar
      PARAMETER (bigar=100)
      DOUBLE PRECISION two
      PARAMETER (two=2.0D0)
      REAL ablow
      PARAMETER (ablow=1.0E-3)
      REAL abhigh
      PARAMETER (abhigh=1.0E10)
      REAL wlow
      PARAMETER (wlow=-1000.0E0)
      REAL whigh
      PARAMETER (whigh=1000.0E0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION dfd,dfn,a,b,w,worig,tau,dum,l,dldw,d2ldw2
      INTEGER i,iw,nw,nab,iab,ia,ib,inun,outun,stdout
      LOGICAL qfile,qscrn,qw,qab
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION warray(bigar)
      DOUBLE PRECISION aarray(bigar),barray(bigar)
      CHARACTER cases(3)*26
C     ..
C     .. External Functions ..
      LOGICAL qyngt
      INTEGER obfugt
      REAL rlgt
      EXTERNAL qyngt,obfugt,rlgt
C     ..
C     .. External Subroutines ..
      EXTERNAL lldrlf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,sqrt,abs,max,log10
C     ..
C     .. Data statements ..
      DATA (cases(i),i=1,3)/'log(f(exp(W)|DFN,DFD))    ',
     +     'log(F(exp(W)|DFN,DFD))    ','log(1 - F(exp(W)|DFN,DFD))'/
C     ..
C     .. Executable Statements ..
C
C     Setup for run
C
C     See if user wants to output to the screen or a file
C
      CALL gtstio(inun,stdout)
      qfile = (qyngt('('' Do you want to output to a file? (y/n)'')'))
      qscrn = .NOT. qfile
      IF (qfile) THEN
          outun = obfugt('('' Enter name of output file.'')',.TRUE.)

      ELSE
          outun = stdout
      END IF
C
C     Ask if user wants to specify A and B or run all possibilities
C
   10 qab = (qyngt('('' Do you want to specify A and B? (y/n)'')'))
C
C     Ask for A and B
C
      IF (qab) THEN
          nab = 1
          aarray(1) = dble(rlgt('('' Enter A.'')',ablow,abhigh))
          barray(1) = dble(rlgt('('' Enter B.'')',ablow,abhigh))

      ELSE
C
C     Otherwise define all A and B possibilities
C
          nab = 0
          DO 20 iab = -3,10,1
              nab = nab + 1
              aarray(nab) = 10.0D0**dble(iab)
              barray(nab) = aarray(nab)
   20     CONTINUE
      END IF
C
C     Ask if user wants to specify W  or run all possibilities
C
      qw = (qyngt('('' Do you want to specify W? (y/n)'')'))
C
C     Ask for W
C
      IF (qw) THEN
          nw = 1
          warray(1) = dble(rlgt(
     +                '('' Enter W (in standard deviations).'')',wlow,
     +                whigh))

      ELSE
C
C     Otherwise define all W possibilities
C
          nw = 0
          DO 30 iw = -1000,-100,100
              nw = nw + 1
              warray(nw) = dble(iw)
   30     CONTINUE
          DO 40 iw = -90,-10,10
              nw = nw + 1
              warray(nw) = dble(iw)
   40     CONTINUE
          DO 50 iw = -9,9,1
              nw = nw + 1
              warray(nw) = dble(iw)
   50     CONTINUE
          DO 60 iw = 10,100,10
              nw = nw + 1
              warray(nw) = dble(iw)
   60     CONTINUE
          DO 70 iw = 200,1000,100
              nw = nw + 1
              warray(nw) = dble(iw)
   70     CONTINUE
      END IF
C
C     loop through all possible combinations of A and B
C
      DO 100 ia = 1,nab
          DO 90 ib = 1,nab
C
C     Caculate tau, the standard deviation of the distribution,
C     and define dfn and dfd
C
              a = aarray(ia)
              b = barray(ib)
              dfn = two*a
              dfd = two*b
              CALL mvlogf(dfn,dfd,dum,tau)
              tau = sqrt(tau)
C
C        output information on a, b and tau
C
              WRITE (outun,9000)
              WRITE (outun,9010)
              WRITE (outun,9020) 'A: ',a
              WRITE (outun,9020) 'B: ',b
              IF ((tau.GE.0.0001) .AND. (tau.LE.99999.0)) THEN
                  WRITE (outun,9030) 'Standard Deviation: ',tau

              ELSE
                  WRITE (outun,9040) 'Standard Deviation: ',tau
              END IF
C
C     loop through all w
C
              DO 80 iw = 1,nw
                  worig = warray(iw)
                  w = worig*tau
C
C          call lldrlf
C
                  WRITE (outun,9050)
                  WRITE (outun,9060) 'W (in standard deviations)',worig
C
C          case = 1
C
                  CALL lldrlf(w,dfn,dfd,1,l,dldw,d2ldw2)
                  WRITE (outun,9070) cases(1),l
                  WRITE (outun,9070) '   First Derivative',dldw
                  WRITE (outun,9070) '   Second Derivative',d2ldw2
C
C          case = 2
C
                  CALL lldrlf(w,dfn,dfd,2,l,dldw,d2ldw2)
                  WRITE (outun,9070) cases(2),l
                  WRITE (outun,9070) '   First Derivative',dldw
                  WRITE (outun,9070) '   Second Derivative',d2ldw2
C
C          case = 3
C
                  CALL lldrlf(w,dfn,dfd,3,l,dldw,d2ldw2)
                  WRITE (outun,9070) cases(3),l
                  WRITE (outun,9070) '   First Derivative',dldw
                  WRITE (outun,9070) '   Second Derivative',d2ldw2
                  IF (qscrn) CALL pause
   80         CONTINUE
   90     CONTINUE
  100 CONTINUE
C
C     Ask if they wish to continue
C
      IF (qyngt('('' Do you want to continue? (y/n)'')')) GO TO 10
C
C     close file and quit
C
      IF (qfile) CLOSE (outun)
C
C     Format Statments
C
 9000 FORMAT (1X,'++++++++++++++++++++++++++++++++++++++++++++++++++')
 9010 FORMAT (1X,'------------------ Values Used -------------------')
 9020 FORMAT (1X,1A,1P,1E10.4)
 9030 FORMAT (1X,1A,1F10.4)
 9040 FORMAT (1X,1A,1P,1E10.4)
 9050 FORMAT (1X,'**************************************************')
 9060 FORMAT (1X,1A30,1F12.4)
 9070 FORMAT (1X,1A30,1P,1E25.16)

      END

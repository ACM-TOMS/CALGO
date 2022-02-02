      PROGRAM DRIVER
*-------------------------------------------------------
*     Update History:
*     July 9,  2004 by G. Howell
*     Nov. 21, 2002 by N. Diaa
*     Nov. 20, 2002 by N. Diaa
*-------------------------------------------------------
C     .. Parameters ..
      INTEGER LDX
      PARAMETER (LDX=500)
C     ..
C     .. Scalars in Common ..
      INTEGER SEED
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONDIT,SIZE2,TOL,XNORM1,XNORM2,XNORM3
      INTEGER AVGNU,I,ICHOIC,IDUM,IFLAG,IGO,IMAT,INFO,INIXSIX,ITEMP,
     +        ITEST,ITESTS,J,M,MAXNU,NUMAX,PRLEVEL,SIZE
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BH(LDX,LDX),DUMV(LDX),MAT0(LDX,LDX),V(LDX,3),
     +                 WMAT(LDX,LDX),XIM(LDX),XRE(LDX)
      INTEGER BHPIVT(LDX),KROW(LDX),KRVECT(LDX),NU(LDX),NUMON(LDX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION RANDOM
      EXTERNAL RANDOM
C     ..
C     .. External Subroutines ..
      EXTERNAL BHAP1,BHAP2,BHAP3,BHAPC,BHBACK,BHESS,BHRPCK,BHUNPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON SEED
C     ..
      PRINT *,'              '
      PRINT *,' BHESS TESTER. '
      PRINT *,' Input an integer on each pause       '
      PRINT *,'              '
      PRINT *,'  BHESS reduces a general square matrix'
      PRINT *,'    to similar small-band Hessenberg form'
      PRINT *,'    using multipliers bounded in root mean'
      PRINT *,'    by a user-specified bound tol .'
      PRINT *,'                  '
      PRINT *,'  Taking tol = 30 typically results in '
      PRINT *,'    an almost tridiagonal matrix from'
      PRINT *,'    which eigenvalues can be extracted'
      PRINT *,'    by BR iteration.'
      PRINT *,'               '
      PRINT *,' DRIVER calls BHESS and some auxiliary'
      PRINT *,'    routines.  It can take input files produced'
      PRINT *,'    by Matlab or Octave and produce files readable'
      PRINT *,'    by those programs, allowing verification of'
      PRINT *,'    correct computation.'
      PRINT *,'                   '
      PRINT *,' In response to system prompts'
      PRINT *,'                  enter integers.'
      PRINT *,' There is additional info in README.'
      PRINT *,'                         '
      PRINT *,'                         '
   10 CONTINUE
      PRINT *,' Input an integer PRINT LEVEL:  '
      PRINT *,'    0 will give no screen output.'
      PRINT *,'    1 will give maximal output.'
      READ *,PRLEVEL
      IF ((PRLEVEL.NE.0) .AND. (PRLEVEL.NE.1)) GO TO 10
   20 CONTINUE
      PRINT *,' Choose:'
      PRINT *,'      1) to randomly generate input matrix '
      PRINT *,'      2) to read input matrix from file BHESS.IN '
      PRINT *,'      3) to quit '

      INIXSIX = 0
      IMAT = 0
      ITEMP = 0

      READ *,ICHOIC

      IF (ICHOIC.EQ.3) GO TO 110
      IF ((ICHOIC.NE.1) .AND. (ICHOIC.NE.2)) GO TO 20

      PRINT *,' Output is to screen'

*----------------------------------------------------------------
*
*     For ICHOIC = 1, an n x n matrix is generated with a given seed.
*     using a linear congruential random number generator.  Entries of
*     the matrix are from a uniform distribution between -1 and 1.
*
*---------------------------------------------------------------

* 1) to randomly generate input matrix
      IF (ICHOIC.EQ.1) THEN
   30     CONTINUE
          PRINT *,'    '
          PRINT *,' Input integer size (2 < size <= ',LDX,' ): '
          READ *,SIZE
          IF ((SIZE.GT.LDX) .OR. (SIZE.LE.2)) GO TO 30
   40     CONTINUE
          PRINT *,'    '
          PRINT *,'Input integer random generator seed > 0:'
          READ *,SEED
          IF (SEED.LE.0) GO TO 40
          DO J = 1,SIZE
              DO I = 1,SIZE
                  MAT0(I,J) = (RANDOM(SEED)-.5)*2.D0
              END DO
          END DO
      END IF

*-----------------------------------------------------------------
*
*     If ICHOIC is 2, then matrix entries are read
*           from the file BHESS.IN
*
*-----------------------------------------------------------------
*
*     To create the file BHESS.IN from matlab,
*       assign a tolerance tol, a size n, and a square matrix A.
*       Use the matlab command
*          save BHESS.IN tol n A -ascii
*
*     To create the file BHESS.IN from the matlab clone octave,
*       assign the same variables and use the comman
*          save BHESS.IN tol n A
*     But you then have to edit the resulting file
*        to remove the nonnumeric lines.
*
*     Typing format long inside octave results in more digits
*        being saved.
*
*  To produce BHESS.IN from matlab
*  assign tolerance tol, size n, square matrix A
*  Use the matlab command
*  Save BHESS.IN tol n A -ascii
*  Enter an arbitary integer to continue
*  Accuracy is
*  constrained by the 8 digits in the ascii
*  format used to write the file BHESS.IN
*  from Matlab
*  In Octave the produced file has extraneous
*  characters which will need to be deleted
*  before BHESS.IN can be used as an input file.
*  Typing
*  format long
*  before saving will give more digits in
*  the saved file.
*--------------------------------------------------------------

* ICHOIC=2 read matrix from BHESS.IN
      IF (ICHOIC.EQ.2) THEN
          PRINT *,'You chose to read input matrix from file BHESS.IN'
          PRINT *,'This option assumes you have a file BHESS.IN'
          PRINT *,'To see how to produce this file see the    '
          PRINT *,'README file or the comments in abh3.m      '

          PRINT *,' Input any integer to continue'
          READ *,IDUM
          OPEN (33,FILE='BHESS.IN',STATUS='UNKNOWN',ERR=100)
          READ (33,FMT=*) TOL
          READ (33,FMT=*) SIZE2
          SIZE = INT(SIZE2)
          IF ((SIZE.GT.LDX) .OR. (SIZE.LE.2) .OR. (TOL.LE.0)) THEN
              PRINT *,' ERROR: Invalid data in BHESS.IN '
              PRINT *,
     +          '       Real number Tol or integer Size in BHESS.IN '
              PRINT *,'        Tol must be real number > 0 '
              PRINT *,'        Size must be integer: 2 < size <= ',LDX,
     +          '  '
              CLOSE (33)
              GO TO 110

          END IF

          DO I = 1,SIZE
              PRINT *,' Reading row ',I,' of ',SIZE
              READ (33,FMT=*) (MAT0(I,J),J=1,SIZE)
          END DO
          CLOSE (33)
      END IF

*--------------------------------------------------------------
*
*       A good general policy is to call a balancing routine
*       (using diagonal similarity transformations)
*       before BHESS (or any reduction to Hessenberg form).
*       A standard routine is balanc.f from EISPACK.
*
*--------------------------------------------------------------
*     Now we have the matrix MAT0,
*     either read from BHESS.IN or randomly generated
*--------------------------------------------------------------
   50 CONTINUE
      PRINT *,'        Choose set of tests                  '
      PRINT *,'   '
      PRINT *,' Input 1 to get output that can be compared  '
      PRINT *,' to Matlab and illustrates the actions of the'
      PRINT *,'        included subroutines'
      PRINT *,'   '
      PRINT *,' Input 2 to test the relation of the input '
      PRINT *,' parameter TOL on computed backward error,'
      PRINT *,' resulting bandwidth, and estimated condition'
      PRINT *,' number of the similarity transformation'
      READ *,ITESTS
      IF ((ITESTS.NE.1) .AND. (ITESTS.NE.2)) GO TO 50

      IF (ITESTS.EQ.2) THEN
   60     CONTINUE
          PRINT *,' Input TOL as a parameter for calling BHESS'
          PRINT *,' TOL is the maximal allowable root mean square'
          PRINT *,' for vectors of paired row-column multipliers'
          PRINT *,' TOL should be a real number bigger than zero'
          READ *,TOL
          IF (TOL.LE.0) GO TO 60

          CALL BHBACK(LDX,SIZE,MAT0,BH,WMAT,V,TOL,NU,BHPIVT,KRVECT,KROW)

          PRINT *,' '
          PRINT *,' '
          PRINT *,'  Input 1 to do another round of tests with the'
          PRINT *,'       same matrix and a new input TOL,       '
          PRINT *,'        Else other integers quit.            '
          READ *,ITESTS
          IF (ITESTS.EQ.1) GO TO 60
          GO TO 110

      END IF

* Choice of randomly generated matrix with some explicit tests of
* routines. This choice may not be very useful.
      IF (ICHOIC.EQ.1) THEN
          PRINT *,' This choice allows you to get a feel for '
          PRINT *,' available subroutines, but does not allow'
          PRINT *,' easy tests of validity.  '
          PRINT *,' For a test of BHESS validity, either use matrices'
          PRINT *,' exported from matlab  and compare the results here'
          PRINT *,' to those from running abh3.m '
          PRINT *,' OR '
          PRINT *,' start over and choose Test set 2'
          PRINT *,' on the step just before this message.'
          PRINT *,'    '

   70     CONTINUE
          PRINT *,' Input Tolerance real number greater than zero  '
          PRINT *,' TOL bounds the product of the root mean squares'
          PRINT *,'   of vectors of column and row multipliers     '
          READ *,TOL
          IF (TOL.LE.0) GO TO 70

   80     CONTINUE
          PRINT *,'    '
          PRINT *,' Input choice: '
          PRINT *,'  1 = To create a MATLAB file called mat.m '
          PRINT *,'      that will contain a copy of your n by n '
          PRINT *,'      input matrix and the size and tol you are '
          PRINT *,'      using. NOTE: If you already have a file, '
          PRINT *,'      mat.m, this program will not overwrite'
          PRINT *,'      it. '
          PRINT *,'      If you run BHAP1, BHAP2, and BHAP3      '
          PRINT *,'      from the ensuing menu, their results   '
          PRINT *,'      will also be stored in mat.m            '
          PRINT *,' '
          PRINT *,'  0 = Do not create a file mat.m'
          PRINT *,' '
          READ *,IMAT
          IF ((IMAT.NE.1) .AND. (IMAT.NE.0)) GO TO 80

* IMAT=1 means write to mat.m
          IF (IMAT.EQ.1) THEN
              OPEN (8,FILE='mat.m',STATUS='NEW',ERR=100)
              WRITE (8,FMT=*) ' n = ',SIZE,';'
              WRITE (8,FMT=*) ' tol = ',TOL,';'
              WRITE (8,FMT=*) ' a = eye(n) ;'
              DO I = 1,SIZE
                  DO J = 1,SIZE
                      WRITE (8,FMT=*) 'a(',I,',',J,')=',MAT0(I,J),';'
                  END DO
              END DO
          END IF

      END IF

*----------------------------------------------------------------
*
*   WARNING!! BHESS requires MAT0 to have SIZE + 1 Columns
*           though the matrix to be reduced has only SIZE Columns
*
*----------------------------------------------------------------

      CALL BHESS(LDX,MAT0,SIZE,BHPIVT,NU,KRVECT,KROW,DUMV,TOL,INFO)

      IF (IMAT.EQ.1) THEN
          DO I = 1,SIZE
              DO J = 1,SIZE
                  WMAT(I,J) = MAT0(I,J)
              END DO
          END DO

          CALL BHRPCK(LDX,SIZE,WMAT,NU)

          PRINT *,' The result of BHESS (with multipliers zeroed '
          PRINT *,' by a call to BHRPCK) will be stored mat.m    '
          PRINT *,' as the matrix B                              '
          PRINT *,'     '
          WRITE (8,FMT=*) ' B = eye(n) ;'
          WRITE (8,FMT=*) ' tol = ',TOL,' ;'
          DO I = 1,SIZE
              DO J = 1,SIZE
                  WRITE (8,FMT=*) 'B(',I,',',J,')=',WMAT(I,J),';'
              END DO
          END DO
      END IF

*--------------------------------------------------------------
*
*       The following continue statement (300) is returned to from
*       the bottom of the routine, allowing the user to
*       retry auxiliary routines
*
*---------------------------------------------------------------

   90 CONTINUE

*--------------------------------------------------------------
*
*     The average bandwidth is computed as AVGNU
*
*--------------------------------------------------------------

      AVGNU = 0
      MAXNU = 0
      IF (PRLEVEL.EQ.1) THEN
          PRINT *,' NU(I) is the number of nonzero entries to the right'
     +      ,' of the diagonal in the ith row      '
          PRINT *,'    '
          PRINT *,' Input any integer to contiune'
          READ *,IDUM
      END IF

      DO I = 1,SIZE
          IF (PRLEVEL.EQ.1) THEN
              PRINT *,'NU(',I,')=',NU(I)
              IF (I.EQ. (I/20)*20) THEN
                  PRINT *,' Input any integer to continue'
                  READ *,IDUM
              END IF

          END IF

          AVGNU = AVGNU + NU(I)
          MAXNU = MAX(MAXNU,NU(I))
      END DO
      AVGNU = (AVGNU+SIZE/2)/SIZE
      PRINT *,' Maximal upper bandwidth of returned matrix=',MAXNU
      PRINT *,'                 '
      IF (PRLEVEL.EQ.1) THEN
          PRINT *,'    '
          PRINT *,'    '
          PRINT *,' Input 1 to see the nonzero entries of the'
          PRINT *,' returned Hessenberg matrix.'
          PRINT *,' Input any other integer to continue.'
          READ *,IGO
          IF (IGO.EQ.1) THEN
              PRINT *,'Row 1', (MAT0(1,J),J=1,1+NU(1))
              DO I = 2,SIZE
                  PRINT *,'Row',I, (MAT0(I,J),J=I-1,I+NU(I))
                  IF (I.EQ. (I/20)*20) THEN
                      PRINT *,' Input any integer to continue'
                      READ *,IDUM
                  END IF

              END DO
          END IF

      END IF
*----------------------------------------------------------------
*
*      The next parts of the code illustrate auxiliary returns
*        for working with the returned matrix.
*
*----------------------------------------------------------------


*----------------------------------------------------------------
*
*      Test BHAP1
*
*----------------------------------------------------------------
      PRINT *,'                      '
      PRINT *,'    '
      PRINT *,' Input 1 to test BHAP1 by multiplying '
      PRINT *,'   by a vector of ones.'
      PRINT *,'  Any other integer will continue to '
      PRINT *,'    the next option. '
      PRINT *,'  BHAP1 applies inv(Z) to a row vector where'
      PRINT *,'    inv(Z)*A*Z= H, H the near tridiagonal matrix'
      PRINT *,'  BHAP1 returns a row vector times inv(Z) '
      PRINT *,'  This operation converts a left eigenvector of'
      PRINT *,'    H to a left eigenvector of A'
      PRINT *,'      The output real and complex parts overwrite the'
      PRINT *,'        input vector'
      PRINT *,'                '

      READ *,ITEST

      IF ((ITEST.EQ.1) .AND. (INIXSIX.EQ.0)) THEN
          DO I = 1,SIZE
              XRE(I) = 1.D0
              XIM(I) = 0.D0
          END DO
          CALL BHAP1(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,IFLAG)
          XNORM1 = 0.D0
          IF (PRLEVEL.EQ.0) THEN
              PRINT *,'                            '
              PRINT *,' Input 1 to see output vector, '
              PRINT *,' else input 0.'
              PRINT *,'               '
              READ *,ITEMP
          END IF

          IF ((ITEMP.EQ.1) .OR. (PRLEVEL.EQ.1)) THEN
              PRINT *,' XRE(I)   I BHPIVT(I)   KROW(I) '
              DO I = 1,SIZE
                  PRINT *,XRE(I),I,BHPIVT(I),KROW(I)
                  IF (I.EQ. (I/20)*20) THEN
                      PRINT *,' Enter any integer to continue'
                      READ *,IDUM
                  END IF

                  XNORM1 = MAX(XNORM1,ABS(XRE(I)))
              END DO
              PRINT *,' The Max Norm is = ',XNORM1
              PRINT *,' XRE(2)= ',XRE(2)
              PRINT *,' XRE(SIZE/2)= ',XRE(SIZE/2)
              PRINT *,'                 '
          END IF

          IF (IMAT.EQ.1) THEN
              PRINT *,'  The resulting vector will be stored as  '
              PRINT *,'  the vector oneinvz '
              PRINT *,'  in the matlab file mat.m                     '
              PRINT *,'               '
              DO J = 1,SIZE
                  WRITE (8,FMT=*) 'oneinvz(1,',J,')=',XRE(J),';'
              END DO
          END IF

      ELSE
          IF ((ITEST.EQ.1) .AND. (INIXSIX.NE.0)) THEN
              PRINT *,
     +          ' Multipliers have been zeroed so multiplication by'
              PRINT *,' Z is infeasible'
          END IF

      END IF

*----------------------------------------------------------------
*
*      Test BHAP2
*
*----------------------------------------------------------------
      PRINT *,'           '
      PRINT *,'           '
      PRINT *,' Input 2 to test BHAP2 by multiplying '
      PRINT *,'   Z by a vector of ones.'
      PRINT *,'           '
      PRINT *,' BHAP2 multiplies an input row vector X'
      PRINT *,'   by Z where '
      PRINT *,'   inv(Z)*A*Z= H, H the near tridiagonal matrix'
      PRINT *,' BHAP2 returns xt = xt*Z '
      PRINT *,'                  '

      READ *,ITEST
      IF ((ITEST.EQ.2) .AND. (INIXSIX.EQ.0)) THEN
          DO I = 1,SIZE
              XRE(I) = 1.D0
          END DO
          CALL BHAP2(SIZE,MAT0,LDX,XRE,BHPIVT,KROW,IFLAG)
          IF (PRLEVEL.EQ.0) THEN
              PRINT *,'                            '
              PRINT *,' Input 1 to see output vector, '
              PRINT *,' else input 0 '
              PRINT *,'               '
              READ *,ITEMP
          END IF

          XNORM2 = 0.D0
          IF ((ITEMP.EQ.1) .OR. (PRLEVEL.EQ.1)) THEN
              PRINT *,' X(I)   I BHPIVT(I)   KROW(I) '
              DO I = 1,SIZE
                  PRINT *,XRE(I),I,BHPIVT(I),KROW(I)
                  IF (I.EQ. (I/20)*20) THEN
                      PRINT *,' Enter any integer to continue'
                      READ *,IDUM
                  END IF

                  XNORM2 = MAX(XNORM2,ABS(XRE(I)))
              END DO
              PRINT *,'                 '
              PRINT *,' The Max Norm is  ',XNORM2
              PRINT *,'               '
              PRINT *,'               '
          END IF

          IF (IMAT.EQ.1) THEN
              PRINT *,' The resulting vector will be stored in '
              PRINT *,
     +          ' the file mat.m as onez                            '
              PRINT *,'               '
              DO J = 1,SIZE
                  WRITE (8,FMT=*) 'onez(1,',J,')=',XRE(J),';'
              END DO
          END IF

      ELSE
          IF ((ITEST.EQ.2) .AND. (INIXSIX.NE.0)) THEN
              PRINT *,'                 '
              PRINT *,
     +          ' Multipliers have been zeroed so multiplication by'
              PRINT *,' Z is infeasible'
          END IF

      END IF

*----------------------------------------------------------------
*
*      Test BHAP3
*
*----------------------------------------------------------------

      PRINT *,'            '
      PRINT *,'            '
      PRINT *,' Input 3 to test BHAP3 by multiplying '
      PRINT *,'   by a vector of ones.'
      PRINT *,' BHAP3 multiplies by Z where'
      PRINT *,'   inv(Z)*A*Z= H, H the near tridiagonal matrix,'
      PRINT *,' The output real and complex parts overwrite the'
      PRINT *,'   input vector.  If the input vector is an'
      PRINT *,'   eigenvector of H, the ouput vector is an '
      PRINT *,'   eigenvector A. '

      READ *,ITEST
      IF ((ITEST.EQ.3) .AND. (INIXSIX.EQ.0)) THEN
          DO I = 1,SIZE
              XRE(I) = 1.D0
              XIM(I) = 0.D0
          END DO
          CALL BHAP3(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,IFLAG)
          IF (PRLEVEL.EQ.0) THEN
              PRINT *,'                            '
              PRINT *,' Input 1 to see output vector, '
              PRINT *,' else input 0 '
              PRINT *,'               '
              READ *,ITEMP
          END IF

          XNORM3 = 0.D0
          IF ((ITEMP.EQ.1) .OR. (PRLEVEL.EQ.1)) THEN
              PRINT *,' XRE(I)   I BHPIVT(I)   KROW(I) '
              DO I = 1,SIZE
                  PRINT *,XRE(I),I,BHPIVT(I),KROW(I)
                  IF (I.EQ. (I/20)*20) THEN
                      PRINT *,' Enter any integer to continue'
                      READ *,IDUM
                  END IF

                  XNORM3 = MAX(XNORM3,ABS(XRE(I)))
              END DO
              PRINT *,'The Max Norm is ',XNORM3
              PRINT *,' XRE(2)= ',XRE(2)
              PRINT *,' XRE(SIZE/2)= ',XRE(SIZE/2)
              PRINT *,'               '
              PRINT *,'  The result given here can be tested '
              PRINT *,'     by comparing to the Matlab-Octave version. '
              PRINT *,'     See comments in the script file abh3.m '
              PRINT *,'               '
              PRINT *,'               '
          END IF

          IF (IMAT.EQ.1) THEN
              PRINT *,' The resulting vector will be stored in '
              PRINT *,
     +          ' the file mat.m as z1                            '
              PRINT *,'               '
              DO J = 1,SIZE
                  WRITE (8,FMT=*) 'z1(',J,',1)=',XRE(J),';'
              END DO
          END IF

      ELSE
          IF ((ITEST.EQ.3) .AND. (INIXSIX.NE.0)) THEN
              PRINT *,'                 '
              PRINT *,
     +          ' Multipliers have been zeroed so multiplication by'
              PRINT *,' Z is infeasible'
          END IF

      END IF

*     ---------------------------------------------------------------
*
*     The next statements call BHAPC -- a LINPACK style
*       estimator of the infinity norm of inv(Z)
*
*     ---------------------------------------------------------------

      PRINT *,'               '
      PRINT *,'               '
      PRINT *,' Input 4 to estimate COND(Z) by multiplying '
      PRINT *,'   by a vector W of ones and minus ones.'
      PRINT *,' BHAPC multiplies by Z where'
      PRINT *,'   inv(Z)*A*Z= H, H the near tridiagonal matrix'
      READ *,ITEST
      IF ((ITEST.EQ.4) .AND. (INIXSIX.EQ.0)) THEN
          CALL BHAPC(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,CONDIT,IFLAG)
          IF (PRLEVEL.EQ.1) THEN
              PRINT *,' Y'
              DO I = 1,SIZE
                  PRINT *,'Y(',I,')=',XRE(I)
              END DO
          END IF

          PRINT *,'                 '
          PRINT *,' ESTIMATED CONDITION OF ACCUMULATED'
          PRINT *,' TRANSFORMATIONS =',CONDIT
          PRINT *,'                 '

      ELSE
          IF ((ITEST.EQ.4) .AND. (INIXSIX.NE.0)) THEN
              PRINT *,
     +          ' Multipliers have been zeroed so multiplication by'
              PRINT *,' Z is infeasible'
          END IF

      END IF

*     ---------------------------------------------------------------
*
*      Test BHUNPK  or  BHRPCK
*
*     ---------------------------------------------------------------
      PRINT *,'                '
      PRINT *,' Input any integer to continue'
      READ *,IDUM
      PRINT *,'             '
      PRINT *,'             '
      PRINT *,' Input 5 to test BHRPCK  '
      PRINT *,'   BHRPCK converts the retuned matrix H to'
      PRINT *,'   column storage with the first column'
      PRINT *,'   the subdiagonal.  This is convenient'
      PRINT *,'   for BR iteration'
      PRINT *,'       '
      PRINT *,'             '
      PRINT *,' Input 6 TO test BHUNPCK '
      PRINT *,'   BHRPCK zeros the multipliers so further tests of'
      PRINT *,'   multiplication of Z*x or xtr*inv(Z) and estimation'
      PRINT *,'   of cond(Z) will not be possible'
      PRINT *,' '
      PRINT *,' Input any other integer to go on.'
      PRINT *,'       '
      READ *,ITEST
      IF (ITEST.EQ.6) INIXSIX = 1
      IF (ITEST.EQ.5) THEN
          CALL BHUNPK(LDX,MAT0,BH,SIZE,NU,NUMON,NUMAX)
          IF (PRLEVEL.EQ.1) THEN
              DO I = 1,SIZE
                  PRINT *, (I,BH(I,J),J=1,NU(I)+2)
              END DO
          END IF

      END IF


      IF (ITEST.EQ.6) THEN

*        --------------------------------------------
*         zero out multipliers in atobhes matrix MAT0
*        --------------------------------------------

          CALL BHRPCK(LDX,SIZE,MAT0,NU)
          IF (PRLEVEL.EQ.1) THEN
              DO I = 1,SIZE
                  M = MIN(SIZE,NU(I)+2)
                  DO J = 1,M
                      IF (I+J-2.NE.0) THEN
                          PRINT *,' MAT0(',I,',',I + J - 2,')=',
     +                      MAT0(I,I+J-2)
                      END IF

                  END DO
              END DO
          END IF

      END IF

      PRINT *,' Input 1 to reloop through auxiliary routines'
      PRINT *,'       with the matrix returned from BHESS.    '
      PRINT *,'   '
      PRINT *,' Input any other integer to quit.          '
      IF (IMAT.EQ.1) THEN
          PRINT *,' ON SOME MACHINES, THE FILE mat.m WILL NOT EXIST '
          PRINT *,'      UNTIL YOU FINISH EXECUTING THE DRIVER     '
      END IF

      READ *,IGO
      IF (IGO.EQ.1) GO TO 90
      IF (IMAT.EQ.1) CLOSE (8)

      GO TO 110

  100 CONTINUE
      PRINT *,' ERROR: One of the following file open errors occured: '
      PRINT *,' -  error trying to create new file mat.m.'
      PRINT *,'    this file may already exist. '
      PRINT *,' -  error trying to open file BHESS.IN. '


  110 CONTINUE
      PRINT *,'  '
      PRINT *,'  '
      PRINT *,' END OF DRIVER PROGRAM. '
      STOP

      END

*
*    ----------------------------------------------

      DOUBLE PRECISION FUNCTION RANDOM(SEED)

*     ----------------------------------------------
*     Returns a pseudo-random number between 0-1.
*     ----------------------------------------------

*     COMMON SEED
C     .. Scalar Arguments ..
      INTEGER SEED
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FMD
      INTEGER I,M,MD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Data statements ..

      DATA M/25173/,I/13849/,MD/65536/,FMD/65536.D0/
C     ..

      SEED = MOD(M*SEED+I,MD)
      RANDOM = SEED/FMD
      RETURN

      END

*     -----------------------------------------------
*
*====================================================
      SUBROUTINE BHBACK(LDX,N,A,B,W,V,TOL,NU,PIV,KRVECT,KROW)

*  This subroutine tests backward error of BHESS
*    by running BHESS to get a Hessenberg matrix
*    then using the stored multipliers to regenerate
*    the original matrix.
*
*  Outputs are the difference between the original matrix
*    and the reconstruction (computed backward error)
*
*  The upper bandwidth of the returned matrix and the estimated
*    condition number of the similariy transformation are
*    also returned.
*
*   Arguments
*-------------------
*
*     N -- Integer  Number of rows and columns of input matrices
*
*     LDX -- Integer Leading Dimension of Variable dimension arrays
*
*     A -- Double Precision Input matrix unchanged
*
*     B and W -- Double Precision Matrices, B must of at
*       least N+1 columns.   On return B is an approximation
*       of the input matrix A (to backward error) and W
*       holds the output from BHESS.
*
*     TOL -- Double precision parameter to control multiplier
*            size in BHESS.  Unchanged.
*
*     V  -- Double precision work vector of three columns.


*
*   Local Variables
*-------------------------------

*
*   Copy A to B
*
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL
      INTEGER LDX,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDX,*),B(LDX,*),V(LDX,3),W(LDX,*)
      INTEGER KROW(*),KRVECT(*),NU(*),PIV(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONDIT,FNORM,RNORM,SUM,SUM2
      INTEGER I,IFLAG,INFO,J,MXBAND
C     ..
C     .. External Subroutines ..
      EXTERNAL BHAP1,BHAP3,BHAPC,BHESS,BHRPCK,DCOPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT,MAX
C     ..
      DO J = 1,N
          DO I = 1,N
              B(I,J) = A(I,J)
          END DO
      END DO
*--------------------------------------------------------
*
*   Call BHESS with B  Warning B requires N + 1 Columns
*
*--------------------------------------------------------

      CALL BHESS(LDX,B,N,PIV,NU,KRVECT,KROW,V,TOL,INFO)
*
*   Copy B to W
*
      DO J = 1,N
          DO I = 1,N
              W(I,J) = B(I,J)
          END DO
      END DO

*
*   Zero out all but upper Hessenberg in B
*
      CALL BHRPCK(LDX,N,B,NU)
*
*   Calculate backward A as B = Z * B * inv(Z)

*
*      First B <-- Z * B
      DO I = 1,N - 1,2
          CALL BHAP3(N,W,LDX,B(1,I),B(1,I+1),PIV,KROW,IFLAG)
      END DO
      IF (N.NE.2* (N/2)) THEN
          CALL BHAP3(N,W,LDX,B(1,N),V(1,2),PIV,KROW,IFLAG)
      END IF
*
*   Then B <-- B * inv(Z)       So B approximates A
*
      DO I = 1,N - 1,2
          CALL DCOPY(N,B(I,1),LDX,V(1,1),1)
          CALL DCOPY(N,B(I+1,1),LDX,V(1,2),1)
          CALL BHAP1(N,W,LDX,V(1,1),V(1,2),PIV,KROW,IFLAG)
          CALL DCOPY(N,V(1,1),1,B(I,1),LDX)
          CALL DCOPY(N,V(1,2),1,B(I+1,1),LDX)
      END DO
      IF (N.NE.2* (N/2)) THEN
          CALL DCOPY(N,B(N,1),LDX,V(1,1),1)
          CALL BHAP1(N,W,LDX,V(1,1),V(1,2),PIV,KROW,IFLAG)
          CALL DCOPY(N,V(1,1),1,B(N,1),LDX)
      END IF
*
*   Calculate || A - B ||
*
      SUM = 0.D0
      SUM2 = 0.D0
      DO I = 1,N
          DO J = 1,N
              SUM = SUM + (A(I,J)-B(I,J))**2
              SUM2 = SUM2 + A(I,J)**2
          END DO
      END DO
* Frobenius norm normalized so that norm(eye(n)) = 1
      FNORM = DSQRT(SUM/N)
* Relative Frobenius norm
      RNORM = DSQRT(SUM/SUM2)
*
*   Estimate conditioning of Z
*
      CALL BHAPC(N,W,LDX,V(1,1),V(1,2),PIV,KROW,CONDIT,IFLAG)
*
*      Calculate Maximal Bandwidth
*
      MXBAND = 1
      DO I = 1,N
          MXBAND = MAX(NU(I),MXBAND)
      END DO
      PRINT *,'    '
      PRINT *,'    '
      PRINT *,'                             n =  ',N
      PRINT *,'                           tol = ',TOL
      PRINT *,' Frobenius Backward Error Norm = ',FNORM
      PRINT *,'  Relative Backward Error Norm = ',RNORM
      PRINT *,'    Estimated Condition Number = ',CONDIT
      PRINT *,'                     Bandwidth =  ',MXBAND

      RETURN

      END

      PROGRAM MAKEPQW
C
C     THIS PROGRAM IS DATED AUGUST 11, 1995, AND GENERATES THE FORTRAN
C     COEFFICIENT FUNCTIONS P(X), Q(X), W(X), AND SUBROUTINE UV WHICH
C     DEFINES THE BOUNDARY CONDITION FUNCTIONS U(X), V(X) FOR SLEIGN2.
C
C     THE DIFFERENTIAL EQUATION IS OF THE FORM
C
C           -(p*y')' + q*y = lambda*w*y
C
C     .. Local Scalars ..
      CHARACTER*1 HQ
      CHARACTER*16 CHANS,TAPE1
      CHARACTER*62 STR
      DOUBLE PRECISION C
C     ..
C     .. External Subroutines ..
      EXTERNAL LC
C     ..
      WRITE(*,*)
      WRITE(*,*) ' HELP may be called at any point where the program   '
      WRITE(*,*) '    halts and displays (h?) by pressing "h <ENTER>". '
      WRITE(*,*) '    To RETURN from HELP, press "r <ENTER>".          '
      WRITE(*,*) '    To QUIT at any program halt, press "q <ENTER>".  '
      WRITE(*,*) ' WOULD YOU LIKE AN OVERVIEW OF HELP ? (Y/N) (h?)     '
      READ(*,1) HQ
      IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y' .OR.
     1            HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(1)
            END IF
C
  100 CONTINUE
         WRITE(*,*) ' SPECIFY OUTPUT FILE NAME (h?)  '
         READ(*,16) CHANS
         IF (CHANS.EQ.'q' .OR. CHANS.EQ.'Q') THEN
            STOP
         ELSE IF (CHANS.EQ.'h' .OR. CHANS.EQ.'H') THEN
            CALL HELP(2)
            GO TO 100
         ELSE
            TAPE1 = CHANS
            END IF
         OPEN(1,FILE=TAPE1,STATUS='NEW')
         WRITE(1,'(A)') 'C'
         WRITE(1,'(A)') 'C               ' // TAPE1
C
      WRITE(*,*) ' THE DIFFERENTIAL EQUATION IS OF THE FORM: '
      WRITE(*,*) '        -(p*y'')'' + q*y = lambda*w*y      '
      WRITE(*,*)
C
  200 CONTINUE
         WRITE(*,*) ' INPUT (h?) p = '
         READ(*,62) STR
         HQ = STR
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 200
         ELSE
            WRITE(1,'(A)') 'C'
            WRITE(1,'(A)') '      FUNCTION P(X)'
            WRITE(1,'(A)') '      DOUBLE PRECISION P,X'
            WRITE(1,'(A)') '      P = ' // STR
            WRITE(1,'(A)') '      RETURN'
            WRITE(1,'(A)') '      END'
         END IF
C
  300 CONTINUE
         WRITE(*,*) ' INPUT (h?) q = '
         READ(*,62) STR
         HQ = STR
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 300
         ELSE
            WRITE(1,'(A)') 'C'
            WRITE(1,'(A)') '      FUNCTION Q(X)'
            WRITE(1,'(A)') '      DOUBLE PRECISION Q,X'
            WRITE(1,'(A)') '      Q = ' // STR
            WRITE(1,'(A)') '      RETURN'
            WRITE(1,'(A)') '      END'
            END IF
C
  400 CONTINUE
         WRITE(*,*) ' INPUT (h?) w = '
         READ(*,62) STR
         HQ = STR
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(3)
            GO TO 400
         ELSE
            WRITE(1,'(A)') 'C'
            WRITE(1,'(A)') '      FUNCTION W(X)'
            WRITE(1,'(A)') '      DOUBLE PRECISION W,X'
            WRITE(1,'(A)') '      W = ' // STR
            WRITE(1,'(A)') '      RETURN'
            WRITE(1,'(A)') '      END'
            END IF
C
  500 CONTINUE
         WRITE(*,*) ' DO YOU REQUIRE INFORMATION ON END-POINT '
         WRITE(*,*) '    CLASSIFICATION ? (Y/N) (h?)  '
         READ(*,1) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y' .OR.
     1            HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(4)
            GO TO 500
         ELSE IF (HQ.EQ.'e' .OR. HQ.EQ.'E') THEN
            GO TO 800
         ELSE
            END IF
C
  600 CONTINUE
         WRITE(*,*) ' DO YOU REQUIRE INFORMATION ON DEFAULT CLASSIFI-'
         WRITE(*,*) '    CATION AND BOUNDARY CONDITIONS ? (Y/N) (h?)  '
         READ(*,1) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y' .OR.
     1            HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(5)
            GO TO 600
         ELSE
            END IF
C
  700 CONTINUE
         WRITE(*,*) ' DO YOU REQUIRE INFORMATION ON LIMIT CIRCLE '
         WRITE(*,*) '    BOUNDARY CONDITIONS ? (Y/N) (h?)  '
         READ(*,1) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y' .OR.
     1            HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(6)
            GO TO 700
         ELSE
            END IF
C
  800 CONTINUE
         WRITE(*,*) ' DO YOU WANT TO USE A LIMIT CIRCLE  '
         WRITE(*,*) '    BOUNDARY CONDITION ? (Y/N) (h?)  '
         READ(*,1) HQ
         IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
            STOP
         ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
            CALL HELP(6)
            GO TO 800
         ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y') THEN
            WRITE(1,'(A)') 'C'
            WRITE(1,'(A)') '      SUBROUTINE UV(X,U,PUP,V,PVP,HU,HV)'
            WRITE(1,'(A)') '      DOUBLE PRECISION X,U,PUP,V,PVP,HU,HV'
  900       CONTINUE
               WRITE(*,*) ' DO YOU WANT TO USE TWO DIFFERENT PAIRS OF '
               WRITE(*,*) '    FUNCTIONS U(X),V(X) ? (Y/N) (h?) '
               READ(*,1) HQ
               IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
                  STOP
               ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                  CALL HELP(6)
                  GO TO 900
               ELSE IF (HQ.EQ.'y' .OR. HQ.EQ.'Y') THEN
 1000             CONTINUE
                     WRITE(*,*) ' ASSUMING THAT ONE PAIR OF FUNCTIONS  '
                     WRITE(*,*) '    U(X),V(X) IS FOR a < X < c, AND   '
                     WRITE(*,*) '    THE OTHER PAIR IS FOR c <= X < b, '
                     WRITE(*,*) '    WHAT IS THE VALUE OF c ? (h?)     '
                     WRITE(*,*)
                     WRITE(*,*) ' c = '
                     READ(*,16) CHANS
                     HQ = CHANS
                     IF (HQ.EQ.'q' .OR. HQ.EQ.'Q') THEN
                        STOP
                     ELSE IF (HQ.EQ.'h' .OR. HQ.EQ.'H') THEN
                        CALL HELP(6)
                        GO TO 1000
                     ELSE
                        READ(CHANS,'(F16.0)') C
                        END IF
                  WRITE(*,*)
                  WRITE(*,*) ' FOR a < X < c :'
                  WRITE(*,*)
                  WRITE(1,'(A,1PE12.5,A)') '      IF (X.LT.',C,') THEN'
                  CALL LC(9)
                  WRITE(*,*)
                  WRITE(*,*) ' FOR c <= X < b :'
                  WRITE(*,*)
                  WRITE(1,'(A)') '      ELSE'
                  CALL LC(9)
                  WRITE(1,'(A)') '         END IF'
               ELSE
                  WRITE(*,*)
                  CALL LC(6)
                  END IF
         ELSE
            WRITE(1,'(A)') 'C'
            WRITE(1,'(A)') '      SUBROUTINE UV'
            END IF
      WRITE(1,'(A)') '      RETURN'
      WRITE(1,'(A)') '      END'
C
      WRITE(1,'(A)') 'C'
      WRITE(1,'(A)') '      SUBROUTINE EXAMP'
      WRITE(1,'(A)') '      RETURN'
      WRITE(1,'(A)') '      END'
      CLOSE(1)
      STOP
    1 FORMAT(A1)
   16 FORMAT(A16)
   62 FORMAT(A62)
      END
C
      SUBROUTINE LC(INDENT)
      INTEGER INDENT
C     .. Local Scalars ..
      CHARACTER*57 STR
C     ..
      IF (INDENT.EQ.6) THEN
         WRITE(*,*) ' INPUT u = '
         READ(*,57) STR
         WRITE(1,'(A)') '      U = ' // STR
         WRITE(*,*) ' INPUT v = '
         READ(*,57) STR
         WRITE(1,'(A)') '      V = ' // STR
         WRITE(*,*) ' INPUT pu'' = '
         READ(*,57) STR
         WRITE(1,'(A)') '      PUP = ' // STR
         WRITE(*,*) ' INPUT pv'' = '
         READ(*,57) STR
         WRITE(1,'(A)') '      PVP = ' // STR
         WRITE(*,*) ' INPUT -(pu'')'' + q*u = '
         READ(*,57) STR
         WRITE(1,'(A)') '      HU = ' // STR
         WRITE(*,*) ' INPUT -(pv'')'' + q*v = '
         READ(*,57) STR
         WRITE(1,'(A)') '      HV = ' // STR
      ELSE
         WRITE(*,*) ' INPUT u = '
         READ(*,57) STR
         WRITE(1,'(A)') '         U = ' // STR
         WRITE(*,*) ' INPUT v = '
         READ(*,57) STR
         WRITE(1,'(A)') '         V = ' // STR
         WRITE(*,*) ' INPUT pu'' = '
         READ(*,57) STR
         WRITE(1,'(A)') '         PUP = ' // STR
         WRITE(*,*) ' INPUT pv'' = '
         READ(*,57) STR
         WRITE(1,'(A)') '         PVP = ' // STR
         WRITE(*,*) ' INPUT -(pu'')'' + q*u = '
         READ(*,57) STR
         WRITE(1,'(A)') '         HU = ' // STR
         WRITE(*,*) ' INPUT -(pv'')'' + q*v = '
         READ(*,57) STR
         WRITE(1,'(A)') '         HV = ' // STR
         END IF
      RETURN
   57 FORMAT(A57)
      END
c
      subroutine help(nh)
      integer i,n,nh
      character*36 x(23),y(23)
      character*1 ans
c
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),NH
c
    1 CONTINUE
      write(*,*) 'H1:  Overview of HELP.'
      x(1)='  This ASCII text file is supplied a'
      y(1)='s a separate file with the SLEIGN2  '
      x(2)='package; it can be accessed on-line '
      y(2)='in both MAKEPQW (if used) and DRIVE.'
      x(3)='  HELP contains information to aid t'
      y(3)='he user in entering data on the     '
      x(4)='coefficient functions p,q,w; on the '
      y(4)='limit circle boundary condition     '
      x(5)='functions u,v; on the end-point clas'
      y(5)='sifications of the differential     '
      x(6)='equation; on DEFAULT entry; on eigen'
      y(6)='value indexes; on IFLAG information;'
      x(7)='and on the general use of the progra'
      y(7)='m SLEIGN2.                          '
      x(8)='  The 17 sections of HELP are:      '
      y(8)=' '
      x(9)=' '
      y(9)=' '
      x(10)='    H1: Overview of HELP.           '
      y(10)=' '
      x(11)='    H2: File name entry.            '
      y(11)=' '
      x(12)='    H3: The differential equation.  '
      y(12)=' '
      x(13)='    H4: End-point classification.   '
      y(13)=' '
      x(14)='    H5: DEFAULT entry.              '
      y(14)=' '
      x(15)='    H6: Limit-circle boundary condit'
      y(15)='ions.                               '
      do 101 i = 1,15
        write(*,*) x(i),y(i)
  101   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='    H7: General boundary conditions.'
      y(1)=' '
      x(2)='    H8: Recording the results.      '
      y(2)=' '
      x(3)='    H9: Type and choice of interval.'
      y(3)=' '
      x(4)='   H10: Entry of end-points.        '
      y(4)=' '
      x(5)='   H11: End-point values of p,q,w.  '
      y(5)=' '
      x(6)='   H12: Initial value problems.     '
      y(6)=' '
      x(7)='   H13: Indexing of eigenvalues.    '
      y(7)=' '
      x(8)='   H14: Entry of eigenvalue index, i'
      y(8)='nitial guess, and tolerance.     '
      x(9)='   H15: IFLAG information.          '
      y(9)=' '
      x(10)='   H16: Plotting.                   '
      y(10)=' '
      x(11)='   H17: Indexing of eigenvalues for '
      y(11)='periodic-type problems.             '
      x(12)=' '
      y(12)=' '
      x(13)='  HELP can be accessed at each point'
      y(13)=' in MAKEPQW and DRIVE where the user'
      x(14)='is asked for input, by pressing "h <'
      y(14)='ENTER>"; this places the user at the'
      x(15)='appropriate HELP section.  Once in H'
      y(15)='ELP, the user can scroll the further'
      x(16)='HELP sections by repeatedly pressing'
      y(16)=' "h <ENTER>", or jump to a specific '
      x(17)='HELP section Hn (n=1,2,...17) by typ'
      y(17)='ing "Hn <ENTER>"; to return to the  '
      x(18)='place in the program from which HELP'
      y(18)=' is called, press "r <ENTER>".      '
      do 102 i = 1,18
        write(*,*) x(i),y(i)
  102   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    2 CONTINUE
      write(*,*) 'H2:  File name entry.'
      x(1)='  MAKEPQW is used to create a FORTRA'
      y(1)='N file containing the coefficients  '
      x(2)='p(x),q(x),w(x), defining the differe'
      y(2)='ntial equation, and the boundary    '
      x(3)='condition functions u(x),v(x) if req'
      y(3)='uired.  The file must be given a NEW'
      x(4)='filename which is acceptable to your'
      y(4)=' FORTRAN compiler.  For example, it '
      x(5)='might be called bessel.f or bessel.f'
      y(5)='or depending upon your compiler.    '
      x(6)='  The same naming considerations app'
      y(6)='ly if the FORTRAN file is prepared  '
      x(7)='other than with the use of MAKEPQW. '
      y(7)=' '
      do 201 i = 1,7
        write(*,*) x(i),y(i)
  201   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    3 CONTINUE
      write(*,*) 'H3:  The differential equation.'
      x(1)='  The prompt "Input p (or q or w) ="'
      y(1)=' requests you to type in a FORTRAN  '
      x(2)='expression defining the function p(x'
      y(2)='), which is one of the three coeffi-'
      x(3)='cient functions defining the Sturm-L'
      y(3)='iouville differential equation      '
      x(4)=' '
      y(4)=' '
      x(5)='                   -(p*y'')'' + q*y = '
      y(5)=' lambda*w*y              (*)        '
      x(6)=' '
      y(6)=' '
      x(7)='to be considered on some interval (a'
      y(7)=',b) of the real line.  The actual   '
      x(8)='interval used in a particular proble'
      y(8)='m can be chosen later, and may be   '
      x(9)='either the whole interval (a,b) wher'
      y(9)='e the coefficient functions p,q,w,  '
      x(10)='etc. are defined or any sub-interval'
      y(10)=' (a'',b'') of (a,b); a = -infinity  '
      x(11)='and/or b = +infinity are allowable c'
      y(11)='hoices for the end-points.          '
      x(12)='  The coefficient functions p,q,w of'
      y(12)=' the differential equation may be   '
      x(13)='chosen arbitrarily but must satisfy '
      y(13)='the following conditions:           '
      x(14)='  (1) p,q,w are real-valued througho'
      y(14)='ut (a,b).                           '
      x(15)='  (2) p,q,w are piece-wise continuou'
      y(15)='s and defined throughout the        '
      x(16)='      interior of the interval (a,b)'
      y(16)='.                                   '
      x(17)='  (3) p and w are strictly positive '
      y(17)='in (a,b).                           '
      x(18)='  For better error analysis in the n'
      y(18)='umerical procedures, condition      '
      x(19)='(2) above is often replaced with    '
      y(19)=' '
      x(20)='  (2'') p,q,w are four times continuo'
      y(20)='usly differentiable on (a,b).       '
      x(21)='  The behavior of p,q,w near the end'
      y(21)='-points a and b is critical to the  '
      x(22)='classification of the differential e'
      y(22)='quation (see H4 and H11).           '
      do 301 i = 1,22
        write(*,*) x(i),y(i)
  301   continue
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    4 CONTINUE
      write(*,*) 'H4:  End-point classification.'
      x(1)='  The correct classification of the '
      y(1)='end-points a and b is essential to  '
      x(2)='the working of the SLEIGN2 program. '
      y(2)=' To classify the end-points, it is  '
      x(3)='convenient to choose a point c in (a'
      y(3)=',b); i.e.,  a < c < b.  Subject to  '
      x(4)='the general conditions on the coeffi'
      y(4)='cient functions p,q,w (see H3):     '
      x(5)='  (1) a is REGULAR (say R) if -infin'
      y(5)='ity < a, p,q,w are piece-wise       '
      x(6)='      continuous on [a,c], and p(a) '
      y(6)='> 0 and w(a) > 0.                   '
      x(7)='  (2) a is WEAKLY REGULAR (say WR) i'
      y(7)='f -infinity < a, a is not R, and    '
      x(8)='               |c                   '
      y(8)=' '
      x(9)='      integral | {1/p+|q|+w} < +infi'
      y(9)='nity.                               '
      x(10)='               |a                   '
      y(10)=' '
      x(11)=' '
      y(11)=' '
      x(12)='    If end-point a is neither R nor '
      y(12)='WR, then a is SINGULAR; that is,    '
      x(13)='  either -infinity = a, or -infinity'
      y(13)=' < a and                            '
      x(14)='               |c                   '
      y(14)=' '
      x(15)='      integral | {1/p+|q|+w} = +infi'
      y(15)='nity.                               '
      x(16)='               |a                   '
      y(16)=' '
      x(17)='  (3) The SINGULAR end-point a is LI'
      y(17)='MIT-CIRCLE NON-OSCILLATORY (say     '
      x(18)='      LCNO) if for some real lambda '
      y(18)='ALL real-valued solutions y of the  '
      x(19)='      differential equation         '
      y(19)=' '
      x(20)=' '
      y(20)=' '
      x(21)='                   -(p*y'')'' + q*y = '
      y(21)=' lambda*w*y on (a,c]     (*)        '
      x(22)=' '
      y(22)=' '
      x(23)='      satisfy the conditions:       '
      y(23)=' '
      do 401 i = 1,23
        write(*,*) x(i),y(i)
  401   continue
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='               |c                   '
      y(1)=' '
      x(2)='      integral | { w*y*y } < +infini'
      y(2)='ty, and                             '
      x(3)='               |a                   '
      y(3)=' '
      x(4)='      y has at most a finite number '
      y(4)='of zeros in (a,c].                  '
      x(5)='  (4) The SINGULAR end-point a is LI'
      y(5)='MIT-CIRCLE OSCILLATORY (say LCO) if '
      x(6)='for some real lambda ALL real-valued'
      y(6)=' solutions of the differential      '
      x(7)='equation (*) satisfy the conditions:'
      y(7)=' '
      x(8)='               |c                   '
      y(8)=' '
      x(9)='      integral | { w*y*y } < +infini'
      y(9)='ty, and                             '
      x(10)='               |a                   '
      y(10)=' '
      x(11)='      y has an infinite number of ze'
      y(11)='ros in (a,c].                       '
      x(12)='  (5) The SINGULAR end-point a is LI'
      y(12)='MIT POINT (say LP) if for some real '
      x(13)='lambda at least one solution of the '
      y(13)='differential equation (*) satisfies '
      x(14)='the condition:                      '
      y(14)=' '
      x(15)='               |c                   '
      y(15)=' '
      x(16)='      integral | {w*y*y} = +infinity'
      y(16)='.                                   '
      x(17)='               |a                   '
      y(17)=' '
      x(18)='  There is a similar classification '
      y(18)='of the end-point b into one of the  '
      x(19)='five distinct cases R, WR, LCNO, LCO'
      y(19)=', LP.                               '
      do 402 i = 1,19
        write(*,*) x(i),y(i)
  402   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  Although the classification of sin'
      y(1)='gular end-points invokes a real     '
      x(2)='value of the parameter lambda, this '
      y(2)='classification is invariant in      '
      x(3)='lambda; all real choices of lambda l'
      y(3)='ead to the same classification.     '
      x(4)='  In determining the classification '
      y(4)='of singular end-points for the      '
      x(5)='differential equation (*), it is oft'
      y(5)='en convenient to start with the     '
      x(6)='choice lambda = 0 in attempting to f'
      y(6)='ind solutions (particularly when    '
      x(7)='q = 0 on (a,b)); however, see exampl'
      y(7)='e 7 below.                          '
      x(8)='  See H6 on the use of maximal domai'
      y(8)='n functions to determine the        '
      x(9)='classification at singular end-point'
      y(9)='s.                                  '
      do 403 i = 1,9
        write(*,*) x(i),y(i)
  403   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      write(*,*) ' EXAMPLES: '
      x(1)='  1. -y'''' = lambda*y is R at both en'
      y(1)='d-points of (a,b) when a and b are  '
      x(2)='     finite.                        '
      y(2)=' '
      x(3)='  2. -y'''' = lambda*y on (-infinity,i'
      y(3)='nfinity) is LP at both end-points.  '
      x(4)='  3. -(sqrt(x)*y''(x))'' = lambda*(1./'
      y(4)='sqrt(x))*y(x) on (0,infinity) is    '
      x(5)='     WR at 0 and LP at +infinity (ta'
      y(5)='ke lambda = 0 in (*)).  See         '
      x(6)='     examples.f, #10 (Weakly Regular'
      y(6)=').                                  '
      x(7)='  4. -((1-x*x)*y''(x))'' = lambda*y(x)'
      y(7)=' on (-1,1) is LCNO at both ends     '
      x(8)='     (take lambda = 0 in (*)).  See '
      y(8)='xamples.f, #1 (Legendre).           '
      x(9)='  5. -y''''(x) + C*(1/(x*x))*y(x) = la'
      y(9)='mbda*y(x) on (0,infinity) is LP at  '
      x(10)='     infinity and at 0 is (take lamb'
      y(10)='da = 0 in (*)):                     '
      x(11)='       LP for C .ge. 3/4 ;          '
      y(11)=' '
      x(12)='       LCNO for -1/4 .le. C .lt. 3/4'
      y(12)=' (but C .ne. 0);                    '
      x(13)='       LCO for C .lt. -1/4.         '
      y(13)=' '
      x(14)='  6. -(x*y''(x))'' - (1/x)*y(x) = lamb'
      y(14)='da*y(x) on (0,infinity) is LCO at 0 '
      x(15)='     and LP at +infinity (take lambd'
      y(15)='a = 0 in (*) with solutions         '
      x(16)='     cos(ln(x)) and sin(ln(x))).  Se'
      y(16)='e xamples.f, #7 (BEZ).              '
      x(17)='  7. -(x*y''(x))'' - x*y(x) = lambda*('
      y(17)='1/x)*y(x) on (0,infinity) is LP at 0'
      x(18)='     and LCO at infinity (take lambd'
      y(18)='a = -1/4 in (*) with solutions      '
      x(19)='     cos(x)/sqrt(x) and sin(x)/sqrt('
      y(19)='x)).  See xamples.f,                '
      x(20)='     #6 (Sears-Titchmarsh).         '
      y(20)=' '
      do 404 i = 1,20
        write(*,*) x(i),y(i)
  404   continue
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    5 CONTINUE
      write(*,*) 'H5:  DEFAULT entry.'
      x(1)='  The complete range of problems for'
      y(1)=' which SLEIGN2 is applicable can    '
      x(2)='only be reached by appropriate entri'
      y(2)='es under end-point classification   '
      x(3)='and boundary conditions.  However, t'
      y(3)='here is a DEFAULT application which '
      x(4)='requires no detailed entry of end-po'
      y(4)='int classification or boundary      '
      x(5)='conditions, subject to:             '
      y(5)=' '
      x(6)='  1) The DEFAULT application CANNOT '
      y(6)='be used at a LCO end-point.         '
      x(7)='  2) If an end-point a is R, then th'
      y(7)='e Dirichlet boundary condition      '
      x(8)='     y(a) = 0 is automatically used.'
      y(8)=' '
      x(9)='  3) If an end-point a is WR, then t'
      y(9)='he following boundary condition     '
      x(10)='     is automatically applied:      '
      y(10)=' '
      x(11)='       if p(a) = 0, and both q(a),w('
      y(11)='a) are bounded, then the Dirichlet  '
      x(12)='       boundary condition y(a) = 0 i'
      y(12)='s used, or                          '
      x(13)='       if p(a) > 0, and q(a) and/or '
      y(13)='w(a)) are not bounded, then the     '
      x(14)='       Neumann boundary condition (p'
      y(14)='y'')(a) = 0 is used.                 '
      x(15)='     If p(a) = 0, and q(a) and/or w('
      y(15)='a) are not bounded, then no reliable'
      x(16)='     information can be given on the'
      y(16)=' DEFAULT boundary condition.        '
      x(17)='  4) If an end-point is LCNO, then i'
      y(17)='n most cases the principal or       '
      x(18)='     Friedrichs boundary condition i'
      y(18)='s applied (see H6).                 '
      x(19)='  5) If an end-point is LP, then the'
      y(19)=' normal LP procedure is applied     '
      x(20)='     (see H7(1.)).                  '
      y(20)=' '
      x(21)='If you choose the DEFAULT condition,'
      y(21)=' then no entry is required for the  '
      x(22)='u,v boundary condition functions.   '
      y(22)=' '
      do 501 i = 1,22
        write(*,*) x(i),y(i)
  501   continue
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    6 CONTINUE
      write(*,*) 'H6:  Limit-circle boundary conditions.'
      x(1)='  At an end-point a, the limit-circl'
      y(1)='e type separated boundary condition '
      x(2)='is of the form (similar remarks thro'
      y(2)='ughout apply to the end-point b)    '
      x(3)=' '
      y(3)=' '
      x(4)='       A1*[y,u](a) + A2*[y,v](a) = 0'
      y(4)=',                          (**)     '
      x(5)=' '
      y(5)=' '
      x(6)='where y is a solution of the differe'
      y(6)='ntial equation                      '
      x(7)=' '
      y(7)=' '
      x(8)='     -(p*y'')'' + q*y = lambda*w*y  on'
      y(8)=' (a,b).                     (*)     '
      x(9)=' '
      y(9)=' '
      x(10)='Here A1, A2 are real numbers; u and '
      y(10)='v are boundary condition functions; '
      x(11)='and for real-valued y and u the form'
      y(11)=' [y,u] is defined by                '
      x(12)=' '
      y(12)=' '
      x(13)='      [y,u](x) = y(x)*(pu'')(x) - u(x'
      y(13)=')*(py'')(x)   for x in (a,b).        '
      x(14)=' '
      y(14)=' '
      x(15)='  The object of this section is to p'
      y(15)='rovide help in choosing appropriate '
      x(16)='functions u and v in (**), given the'
      y(16)=' differential equation (*).  Full   '
      x(17)='details of the boundary conditions f'
      y(17)='or (*) are discussed in H7; here it '
      x(18)='is sufficient to say that the limit-'
      y(18)='circle type boundary condition (**) '
      x(19)='can be applied at any end-point in t'
      y(19)='he LCNO, LCO classification, but    '
      x(20)='also in the R, WR classification sub'
      y(20)='ject to the appropriate choice of   '
      x(21)='u and v.                            '
      y(21)=' '
      do 601 i = 1,21
        write(*,*) x(i),y(i)
  601   continue
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  Let (*) be R, WR, LCNO, or LCO at '
      y(1)='end-point a and choose c in (a,b).  '
      x(2)='Then either                         '
      y(2)=' '
      x(3)='    u and v are a pair of linearly i'
      y(3)='ndependent solutions of (*) on (a,c]'
      x(4)='  for any chosen real values of lamb'
      y(4)='da, or                              '
      x(5)='    u and v are a pair of real-value'
      y(5)='d maximal domain functions defined  '
      x(6)='  on (a,c] satisfying [u,v](a) .ne. '
      y(6)='0.  The maximal domain D(a,c] is    '
      x(7)='  defined by                        '
      y(7)=' '
      x(8)=' '
      y(8)=' '
      x(9)='       D(a,c] = {f:(a,c]->R:: f,pf'' '
      y(9)='in AC(a,c];                         '
      x(10)='                   f, ((-pf'')''+qf)/w'
      y(10)=' in L2((a,c;w)}                     '
      x(11)=' '
      y(11)=' '
      x(12)='  It is known that for all f,g in D('
      y(12)='a,c] the limit                      '
      x(13)=' '
      y(13)=' '
      x(14)='       [f,g](a) = lim[f,g](x) as x->'
      y(14)='a                                   '
      x(15)=' '
      y(15)=' '
      x(16)='  exists and is finite.  If (*) is L'
      y(16)='CNO or LCO at a, then all solutions '
      x(17)='  of (*) belong to D(a,c] for all va'
      y(17)='lues of lambda.                     '
      do 602 i = 1,17
        write(*,*) x(i),y(i)
  602   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  The boundary condition (**) is ess'
      y(1)='ential in the LCNO and LCO cases but'
      x(2)='can also be used with advantage in s'
      y(2)='ome R and WR cases.  In the R, WR, '
      x(3)='and LCNO cases, but not in the LCO c'
      y(3)='ase, the boundary condition         '
      x(4)='functions can always be chosen so th'
      y(4)='at                                  '
      x(5)='        lim u(x)/v(x) = 0 as x->a,  '
      y(5)=' '
      x(6)='and it is recommended that this norm'
      y(6)='alisation be effected; this has been'
      x(7)='done in the examples given below. In'
      y(7)=' this case, the boundary condition  '
      x(8)='[y,u](a) = 0 (i.e., A1 = 1, A2 = 0 i'
      y(8)='n (**)) is called the principal or  '
      x(9)='Friedrichs boundary condition.      '
      y(9)=' '
      x(10)=' '
      y(10)=' '
      x(11)='  In the case when end-points a and '
      y(11)='b are, independently, in R, WR,     '
      x(12)='LCNO, or LCO classification, it may '
      y(12)='be that symmetry or other reasons   '
      x(13)='permit one set of boundary condition'
      y(13)=' functions to be used at both end-  '
      x(14)='points (see xamples.f, #1 (Legendre)'
      y(14)=').  In other cases, different pairs '
      x(15)='must be chosen for each end-point (s'
      y(15)='ee xamples.f: #16 (Jacobi),         '
      x(16)='#18 (Dunsch), and #19 (Donsch)).    '
      y(16)=' '
      x(17)=' '
      y(17)=' '
      x(18)='  Note that a solution pair u,v is a'
      y(18)='lways a maximal domain pair, but not'
      x(19)='necessarily vice versa.             '
      y(19)=' '
      do 603 i = 1,19
        write(*,*) x(i),y(i)
  603   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)=' EXAMPLES:                          '
      y(1)=' '
      x(2)='1. -y''''(x) = lambda*y(x) on [0,pi] i'
      y(2)='s R at 0 and R at pi.               '
      x(3)='   At 0, with lambda = 0, a solution'
      y(3)=' pair is u(x) = x, v(x) = 1.        '
      x(4)='   At pi, with lambda = 1, a solutio'
      y(4)='n pair is                           '
      x(5)='     u(x) = sin(x), v(x) = cos(x).  '
      y(5) =' '
      x(6)='2. -(sqrt(x)*y''(x))'' = lambda*y(x)/s'
      y(6)='qrt(x) on (0,1] is                  '
      x(7)='     WR at 0 and R at 1.            '
      y(7)=' '
      x(8)='   (The general solutions of this eq'
      y(8)='uation are                          '
      x(9)='     u(x) = cos(2*sqrt(x*lambda)), v'
      y(9)='(x) = sin(2*sqrt(x*lambda)).)       '
      x(10)='   At 0, with lambda = 0, a solution'
      y(10)=' pair is                            '
      x(11)='     u(x) = 2*sqrt(x), v(x) = 1.    '
      y(11)=' '
      x(12)='   At 1, with lambda = pi*pi/4, a so'
      y(12)='lution pair is                      '
      x(13)='     u(x) = sin(pi*sqrt(x)), v(x) = '
      y(13)='cos(pi*sqrt(x)).                    '
      x(14)='   At 1, with lambda = 0, a solution'
      y(14)=' pair is                            '
      x(15)='     u(x) = 2*(1-sqrt(x)), v(x) = 1.'
      y(15)=' '
      x(16)='   See also xamples.f, #10 (Weakly R'
      y(16)='egular).                            '
      x(17)='3. -((1-x*x)*y''(x))'' = lambda*y(x) o'
      y(17)='n (-1,1) is LCNO at both ends.      '
      x(18)='   At +-1, with lambda = 0, a soluti'
      y(18)='on pair is                          '
      x(19)='     u(x) = 1, v(x) = 0.5*log((1+x)/'
      y(19)='(1-x)).                             '
      x(20)='   At 1, a maximal domain pair is u('
      y(20)='x) = 1, v(x) = log(1-x)             '
      x(21)='   At -1, a maximal domain pair is u'
      y(21)='(x) = 1, v(x) = log(1+x).           '
      x(22)='   See also xamples.f, #1 (Legendre)'
      y(22)='.                                   '
      do 604 i = 1,22
        write(*,*) x(i),y(i)
  604   continue
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='4. -y''''(x) - (1/(4x*x))*y(x) = lambd'
      y(1)='a*y(x) on (0,infinity) is           '
      x(2)='     LCNO at 0 and LP at +infinity. '
      y(2)=' '
      x(3)='   At 0, a maximal domain pair is   '
      y(3)=' '
      x(4)='     u(x) = sqrt(x), v(x) = sqrt(x)*'
      y(4)='log(x).                             '
      x(5)='   See also xamples.f, #2 (Bessel). '
      y(5)=' '
      x(6)='5. -y''''(x) - 5*(1/(4*x*x))*y(x) = la'
      y(6)='mbda*y(x) on (0,infinity) is        '
      x(7)='     LCO at 0 and LP at +infinity.  '
      y(7)=' '
      x(8)='   At 0, with lambda = 0, a solution'
      y(8)=' pair is                            '
      x(9)='     u(x) = sqrt(x)*cos(log(x)), v(x'
      y(9)=') = sqrt(x)*sin(log(x))             '
      x(10)='   See also xamples.f, #20 (Krall). '
      y(10)=' '
      x(11)='6. -y''''(x) - (1/x)*y(x) = lambda*y(x'
      y(11)=') on (0,infinity) is               '
      x(12)='     LCNO at 0 and LP at +infinity.'
      y(12)=' '
      x(13)='   At 0, a maximal domain pair is   '
      y(13)=' '
      x(14)='     u(x) = x, v(x) = 1 -x*log(x).  '
      y(14)=' '
      x(15)='   See also xamples.f, #4(Boyd).    '
      y(15)=' '
      x(16)='7. -((1/x)*y''(x))'' + (k/(x*x) + k*k/'
      y(16)='x)*y(x) = lambda*y(x) on (0,1],     '
      x(17)='     with k real and .ne. 0, is LCNO'
      y(17)=' at 0 and R at 1.                   '
      x(18)='   At 0, a maximal domain pair is   '
      y(18)=' '
      x(19)='     u(x) = x*x, v(x) = x - 1/k.    '
      y(19)=' '
      x(20)='   See also xamples.f, #8 (Laplace T'
      y(20)='idal Wave).                         '
      do 605 i = 1,20
        write(*,*) x(i),y(i)
  605   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    7 CONTINUE
      write(*,*) 'H7:  General boundary conditions.'
      x(1)='  Boundary conditions for Sturm-Liou'
      y(1)='ville boundary value problems       '
      x(2)=' '
      y(2)=' '
      x(3)='                   -(p*y'')'' + q*y = '
      y(3)=' lambda*w*y              (*)        '
      x(4)=' '
      y(4)=' '
      x(5)='on an interval (a,b) are either     '
      y(5)=' '
      x(6)='    SEPARATED, with at most one cond'
      y(6)='ition at end-point a and at most    '
      x(7)='  one condition at end-point b, or  '
      y(7)=' '
      x(8)='    COUPLED, when both a and b are, '
      y(8)='independently, in one of the end-   '
      x(9)='  point classifications R, WR, LCNO,'
      y(9)=' LCO, in which case two independent '
      x(10)='  ent boundary conditions are requir'
      y(10)='ed which link the solution values   '
      x(11)='  near a to those near b.           '
      y(11)=' '
      x(12)='The SLEIGN2 program allows for all s'
      y(12)='eparated conditions; and special    '
      x(13)='cases of the coupled conditions -- t'
      y(13)='he so-called periodic boundary      '
      x(14)='conditions applicable only when the '
      y(14)='interval (a,b) is finite and both   '
      x(15)='a and b are R.                      '
      y(15)=' '
      do 701 i = 1,15
        write(*,*) x(i),y(i)
  701   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      write(*,*) '    Separated Conditions:      '
      write(*,*) '    ---------------------      '
      x(1)='  The boundary conditions to be sele'
      y(1)='cted depend upon the classification '
      x(2)='of the differential equation at the '
      y(2)='end-point, say, a:                  '
      x(3)='    1. If the end-point a is LP, the'
      y(3)='n no boundary condition is required '
      x(4)='  or allowed.                       '
      y(4)=' '
      x(5)='    2. If the end-point a is R or WR'
      y(5)=', then a separated boundary         '
      x(6)='  condition is of the form          '
      y(6)=' '
      x(7)='         A1*y(a) + A2*(py'')(a) = 0,'
      y(7)=' '
      x(8)='  where A1, A2 are real constants yo'
      y(8)='u must choose, not both zero.       '
      x(9)='    3. If the end-point a is LCNO or'
      y(9)=' LCO, then a separated boundary     '
      x(10)='  condition is of the form          '
      y(10)=' '
      x(11)='         A1*[y,u](a) + A2*[y,v](a) ='
      y(11)=' 0,                                 '
      x(12)='  where A1, A2 are real constants yo'
      y(12)='u must choose, not both zero;       '
      x(13)='  here, u,v are the pair of boundary'
      y(13)=' condition functions you have       '
      x(14)='  previously selected when the input'
      y(14)=' FORTRAN file was being prepared.   '
      x(15)='    4. If the end-point a is LCNO an'
      y(15)='d the boundary condition pair       '
      x(16)='  u,v has been chosen so that       '
      y(16)=' '
      x(17)='         lim u(x)/v(x) = 0  as x->a '
      y(17)=' '
      x(18)='  (which is always possible), then A'
      y(18)='1 = 1, A2 = 0 (i.e., [y,u](a) = 0)  '
      x(19)='  gives the principal (Friedrichs) b'
      y(19)='oundary condition at a.             '
      do 702 i = 1,19
        write(*,*) x(i),y(i)
  702   continue
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='    5. If a is R or WR and boundary '
      y(1)='condition functions u,v have been   '
      x(2)='  entered in the FORTRAN input file,'
      y(2)=' then (3.,4.) above apply to        '
      x(3)='  entering separated boundary condit'
      y(3)='ions at such an end-point; the      '
      x(4)='  boundary conditions in this form a'
      y(4)='re equivalent to the point-wise     '
      x(5)='  conditions in (2.) (subject to car'
      y(5)='e in choosing A1, A2).  This        '
      x(6)='  singular form of a regular boundar'
      y(6)='y condition may be particularly     '
      x(7)='  effective in the WR case if the bo'
      y(7)='undary condition form in (2.) leads '
      x(8)='  to numerical difficulties.        '
      y(8)=' '
      x(9)=' '
      y(9)=' '
      x(10)='  Conditions (2.,3.,4.,5.) apply sim'
      y(10)='ilarly at end-point b.              '
      x(11)=' '
      y(11)=' '
      x(12)='    6. If a is R, WR, LCNO, or LCO a'
      y(12)='nd b is LP, then only a separated   '
      x(13)='  condition at a is required and all'
      y(13)='owed (or instead at b if a and b    '
      x(14)='  are interchanged).                '
      y(14)=' '
      x(15)='    7. If both end-points a and b ar'
      y(15)='e LP, then no boundary conditions   '
      x(16)='  are required or allowed.          '
      y(16)=' '
      x(17)=' '
      y(17)=' '
      x(18)='  The indexing of eigenvalues for bo'
      y(18)='undary value problems with separated'
      x(19)='  conditions is discussed in H13.   '
      y(19)=' '
      do 703 i = 1,19
        write(*,*) x(i),y(i)
  703   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),n
      write(*,*) '    Coupled Conditions:        '
      write(*,*) '    -------------------        '
      x(1)='    8. Periodic-type boundary condit'
      y(1)='ions on (a,b) apply only when       '
      x(2)='  both end-points a and b are R; the'
      y(2)='se conditions are of the form       '
      x(3)=' '
      y(3)=' '
      x(4)='        y(b) = c*y(a),  (py'')(b) = ('
      y(4)='py'')(a)/c,                          '
      x(5)=' '
      y(5)=' '
      x(6)='  where c may be chosen to be any re'
      y(6)='al number not equal to 0.  The case '
      x(7)='  c = 1 is called periodic, the case'
      y(7)=' c = -1 is called semi-periodic.    '
      x(8)=' '
      y(8)=' '
      x(9)='  The indexing of eigenvalues for pe'
      y(9)='riodic-type boundary conditions is  '
      x(10)='  discussed in H17.                 '
      y(10)=' '
      do 704 i = 1,10
        write(*,*) x(i),y(i)
  704   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    8 CONTINUE
      write(*,*) 'H8:  Recording the results.'
      x(1)='  If you choose to have a record kep'
      y(1)='t of the results, then the following'
      x(2)='information is stored in a file with'
      y(2)=' the name you select:               '
      x(3)=' '
      y(3)=' '
      x(4)='  1. The file name.                 '
      y(4)=' '
      x(5)='  2. The header line prompted for (u'
      y(5)='p to 32 characters of your choice). '
      x(6)='  3. The interval (a,b) which was us'
      y(6)='ed.                                 '
      x(7)=' '
      y(7)=' '
      x(8)='  For SEPARATED boundary conditions:'
      y(8)=' '
      x(9)='  4. The end-point classification.  '
      y(9)=' '
      x(10)='  5. A summary of coefficient inform'
      y(10)='ation at WR, LCNO, LCO end-points.  '
      x(11)='  6. The boundary condition constant'
      y(11)='s (A1,A2), (B1,B2) if entered.      '
      x(12)='  7. (NUMEIG,EIG,TOL) or (NUMEIG1,NU'
      y(12)='MEIG2,TOL), as entered.             '
      x(13)=' '
      y(13)=' '
      x(14)='  For COUPLED boundary conditions:  '
      y(14)=' '
      x(15)='  8. The boundary condition paramete'
      y(15)='r, c.                               '
      x(16)=' '
      y(16)=' '
      x(17)='  For ALL boundary conditions:      '
      y(17)=' '
      x(18)='  9. The computed eigenvalue, EIG, a'
      y(18)='nd its estimated accuracy, TOL.     '
      x(19)=' 10. IFLAG reported (see H15).      '
      y(19)=' '
      do 801 i = 1,19
        write(*,*) x(i),y(i)
  801   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
    9 CONTINUE
      write(*,*) 'H9:  Type and choice of interval.'
      x(1)='  You may enter any interval (a,b) f'
      y(1)='or which the coefficients p,q,w are '
      x(2)='well defined by your FORTRAN stateme'
      y(2)='nts in the input file, provided that'
      x(3)='(a,b) contains no interior singulari'
      y(3)='ities.                              '
      do 901 i = 1,3
        write(*,*) x(i),y(i)
  901   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   10 CONTINUE
      write(*,*) 'H10:  Entry of end-points.'
      x(1)='  End-points a and b must be entered'
      y(1)=' as real numbers.  There is no      '
      x(2)='symbolic entry; e.g., pi must be ent'
      y(2)='ered as 3.14159... to an appropriate'
      x(3)='number of decimal places.           '
      y(3)=' '
      do 1001 i = 1,3
        write(*,*) x(i),y(i)
 1001   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   11 CONTINUE
      write(*,*) 'H11:  End-point values of p,q,w.'
      x(1)='  The program SLEIGN2 needs to know '
      y(1)='whether the coefficient functions   '
      x(2)='p(x),q(x),w(x) defined by the FORTRA'
      y(2)='N expressions entered in the input  '
      x(3)='file can be evaluated numerically wi'
      y(3)='thout running into difficulty.  If, '
      x(4)='for example, either q or w is unboun'
      y(4)='ded at a, or p(a) is 0, then SLEIGN2'
      x(5)='needs to know this so that a is not '
      y(5)='chosen for functional evaluation.   '
      do 1101 i = 1,5
        write(*,*) x(i),y(i)
 1101   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   12 CONTINUE
      write(*,*) 'H12:  Initial value problems.'
      x(1)='  The initial value problem facility'
      y(1)=' for Sturm-Liouville problems       '
      x(2)=' '
      y(2)=' '
      x(3)='                   -(p*y'')'' + q*y = '
      y(3)=' lambda*w*y              (*)        '
      x(4)=' '
      y(4)=' '
      x(5)='allows for the computation of a solu'
      y(5)='tion of (*) with a user-chosen      '
      x(6)='value lambda and any one of the foll'
      y(6)='owing initial conditions:           '
      x(7)='  1. From end-point a of any classif'
      y(7)='ication except LP towards           '
      x(8)='end-point b of any classification,  '
      y(8)=' '
      x(9)='  2. From end-point b of any classif'
      y(9)='ication except LP back towards      '
      x(10)='end-point a of any classification,  '
      y(10)=' '
      x(11)='  3. From end-points a and b of any '
      y(11)='classifications except LP towards an'
      x(12)='interior point of (a,b) selected by '
      y(12)='the program.                        '
      x(13)=' '
      y(13)=' '
      x(14)='  Initial values at a are of the for'
      y(14)='m y(a) = alpha1, (p*y'')a = alpha2, '
      x(15)='when a is R or WR; and [y,u](a) = al'
      y(15)='pha1, [y,v](a) = alpha2, when a is  '
      x(16)='LCNO or LCO.                        '
      y(16)=' '
      x(17)='  Initial values at b are of the for'
      y(17)='m y(b) = beta1, (p*y'')b = beta2,   '
      x(18)='when b is R or WR; and [y,u](b) = be'
      y(18)='ta1, [y,v](b) = beta2, when b is  '
      x(19)='LCNO or LCO.                        '
      y(19)=' '
      x(20)='  In (*), lambda is a user-chosen re'
      y(20)='al number; while in the above       '
      x(21)='initial values, (alpha1,alpha2) and '
      y(21)='(beta1,beta2) are user-chosen pairs '
      x(22)='of real numbers not both zero.      '
      y(22)=' '
      do 1201 i = 1,22
        write(*,*) x(i),y(i)
 1201   continue
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  In the initial value case (3.) abo'
      y(1)='ve when the interval (a,b) is       '
      x(2)='finite, the interior point selected '
      y(2)='by the program is generally near the'
      x(3)='midpoint of (a,b); when (a,b) is inf'
      y(3)='inite, no general rule can be given.'
      x(4)='Also if, given (alpha1,alpha2) and ('
      y(4)='beta1,beta2), the lambda chosen is  '
      x(5)='an eigenvalue of the associated boun'
      y(5)='dary value problem, the computed    '
      x(6)='solution may not be the correspondin'
      y(6)='g eigenfunction -- the signs of the '
      x(7)='computed solutions on either side of'
      y(7)=' the interior point may be opposite.'
      x(8)='  The output for a solution of an in'
      y(8)='itial value problem is in the form  '
      x(9)='of stored numerical data which can b'
      y(9)='e plotted on the screen (see H16),  '
      x(10)='or printed out in graphical form if '
      y(10)='graphics software is available.     '
      do 1202 i = 1,10
        write(*,*) x(i),y(i)
 1202   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   13 continue
      write(*,*) 'H13:  Indexing of eigenvalues.'
      x(1)='  The indexing of eigenvalues is an '
      y(1)='automatic facility in SLEIGN2.  The '
      x(2)='following general results hold for t'
      y(2)='he separated boundary condition     '
      x(3)='problem (see H7):                   '
      y(3)=' '
      x(4)='  1. If neither end-point a or b is '
      y(4)='LP or LCO, then the spectrum of the '
      x(5)='eigenvalue problem is discrete (eige'
      y(5)='nvalues only), simple (eigenvalues  '
      x(6)='all of multiplicity 1), and bounded '
      y(6)='below with a single cluster point at'
      x(7)='+infinity.  The eigenvalues are inde'
      y(7)='xed as {lambda(n): n=0,1,2,...},    '
      x(8)='where lambda(n) < lambda(n+1) (n=0,1'
      y(8)=',2,...), lim lambda(n) -> +infinity;'
      x(9)='and if {psi(n): n=0,1,2,...} are the'
      y(9)=' corresponding eigenfunctions, then '
      x(10)='psi(n) has exactly n zeros in the op'
      y(10)='en interval (a,b).                 '
      x(11)='  2. If neither end-point a or b is '
      y(11)='LP but at least one end-point is    '
      x(12)='LCO, then the spectrum is discrete a'
      y(12)='nd simple as for (1.), but with     '
      x(13)='cluster points at both +infinity and'
      y(13)=' -infinity.  The eigenvalues are    '
      x(14)='indexed as {lambda(n): n=0,1,-1,2,-2'
      y(14)=',...}, where                        '
      x(15)='lambda(n) < lambda(n+1) (n=...-2,-1,'
      y(15)=',0,1,2,...) with lambda(0) the      '
      x(16)='smallest non-negative eigenvalue and'
      y(16)=' lim lambda(n) -> +infinity or      '
      x(17)='-> -infinity with n; and if {psi(n):'
      y(17)=' n=0,1,-1,2,-2,...} are the         '
      x(18)='corresponding eigenfunctions, then e'
      y(18)='very psi(n) has infinitely many     '
      x(19)='zeros in (a,b).                     '
      y(19)=' '
      do 1301 i = 1,19
        write(*,*) x(i),y(i)
 1301   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  3. If one or both end-points is LP'
      y(1)=', then there can be one or more     '
      x(2)='intervals of continuous spectrum for'
      y(2)=' the boundary value problem in      '
      x(3)='addition to some (necessarily simple'
      y(3)=') eigenvalues.  For these           '
      x(4)='essentially more difficult problems,'
      y(4)=' SLEIGN2 can be used as an          '
      x(5)='investigative tool to give qualitati'
      y(5)='ve and possibly quantitative        '
      x(6)='information on the spectrum.        '
      y(6)=' '
      x(7)='     For example, if a problem has a'
      y(7)=' single interval of continuous      '
      x(8)='spectrum bounded below by K, then th'
      y(8)='ere may be any number of eigenvalues'
      x(9)='below K.  In some cases, SLEIGN2 can'
      y(9)=' compute K, and determine the number'
      x(10)='of these eigenvalues and compute the'
      y(10)='m.  In this respect, see xamples.f: '
      x(11)='#13 (Hydrogen Atom), #17 (Morse Osci'
      y(11)='llator), #21 (Fourier), and         '
      x(12)='#27 (Joergens) as examples of succes'
      y(12)='s; and #2 (Mathieu), #14 (Marletta),'
      x(13)='and #28 (Behnke-Goerisch) as example'
      y(13)='s of failure.                       '
      x(14)='     The problem need not have a con'
      y(14)='tinuous spectrum, in which case if  '
      x(15)='its discrete spectrum is bounded bel'
      y(15)='ow, then the eigenvalues are indexed'
      x(16)='and the eigenfunctions have zero cou'
      y(16)='nts as in (1.).  If, on the other   '
      x(17)='hand, the discrete spectrum is unbou'
      y(17)='nded below, then all the            '
      x(18)='eigenfunctions have infinitely many '
      y(18)='zeros in the interval.              '
      do 1302 i = 1,18
        write(*,*) x(i),y(i)
 1302   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  In respect to the three classes of'
      y(1)=' end-points, the following          '
      x(2)='identified examples from xamples.f i'
      y(2)='llustrate the spectral property     '
      x(3)='of these boundary value problems:   '
      y(3)=' '
      x(4)='  1. Neither end-point is LP or LCO.'
      y(4)=' '
      x(5)='       #1 (Legendre)                '
      y(5)=' '
      x(6)='       #2 (Bessel) with -1/4 < c < 3'
      y(6)='/4                                  '
      x(7)='       #4 (Boyd)                    '
      y(7)=' '
      x(8)='       #5 (Latzko)                  '
      y(8)=' '
      x(9)='  2. Neither end-point is LP, but at'
      y(9)=' least one is LCO.                  '
      x(10)='       #6 (Sears-Titchmarsh)        '
      y(10)=' '
      x(11)='       #7 (BEZ)                     '
      y(11)=' '
      x(12)='      #19 (Donsch)                  '
      y(12)=' '
      x(13)='  3. At least one end-point is LP.  '
      y(13)=' '
      x(14)='      #13 (Hydrogen Atom)           '
      y(14)=' '
      x(15)='      #14 (Marletta)                '
      y(15)=' '
      x(16)='      #20 (Krall)                   '
      y(16)=' '
      x(17)='      #21 (Fourier) on [0,infinity) '
      y(17)=' '
      do 1303 i = 1,17
        write(*,*) x(i),y(i)
 1303   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   14 CONTINUE
      write(*,*) 'H14:  Entry of eigenvalue index, initial guess,'//
     1           ' and tolerance.'
      x(1)='  For SEPARATED boundary condition p'
      y(1)='roblems (see H7), SLEIGN2 calls for '
      x(2)='input information options to compute'
      y(2)=' either                             '
      x(3)='  1. a single eigenvalue, or        '
      y(3)=' '
      x(4)='  2. a series of eigenvalues.       '
      y(4)=' '
      x(5)='In each case indexing of eigenvalues'
      y(5)=' is called for (see H13).           '
      x(6)='  (1.) above asks for data triples N'
      y(6)='UMEIG, EIG, TOL separated by commas.'
      x(7)='Here NUMEIG is the integer index of '
      y(7)='the desired eigenvalue; NUMEIG can  '
      x(8)='be negative only when the problem is'
      y(8)=' LCO at one or both end-points.     '
      x(9)='EIG allows for the entry of an initi'
      y(9)='al guess for the requested          '
      x(10)='eigenvalue (if an especially good on'
      y(10)='e is available), or can be set to 0 '
      x(11)='in which case an initial guess is ge'
      y(11)='nerated by SLEIGN2 itself.          '
      x(12)='TOL is the desired accuracy of the c'
      y(12)='omputed eigenvalue.  It is an       '
      x(13)='absolute accuracy if the magnitude o'
      y(13)='f the eigenvalue is 1 or less, and  '
      x(14)='is a relative accuracy otherwise.  T'
      y(14)='ypical values might be .001 for     '
      x(15)='moderate accuracy and .0000001 for h'
      y(15)='igh accuracy in single precision.   '
      x(16)='If TOL is set to 0, the maximum achi'
      y(16)='evable accuracy is requested.       '
      x(17)='  If the input data list is truncate'
      y(17)='d with a "/" after NUMEIG or EIG,   '
      x(18)='then the remaining elements default '
      y(18)='to 0.                               '
      do 1401 i = 1,18
        write(*,*) x(i),y(i)
 1401   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      x(1)='  (2.) above asks for data triples N'
      y(1)='UMEIG1, NUMEIG2, TOL separated by   '
      x(2)='commas.  Here NUMEIG1 and NUMEIG2 ar'
      y(2)='e the first and last integer indices'
      x(3)='of the sequence of desired eigenvalu'
      y(3)='es, NUMEIG1 < NUMEIG2; they can be  '
      x(4)='negative only when the problem is LC'
      y(4)='O at one or both end-points.        '
      x(5)='TOL is the desired accuracy of the c'
      y(5)='omputed eigenvalues.  It is an      '
      x(6)='absolute accuracy if the magnitude o'
      y(6)='f an eigenvalue is 1 or less, and   '
      x(7)='is a relative accuracy otherwise.  T'
      y(7)='ypical values might be .001 for     '
      x(8)='moderate accuracy and .0000001 for h'
      y(8)='igh accuracy in single precision.   '
      x(9)='If TOL is set to 0, the maximum achi'
      y(9)='evable accuracy is requested.       '
      x(10)='  If the input data list is truncate'
      y(10)='d with a "/" after NUMEIG2, then TOL'
      x(11)='defaults to 0.                      '
      y(11)=' '
      x(12)=' '
      y(12)=' '
      x(13)='  For COUPLED periodic-type boundary'
      y(13)=' condition problems (see H7 and     '
      x(14)='H17), SLEIGN2 asks only for NUMEIG, '
      y(14)='the non-negative integer index of   '
      x(15)='the desired eigenvalue; TOL is set i'
      y(15)='nternally.                          '
      do 1402 i = 1,15
        write(*,*) x(i),y(i)
 1402   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   15 CONTINUE
      write(*,*) 'H15:  IFLAG information.'
      x(1)='  All results are reported by SLEIGN'
      y(1)='2 with a flag identification.  There'
      x(2)='are four values of IFLAG:           '
      y(2)=' '
      x(3)=' '
      y(3)=' '
      x(4)='  1 - The computed eigenvalue has an'
      y(4)=' estimated accuracy within the      '
      x(5)='      tolerance requested.          '
      y(5)=' '
      x(6)=' '
      y(6)=' '
      x(7)='  2 - The computed eigenvalue does n'
      y(7)='ot have an estimated accuracy within'
      x(8)='      the tolerance requested, but i'
      y(8)='s the best the program could obtain.'
      x(9)=' '
      y(9)=' '
      x(10)='  3 - There seems to be no eigenvalu'
      y(10)='e of index equal to NUMEIG.         '
      x(11)=' '
      y(11)=' '
      x(12)='  4 - The program has been unable to'
      y(12)=' compute the requested eigenvalue.  '
      do 1501 i = 1,12
        write(*,*) x(i),y(i)
 1501   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   16 CONTINUE
      write(*,*) 'H16:  Plotting.'
      x(1)='  After computing a single eigenvalu'
      y(1)='e (see H14(1.)), but not a sequence '
      x(2)='of eigenvalues (see H14(2.)), the ei'
      y(2)='genfunction can be plotted.  If this'
      x(3)='is desired, respond "y" when asked s'
      y(3)='o that SLEIGN2 will compute some    '
      x(4)='eigenfunction data and store them.  '
      y(4)=' '
      x(5)='  One can ask that the eigenfunction'
      y(5)=' data be in the form of either      '
      x(6)='points (x,y) for x in (a,b), or poin'
      y(6)='ts (t,y) for t in the standardized  '
      x(7)='interval (-1,1) mapped onto from (a,'
      y(7)='b); the t- choice can be especially '
      x(8)='helpful when the original interval i'
      y(8)='s infinite.  Additionally, one can  '
      x(9)='ask for a plot of the so-called Pruf'
      y(9)='er angle, in x- or t- variables.    '
      x(10)='  In both forms, once the choice has'
      y(10)=' been made of the function to be    '
      x(11)='plotted, a crude plot is displayed o'
      y(11)='n the monitor screen and you are    '
      x(12)='asked whether you wish to save the c'
      y(12)='omputed plot points in a file. '
      do 1601 i = 1,12
        write(*,*) x(i),y(i)
 1601   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
   17 CONTINUE
      write(*,*) 'H17:  Indexing of eigenvalues for'//
     1           ' periodic-type problems.'
      x(1)='  The indexing of eigenvalues is an '
      y(1)='automatic facility in SLEIGN2.  The '
      x(2)='following general result holds for t'
      y(2)='he periodic-type boundary condition '
      x(3)='problem (see H7):                   '
      y(3)=' '
      x(4)='  The spectrum of the eigenvalue pro'
      y(4)='blem is discrete (eigenvalues only),'
      x(5)='and bounded below with a single clus'
      y(5)='ter point at +infinity.  In general,'
      x(6)='the spectrum is not simple, but no e'
      y(6)='igenvalue exceeds multiplicity 2.   '
      x(7)='The eigenvalues are indexed as {lamb'
      y(7)='da(n): n=0,1,2,...}, where          '
      x(8)='lambda(n) .le. lambda(n+1) (n=0,1,2,'
      y(8)='...), lim lambda(n) -> +infinity.   '
      x(9)='  The connection between the index n'
      y(9)=' and the number of zeros of the     '
      x(10)='corresponding eigenfunction psi(x,n)'
      y(10)=' is not as simple as for separated  '
      x(11)='conditions; that is, psi(x,n) need n'
      y(11)='ot have exactly n zeros in (a,b).   '
      x(12)='  The following identified examples '
      y(12)='from xamples.f are of special       '
      x(13)='interest:                           '
      y(13)=' '
      x(14)='    #11 (Plum) on [0,pi]            '
      y(14)=' '
      x(15)='    #21 (Fourier) on [0,pi]         '
      y(15)=' '
      x(16)='    #25 (Meissner) on [-0.5,0.5]    '
      y(16)=' '
      do 1701 i = 1,16
        write(*,*) x(i),y(i)
 1701   continue
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) '-----------------------------------------------'
      read(*,999) ans,n
      if (ans.eq.'r' .or. ans.eq.'R') return
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),N
      go to 1
  999 FORMAT(A1,I2)
      end

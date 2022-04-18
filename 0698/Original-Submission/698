#!/bin/sh
# This is a shell archive (produced by GNU sharutils 4.1).
# To extract the files from this archive, save it to some FILE, remove
# everything before the `!/bin/sh' line above, then type `sh FILE'.
#
# Existing files will *not* be overwritten unless `-c' is specified.
#
# This shar contains:
# length mode       name
# ------ ---------- ------------------------------------------
#   1456 -rw-r--r-- READ.ME
#   1977 -rw-r--r-- makefile
#   1559 -rw-r--r-- dtest1.f
#   8207 -rw-r--r-- dtest2.f
#  17921 -rw-r--r-- dcuhre.f
#   8328 -rw-r--r-- dchhre.f
#  16631 -rw-r--r-- dadhre.f
#   5933 -rw-r--r-- dtrhre.f
#   4084 -rw-r--r-- dinhre.f
#   5952 -rw-r--r-- d132re.f
#   6203 -rw-r--r-- d113re.f
#   7863 -rw-r--r-- d09hre.f
#   4921 -rw-r--r-- d07hre.f
#   7938 -rw-r--r-- drlhre.f
#   3846 -rw-r--r-- dfshre.f
#   1524 -rw-r--r-- stest1.f
#   7616 -rw-r--r-- stest2.f
#  17849 -rw-r--r-- scuhre.f
#   8316 -rw-r--r-- schhre.f
#  16487 -rw-r--r-- sadhre.f
#   5837 -rw-r--r-- strhre.f
#   4036 -rw-r--r-- sinhre.f
#   5904 -rw-r--r-- s132re.f
#   6155 -rw-r--r-- s113re.f
#   7827 -rw-r--r-- s09hre.f
#   4885 -rw-r--r-- s07hre.f
#   7902 -rw-r--r-- srlhre.f
#   3822 -rw-r--r-- sfshre.f
#
touch -am 1231235999 $$.touch >/dev/null 2>&1
if test ! -f 1231235999 && test -f $$.touch; then
  shar_touch=touch
else
  shar_touch=:
  echo
  echo 'WARNING: not restoring timestamps.  Consider getting and'
  echo "installing GNU \`touch', distributed in GNU File Utilities..."
  echo
fi
rm -f 1231235999 $$.touch
#
# ============= READ.ME ==============
if test -f 'READ.ME' && test X"$1" != X"-c"; then
  echo 'x - skipping READ.ME (file already exists)'
else
  echo 'x - extracting READ.ME (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'READ.ME' &&
C      ALGORITHM 698 , COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 17, NO. 4, DECEMBER, 1991, PP. 452-456.
C
C     This file summarizes the DOUBLE PRECISION routines
C     in this directory. SINGLE PRECISION names are obtained
C     by replacing the leading D with S in all names.
C
C     DTEST1   Simple test driver for DCUHRE.
C              Sample output from a SUN 3/50 is included.
C     DTEST2   More detailed test driver for DCUHRE.
C              Sample output from a SUN 3/50 is included.
C     DCUHRE   Main Integrator.DCUHRE calls DCHHRE and DADHRE.
C     DCHHRE   Checks the input to DCUHRE.
C     DADHRE   The adaptive integration routine.
C              DADHRE calls DTRHRE, DINHRE and DRLHRE.
C     DTRHRE   Maintaines the heap of subtriangles.
C     DINHRE   Computes weights and abscissas of the integration
C              rule. DINHRE calls D132RE, D112RE, D09HRE and D07HRE.
C     D132RE   Computes weights and abscissas for a 2-dimensional
C              rule of degree 13.
C     D113RE   Computes weights and abscissas for a 3-dimensional
C              rule of degree 11.
C     D09HRE   Computes weights and abscissas for a degree 9 rule.
C     D07HRE   Computes weights and abscissas for a degree 7 rule.
C     DRLHRE   Computes estimates of integral and error over
C              subtriangles.
C     DFSHRE   Computes fully symmetric sums of function evaluations.
SHAR_EOF
  $shar_touch -am 0531091995 'READ.ME' &&
  chmod 0644 'READ.ME' ||
  echo 'restore of READ.ME failed'
  shar_count="`wc -c < 'READ.ME'`"
  test 1456 -eq "$shar_count" ||
    echo "READ.ME: original size 1456, current size $shar_count"
fi
# ============= makefile ==============
if test -f 'makefile' && test X"$1" != X"-c"; then
  echo 'x - skipping makefile (file already exists)'
else
  echo 'x - extracting makefile (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'makefile' &&
dtest1: dtest1.o dcuhre.o d07hre.o d09hre.o d113re.o d132re.o \
X        dadhre.o dchhre.o dfshre.o dinhre.o drlhre.o dtrhre.o
X	f77 -o dtest1 dtest1.o dcuhre.o d07hre.o d09hre.o \
X        d113re.o d132re.o \
X        dadhre.o dchhre.o dfshre.o dinhre.o drlhre.o dtrhre.o
X 
dtest2: dtest2.o dcuhre.o d07hre.o d09hre.o d113re.o d132re.o \
X        dadhre.o dchhre.o dfshre.o dinhre.o drlhre.o dtrhre.o
X	f77 -o dtest2 dtest2.o dcuhre.o d07hre.o d09hre.o \
X        d113re.o d132re.o \
X        dadhre.o dchhre.o dfshre.o dinhre.o drlhre.o dtrhre.o
X 
dtest1.o: dtest1.f
X	f77 -c dtest1.f
dtest2.o: dtest2.f
X	f77 -c dtest2.f
dcuhre.o: dcuhre.f
X	f77 -c dcuhre.f
d07hre.o: d07hre.f
X	f77 -c d07hre.f
d09hre.o: d09hre.f
X	f77 -c d09hre.f
d113re.o: d113re.f
X	f77 -c d113re.f
d132re.o: d132re.f
X	f77 -c d132re.f
dadhre.o: dadhre.f
X	f77 -c dadhre.f
dchhre.o: dchhre.f
X	f77 -c dchhre.f
dfshre.o: dfshre.f
X	f77 -c dfshre.f
dinhre.o: dinhre.f
X	f77 -c dinhre.f
drlhre.o: drlhre.f
X	f77 -c drlhre.f
dtrhre.o: dtrhre.f
X	f77 -c dtrhre.f
X
stest1: stest1.o scuhre.o s07hre.o s09hre.o s113re.o s132re.o \
X	sadhre.o schhre.o sfshre.o sinhre.o srlhre.o strhre.o
X	f77 -o stest1 stest1.o scuhre.o s07hre.o s09hre.o \
X	s113re.o s132re.o \
X	sadhre.o schhre.o sfshre.o sinhre.o srlhre.o strhre.o
X
stest2: stest2.o scuhre.o s07hre.o s09hre.o s113re.o s132re.o \
X	sadhre.o schhre.o sfshre.o sinhre.o srlhre.o strhre.o
X	f77 -o stest2 stest2.o scuhre.o s07hre.o s09hre.o \
X	s113re.o s132re.o \
X	sadhre.o schhre.o sfshre.o sinhre.o srlhre.o strhre.o
X
stest1.o: stest1.f
X	f77 -c stest1.f
stest2.o: stest2.f
X	f77 -c stest2.f
scuhre.o: scuhre.f
X	f77 -c scuhre.f
s07hre.o: s07hre.f
X	f77 -c s07hre.f
s09hre.o: s09hre.f
X	f77 -c s09hre.f
s113re.o: s113re.f
X	f77 -c s113re.f
s132re.o: s132re.f
X	f77 -c s132re.f
sadhre.o: sadhre.f
X	f77 -c sadhre.f
schhre.o: schhre.f
X	f77 -c schhre.f
sfshre.o: sfshre.f
X	f77 -c sfshre.f
sinhre.o: sinhre.f
X	f77 -c sinhre.f
srlhre.o: srlhre.f
X	f77 -c srlhre.f
strhre.o: strhre.f
X	f77 -c strhre.f
SHAR_EOF
  $shar_touch -am 0531091995 'makefile' &&
  chmod 0644 'makefile' ||
  echo 'restore of makefile failed'
  shar_count="`wc -c < 'makefile'`"
  test 1977 -eq "$shar_count" ||
    echo "makefile: original size 1010, current size $shar_count"
fi
# ============= dtest1.f ==============
if test -f 'dtest1.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dtest1.f (file already exists)'
else
  echo 'x - extracting dtest1.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dtest1.f' &&
C
C   DTEST1 is a simple test driver for DCUHRE.
C
C   Output produced on a SUN 3/50.
c
C       DCUHRE TEST RESULTS
C
C    FTEST CALLS = 3549, IFAIL =  0
C   N   ESTIMATED ERROR   INTEGRAL
C   1     0.00000010     0.13850818
C   2     0.00000013     0.06369469
C   3     0.00000874     0.05861748
C   4     0.00000021     0.05407034
C   5     0.00000019     0.05005614
C   6     0.00000009     0.04654608
C
X      PROGRAM DTEST1
X      EXTERNAL FTEST
X      INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
X      PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
X      DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
X      DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
X      DO 10 N = 1,NDIM
X         A(N) = 0
X         B(N) = 1
X   10 CONTINUE
X      MINCLS = 0
X      MAXCLS = 10000
X      KEY = 0
X      ABSREQ = 0
X      RELREQ = 1E-3
X      CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
X     * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
X      PRINT 9999, NEVAL, IFAIL
X 9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
X     * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL')
X      DO 20 N = 1,NF
X         PRINT 9998, N, ABSEST(N), FINEST(N)
X 9998    FORMAT (3X, I2, 2F15.8)
X   20 CONTINUE
X      END
X      SUBROUTINE FTEST(NDIM, Z, NFUN, F)
X      INTEGER N, NDIM, NFUN
X      DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
X      SUM = 0
X      DO 10 N = 1,NDIM
X         SUM = SUM + N*Z(N)**2
X   10 CONTINUE
X      F(1) = EXP(-SUM/2)
X      DO 20 N = 1,NDIM
X         F(N+1) = Z(N)*F(1)
X   20 CONTINUE
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dtest1.f' &&
  chmod 0644 'dtest1.f' ||
  echo 'restore of dtest1.f failed'
  shar_count="`wc -c < 'dtest1.f'`"
  test 1559 -eq "$shar_count" ||
    echo "dtest1.f: original size 1559, current size $shar_count"
fi
# ============= dtest2.f ==============
if test -f 'dtest2.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dtest2.f (file already exists)'
else
  echo 'x - extracting dtest2.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dtest2.f' &&
X      PROGRAM DTEST2
C
C   DTEST2 tests some of the features of DCUHRE.
C   DTEST2 checks that DCUHRE integrates to machine
C   precision some of the monomials that DCUHRE is
C   supposed to integrate to machine precision.
C   DTEST2 checks that the restart feature of DCUHRE works.
C   DTEST2 runs small tests in dimensions 2, 3, 5, 7 and 10.
C
C   Output produced on a SUN 3/50.
C
C
C
C    DCUHRE TEST WITH NDIM =   2, KEY =  1
C    SUBROUTINE CALLS =    195, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.00e+00   1.0000000000000000
C      2      0.22e-15   1.0000000000000002
C      3      0.22e-15   1.0000000000000002
C      4      0.22e-15   1.0000000000000002
C      5      0.22e-15   1.0000000000000002
C      6      0.00e+00   1.0000000000000000
C      7      0.22e-15   1.0000000000000002
C      8      0.00e+00   1.0000000000000000
C      9      0.00e+00   1.0000000000000000
C     10      0.22e-15   1.0000000000000002
C     11      0.00e+00   1.0000000000000000
C     12      0.22e-15   1.0000000000000002
C     13      0.22e-15   1.0000000000000002
C     14      0.22e-15   1.0000000000000002
C     15      0.22e-15   1.0000000000000002
C     16      0.22e-15   1.0000000000000002
C
C    DCUHRE TEST WITH NDIM =   3, KEY =  2
C    SUBROUTINE CALLS =    381, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.00e+00   1.0000000000000000
C      2      0.00e+00   1.0000000000000000
C      3      0.00e+00   1.0000000000000000
C      4      0.00e+00   1.0000000000000000
C      5      0.11e-15   0.9999999999999999
C      6      0.00e+00   1.0000000000000000
C      7      0.00e+00   1.0000000000000000
C      8      0.11e-15   0.9999999999999999
C      9      0.00e+00   1.0000000000000000
C     10      0.00e+00   1.0000000000000000
C     11      0.00e+00   1.0000000000000000
C     12      0.00e+00   1.0000000000000000
C     13      0.11e-15   0.9999999999999999
C     14      0.00e+00   1.0000000000000000
C     15      0.00e+00   1.0000000000000000
C     16      0.87e-09   1.0000000008661680
C
C    DCUHRE TEST WITH NDIM =   4, KEY =  3
C    SUBROUTINE CALLS =    459, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.44e-15   1.0000000000000004
C      2      0.00e+00   1.0000000000000000
C      3      0.44e-15   1.0000000000000004
C      4      0.00e+00   1.0000000000000000
C      5      0.00e+00   1.0000000000000000
C      6      0.44e-15   1.0000000000000004
C      7      0.22e-15   1.0000000000000002
C      8      0.22e-15   1.0000000000000002
C      9      0.11e-15   0.9999999999999999
C     10      0.22e-15   1.0000000000000002
C     11      0.44e-15   1.0000000000000004
C     12      0.44e-15   1.0000000000000004
C     13      0.00e+00   1.0000000000000000
C     14      0.11e-15   0.9999999999999999
C     15      0.24e-07   0.9999999758753315
C     16      0.33e-06   0.9999996739089504
C
C    DCUHRE TEST WITH NDIM =   5, KEY =  4
C    SUBROUTINE CALLS =    309, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.00e+00   1.0000000000000000
C      2      0.00e+00   1.0000000000000000
C      3      0.00e+00   1.0000000000000000
C      4      0.00e+00   1.0000000000000000
C      5      0.22e-15   0.9999999999999998
C      6      0.00e+00   1.0000000000000000
C      7      0.11e-15   0.9999999999999999
C      8      0.00e+00   1.0000000000000000
C      9      0.22e-15   1.0000000000000002
C     10      0.22e-15   1.0000000000000002
C     11      0.22e-15   1.0000000000000002
C     12      0.22e-15   1.0000000000000002
C     13      0.22e-15   1.0000000000000002
C     14      0.15e-05   1.0000015013499370
C     15      0.15e-04   1.0000147498974314
C     16      0.76e-04   1.0000756306700240
C
C    DCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =   2737, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.23e-05   1.0000023460841487
C
C    DCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =   5957, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.11e-05   1.0000010801329045
C
C    DCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =  11753, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.62e-06   1.0000006230062635
C
C    DCUHRE TEST WITH NDIM =   2, KEY =  1
C    SUBROUTINE CALLS =    455, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.19e-09   0.9999999998085931
C
C    DCUHRE TEST WITH NDIM =   3, KEY =  2
C    SUBROUTINE CALLS =   1397, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.57e-07   1.0000000565042295
C
C    DCUHRE TEST WITH NDIM =   5, KEY =  3
C    SUBROUTINE CALLS =   4641, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.13e-05   0.9999987098089600
C
C    DCUHRE TEST WITH NDIM =   7, KEY =  4
C    SUBROUTINE CALLS =   9945, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.35e-06   1.0000003517922749
C
C    DCUHRE TEST WITH NDIM =  10, KEY =  4
C    SUBROUTINE CALLS =  18975, IFAIL =  1
C      N   ABSOLUTE ERROR    INTEGRAL
C      1      0.39e-04   0.9999614789438658
C
X      EXTERNAL FTESTP,FTESTO,FTESTX
X      INTEGER N,NW
X      PARAMETER (NW = 5000)
X      DOUBLE PRECISION A(10),B(10),WRKSTR(NW),ABSERR
X      DO 10 N = 1,10
X        A(N) = 0
X        B(N) = 1
X   10 CONTINUE
X      ABSERR = 1E-10
C
C    TEST FOR INTEGRATING POLYNOMIALS
C     Selected monomials, degrees 0-13
C
C         Degree 13 rule
X      CALL ATEST(2,A,B,195,16,FTESTP,ABSERR,1,NW,0,WRKSTR)
C         Degree 11 rule
X      CALL ATEST(3,A,B,381,16,FTESTP,ABSERR,2,NW,0,WRKSTR)
C         Degree  9 rule
X      CALL ATEST(4,A,B,459,16,FTESTP,ABSERR,3,NW,0,WRKSTR)
C         Degree  7 rule
X      CALL ATEST(5,A,B,309,16,FTESTP,ABSERR,4,NW,0,WRKSTR)
C
C    TEST RESTART
C
X      CALL ATEST(6,A,B,3000,1,FTESTO,ABSERR,4,NW,0,WRKSTR)
X      CALL ATEST(6,A,B,6000,1,FTESTO,ABSERR,4,NW,1,WRKSTR)
X      CALL ATEST(6,A,B,12000,1,FTESTO,ABSERR,4,NW,1,WRKSTR)
C
C    TEST WITH NDIM = 2, 3, 5, 7, 10
C
X      CALL ATEST(2,A,B,500,1,FTESTX,ABSERR,1,NW,0,WRKSTR)
X      CALL ATEST(3,A,B,1500,1,FTESTX,ABSERR,2,NW,0,WRKSTR)
X      CALL ATEST(5,A,B,5000,1,FTESTX,ABSERR,3,NW,0,WRKSTR)
X      CALL ATEST(7,A,B,10000,1,FTESTX,ABSERR,4,NW,0,WRKSTR)
X      CALL ATEST(10,A,B,20000,1,FTESTX,ABSERR,4,NW,0,WRKSTR)
C
X      END
X      SUBROUTINE ATEST(NDIM, A, B, MAXCLS, NFUN, TSTSUB,
X     * ABSERR, KEY, LENWRK, IREST, WRKSTR)
X      EXTERNAL TSTSUB
X      INTEGER NDIM, LENWRK, KEY, IREST, NEVAL
X      DOUBLE PRECISION A(NDIM), B(NDIM), ABSEST(20), FINEST(20),
X     * WRKSTR(LENWRK), ABSERR, REL
X      SAVE NEVAL, ABSEST, FINEST
X      INTEGER N, MAXCLS, NFUN, IFAIL
X      REL = 0
X      CALL DCUHRE(NDIM, NFUN, A, B, 0, MAXCLS, TSTSUB, ABSERR, REL,
X     * KEY, LENWRK, IREST, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
X      WRITE (*,99999) NDIM, KEY
99999 FORMAT (/5X,'DCUHRE TEST WITH NDIM = ',I3,', KEY = ',I2)
X      WRITE (*,99998) NEVAL, IFAIL
99998 FORMAT (5X, 'SUBROUTINE CALLS = ', I6, ', IFAIL = ', I2)
X      WRITE (*,99997)
99997 FORMAT (7X, 'N   ABSOLUTE ERROR    INTEGRAL')
X      DO 10 N = 1,NFUN
X        WRITE (*,99996) N, ABS(FINEST(N)-1), FINEST(N)
99996   FORMAT (6X, I2, E14.2, F21.16)
X   10 CONTINUE
X      END
X      SUBROUTINE FTESTP(NDIM, Z, NFUN, F)
C
C       Selected monomials, degree 0-13
C
X      INTEGER NDIM, NFUN
X      DOUBLE PRECISION Z(NDIM), F(NFUN)
X      F(1) = 1
X      F(2) = 2*Z(1)
X      F(3) = 3*Z(1)**2
X      F(4) = F(2)*2*Z(2)
X      F(5) = 4*Z(1)**3
X      F(6) = F(3)*2*Z(2)
X      F(7) = 5*Z(1)**4
X      F(8) = F(5)*2*Z(2)
X      F(9) = F(3)*3*Z(2)**2
X      F(10) = 6*Z(1)**5
X      F(11) = F(7)*2*Z(2)
X      F(12) = F(5)*3*Z(2)**2
X      F(13) = 8*Z(1)**7
X      F(14) = 10*Z(1)**9
X      F(15) = 12*Z(1)**11
X      F(16) = 14*Z(1)**13
X      END
X      SUBROUTINE FTESTO(NDIM, Z, NFUN, F)
C
C     Corner Peak
C
X      INTEGER NDIM, NFUN
X      DOUBLE PRECISION Z(NDIM), F(NFUN)
X      F(1) = 10/(1+0.1*Z(1)+0.2*Z(2)+0.3*Z(3)+0.4*Z(4)+0.5*Z(5)+0.6*
X     * Z(6))**6/0.2057746
X      END
X      SUBROUTINE FTESTX(NDIM, Z, NFUN, F)
C
C     Sum of Cosines
C
X      INTEGER N, NDIM, NFUN
X      DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
X      SUM = 0
X      DO 10 N = 1,2
X        SUM = SUM - COS(10*Z(N))/0.0544021110889370
X   10 CONTINUE
X      F(1) = SUM/2
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dtest2.f' &&
  chmod 0644 'dtest2.f' ||
  echo 'restore of dtest2.f failed'
  shar_count="`wc -c < 'dtest2.f'`"
  test 8207 -eq "$shar_count" ||
    echo "dtest2.f: original size 8207, current size $shar_count"
fi
# ============= dcuhre.f ==============
if test -f 'dcuhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dcuhre.f (file already exists)'
else
  echo 'x - extracting dcuhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dcuhre.f' &&
X      SUBROUTINE DCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,
X     +                  EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,
X     +                  WORK)
C***BEGIN PROLOGUE DCUHRE
C***DATE WRITTEN   900116   (YYMMDD)
C***REVISION DATE  900116   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C            Alan Genz, Computer Science Department, Washington State
C            University, Pullman, WA 99163-1210, USA
C            Phone.. 509-335-2131
C            Email..  acg@cs2.cs.wsu.edu
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals
C
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
C              I   I  1  2      NDIM
C
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DCUHRE is a driver for the integration routine
C            DADHRE, which repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with greatest
C            estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            For NDIM = 2 the default integration rule is of
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C
C            When DCUHRE computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DCUHRE to all integrands in one call.
C            For integrals that vary continuously as functions of
C            some parameter, the estimates produced by DCUHRE will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in DCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for DCUHRE that may
C            be changed (with respect to the previous call of DCUHRE)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUHRE.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for DCUHRE
C              to obtain the required accuracy. In this case DCUHRE
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN is less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2)
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
C              WRKSUB) contain the centers of the subregions.
C            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
C              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths.
C            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
C              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of
C              subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)*
C              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary
C              storage in DADHRE.
C
C
C        DCUHRE Example Test Program
C
C
C   DTEST1 is a simple test driver for DCUHRE.
C
C   Output produced on a SUN 3/50.
c
C       DCUHRE TEST RESULTS
C
C    FTEST CALLS = 3549, IFAIL =  0
C   N   ESTIMATED ERROR   INTEGRAL
C   1     0.00000010     0.13850818
C   2     0.00000013     0.06369469
C   3     0.00000874     0.05861748
C   4     0.00000021     0.05407034
C   5     0.00000019     0.05005614
C   6     0.00000009     0.04654608
C
C     PROGRAM DTEST1
C     EXTERNAL FTEST
C     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
C     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
C     DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
C     DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
C     DO 10 N = 1,NDIM
C        A(N) = 0
C        B(N) = 1
C  10 CONTINUE
C     MINCLS = 0
C     MAXCLS = 10000
C     KEY = 0
C     ABSREQ = 0
C     RELREQ = 1E-3
C     CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
C    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
C     PRINT 9999, NEVAL, IFAIL
C9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
C    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL')
C     DO 20 N = 1,NF
C        PRINT 9998, N, ABSEST(N), FINEST(N)
C9998    FORMAT (3X, I2, 2F15.8)
C  20 CONTINUE
C     END
C     SUBROUTINE FTEST(NDIM, Z, NFUN, F)
C     INTEGER N, NDIM, NFUN
C     DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
C     SUM = 0
C     DO 10 N = 1,NDIM
C        SUM = SUM + N*Z(N)**2
C  10 CONTINUE
C     F(1) = EXP(-SUM/2)
C     DO 20 N = 1,NDIM
C        F(N+1) = Z(N)*F(1)
C  20 CONTINUE
C     END
C
C***LONG DESCRIPTION
C
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to DADHRE, WORK is split into eight
C   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   OLDRES and WORK.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   CENTRS contains the centers of the subregions.
C   HWIDTS contains the half widths of the subregions.
C   GREATE contains the greatest estimated error for each subregion.
C   DIR    contains the directions for further subdivision.
C   OLDRES and WORK are used as work arrays in DADHRE.
C
C   The data structures for the subregions are in DADHRE organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   DTRHRE.
C
C   The subroutine DADHRE is written for efficient execution on shared
C   memory parallel computer. On a computer with NPROC processors we wil
C   in each subdivision step divide MDIV regions, where MDIV is
C   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new region
C   Each processor will then compute estimates of the integrals and erro
C   over 2*MDIV/NPROC subregions in each subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, DRLHRE, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of subregions, compiler
C   directives should be placed in front of the DO 200
C   loop in DADHRE on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   To be published.
C
C   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, To be published.
C
C***ROUTINES CALLED DCHHRE,DADHRE
C***END PROLOGUE DCUHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR
X      INTEGER NEVAL,IFAIL
X      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
X      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MDIV   Integer.
C          MDIV is the number of subregions that are divided in
C          each subdivision step in DADHRE.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MDIV,MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS, KEY and NDIM.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN, NDIM and MDIV. This determines the length
C          of the main work arrays.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C
X      INTEGER MDIV,MAXWT,WTLENG,MAXDIM,LENW2,MAXSUB,MINSUB
X      INTEGER NUM,NSUB,LENW,KEYF
X      PARAMETER (MDIV=1)
X      PARAMETER (MAXDIM=15)
X      PARAMETER (MAXWT=14)
X      PARAMETER (LENW2=2*MDIV*MAXDIM* (MAXWT+1)+12*MAXWT+2*MAXDIM)
X      INTEGER WRKSUB,I1,I2,I3,I4,I5,I6,I7,I8,K1,K2,K3,K4,K5,K6,K7,K8
X      DOUBLE PRECISION WORK2(LENW2)
C
C***FIRST EXECUTABLE STATEMENT DCUHRE
C
C   Compute NUM, WTLENG, MAXSUB and MINSUB,
C   and check the input parameters.
C
X      CALL DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,EPSABS,
X     +            EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,KEYF,
X     +            IFAIL,WTLENG)
X      WRKSUB = (NW - 1 - 17*MDIV*NUMFUN)/(2*NDIM + 2*NUMFUN + 2)
X      IF (IFAIL.NE.0) THEN
X          GO TO 999
X      END IF
C
C   Split up the work space.
C
X      I1 = 1
X      I2 = I1 + WRKSUB*NUMFUN
X      I3 = I2 + WRKSUB*NUMFUN
X      I4 = I3 + WRKSUB*NDIM
X      I5 = I4 + WRKSUB*NDIM
X      I6 = I5 + WRKSUB
X      I7 = I6 + WRKSUB
X      I8 = I7 + NUMFUN*MDIV
X      K1 = 1
X      K2 = K1 + 2*MDIV*WTLENG*NDIM
X      K3 = K2 + WTLENG*5
X      K4 = K3 + WTLENG
X      K5 = K4 + NDIM
X      K6 = K5 + NDIM
X      K7 = K6 + 2*MDIV*NDIM
X      K8 = K7 + 3*WTLENG
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
X      IF (RESTAR.EQ.1) THEN
X          NSUB = WORK(NW)
X      END IF
C
C   Compute the size of the temporary work space needed in DADHRE.
C
X      LENW = 16*MDIV*NUMFUN
X      CALL DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,EPSABS,
X     +            EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,
X     +            RESULT,ABSERR,NEVAL,NSUB,IFAIL,WORK(I1),WORK(I2),
X     +            WORK(I3),WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),
X     +            WORK2(K1),WORK2(K2),WORK2(K3),WORK2(K4),WORK2(K5),
X     +            WORK2(K6),WORK2(K7),WORK2(K8))
X      WORK(NW) = NSUB
999   RETURN
C
C***END DCUHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dcuhre.f' &&
  chmod 0644 'dcuhre.f' ||
  echo 'restore of dcuhre.f failed'
  shar_count="`wc -c < 'dcuhre.f'`"
  test 17921 -eq "$shar_count" ||
    echo "dcuhre.f: original size 17921, current size $shar_count"
fi
# ============= dchhre.f ==============
if test -f 'dchhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dchhre.f (file already exists)'
else
  echo 'x - extracting dchhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dchhre.f' &&
X      SUBROUTINE DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,
X     +                  EPSABS,EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,
X     +                  KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE DCHHRE
C***PURPOSE  DCHHRE checks the validity of the
C            input parameters to DCUHRE.
C***DESCRIPTION
C            DCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and
C            IFAIL as functions of the input parameters to DCUHRE,
C            and checks the validity of the input parameters to DCUHRE.
C
C   ON ENTRY
C
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of subregions that are divided in
C            each subdivision step in DADHRE.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in DCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C
C   ON RETURN
C
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS, KEY and NDIM.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS, KEY and NDIM.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than
C                      MAXDIM.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     KEYF   Integer.
C            Key to selected integration rule.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DCHHRE
C
C   Global variables.
C
X      INTEGER NDIM,NUMFUN,MDIV,MINPTS,MAXPTS,KEY,NW,MINSUB,MAXSUB
X      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG
X      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
C
C   Local variables.
C
X      INTEGER LIMIT,J
C
C***FIRST EXECUTABLE STATEMENT DCHHRE
C
X      IFAIL = 0
C
C   Check on legal KEY.
C
X      IF (KEY.LT.0 .OR. KEY.GT.4) THEN
X          IFAIL = 2
X          GO TO 999
X      END IF
C
C   Check on legal NDIM.
C
X      IF (NDIM.LT.2 .OR. NDIM.GT.MAXDIM) THEN
X          IFAIL = 3
X          GO TO 999
X      END IF
C
C   For KEY = 1, NDIM must be equal to 2.
C
X      IF (KEY.EQ.1 .AND. NDIM.NE.2) THEN
X          IFAIL = 4
X          GO TO 999
X      END IF
C
C   For KEY = 2, NDIM must be equal to 3.
C
X      IF (KEY.EQ.2 .AND. NDIM.NE.3) THEN
X          IFAIL = 5
X          GO TO 999
X      END IF
C
C   For KEY = 0, we point at the selected integration rule.
C
X      IF (KEY.EQ.0) THEN
X          IF (NDIM.EQ.2) THEN
X              KEYF = 1
X          ELSE IF (NDIM.EQ.3) THEN
X              KEYF = 2
X          ELSE
X              KEYF = 3
X          ENDIF
X      ELSE
X          KEYF = KEY
X      ENDIF
C
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C
X      IF (KEYF.EQ.1) THEN
X          NUM = 65
X          WTLENG = 14
X      ELSE IF (KEYF.EQ.2) THEN
X          NUM = 127
X          WTLENG = 13
X      ELSE IF (KEYF.EQ.3) THEN
X          NUM = 1 + 4*2*NDIM + 2*NDIM* (NDIM-1) + 4*NDIM* (NDIM-1) +
X     +          4*NDIM* (NDIM-1)* (NDIM-2)/3 + 2**NDIM
X          WTLENG = 9
X          IF (NDIM.EQ.2) WTLENG = 8
X      ELSE IF (KEYF.EQ.4) THEN
X          NUM = 1 + 3*2*NDIM + 2*NDIM* (NDIM-1) + 2**NDIM
X          WTLENG = 6
X      END IF
C
C   Compute MAXSUB.
C
X      MAXSUB = (MAXPTS-NUM)/ (2*NUM) + 1
C
C   Compute MINSUB.
C
X      MINSUB = (MINPTS-NUM)/ (2*NUM) + 1
X      IF (MOD(MINPTS-NUM,2*NUM).NE.0) THEN
X          MINSUB = MINSUB + 1
X      END IF
X      MINSUB = MAX(2,MINSUB)
C
C   Check on positive NUMFUN.
C
X      IF (NUMFUN.LT.1) THEN
X          IFAIL = 6
X          GO TO 999
X      END IF
C
C   Check on legal upper and lower limits of integration.
C
X      DO 10 J = 1,NDIM
X          IF (A(J)-B(J).EQ.0) THEN
X              IFAIL = 7
X              GO TO 999
X          END IF
10    CONTINUE
C
C   Check on MAXPTS < 3*NUM.
C
X      IF (MAXPTS.LT.3*NUM) THEN
X          IFAIL = 8
X          GO TO 999
X      END IF
C
C   Check on MAXPTS >= MINPTS.
C
X      IF (MAXPTS.LT.MINPTS) THEN
X          IFAIL = 9
X          GO TO 999
X      END IF
C
C   Check on legal accuracy requests.
C
X      IF (EPSABS.LT.0 .AND. EPSREL.LT.0) THEN
X          IFAIL = 10
X          GO TO 999
X      END IF
C
C   Check on big enough double precision workspace.
C
X      LIMIT = MAXSUB* (2*NDIM+2*NUMFUN+2) + 17*MDIV*NUMFUN + 1
X      IF (NW.LT.LIMIT) THEN
X          IFAIL = 11
X          GO TO 999
X      END IF
C
C    Check on legal RESTAR.
C
X      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
X          IFAIL = 12
X          GO TO 999
X      END IF
999   RETURN
C
C***END DCHHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dchhre.f' &&
  chmod 0644 'dchhre.f' ||
  echo 'restore of dchhre.f failed'
  shar_count="`wc -c < 'dchhre.f'`"
  test 8328 -eq "$shar_count" ||
    echo "dchhre.f: original size 8328, current size $shar_count"
fi
# ============= dadhre.f ==============
if test -f 'dadhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dadhre.f (file already exists)'
else
  echo 'x - extracting dadhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dadhre.f' &&
X      SUBROUTINE DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,
X     +                  EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,
X     +                  RESULT,ABSERR,NEVAL,NSUB,IFAIL,VALUES,
X     +                  ERRORS,CENTRS,HWIDTS,GREATE,DIR,OLDRES,WORK,G,W,
X     +                  RULPTS,CENTER,HWIDTH,X,SCALES,NORMS)
C***BEGIN PROLOGUE DADHRE
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals, I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DADHRE repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with  greatest
C            estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The regions are divided in two equally sized parts along
C            the direction with greatest absolute fourth divided
C            difference.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            Defines the number of new subregions that are divided
C            in each subdivision step.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINSUB Integer.
C            The computations proceed until there are at least
C            MINSUB subregions in the data structure.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from DADHRE
C            must not be changed since the last
C            exit from DADHRE.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*MDIV*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to DADHRE.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by DADHRE.
C     NSUB   Integer.
C            Number of stored subregions.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for DADHRE
C              to obtain the required accuracy. In this case DADHRE
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB.
C            DIR is used to store the directions for
C            further subdivision.
C     OLDRES Real array of dimension (NUMFUN,MDIV).
C            Used to store old estimates of the integrals over the
C            subregions.
C     WORK   Real array of dimension LENW.
C            Used  in DRLHRE and DTRHRE.
C     G      Real array of dimension (NDIM,WTLENG,2*MDIV).
C            The fully symmetric sum generators for the rules.
C            G(1,J,1),...,G(NDIM,J,1) are the generators for the
C            points associated with the Jth weights.
C            When MDIV subregions are divided in 2*MDIV
C            subregions, the subregions may be processed on different
C            processors and we must make a copy of the generators
C            for each processor.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in DINHRE.
C     CENTER Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     X      Real array of dimension (NDIM,2*MDIV).
C            Work array used in DRLHRE.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C
C***REFERENCES
C
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
C
C***ROUTINES CALLED DTRHRE,DINHRE,DRLHRE
C***END PROLOGUE DADHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN,MDIV,MINSUB,MAXSUB,KEY,LENW,RESTAR
X      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG
X      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
X      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN)
X      DOUBLE PRECISION VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB)
X      DOUBLE PRECISION CENTRS(NDIM,MAXSUB)
X      DOUBLE PRECISION HWIDTS(NDIM,MAXSUB)
X      DOUBLE PRECISION GREATE(MAXSUB),DIR(MAXSUB)
X      DOUBLE PRECISION OLDRES(NUMFUN,MDIV)
X      DOUBLE PRECISION WORK(LENW),RULPTS(WTLENG)
X      DOUBLE PRECISION G(NDIM,WTLENG,2*MDIV),W(5,WTLENG)
X      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM,2*MDIV)
X      DOUBLE PRECISION SCALES(3,WTLENG),NORMS(3,WTLENG)
C
C   Local variables.
C
C   INTSGN is used to get correct sign on the integral.
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coeff. defined in DINHRE and used by DRLHRE
C          and DADHRE.
C
X      INTEGER I,J,K
X      INTEGER INTSGN,SBRGNS
X      INTEGER L1
X      INTEGER NDIV,POINTR,DIRECT,INDEX
X      DOUBLE PRECISION OLDCEN,EST1,EST2,ERRCOF(6)
C
C***FIRST EXECUTABLE STATEMENT DADHRE
C
C   Get the correct sign on the integral.
C
X      INTSGN = 1
X      DO 10 J = 1,NDIM
X          IF (B(J).LT.A(J)) THEN
X              INTSGN = - INTSGN
X          END IF
10    CONTINUE
C
C   Call DINHRE to compute the weights and abscissas of
C   the function evaluation points.
C
X      CALL DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C
C   If RESTAR = 1, then this is a restart run.
C
X      IF (RESTAR.EQ.1) THEN
X          SBRGNS = NSUB
X          GO TO 110
X      END IF
C
C   Initialize the SBRGNS, CENTRS and HWIDTS.
C
X      SBRGNS = 1
X      DO 15 J = 1,NDIM
X          CENTRS(J,1) = (A(J)+B(J))/2
X          HWIDTS(J,1) = ABS(B(J)-A(J))/2
15    CONTINUE
C
C   Initialize RESULT, ABSERR and NEVAL.
C
X      DO 20 J = 1,NUMFUN
X          RESULT(J) = 0
X          ABSERR(J) = 0
20    CONTINUE
X      NEVAL = 0
C
C   Apply DRLHRE over the whole region.
C
X      CALL DRLHRE(NDIM,CENTRS(1,1),HWIDTS(1,1),WTLENG,G,W,ERRCOF,NUMFUN,
X     +            FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,1),ERRORS(1,1),
X     +            DIR(1))
X      NEVAL = NEVAL + NUM
C
C   Add the computed values to RESULT and ABSERR.
C
X      DO 55 J = 1,NUMFUN
X          RESULT(J) = RESULT(J) + VALUES(J,1)
55    CONTINUE
X      DO 65 J = 1,NUMFUN
X          ABSERR(J) = ABSERR(J) + ERRORS(J,1)
65    CONTINUE
C
C   Store results in heap.
C
X      INDEX = 1
X      CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,HWIDTS,
X     +            GREATE,WORK(1),WORK(NUMFUN+1),CENTER,HWIDTH,DIR)
C
C***End initialisation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+1 is less than MAXSUB.
C
110   IF (SBRGNS+1.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over each half of the
C   NDIV subregions with greatest errors.
C   If MAXSUB is great enough, NDIV = MDIV
C
X          IF (MDIV.GT.1) THEN
X              NDIV = MAXSUB - SBRGNS
X              NDIV = MIN(NDIV,MDIV,SBRGNS)
X          ELSE
X              NDIV = 1
X          END IF
C
C   Divide the NDIV subregions in two halves, and compute
C   integral and error over each half.
C
X          DO 150 I = 1,NDIV
X              POINTR = SBRGNS + NDIV + 1 - I
C
C   Adjust RESULT and ABSERR.
C
X              DO 115 J = 1,NUMFUN
X                  RESULT(J) = RESULT(J) - VALUES(J,1)
X                  ABSERR(J) = ABSERR(J) - ERRORS(J,1)
115           CONTINUE
C
C   Compute first half region.
C
X              DO 120 J = 1,NDIM
X                  CENTRS(J,POINTR) = CENTRS(J,1)
X                  HWIDTS(J,POINTR) = HWIDTS(J,1)
120           CONTINUE
X              DIRECT = DIR(1)
X              DIR(POINTR) = DIRECT
X              HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,1)/2
X              OLDCEN = CENTRS(DIRECT,1)
X              CENTRS(DIRECT,POINTR) = OLDCEN - HWIDTS(DIRECT,POINTR)
C
C   Save the computed values of the integrals.
C
X              DO 125 J = 1,NUMFUN
X                  OLDRES(J,NDIV-I+1) = VALUES(J,1)
125           CONTINUE
C
C   Adjust the heap.
C
X              CALL DTRHRE(1,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
X     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
X     +                    HWIDTH,DIR)
C
C   Compute second half region.
C
X              DO 130 J = 1,NDIM
X                  CENTRS(J,POINTR-1) = CENTRS(J,POINTR)
X                  HWIDTS(J,POINTR-1) = HWIDTS(J,POINTR)
130           CONTINUE
X              CENTRS(DIRECT,POINTR-1) = OLDCEN + HWIDTS(DIRECT,POINTR)
X              HWIDTS(DIRECT,POINTR-1) = HWIDTS(DIRECT,POINTR)
X              DIR(POINTR-1) = DIRECT
150       CONTINUE
C
C   Make copies of the generators for each processor.
C
X          DO 190 I = 2,2*NDIV
X              DO 190 J = 1,NDIM
X                  DO 190 K = 1,WTLENG
X                      G(J,K,I) = G(J,K,1)
190       CONTINUE
C
C   Apply basic rule.
C
Cvd$l cncall
X          DO 200 I = 1,2*NDIV
X              INDEX = SBRGNS + I
X              L1 = 1 + (I-1)*8*NUMFUN
X              CALL DRLHRE(NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,
X     +                    G(1,1,I),W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,
X     +                    X(1,I),WORK(L1),VALUES(1,INDEX),
X     +                    ERRORS(1,INDEX),DIR(INDEX))
200       CONTINUE
X          NEVAL = NEVAL + 2*NDIV*NUM
C
C   Add new contributions to RESULT.
C
X          DO 220 I = 1,2*NDIV
X              DO 210 J = 1,NUMFUN
X                  RESULT(J) = RESULT(J) + VALUES(J,SBRGNS+I)
210           CONTINUE
220       CONTINUE
C
C   Check consistency of results and if necessary adjust
C   the estimated errors.
C
X          DO 240 I = 1,NDIV
X              GREATE(SBRGNS+2*I-1) = 0
X              GREATE(SBRGNS+2*I) = 0
X              DO 230 J = 1,NUMFUN
X                  EST1 = ABS(OLDRES(J,I)- (VALUES(J,
X     +                   SBRGNS+2*I-1)+VALUES(J,SBRGNS+2*I)))
X                  EST2 = ERRORS(J,SBRGNS+2*I-1) + ERRORS(J,SBRGNS+2*I)
X                  IF (EST2.GT.0) THEN
X                      ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)*
X     +                  (1+ERRCOF(5)*EST1/EST2)
X                      ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)*
X     +                                       (1+ERRCOF(5)*EST1/EST2)
X                  END IF
X                  ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1) +
X     +                                     ERRCOF(6)*EST1
X                  ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I) +
X     +                                   ERRCOF(6)*EST1
X                  IF (ERRORS(J,SBRGNS+2*I-1).GT.
X     +                GREATE(SBRGNS+2*I-1)) THEN
X                      GREATE(SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)
X                  END IF
X                  IF (ERRORS(J,SBRGNS+2*I).GT.GREATE(SBRGNS+2*I)) THEN
X                      GREATE(SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)
X                  END IF
X                  ABSERR(J) = ABSERR(J) + ERRORS(J,SBRGNS+2*I-1) +
X     +                        ERRORS(J,SBRGNS+2*I)
230           CONTINUE
240       CONTINUE
C
C   Store results in heap.
C
X          DO 250 I = 1,2*NDIV
X              INDEX = SBRGNS + I
X              CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,
X     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
X     +                    HWIDTH,DIR)
250       CONTINUE
X          SBRGNS = SBRGNS + 2*NDIV
C
C   Check for termination.
C
X          IF (SBRGNS.LT.MINSUB) THEN
X              GO TO 110
X          END IF
X          DO 255 J = 1,NUMFUN
X              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
X     +            ABSERR(J).GT.EPSABS) THEN
X                  GO TO 110
X              END IF
255       CONTINUE
X          IFAIL = 0
X          GO TO 499
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
X      ELSE
X          IFAIL = 1
X      END IF
C
C   Compute more accurate values of RESULT and ABSERR.
C
499   CONTINUE
X      DO 500 J = 1,NUMFUN
X          RESULT(J) = 0
X          ABSERR(J) = 0
500   CONTINUE
X      DO 510 I = 1,SBRGNS
X          DO 505 J = 1,NUMFUN
X              RESULT(J) = RESULT(J) + VALUES(J,I)
X              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
505       CONTINUE
510   CONTINUE
C
C   Compute correct sign on the integral.
C
X      DO 600 J = 1,NUMFUN
X          RESULT(J) = RESULT(J)*INTSGN
600   CONTINUE
X      NSUB = SBRGNS
X      RETURN
C
C***END DADHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dadhre.f' &&
  chmod 0644 'dadhre.f' ||
  echo 'restore of dadhre.f failed'
  shar_count="`wc -c < 'dadhre.f'`"
  test 16631 -eq "$shar_count" ||
    echo "dadhre.f: original size 16631, current size $shar_count"
fi
# ============= dtrhre.f ==============
if test -f 'dtrhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dtrhre.f (file already exists)'
else
  echo 'x - extracting dtrhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dtrhre.f' &&
X      SUBROUTINE DTRHRE(DVFLAG,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
X     +                  HWIDTS,GREATE,ERROR,VALUE,CENTER,HWIDTH,DIR)
C***BEGIN PROLOGUE DTRHRE
C***PURPOSE DTRHRE maintains a heap of subregions.
C***DESCRIPTION DTRHRE maintains a heap of subregions.
C            The subregions are ordered according to the size
C            of the greatest error estimates of each subregion(GREATE).
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     NDIM   Integer.
C            Number of variables.
C     NUMFUN Integer.
C            Number of components of the integral.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     VALUES Real array of dimension (NUMFUN,SBRGNS).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,SBRGNS).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,SBRGNS).
C            Used to store the center limits of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,SBRGNS).
C            Used to store the hwidth limits of the stored subregions.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     ERROR  Real array of dimension NUMFUN.
C            Used as intermediate storage for the error of a subregion.
C     VALUE  Real array of dimension NUMFUN.
C            Used as intermediate storage for the estimate
C            of the integral over a subregion.
C     CENTER Real array of dimension NDIM.
C            Used as intermediate storage for the center of
C            the subregion.
C     HWIDTH Real array of dimension NDIM.
C            Used as intermediate storage for the half width of
C            the subregion.
C     DIR    Integer array of dimension SBRGNS.
C            DIR is used to store the directions for
C            further subdivision.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DTRHRE
C
C   Global variables.
C
X      INTEGER DVFLAG,NDIM,NUMFUN,SBRGNS
X      DOUBLE PRECISION VALUES(NUMFUN,*),ERRORS(NUMFUN,*)
X      DOUBLE PRECISION CENTRS(NDIM,*)
X      DOUBLE PRECISION HWIDTS(NDIM,*)
X      DOUBLE PRECISION GREATE(*)
X      DOUBLE PRECISION ERROR(NUMFUN),VALUE(NUMFUN)
X      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM)
X      DOUBLE PRECISION DIR(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   DIRECT is used as intermediate storage for the direction of further
C          subdivision.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
C
X      INTEGER J,SUBRGN,SUBTMP
X      DOUBLE PRECISION GREAT,DIRECT
C
C***FIRST EXECUTABLE STATEMENT DTRHRE
C
C   Save values to be stored in their correct place in the heap.
C
X      GREAT = GREATE(SBRGNS)
X      DIRECT = DIR(SBRGNS)
X      DO 5 J = 1,NUMFUN
X          ERROR(J) = ERRORS(J,SBRGNS)
X          VALUE(J) = VALUES(J,SBRGNS)
5     CONTINUE
X      DO 10 J = 1,NDIM
X          CENTER(J) = CENTRS(J,SBRGNS)
X          HWIDTH(J) = HWIDTS(J,SBRGNS)
10    CONTINUE
C
C    If DVFLAG = 1, we will remove the region
C    with greatest estimated error from the heap.
C
X      IF (DVFLAG.EQ.1) THEN
X          SBRGNS = SBRGNS - 1
X          SUBRGN = 1
20        SUBTMP = 2*SUBRGN
X          IF (SUBTMP.LE.SBRGNS) THEN
X              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
X                  IF (GREATE(SUBTMP).LT.GREATE(SUBTMP+1)) THEN
X                      SUBTMP = SUBTMP + 1
X                  END IF
X              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
X              IF (GREAT.LT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp up the heap.
C
X                  GREATE(SUBRGN) = GREATE(SUBTMP)
X                  DO 25 J = 1,NUMFUN
X                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
X                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
25                CONTINUE
X                  DIR(SUBRGN) = DIR(SUBTMP)
X                  DO 30 J = 1,NDIM
X                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
X                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
30                CONTINUE
X                  SUBRGN = SUBTMP
X                  GO TO 20
X              END IF
X          END IF
X      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, then insert new region in the heap.
C
X          SUBRGN = SBRGNS
40        SUBTMP = SUBRGN/2
X          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
X              IF (GREAT.GT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp down the heap.
C
X                  GREATE(SUBRGN) = GREATE(SUBTMP)
X                  DO 45 J = 1,NUMFUN
X                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
X                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
45                CONTINUE
X                  DIR(SUBRGN) = DIR(SUBTMP)
X                  DO 50 J = 1,NDIM
X                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
X                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
50                CONTINUE
X                  SUBRGN = SUBTMP
X                  GO TO 40
X              END IF
X          END IF
X      END IF
C
C    Insert the saved values in their correct places.
C
X      IF (SBRGNS.GT.0) THEN
X          GREATE(SUBRGN) = GREAT
X          DO 55 J = 1,NUMFUN
X              ERRORS(J,SUBRGN) = ERROR(J)
X              VALUES(J,SUBRGN) = VALUE(J)
55        CONTINUE
X          DIR(SUBRGN) = DIRECT
X          DO 60 J = 1,NDIM
X              CENTRS(J,SUBRGN) = CENTER(J)
X              HWIDTS(J,SUBRGN) = HWIDTH(J)
60        CONTINUE
X      END IF
C
C***END DTRHRE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dtrhre.f' &&
  chmod 0644 'dtrhre.f' ||
  echo 'restore of dtrhre.f failed'
  shar_count="`wc -c < 'dtrhre.f'`"
  test 5933 -eq "$shar_count" ||
    echo "dtrhre.f: original size 5933, current size $shar_count"
fi
# ============= dinhre.f ==============
if test -f 'dinhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dinhre.f (file already exists)'
else
  echo 'x - extracting dinhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dinhre.f' &&
X      SUBROUTINE DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C***BEGIN PROLOGUE DINHRE
C***PURPOSE DINHRE computes abscissas and weights of the integration
C            rule and the null rules to be used in error estimation.
C            These are computed as functions of NDIM and KEY.
C***DESCRIPTION DINHRE will for given values of NDIM and KEY compute or
C            select the correct values of the abscissas and
C            corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W.
C            The heuristic error coefficients ERRCOF
C            will be computed as a function of KEY.
C            Scaling factors SCALES and normalization factors NORMS
C            used in the error estimation are computed.
C
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables.
C     KEY    Integer.
C            Key to selected local integration rule.
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C            It is assumed that the error is computed using:
C             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C               THEN ERROR = ERRCOF(3)*N1
C               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C            where N1-N3 are the null rules, EP is the error for
C            the parent
C            subregion and ES is the error for the sibling subregion.
C     RULPTS Real array of dimension WTLENG.
C            A work array containing the number of points produced by
C            each generator of the selected rule.
C     SCALES Real array of dimension (3,WTLENG).
C            Scaling factors used to construct new null rules,
C            N1, N2 and N3,
C            based on a linear combination of two successive null rules
C            in the sequence of null rules.
C     NORMS  Real array of dimension (3,WTLENG).
C            2**NDIM/(1-norm of the null rule constructed by each of
C            the scaling factors.)
C
C***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE
C***END PROLOGUE DINHRE
C
C   Global variables.
C
X      INTEGER NDIM,KEY,WTLENG
X      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
X      DOUBLE PRECISION RULPTS(WTLENG),SCALES(3,WTLENG)
X      DOUBLE PRECISION NORMS(3,WTLENG)
C
C   Local variables.
C
X      INTEGER I,J,K
X      DOUBLE PRECISION WE(14)
C
C***FIRST EXECUTABLE STATEMENT DINHRE
C
C   Compute W, G and ERRCOF.
C
X      IF (KEY.EQ.1) THEN
X          CALL D132RE(WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.2) THEN
X          CALL D113RE(WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.3) THEN
X          CALL D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.4) THEN
X          CALL D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
X      END IF
C
C   Compute SCALES and NORMS.
C
X      DO 100 K = 1,3
X          DO 50 I = 1,WTLENG
X              IF (W(K+1,I).NE.0) THEN
X                  SCALES(K,I) = - W(K+2,I)/W(K+1,I)
X              ELSE
X                  SCALES(K,I) = 100
X              END IF
X              DO 30 J = 1,WTLENG
X                  WE(J) = W(K+2,J) + SCALES(K,I)*W(K+1,J)
30            CONTINUE
X              NORMS(K,I) = 0
X              DO 40 J = 1,WTLENG
X                  NORMS(K,I) = NORMS(K,I) + RULPTS(J)*ABS(WE(J))
40            CONTINUE
X              NORMS(K,I) = 2**NDIM/NORMS(K,I)
50        CONTINUE
100   CONTINUE
X      RETURN
C
C***END DINHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dinhre.f' &&
  chmod 0644 'dinhre.f' ||
  echo 'restore of dinhre.f failed'
  shar_count="`wc -c < 'dinhre.f'`"
  test 4084 -eq "$shar_count" ||
    echo "dinhre.f: original size 4084, current size $shar_count"
fi
# ============= d132re.f ==============
if test -f 'd132re.f' && test X"$1" != X"-c"; then
  echo 'x - skipping d132re.f (file already exists)'
else
  echo 'x - extracting d132re.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'd132re.f' &&
X      SUBROUTINE D132RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D132RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D132RE computes abscissas and weights of a 2 dimensional
C            integration rule of degree 13.
C            Two null rules of degree 11, one null rule of degree 9
C            and one null rule of degree 7 to be used in error
C            estimation are also computed.
C ***DESCRIPTION D132RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W. The heuristic error coefficients ERRCOF
C            will also be assigned.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points produced by each generator.
C***REFERENCES S.Eriksen,
C              Thesis of the degree cand.scient, Dept. of Informatics,
C              Univ. of Bergen,Norway, 1984.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE D132RE
C
C   Global variables
C
X      INTEGER WTLENG
X      DOUBLE PRECISION W(5,WTLENG),G(2,WTLENG),ERRCOF(6)
X      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
X      INTEGER I,J
X      DOUBLE PRECISION DIM2G(16)
X      DOUBLE PRECISION DIM2W(14,5)
C
X      DATA (DIM2G(I),I=1,16)/0.2517129343453109D+00,
X     +     0.7013933644534266D+00,0.9590960631619962D+00,
X     +     0.9956010478552127D+00,0.5000000000000000D+00,
X     +     0.1594544658297559D+00,0.3808991135940188D+00,
X     +     0.6582769255267192D+00,0.8761473165029315D+00,
X     +     0.9982431840531980D+00,0.9790222658168462D+00,
X     +     0.6492284325645389D+00,0.8727421201131239D+00,
X     +     0.3582614645881228D+00,0.5666666666666666D+00,
X     +     0.2077777777777778D+00/
C
X      DATA (DIM2W(I,1),I=1,14)/0.3379692360134460D-01,
X     +     0.9508589607597761D-01,0.1176006468056962D+00,
X     +     0.2657774586326950D-01,0.1701441770200640D-01,
X     +     0.0000000000000000D+00,0.1626593098637410D-01,
X     +     0.1344892658526199D+00,0.1328032165460149D+00,
X     +     0.5637474769991870D-01,0.3908279081310500D-02,
X     +     0.3012798777432150D-01,0.1030873234689166D+00,
X     +     0.6250000000000000D-01/
C
X      DATA (DIM2W(I,2),I=1,14)/0.3213775489050763D+00,
X     +     - .1767341636743844D+00,0.7347600537466072D-01,
X     +     - .3638022004364754D-01,0.2125297922098712D-01,
X     +     0.1460984204026913D+00,0.1747613286152099D-01,
X     +     0.1444954045641582D+00,0.1307687976001325D-03,
X     +     0.5380992313941161D-03,0.1042259576889814D-03,
X     +     - .1401152865045733D-02,0.8041788181514763D-02,
X     +     - .1420416552759383D+00/
C
X      DATA (DIM2W(I,3),I=1,14)/0.3372900883288987D+00,
X     +     - .1644903060344491D+00,0.7707849911634622D-01,
X     +     - .3804478358506310D-01,0.2223559940380806D-01,
X     +     0.1480693879765931D+00,0.4467143702185814D-05,
X     +     0.1508944767074130D+00,0.3647200107516215D-04,
X     +     0.5777198999013880D-03,0.1041757313688177D-03,
X     +     - .1452822267047819D-02,0.8338339968783705D-02,
X     +     - .1472796329231960D+00/
C
X      DATA (DIM2W(I,4),I=1,14)/ - .8264123822525677D+00,
X     +     0.3065838614094360D+00,0.2389292538329435D-02,
X     +     - .1343024157997222D+00,0.8833366840533900D-01,
X     +     0.0000000000000000D+00,0.9786283074168292D-03,
X     +     - .1319227889147519D+00,0.7990012200150630D-02,
X     +     0.3391747079760626D-02,0.2294915718283264D-02,
X     +     - .1358584986119197D-01,0.4025866859057809D-01,
X     +     0.3760268580063992D-02/
C
X      DATA (DIM2W(I,5),I=1,14)/0.6539094339575232D+00,
X     +     - .2041614154424632D+00, - .1746981515794990D+00,
X     +     0.3937939671417803D-01,0.6974520545933992D-02,
X     +     0.0000000000000000D+00,0.6667702171778258D-02,
X     +     0.5512960621544304D-01,0.5443846381278607D-01,
X     +     0.2310903863953934D-01,0.1506937747477189D-01,
X     +     - .6057021648901890D-01,0.4225737654686337D-01,
X     +     0.2561989142123099D-01/
C
C***FIRST EXECUTABLE STATEMENT D132RE
C
C   Assign values to W.
C
X      DO 10 I = 1,14
X          DO 10 J = 1,5
X              W(J,I) = DIM2W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
X      DO 20 I = 1,2
X          DO 20 J = 1,14
X              G(I,J) = 0
20    CONTINUE
X      G(1,2) = DIM2G(1)
X      G(1,3) = DIM2G(2)
X      G(1,4) = DIM2G(3)
X      G(1,5) = DIM2G(4)
X      G(1,6) = DIM2G(5)
X      G(1,7) = DIM2G(6)
X      G(2,7) = G(1,7)
X      G(1,8) = DIM2G(7)
X      G(2,8) = G(1,8)
X      G(1,9) = DIM2G(8)
X      G(2,9) = G(1,9)
X      G(1,10) = DIM2G(9)
X      G(2,10) = G(1,10)
X      G(1,11) = DIM2G(10)
X      G(2,11) = G(1,11)
X      G(1,12) = DIM2G(11)
X      G(2,12) = DIM2G(12)
X      G(1,13) = DIM2G(13)
X      G(2,13) = DIM2G(14)
X      G(1,14) = DIM2G(15)
X      G(2,14) = DIM2G(16)
C
C   Assign values to RULPTS.
C
X      RULPTS(1) = 1
X      DO 30 I = 2,11
X          RULPTS(I) = 4
30    CONTINUE
X      RULPTS(12) = 8
X      RULPTS(13) = 8
X      RULPTS(14) = 8
C
C   Assign values to ERRCOF.
C
X      ERRCOF(1) = 10
X      ERRCOF(2) = 10
X      ERRCOF(3) = 1.
X      ERRCOF(4) = 5.
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END D132RE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'd132re.f' &&
  chmod 0644 'd132re.f' ||
  echo 'restore of d132re.f failed'
  shar_count="`wc -c < 'd132re.f'`"
  test 5952 -eq "$shar_count" ||
    echo "d132re.f: original size 5952, current size $shar_count"
fi
# ============= d113re.f ==============
if test -f 'd113re.f' && test X"$1" != X"-c"; then
  echo 'x - skipping d113re.f (file already exists)'
else
  echo 'x - extracting d113re.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'd113re.f' &&
X      SUBROUTINE D113RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D113RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D113RE computes abscissas and weights of a 3 dimensional
C            integration rule of degree 11.
C            Two null rules of degree 9, one null rule of degree 7
C            and one null rule of degree 5 to be used in error
C            estimation are also computed.
C***DESCRIPTION D113RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to G
C            and W.
C            The heuristic error coefficients ERRCOF
C            will also be computed.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points used by each generator.
C
C***REFERENCES  J.Berntsen, Cautious adaptive numerical integration
C               over the 3-cube, Reports in Informatics 17, Dept. of
C               Inf.,Univ. of Bergen, Norway, 1985.
C               J.Berntsen and T.O.Espelid, On the construction of
C               higher degree three-dimensional embedded integration
C               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234,
C               1988.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D113RE
C
C   Global variables.
C
X      INTEGER WTLENG
X      DOUBLE PRECISION W(5,WTLENG),G(3,WTLENG),ERRCOF(6)
X      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
X      INTEGER I,J
X      DOUBLE PRECISION DIM3G(14)
X      DOUBLE PRECISION DIM3W(13,5)
C
X      DATA (DIM3G(I),I=1,14)/0.1900000000000000D+00,
X     +     0.5000000000000000D+00,0.7500000000000000D+00,
X     +     0.8000000000000000D+00,0.9949999999999999D+00,
X     +     0.9987344998351400D+00,0.7793703685672423D+00,
X     +     0.9999698993088767D+00,0.7902637224771788D+00,
X     +     0.4403396687650737D+00,0.4378478459006862D+00,
X     +     0.9549373822794593D+00,0.9661093133630748D+00,
X     +     0.4577105877763134D+00/
C
X      DATA (DIM3W(I,1),I=1,13)/0.7923078151105734D-02,
X     +     0.6797177392788080D-01,0.1086986538805825D-02,
X     +     0.1838633662212829D+00,0.3362119777829031D-01,
X     +     0.1013751123334062D-01,0.1687648683985235D-02,
X     +     0.1346468564512807D+00,0.1750145884600386D-02,
X     +     0.7752336383837454D-01,0.2461864902770251D+00,
X     +     0.6797944868483039D-01,0.1419962823300713D-01/
C
X      DATA (DIM3W(I,2),I=1,13)/0.1715006248224684D+01,
X     +     - .3755893815889209D+00,0.1488632145140549D+00,
X     +     - .2497046640620823D+00,0.1792501419135204D+00,
X     +     0.3446126758973890D-02, - .5140483185555825D-02,
X     +     0.6536017839876425D-02, - .6513454939229700D-03,
X     +     - .6304672433547204D-02,0.1266959399788263D-01,
X     +     - .5454241018647931D-02,0.4826995274768427D-02/
C
X      DATA (DIM3W(I,3),I=1,13)/0.1936014978949526D+01,
X     +     - .3673449403754268D+00,0.2929778657898176D-01,
X     +     - .1151883520260315D+00,0.5086658220872218D-01,
X     +     0.4453911087786469D-01, - .2287828257125900D-01,
X     +     0.2908926216345833D-01, - .2898884350669207D-02,
X     +     - .2805963413307495D-01,0.5638741361145884D-01,
X     +     - .2427469611942451D-01,0.2148307034182882D-01/
C
X      DATA (DIM3W(I,4),I=1,13)/0.5170828195605760D+00,
X     +     0.1445269144914044D-01, - .3601489663995932D+00,
X     +     0.3628307003418485D+00,0.7148802650872729D-02,
X     +     - .9222852896022966D-01,0.1719339732471725D-01,
X     +     - .1021416537460350D+00, - .7504397861080493D-02,
X     +     0.1648362537726711D-01,0.5234610158469334D-01,
X     +     0.1445432331613066D-01,0.3019236275367777D-02/
C
X      DATA (DIM3W(I,5),I=1,13)/0.2054404503818520D+01,
X     +     0.1377759988490120D-01, - .5768062917904410D+00,
X     +     0.3726835047700328D-01,0.6814878939777219D-02,
X     +     0.5723169733851849D-01, - .4493018743811285D-01,
X     +     0.2729236573866348D-01,0.3547473950556990D-03,
X     +     0.1571366799739551D-01,0.4990099219278567D-01,
X     +     0.1377915552666770D-01,0.2878206423099872D-02/
C
C***FIRST EXECUTABLE STATEMENT D113RE
C
C   Assign values to W.
C
X      DO 10 I = 1,13
X          DO 10 J = 1,5
X              W(J,I) = DIM3W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
X      DO 20 I = 1,3
X          DO 20 J = 1,13
X              G(I,J) = 0
20    CONTINUE
X      G(1,2) = DIM3G(1)
X      G(1,3) = DIM3G(2)
X      G(1,4) = DIM3G(3)
X      G(1,5) = DIM3G(4)
X      G(1,6) = DIM3G(5)
X      G(1,7) = DIM3G(6)
X      G(2,7) = G(1,7)
X      G(1,8) = DIM3G(7)
X      G(2,8) = G(1,8)
X      G(1,9) = DIM3G(8)
X      G(2,9) = G(1,9)
X      G(3,9) = G(1,9)
X      G(1,10) = DIM3G(9)
X      G(2,10) = G(1,10)
X      G(3,10) = G(1,10)
X      G(1,11) = DIM3G(10)
X      G(2,11) = G(1,11)
X      G(3,11) = G(1,11)
X      G(1,12) = DIM3G(12)
X      G(2,12) = DIM3G(11)
X      G(3,12) = G(2,12)
X      G(1,13) = DIM3G(13)
X      G(2,13) = G(1,13)
X      G(3,13) = DIM3G(14)
C
C   Assign values to RULPTS.
C
X      RULPTS(1) = 1
X      RULPTS(2) = 6
X      RULPTS(3) = 6
X      RULPTS(4) = 6
X      RULPTS(5) = 6
X      RULPTS(6) = 6
X      RULPTS(7) = 12
X      RULPTS(8) = 12
X      RULPTS(9) = 8
X      RULPTS(10) = 8
X      RULPTS(11) = 8
X      RULPTS(12) = 24
X      RULPTS(13) = 24
C
C   Assign values to ERRCOF.
C
X      ERRCOF(1) = 4
X      ERRCOF(2) = 4.
X      ERRCOF(3) = 0.5
X      ERRCOF(4) = 3.
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END D113RE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'd113re.f' &&
  chmod 0644 'd113re.f' ||
  echo 'restore of d113re.f failed'
  shar_count="`wc -c < 'd113re.f'`"
  test 6203 -eq "$shar_count" ||
    echo "d113re.f: original size 6203, current size $shar_count"
fi
# ============= d09hre.f ==============
if test -f 'd09hre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping d09hre.f (file already exists)'
else
  echo 'x - extracting d09hre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'd09hre.f' &&
X      SUBROUTINE D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D09HRE
C***KEYWORDS basic integration rule, degree 9
C***PURPOSE  To initialize a degree 9 basic rule and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-20
C***DESCRIPTION  D09HRE initializes a degree 9 integration rule,
C            two degree 7 null rules, one degree 5 null rule and one
C            degree 3 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D09HRE
C
C   Global variables
C
X      INTEGER NDIM,WTLENG
X      DOUBLE PRECISION W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
X      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local Variables
C
X      DOUBLE PRECISION RATIO,LAM0,LAM1,LAM2,LAM3,LAMP,TWONDM
X      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT D09HRE
C
C
C     Initialize generators, weights and RULPTS
C
X      DO 30 J = 1,WTLENG
X          DO 10 I = 1,NDIM
X              G(I,J) = 0
10        CONTINUE
X          DO 20 I = 1,5
X              W(I,J) = 0
20        CONTINUE
X          RULPTS(J) = 2*NDIM
30    CONTINUE
X      TWONDM = 2**NDIM
X      RULPTS(WTLENG) = TWONDM
X      IF (NDIM.GT.2) RULPTS(8) = (4*NDIM* (NDIM-1)* (NDIM-2))/3
X      RULPTS(7) = 4*NDIM* (NDIM-1)
X      RULPTS(6) = 2*NDIM* (NDIM-1)
X      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
X      LAM0 = 0.4707
X      LAM1 = 4/ (15-5/LAM0)
X      RATIO = (1-LAM1/LAM0)/27
X      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
X      RATIO = RATIO* (1-LAM2/LAM0)/3
X      LAM3 = (7-9* (LAM2+LAM1)+63*LAM2*LAM1/5-63*RATIO)/
X     +       (9-63* (LAM2+LAM1)/5+21*LAM2*LAM1-63*RATIO/LAM0)
X      LAMP = 0.0625
C
C     Compute degree 9 rule weights
C
X      W(1,WTLENG) = 1/ (3*LAM0)**4/TWONDM
X      IF (NDIM.GT.2) W(1,8) = (1-1/ (3*LAM0))/ (6*LAM1)**3
X      W(1,7) = (1-7* (LAM0+LAM1)/5+7*LAM0*LAM1/3)/
X     +         (84*LAM1*LAM2* (LAM2-LAM0)* (LAM2-LAM1))
X      W(1,6) = (1-7* (LAM0+LAM2)/5+7*LAM0*LAM2/3)/
X     +         (84*LAM1*LAM1* (LAM1-LAM0)* (LAM1-LAM2)) -
X     +         W(1,7)*LAM2/LAM1 - 2* (NDIM-2)*W(1,8)
X      W(1,4) = (1-9* ((LAM0+LAM1+LAM2)/7- (LAM0*LAM1+LAM0*LAM2+
X     +         LAM1*LAM2)/5)-3*LAM0*LAM1*LAM2)/
X     +         (18*LAM3* (LAM3-LAM0)* (LAM3-LAM1)* (LAM3-LAM2))
X      W(1,3) = (1-9* ((LAM0+LAM1+LAM3)/7- (LAM0*LAM1+LAM0*LAM3+
X     +         LAM1*LAM3)/5)-3*LAM0*LAM1*LAM3)/
X     +         (18*LAM2* (LAM2-LAM0)* (LAM2-LAM1)* (LAM2-LAM3)) -
X     +         2* (NDIM-1)*W(1,7)
X      W(1,2) = (1-9* ((LAM0+LAM2+LAM3)/7- (LAM0*LAM2+LAM0*LAM3+
X     +         LAM2*LAM3)/5)-3*LAM0*LAM2*LAM3)/
X     +         (18*LAM1* (LAM1-LAM0)* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(1,7)+W(1,6)+ (NDIM-2)*W(1,8))
C
C     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
C
X      W(2,WTLENG) = 1/ (108*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(2,8) = (1-27*TWONDM*W(2,9)*LAM0**3)/ (6*LAM1)**3
X      W(2,7) = (1-5*LAM1/3-15*TWONDM*W(2,WTLENG)*LAM0**2* (LAM0-LAM1))/
X     +          (60*LAM1*LAM2* (LAM2-LAM1))
X      W(2,6) = (1-9* (8*LAM1*LAM2*W(2,7)+TWONDM*W(2,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(2,8)* (NDIM-2)
X      W(2,4) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
X      W(2,3) = (1-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(2,7)
X      W(2,2) = (1-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(2,7)+W(2,6)+ (NDIM-2)*W(2,8))
X      W(3,WTLENG) = 5/ (324*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(3,8) = (1-27*TWONDM*W(3,9)*LAM0**3)/ (6*LAM1)**3
X      W(3,7) = (1-5*LAM1/3-15*TWONDM*W(3,WTLENG)*LAM0**2* (LAM0-LAM1))/
X     +          (60*LAM1*LAM2* (LAM2-LAM1))
X      W(3,6) = (1-9* (8*LAM1*LAM2*W(3,7)+TWONDM*W(3,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(3,8)* (NDIM-2)
X      W(3,5) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAMP* (LAMP-LAM1)* (LAMP-LAM2))
X      W(3,3) = (1-7* ((LAM1+LAMP)/5-LAM1*LAMP/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAMP)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAMP)) - 2* (NDIM-1)*W(3,7)
X      W(3,2) = (1-7* ((LAM2+LAMP)/5-LAM2*LAMP/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAMP)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAMP)) -
X     +         2* (NDIM-1)* (W(3,7)+W(3,6)+ (NDIM-2)*W(3,8))
X      W(4,WTLENG) = 2/ (81*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(4,8) = (2-27*TWONDM*W(4,9)*LAM0**3)/ (6*LAM1)**3
X      W(4,7) = (2-15*LAM1/9-15*TWONDM*W(4,WTLENG)*LAM0* (LAM0-LAM1))/
X     +         (60*LAM1*LAM2* (LAM2-LAM1))
X      W(4,6) = (1-9* (8*LAM1*LAM2*W(4,7)+TWONDM*W(4,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(4,8)* (NDIM-2)
X      W(4,4) = (2-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
X      W(4,3) = (2-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(4,7)
X      W(4,2) = (2-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(4,7)+W(4,6)+ (NDIM-2)*W(4,8))
X      W(5,2) = 1/ (6*LAM1)
C
C     Set generator values
C
X      LAM0 = SQRT(LAM0)
X      LAM1 = SQRT(LAM1)
X      LAM2 = SQRT(LAM2)
X      LAM3 = SQRT(LAM3)
X      LAMP = SQRT(LAMP)
X      DO 40 I = 1,NDIM
X          G(I,WTLENG) = LAM0
40    CONTINUE
X      IF (NDIM.GT.2) THEN
X          G(1,8) = LAM1
X          G(2,8) = LAM1
X          G(3,8) = LAM1
X      END IF
X      G(1,7) = LAM1
X      G(2,7) = LAM2
X      G(1,6) = LAM1
X      G(2,6) = LAM1
X      G(1,5) = LAMP
X      G(1,4) = LAM3
X      G(1,3) = LAM2
X      G(1,2) = LAM1
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 9 rule weights and lower degree rule weights.
C
X      W(1,1) = TWONDM
X      DO 70 J = 2,5
X          DO 50 I = 2,WTLENG
X              W(J,I) = W(J,I) - W(1,I)
X              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
X      DO 80 I = 2,WTLENG
X          W(1,I) = TWONDM*W(1,I)
X          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
X      ERRCOF(1) = 5
X      ERRCOF(2) = 5
X      ERRCOF(3) = 1.
X      ERRCOF(4) = 5
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END D09HRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'd09hre.f' &&
  chmod 0644 'd09hre.f' ||
  echo 'restore of d09hre.f failed'
  shar_count="`wc -c < 'd09hre.f'`"
  test 7863 -eq "$shar_count" ||
    echo "d09hre.f: original size 7863, current size $shar_count"
fi
# ============= d07hre.f ==============
if test -f 'd07hre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping d07hre.f (file already exists)'
else
  echo 'x - extracting d07hre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'd07hre.f' &&
X      SUBROUTINE D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D07HRE
C***KEYWORDS basic integration rule, degree 7
C***PURPOSE  To initialize a degree 7 basic rule, and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-31
C***DESCRIPTION  D07HRE initializes a degree 7 integration rule,
C            two degree 5 null rules, one degree 3 null rule and one
C            degree 1 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C          WTLENG MUST be set equal to 6.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D07HRE
C
C   Global variables
C
X      INTEGER NDIM,WTLENG
X      DOUBLE PRECISION W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
X      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local Variables
C
X      DOUBLE PRECISION RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
X      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT D07HRE
C
C
C     Initialize generators, weights and RULPTS
C
X      DO 30 J = 1,WTLENG
X          DO 10 I = 1,NDIM
X              G(I,J) = 0
10        CONTINUE
X          DO 20 I = 1,5
X              W(I,J) = 0
20        CONTINUE
X          RULPTS(J) = 2*NDIM
30    CONTINUE
X      TWONDM = 2**NDIM
X      RULPTS(WTLENG) = TWONDM
X      RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
X      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
X      LAM0 = 0.4707
X      LAMP = 0.5625
X      LAM1 = 4/ (15-5/LAM0)
X      RATIO = (1-LAM1/LAM0)/27
X      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
C
C     Compute degree 7 rule weights
C
X      W(1,6) = 1/ (3*LAM0)**3/TWONDM
X      W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
X      W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/
X     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
X      W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAM2* (LAM2-LAM1))
C
C     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
C
X      W(2,6) = 1/ (36*LAM0**3)/TWONDM
X      W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
X      W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/
X     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
X      W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAM2* (LAM2-LAM1))
X      W(3,6) = 5/ (108*LAM0**3)/TWONDM
X      W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
X      W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/
X     +         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
X      W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAMP* (LAMP-LAM1))
X      W(4,6) = 1/ (54*LAM0**3)/TWONDM
X      W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
X      W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/
X     +         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
X      W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/
X     +         (20*LAM2* (LAM2-LAM1))
C
C     Set generator values
C
X      LAM0 = SQRT(LAM0)
X      LAM1 = SQRT(LAM1)
X      LAM2 = SQRT(LAM2)
X      LAMP = SQRT(LAMP)
X      DO 40 I = 1,NDIM
X          G(I,WTLENG) = LAM0
40    CONTINUE
X      G(1,WTLENG-1) = LAM1
X      G(2,WTLENG-1) = LAM1
X      G(1,WTLENG-4) = LAM2
X      G(1,WTLENG-3) = LAM1
X      G(1,WTLENG-2) = LAMP
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 7 rule weights and lower degree rule weights.
C
X      W(1,1) = TWONDM
X      DO 70 J = 2,5
X          DO 50 I = 2,WTLENG
X              W(J,I) = W(J,I) - W(1,I)
X              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
X      DO 80 I = 2,WTLENG
X          W(1,I) = TWONDM*W(1,I)
X          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
X      ERRCOF(1) = 5
X      ERRCOF(2) = 5
X      ERRCOF(3) = 1
X      ERRCOF(4) = 5
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END D07HRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'd07hre.f' &&
  chmod 0644 'd07hre.f' ||
  echo 'restore of d07hre.f failed'
  shar_count="`wc -c < 'd07hre.f'`"
  test 4921 -eq "$shar_count" ||
    echo "d07hre.f: original size 4921, current size $shar_count"
fi
# ============= drlhre.f ==============
if test -f 'drlhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping drlhre.f (file already exists)'
else
  echo 'x - extracting drlhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'drlhre.f' &&
X      SUBROUTINE DRLHRE(NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,
X     +                  FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT)
C***BEGIN PROLOGUE DRLHRE
C***KEYWORDS basic numerical integration rule
C***PURPOSE  To compute basic integration rule values.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 90-02-06
C***DESCRIPTION DRLHRE computes basic integration rule values for a
C            vector of integrands over a hyper-rectangular region.
C            These are estimates for the integrals. DRLHRE also computes
C            estimates for the errors and determines the coordinate axis
C            where the fourth difference for the integrands is largest.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   WTLENG Integer.
C          The number of weights in the basic integration rule.
C   G      Real array of dimension (NDIM,WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1,J), ..., G(NDIM,J) are the are the generators for the
C          points associated with the Jth weights.
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   ERRCOF Real array of dimension 6.
C          The error coefficients for the rules.
C          It is assumed that the error is computed using:
C           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C             THEN ERROR = ERRCOF(3)*N1
C             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C          where N1-N4 are the null rules, EP is the error
C          for the parent
C          subregion and ES is the error for the sibling subregion.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   SCALES Real array of dimension (3,WTLENG).
C          Scaling factors used to construct new null rules based
C          on a linear combination of two successive null rules
C          in the sequence of null rules.
C   NORMS  Real array of dimension (3,WTLENG).
C          2**NDIM/(1-norm of the null rule constructed by each of the
C          scaling factors.)
C   X      Real Array of dimension NDIM.
C          A work array.
C   NULL   Real array of dimension (NUMFUN, 8)
C          A work array.
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   DIRECT Real.
C          The coordinate axis where the fourth difference of the
C          integrand values is largest.
C
C***REFERENCES
C   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
C   integration over an N-dimensional rectangular region,
C   J.Comp.Appl.Math., 6:295-302, 1980.
C
C   T.O.Espelid, Integration Rules, Null Rules and Error
C   Estimation, Reports in Informatics 33, Dept. of Informatics,
C   Univ. of Bergen, 1988.
C
C***ROUTINES CALLED: DFSHRE, FUNSUB
C
C***END PROLOGUE DRLHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER WTLENG,NUMFUN,NDIM
X      DOUBLE PRECISION CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),
X     +                 RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),
X     +                 G(NDIM,WTLENG),ERRCOF(6),DIRECT,SCALES(3,WTLENG),
X     +                 NORMS(3,WTLENG)
C
C   Local variables.
C
X      DOUBLE PRECISION RGNVOL,DIFSUM,DIFMAX,FRTHDF
X      INTEGER I,J,K,DIVAXN
X      DOUBLE PRECISION SEARCH,RATIO
C
C***FIRST EXECUTABLE STATEMENT DRLHRE
C
C
C       Compute volume of subregion, initialize DIVAXN and rule sums;
C       compute fourth differences and new DIVAXN (RGNERR is used
C       for a work array here). The integrand values used for the
C       fourth divided differences are accumulated in rule arrays.
C
X      RGNVOL = 1
X      DIVAXN = 1
X      DO 10 I = 1,NDIM
X          RGNVOL = RGNVOL*HWIDTH(I)
X          X(I) = CENTER(I)
X          IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN = I
10    CONTINUE
X      CALL FUNSUB(NDIM,X,NUMFUN,RGNERR)
X      DO 30 J = 1,NUMFUN
X          BASVAL(J) = W(1,1)*RGNERR(J)
X          DO 20 K = 1,4
X              NULL(J,K) = W(K+1,1)*RGNERR(J)
20        CONTINUE
30    CONTINUE
X      DIFMAX = 0
X      RATIO = (G(1,3)/G(1,2))**2
X      DO 60 I = 1,NDIM
X          X(I) = CENTER(I) - HWIDTH(I)*G(1,2)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,5))
X          X(I) = CENTER(I) + HWIDTH(I)*G(1,2)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,6))
X          X(I) = CENTER(I) - HWIDTH(I)*G(1,3)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,7))
X          X(I) = CENTER(I) + HWIDTH(I)*G(1,3)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,8))
X          X(I) = CENTER(I)
X          DIFSUM = 0
X          DO 50 J = 1,NUMFUN
X              FRTHDF = 2* (1-RATIO)*RGNERR(J) - (NULL(J,7)+NULL(J,8)) +
X     +                 RATIO* (NULL(J,5)+NULL(J,6))
C
C           Ignore differences below roundoff
C
X              IF (RGNERR(J)+FRTHDF/4.NE.RGNERR(J)) DIFSUM = DIFSUM +
X     +            ABS(FRTHDF)
X              DO 40 K = 1,4
X                  NULL(J,K) = NULL(J,K) + W(K+1,2)*
X     +                        (NULL(J,5)+NULL(J,6)) +
X     +                        W(K+1,3)* (NULL(J,7)+NULL(J,8))
40            CONTINUE
X              BASVAL(J) = BASVAL(J) + W(1,2)* (NULL(J,5)+NULL(J,6)) +
X     +                    W(1,3)* (NULL(J,7)+NULL(J,8))
50        CONTINUE
X          IF (DIFSUM.GT.DIFMAX) THEN
X              DIFMAX = DIFSUM
X              DIVAXN = I
X          END IF
60    CONTINUE
X      DIRECT = DIVAXN
C
C    Finish computing the rule values.
C
X      DO 90 I = 4,WTLENG
X          CALL DFSHRE(NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,
X     +                NULL(1,5))
X          DO 80 J = 1,NUMFUN
X              BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
X              DO 70 K = 1,4
X                  NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
70            CONTINUE
80        CONTINUE
90    CONTINUE
C
C    Compute errors.
C
X      DO 130 J = 1,NUMFUN
C
C    We search for the null rule, in the linear space spanned by two
C    successive null rules in our sequence, which gives the greatest
C    error estimate among all normalized (1-norm) null rules in this
C    space.
C
X          DO 110 I = 1,3
X              SEARCH = 0
X              DO 100 K = 1,WTLENG
X                  SEARCH = MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,
X     +                     K)*NULL(J,I))*NORMS(I,K))
100           CONTINUE
X              NULL(J,I) = SEARCH
110       CONTINUE
X          IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2) .AND.
X     +        ERRCOF(2)*NULL(J,2).LE.NULL(J,3)) THEN
X              RGNERR(J) = ERRCOF(3)*NULL(J,1)
X          ELSE
X              RGNERR(J) = ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
X          END IF
X          RGNERR(J) = RGNVOL*RGNERR(J)
X          BASVAL(J) = RGNVOL*BASVAL(J)
130   CONTINUE
C
C***END DRLHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'drlhre.f' &&
  chmod 0644 'drlhre.f' ||
  echo 'restore of drlhre.f failed'
  shar_count="`wc -c < 'drlhre.f'`"
  test 7938 -eq "$shar_count" ||
    echo "drlhre.f: original size 7938, current size $shar_count"
fi
# ============= dfshre.f ==============
if test -f 'dfshre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping dfshre.f (file already exists)'
else
  echo 'x - extracting dfshre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'dfshre.f' &&
X      SUBROUTINE DFSHRE(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
X     +                  FUNVLS)
C***BEGIN PROLOGUE DFSHRE
C***KEYWORDS fully symmetric sum
C***PURPOSE  To compute fully symmetric basic rule sums
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DFSHRE computes a fully symmetric sum for a vector
C            of integrand values over a hyper-rectangular region.
C            The sum is fully symmetric with respect to the center of
C            the region and is taken over all sign changes and
C            permutations of the generators for the sum.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   X      Real Array of dimension NDIM.
C          A work array.
C   G      Real Array of dimension NDIM.
C          The generators for the fully symmetric sum. These MUST BE
C          non-negative and non-increasing.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   ON RETURN
C
C   FULSMS Real array of dimension NUMFUN.
C          The values for the fully symmetric sums for each component
C          of the integrand.
C   FUNVLS Real array of dimension NUMFUN.
C          A work array.
C
C***ROUTINES CALLED: FUNSUB
C
C***END PROLOGUE DFSHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN
X      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
X     +                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
C
C   Local variables.
C
X      INTEGER IXCHNG,LXCHNG,I,J,L
X      DOUBLE PRECISION GL,GI
C
C***FIRST EXECUTABLE STATEMENT DFSHRE
C
X      DO 10 J = 1,NUMFUN
X          FULSMS(J) = 0
10    CONTINUE
C
C     Compute centrally symmetric sum for permutation of G
C
20    DO 30 I = 1,NDIM
X          X(I) = CENTER(I) + G(I)*HWIDTH(I)
30    CONTINUE
40    CALL FUNSUB(NDIM,X,NUMFUN,FUNVLS)
X      DO 50 J = 1,NUMFUN
X          FULSMS(J) = FULSMS(J) + FUNVLS(J)
50    CONTINUE
X      DO 60 I = 1,NDIM
X          G(I) = - G(I)
X          X(I) = CENTER(I) + G(I)*HWIDTH(I)
X          IF (G(I).LT.0) GO TO 40
60    CONTINUE
C
C       Find next distinct permutation of G and loop back for next sum.
C       Permutations are generated in reverse lexicographic order.
C
X      DO 80 I = 2,NDIM
X          IF (G(I-1).GT.G(I)) THEN
X              GI = G(I)
X              IXCHNG = I - 1
X              DO 70 L = 1, (I-1)/2
X                  GL = G(L)
X                  G(L) = G(I-L)
X                  G(I-L) = GL
X                  IF (GL.LE.GI) IXCHNG = IXCHNG - 1
X                  IF (G(L).GT.GI) LXCHNG = L
70            CONTINUE
X              IF (G(IXCHNG).LE.GI) IXCHNG = LXCHNG
X              G(I) = G(IXCHNG)
X              G(IXCHNG) = GI
X              GO TO 20
X          END IF
80    CONTINUE
C
C     Restore original order to generators
C
X      DO 90 I = 1,NDIM/2
X          GI = G(I)
X          G(I) = G(NDIM-I+1)
X          G(NDIM-I+1) = GI
90    CONTINUE
C
C***END DFSHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'dfshre.f' &&
  chmod 0644 'dfshre.f' ||
  echo 'restore of dfshre.f failed'
  shar_count="`wc -c < 'dfshre.f'`"
  test 3846 -eq "$shar_count" ||
    echo "dfshre.f: original size 3846, current size $shar_count"
fi
# ============= stest1.f ==============
if test -f 'stest1.f' && test X"$1" != X"-c"; then
  echo 'x - skipping stest1.f (file already exists)'
else
  echo 'x - extracting stest1.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'stest1.f' &&
C
C   STEST1 is a simple test driver for SCUHRE.
C
C   Output produced on a SUN 3/50.
c
C       SCUHRE TEST RESULTS
C
C    FTEST CALLS = 3549, IFAIL =  0
C   N   ESTIMATED ERROR    INTEGRAL
C   1     0.00000013     0.13850819
C   2     0.00000015     0.06369469
C   3     0.00000875     0.05861748
C   4     0.00000020     0.05407035
C   5     0.00000020     0.05005614
C   6     0.00000009     0.04654606
C
X      PROGRAM STEST1
X      EXTERNAL FTEST
X      INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
X      PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
X      REAL A(NDIM), B(NDIM), WRKSTR(NW)
X      REAL ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
X      DO 10 N = 1,NDIM
X         A(N) = 0
X         B(N) = 1
X   10 CONTINUE
X      MINCLS = 0
X      MAXCLS = 10000
X      KEY = 0
X      ABSREQ = 0
X      RELREQ = 1E-3
X      CALL SCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
X     * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
X      PRINT 9999, NEVAL, IFAIL
X 9999 FORMAT (8X, 'SCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
X     * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL')
X      DO 20 N = 1,NF
X         PRINT 9998, N, ABSEST(N), FINEST(N)
X 9998    FORMAT (3X, I2, 2F15.8)
X   20 CONTINUE
X      END
X      SUBROUTINE FTEST(NDIM, Z, NFUN, F)
X      INTEGER N, NDIM, NFUN
X      REAL Z(NDIM), F(NFUN), SUM
X      SUM = 0
X      DO 10 N = 1,NDIM
X         SUM = SUM + N*Z(N)**2
X   10 CONTINUE
X      F(1) = EXP(-SUM/2)
X      DO 20 N = 1,NDIM
X         F(N+1) = Z(N)*F(1)
X   20 CONTINUE
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'stest1.f' &&
  chmod 0644 'stest1.f' ||
  echo 'restore of stest1.f failed'
  shar_count="`wc -c < 'stest1.f'`"
  test 1524 -eq "$shar_count" ||
    echo "stest1.f: original size 1524, current size $shar_count"
fi
# ============= stest2.f ==============
if test -f 'stest2.f' && test X"$1" != X"-c"; then
  echo 'x - skipping stest2.f (file already exists)'
else
  echo 'x - extracting stest2.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'stest2.f' &&
X      PROGRAM STEST2
C
C   STEST2 tests some of the features of SCUHRE.
C   STEST2 checks that SCUHRE integrates to machine
C   precision some of the monomials that SCUHRE is
C   supposed to integrate to machine precision.
C   STEST2 checks that the restart feature of SCUHRE works.
C   STEST2 runs small tests in dimensions 2, 3, 5, 7 and 10.
C
C   Output produced on a SUN 3/50.
C
C
C
C    SCUHRE TEST WITH NDIM =   2, KEY =  1
C    SUBROUTINE CALLS =    195, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.00e+00    1.00000000
C      2      0.00e+00    1.00000000
C      3      0.00e+00    1.00000000
C      4      0.00e+00    1.00000000
C      5      0.00e+00    1.00000000
C      6      0.00e+00    1.00000000
C      7      0.00e+00    1.00000000
C      8      0.00e+00    1.00000000
C      9      0.00e+00    1.00000000
C     10      0.00e+00    1.00000000
C     11      0.60e-07    0.99999994
C     12      0.00e+00    1.00000000
C     13      0.00e+00    1.00000000
C     14      0.12e-06    0.99999988
C     15      0.12e-06    1.00000012
C     16      0.12e-06    1.00000012
C
C    SCUHRE TEST WITH NDIM =   3, KEY =  2
C    SUBROUTINE CALLS =    381, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.00e+00    1.00000000
C      2      0.00e+00    1.00000000
C      3      0.00e+00    1.00000000
C      4      0.00e+00    1.00000000
C      5      0.00e+00    1.00000000
C      6      0.60e-07    0.99999994
C      7      0.00e+00    1.00000000
C      8      0.00e+00    1.00000000
C      9      0.60e-07    0.99999994
C     10      0.00e+00    1.00000000
C     11      0.00e+00    1.00000000
C     12      0.00e+00    1.00000000
C     13      0.00e+00    1.00000000
C     14      0.00e+00    1.00000000
C     15      0.60e-07    0.99999994
C     16      0.12e-06    0.99999988
C
C    SCUHRE TEST WITH NDIM =   4, KEY =  3
C    SUBROUTINE CALLS =    459, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.12e-06    0.99999988
C      2      0.00e+00    1.00000000
C      3      0.12e-06    0.99999988
C      4      0.00e+00    1.00000000
C      5      0.30e-06    0.99999970
C      6      0.60e-07    0.99999994
C      7      0.24e-06    0.99999976
C      8      0.24e-06    0.99999976
C      9      0.36e-06    0.99999964
C     10      0.00e+00    1.00000000
C     11      0.12e-06    0.99999988
C     12      0.12e-06    0.99999988
C     13      0.48e-06    0.99999952
C     14      0.36e-06    0.99999964
C     15      0.11e-05    0.99999893
C     16      0.54e-06    0.99999946
C
C    SCUHRE TEST WITH NDIM =   5, KEY =  4
C    SUBROUTINE CALLS =    309, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.00e+00    1.00000000
C      2      0.12e-06    1.00000012
C      3      0.12e-06    1.00000012
C      4      0.12e-06    1.00000012
C      5      0.60e-07    0.99999994
C      6      0.12e-06    1.00000012
C      7      0.00e+00    1.00000000
C      8      0.00e+00    1.00000000
C      9      0.12e-06    1.00000012
C     10      0.00e+00    1.00000000
C     11      0.12e-06    1.00000012
C     12      0.12e-06    1.00000012
C     13      0.12e-06    1.00000012
C     14      0.15e-05    1.00000155
C     15      0.15e-04    1.00001466
C     16      0.76e-04    1.00007582
C
C    SCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =   2737, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.23e-05    1.00000227
C
C    SCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =   5957, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.95e-06    1.00000095
C
C    SCUHRE TEST WITH NDIM =   6, KEY =  4
C    SUBROUTINE CALLS =  11753, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.60e-06    1.00000060
C
C    SCUHRE TEST WITH NDIM =   2, KEY =  1
C    SUBROUTINE CALLS =    455, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.83e-06    1.00000083
C
C    SCUHRE TEST WITH NDIM =   3, KEY =  2
C    SUBROUTINE CALLS =   1397, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.48e-06    1.00000048
C
C    SCUHRE TEST WITH NDIM =   5, KEY =  3
C    SUBROUTINE CALLS =   4641, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.18e-05    0.99999821
C
C    SCUHRE TEST WITH NDIM =   7, KEY =  4
C    SUBROUTINE CALLS =   9945, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.48e-06    1.00000048
C
C    SCUHRE TEST WITH NDIM =  10, KEY =  4
C    SUBROUTINE CALLS =  18975, IFAIL =  1
C      N   ABSOLUTE ERROR  INTEGRAL
C      1      0.32e-04    0.99996787
C
X      EXTERNAL FTESTP,FTESTO,FTESTX
X      INTEGER N,NW
X      PARAMETER (NW = 5000)
X      REAL A(10),B(10),WRKSTR(NW),ABSERR
X      DO 10 N = 1,10
X        A(N) = 0
X        B(N) = 1
X   10 CONTINUE
X      ABSERR = 1E-10
C
C    TEST FOR INTEGRATING POLYNOMIALS
C     Selected monomials, degrees 0-13
C
C         Degree 13 rule
X      CALL ATEST(2,A,B,195,16,FTESTP,ABSERR,1,NW,0,WRKSTR)
C         Degree 11 rule
X      CALL ATEST(3,A,B,381,16,FTESTP,ABSERR,2,NW,0,WRKSTR)
C         Degree  9 rule
X      CALL ATEST(4,A,B,459,16,FTESTP,ABSERR,3,NW,0,WRKSTR)
C         Degree  7 rule
X      CALL ATEST(5,A,B,309,16,FTESTP,ABSERR,4,NW,0,WRKSTR)
C
C    TEST RESTART
C
X      CALL ATEST(6,A,B,3000,1,FTESTO,ABSERR,4,NW,0,WRKSTR)
X      CALL ATEST(6,A,B,6000,1,FTESTO,ABSERR,4,NW,1,WRKSTR)
X      CALL ATEST(6,A,B,12000,1,FTESTO,ABSERR,4,NW,1,WRKSTR)
C
C    TEST WITH NDIM = 2, 3, 5, 7, 10
C
X      CALL ATEST(2,A,B,500,1,FTESTX,ABSERR,1,NW,0,WRKSTR)
X      CALL ATEST(3,A,B,1500,1,FTESTX,ABSERR,2,NW,0,WRKSTR)
X      CALL ATEST(5,A,B,5000,1,FTESTX,ABSERR,3,NW,0,WRKSTR)
X      CALL ATEST(7,A,B,10000,1,FTESTX,ABSERR,4,NW,0,WRKSTR)
X      CALL ATEST(10,A,B,20000,1,FTESTX,ABSERR,4,NW,0,WRKSTR)
C
X      END
X      SUBROUTINE ATEST(NDIM, A, B, MAXCLS, NFUN, TSTSUB,
X     * ABSERR, KEY, LENWRK, IREST, WRKSTR)
X      EXTERNAL TSTSUB
X      INTEGER NDIM, LENWRK, KEY, IREST, NEVAL
X      REAL A(NDIM), B(NDIM), ABSEST(20), FINEST(20),
X     * WRKSTR(LENWRK), ABSERR, REL
X      SAVE NEVAL, ABSEST, FINEST
X      INTEGER N, MAXCLS, NFUN, IFAIL
X      REL = 0
X      CALL SCUHRE(NDIM, NFUN, A, B, 0, MAXCLS, TSTSUB, ABSERR, REL,
X     * KEY, LENWRK, IREST, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
X      WRITE (*,99999) NDIM, KEY
99999 FORMAT (/5X,'SCUHRE TEST WITH NDIM = ',I3,', KEY = ',I2)
X      WRITE (*,99998) NEVAL, IFAIL
99998 FORMAT (5X, 'SUBROUTINE CALLS = ', I6, ', IFAIL = ', I2)
X      WRITE (*,99997)
99997 FORMAT (7X, 'N   ABSOLUTE ERROR  INTEGRAL')
X      DO 10 N = 1,NFUN
X        WRITE (*,99996) N, ABS(FINEST(N)-1), FINEST(N)
99996   FORMAT (6X, I2, E14.2, F14.8)
X   10 CONTINUE
X      END
X      SUBROUTINE FTESTP(NDIM, Z, NFUN, F)
C
C       Selected monomials, degree 0-13
C
X      INTEGER NDIM, NFUN
X      REAL Z(NDIM), F(NFUN)
X      F(1) = 1
X      F(2) = 2*Z(1)
X      F(3) = 3*Z(1)**2
X      F(4) = F(2)*2*Z(2)
X      F(5) = 4*Z(1)**3
X      F(6) = F(3)*2*Z(2)
X      F(7) = 5*Z(1)**4
X      F(8) = F(5)*2*Z(2)
X      F(9) = F(3)*3*Z(2)**2
X      F(10) = 6*Z(1)**5
X      F(11) = F(7)*2*Z(2)
X      F(12) = F(5)*3*Z(2)**2
X      F(13) = 8*Z(1)**7
X      F(14) = 10*Z(1)**9
X      F(15) = 12*Z(1)**11
X      F(16) = 14*Z(1)**13
X      END
X      SUBROUTINE FTESTO(NDIM, Z, NFUN, F)
C
C     Corner Peak
C
X      INTEGER NDIM, NFUN
X      REAL Z(NDIM), F(NFUN)
X      F(1) = 10/(1+0.1*Z(1)+0.2*Z(2)+0.3*Z(3)+0.4*Z(4)+0.5*Z(5)+0.6*
X     * Z(6))**6/0.2057746
X      END
X      SUBROUTINE FTESTX(NDIM, Z, NFUN, F)
C
C     Sum of Cosines
C
X      INTEGER N, NDIM, NFUN
X      REAL Z(NDIM), F(NFUN), SUM
X      SUM = 0
X      DO 10 N = 1,2
X        SUM = SUM - COS(10*Z(N))/0.0544021110889370
X   10 CONTINUE
X      F(1) = SUM/2
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'stest2.f' &&
  chmod 0644 'stest2.f' ||
  echo 'restore of stest2.f failed'
  shar_count="`wc -c < 'stest2.f'`"
  test 7616 -eq "$shar_count" ||
    echo "stest2.f: original size 7616, current size $shar_count"
fi
# ============= scuhre.f ==============
if test -f 'scuhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping scuhre.f (file already exists)'
else
  echo 'x - extracting scuhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'scuhre.f' &&
X      SUBROUTINE SCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,
X     +                  EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,
X     +                  WORK)
C***BEGIN PROLOGUE SCUHRE
C***DATE WRITTEN   900116   (YYMMDD)
C***REVISION DATE  900116   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C            Alan Genz, Computer Science Department, Washington State
C            University, Pullman, WA 99163-1210, USA
C            Phone.. 509-335-2131
C            Email..  acg@cs2.cs.wsu.edu
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals
C
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
C              I   I  1  2      NDIM
C
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            SCUHRE is a driver for the integration routine
C            SADHRE, which repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with greatest
C            estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            For NDIM = 2 the default integration rule is of
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C
C            When SCUHRE computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying SCUHRE to all integrands in one call.
C            For integrals that vary continuously as functions of
C            some parameter, the estimates produced by SCUHRE will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in SCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for SCUHRE that may
C            be changed (with respect to the previous call of SCUHRE)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by SCUHRE.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for SCUHRE
C              to obtain the required accuracy. In this case SCUHRE
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN is less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2)
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
C              WRKSUB) contain the centers of the subregions.
C            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
C              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths.
C            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
C              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of
C              subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)*
C              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary
C              storage in SADHRE.
C
C
C        SCUHRE Example Test Program
C
C
C   STEST1 is a simple test driver for SCUHRE.
C
C   Output produced on a SUN 3/50.
c
C       SCUHRE TEST RESULTS
C
C    FTEST CALLS = 3549, IFAIL =  0
C   N   ESTIMATED ERROR   INTEGRAL
C   1     0.00000013     0.13850819
C   2     0.00000015     0.06369469
C   3     0.00000875     0.05861748
C   4     0.00000020     0.05407035
C   5     0.00000020     0.05005614
C   6     0.00000009     0.04654606
C
C     PROGRAM STEST1
C     EXTERNAL FTEST
C     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
C     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
C     REAL A(NDIM), B(NDIM), WRKSTR(NW)
C     REAL ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
C     DO 10 N = 1,NDIM
C        A(N) = 0
C        B(N) = 1
C  10 CONTINUE
C     MINCLS = 0
C     MAXCLS = 10000
C     KEY = 0
C     ABSREQ = 0
C     RELREQ = 1E-3
C     CALL SCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
C    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
C     PRINT 9999, NEVAL, IFAIL
C9999 FORMAT (8X, 'SCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
C    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL')
C     DO 20 N = 1,NF
C        PRINT 9998, N, ABSEST(N), FINEST(N)
C9998    FORMAT (3X, I2, 2F15.8)
C  20 CONTINUE
C     END
C     SUBROUTINE FTEST(NDIM, Z, NFUN, F)
C     INTEGER N, NDIM, NFUN
C     REAL Z(NDIM), F(NFUN), SUM
C     SUM = 0
C     DO 10 N = 1,NDIM
C        SUM = SUM + N*Z(N)**2
C  10 CONTINUE
C     F(1) = EXP(-SUM/2)
C     DO 20 N = 1,NDIM
C        F(N+1) = Z(N)*F(1)
C  20 CONTINUE
C     END
C
C***LONG DESCRIPTION
C
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to SADHRE, WORK is split into eight
C   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   OLDRES and WORK.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   CENTRS contains the centers of the subregions.
C   HWIDTS contains the half widths of the subregions.
C   GREATE contains the greatest estimated error for each subregion.
C   DIR    contains the directions for further subdivision.
C   OLDRES and WORK are used as work arrays in SADHRE.
C
C   The data structures for the subregions are in SADHRE organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   STRHRE.
C
C   The subroutine SADHRE is written for efficient execution on shared
C   memory parallel computer. On a computer with NPROC processors we wil
C   in each subdivision step divide MDIV regions, where MDIV is
C   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new region
C   Each processor will then compute estimates of the integrals and erro
C   over 2*MDIV/NPROC subregions in each subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, SRLHRE, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of subregions, compiler
C   directives should be placed in front of the DO 200
C   loop in SADHRE on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   To be published.
C
C   J.Berntsen, T.O.Espelid and A.Genz, SCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, To be published.
C
C***ROUTINES CALLED SCHHRE,SADHRE
C***END PROLOGUE SCUHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR
X      INTEGER NEVAL,IFAIL
X      REAL A(NDIM),B(NDIM),EPSABS,EPSREL
X      REAL RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MDIV   Integer.
C          MDIV is the number of subregions that are divided in
C          each subdivision step in SADHRE.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MDIV,MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS, KEY and NDIM.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN, NDIM and MDIV. This determines the length
C          of the main work arrays.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C
X      INTEGER MDIV,MAXWT,WTLENG,MAXDIM,LENW2,MAXSUB,MINSUB
X      INTEGER NUM,NSUB,LENW,KEYF
X      PARAMETER (MDIV=1)
X      PARAMETER (MAXDIM=15)
X      PARAMETER (MAXWT=14)
X      PARAMETER (LENW2=2*MDIV*MAXDIM* (MAXWT+1)+12*MAXWT+2*MAXDIM)
X      INTEGER WRKSUB,I1,I2,I3,I4,I5,I6,I7,I8,K1,K2,K3,K4,K5,K6,K7,K8
X      REAL WORK2(LENW2)
C
C***FIRST EXECUTABLE STATEMENT SCUHRE
C
C   Compute NUM, WTLENG, MAXSUB and MINSUB,
C   and check the input parameters.
C
X      CALL SCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,EPSABS,
X     +            EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,KEYF,
X     +            IFAIL,WTLENG)
X      WRKSUB = (NW - 1 - 17*MDIV*NUMFUN)/(2*NDIM + 2*NUMFUN + 2)
X      IF (IFAIL.NE.0) THEN
X          GO TO 999
X      END IF
C
C   Split up the work space.
C
X      I1 = 1
X      I2 = I1 + WRKSUB*NUMFUN
X      I3 = I2 + WRKSUB*NUMFUN
X      I4 = I3 + WRKSUB*NDIM
X      I5 = I4 + WRKSUB*NDIM
X      I6 = I5 + WRKSUB
X      I7 = I6 + WRKSUB
X      I8 = I7 + NUMFUN*MDIV
X      K1 = 1
X      K2 = K1 + 2*MDIV*WTLENG*NDIM
X      K3 = K2 + WTLENG*5
X      K4 = K3 + WTLENG
X      K5 = K4 + NDIM
X      K6 = K5 + NDIM
X      K7 = K6 + 2*MDIV*NDIM
X      K8 = K7 + 3*WTLENG
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
X      IF (RESTAR.EQ.1) THEN
X          NSUB = WORK(NW)
X      END IF
C
C   Compute the size of the temporary work space needed in SADHRE.
C
X      LENW = 16*MDIV*NUMFUN
X      CALL SADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,EPSABS,
X     +            EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,
X     +            RESULT,ABSERR,NEVAL,NSUB,IFAIL,WORK(I1),WORK(I2),
X     +            WORK(I3),WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),
X     +            WORK2(K1),WORK2(K2),WORK2(K3),WORK2(K4),WORK2(K5),
X     +            WORK2(K6),WORK2(K7),WORK2(K8))
X      WORK(NW) = NSUB
999   RETURN
C
C***END SCUHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'scuhre.f' &&
  chmod 0644 'scuhre.f' ||
  echo 'restore of scuhre.f failed'
  shar_count="`wc -c < 'scuhre.f'`"
  test 17849 -eq "$shar_count" ||
    echo "scuhre.f: original size 17849, current size $shar_count"
fi
# ============= schhre.f ==============
if test -f 'schhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping schhre.f (file already exists)'
else
  echo 'x - extracting schhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'schhre.f' &&
X      SUBROUTINE SCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,
X     +                  EPSABS,EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,
X     +                  KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE SCHHRE
C***PURPOSE  SCHHRE checks the validity of the
C            input parameters to SCUHRE.
C***DESCRIPTION
C            SCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and
C            IFAIL as functions of the input parameters to SCUHRE,
C            and checks the validity of the input parameters to SCUHRE.
C
C   ON ENTRY
C
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of subregions that are divided in
C            each subdivision step in SADHRE.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in SCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C
C   ON RETURN
C
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS, KEY and NDIM.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS, KEY and NDIM.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than
C                      MAXDIM.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     KEYF   Integer.
C            Key to selected integration rule.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE SCHHRE
C
C   Global variables.
C
X      INTEGER NDIM,NUMFUN,MDIV,MINPTS,MAXPTS,KEY,NW,MINSUB,MAXSUB
X      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG
X      REAL A(NDIM),B(NDIM),EPSABS,EPSREL
C
C   Local variables.
C
X      INTEGER LIMIT,J
C
C***FIRST EXECUTABLE STATEMENT SCHHRE
C
X      IFAIL = 0
C
C   Check on legal KEY.
C
X      IF (KEY.LT.0 .OR. KEY.GT.4) THEN
X          IFAIL = 2
X          GO TO 999
X      END IF
C
C   Check on legal NDIM.
C
X      IF (NDIM.LT.2 .OR. NDIM.GT.MAXDIM) THEN
X          IFAIL = 3
X          GO TO 999
X      END IF
C
C   For KEY = 1, NDIM must be equal to 2.
C
X      IF (KEY.EQ.1 .AND. NDIM.NE.2) THEN
X          IFAIL = 4
X          GO TO 999
X      END IF
C
C   For KEY = 2, NDIM must be equal to 3.
C
X      IF (KEY.EQ.2 .AND. NDIM.NE.3) THEN
X          IFAIL = 5
X          GO TO 999
X      END IF
C
C   For KEY = 0, we point at the selected integration rule.
C
X      IF (KEY.EQ.0) THEN
X          IF (NDIM.EQ.2) THEN
X              KEYF = 1
X          ELSE IF (NDIM.EQ.3) THEN
X              KEYF = 2
X          ELSE
X              KEYF = 3
X          ENDIF
X      ELSE
X          KEYF = KEY
X      ENDIF
C
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C
X      IF (KEYF.EQ.1) THEN
X          NUM = 65
X          WTLENG = 14
X      ELSE IF (KEYF.EQ.2) THEN
X          NUM = 127
X          WTLENG = 13
X      ELSE IF (KEYF.EQ.3) THEN
X          NUM = 1 + 4*2*NDIM + 2*NDIM* (NDIM-1) + 4*NDIM* (NDIM-1) +
X     +          4*NDIM* (NDIM-1)* (NDIM-2)/3 + 2**NDIM
X          WTLENG = 9
X          IF (NDIM.EQ.2) WTLENG = 8
X      ELSE IF (KEYF.EQ.4) THEN
X          NUM = 1 + 3*2*NDIM + 2*NDIM* (NDIM-1) + 2**NDIM
X          WTLENG = 6
X      END IF
C
C   Compute MAXSUB.
C
X      MAXSUB = (MAXPTS-NUM)/ (2*NUM) + 1
C
C   Compute MINSUB.
C
X      MINSUB = (MINPTS-NUM)/ (2*NUM) + 1
X      IF (MOD(MINPTS-NUM,2*NUM).NE.0) THEN
X          MINSUB = MINSUB + 1
X      END IF
X      MINSUB = MAX(2,MINSUB)
C
C   Check on positive NUMFUN.
C
X      IF (NUMFUN.LT.1) THEN
X          IFAIL = 6
X          GO TO 999
X      END IF
C
C   Check on legal upper and lower limits of integration.
C
X      DO 10 J = 1,NDIM
X          IF (A(J)-B(J).EQ.0) THEN
X              IFAIL = 7
X              GO TO 999
X          END IF
10    CONTINUE
C
C   Check on MAXPTS < 3*NUM.
C
X      IF (MAXPTS.LT.3*NUM) THEN
X          IFAIL = 8
X          GO TO 999
X      END IF
C
C   Check on MAXPTS >= MINPTS.
C
X      IF (MAXPTS.LT.MINPTS) THEN
X          IFAIL = 9
X          GO TO 999
X      END IF
C
C   Check on legal accuracy requests.
C
X      IF (EPSABS.LT.0 .AND. EPSREL.LT.0) THEN
X          IFAIL = 10
X          GO TO 999
X      END IF
C
C   Check on big enough double precision workspace.
C
X      LIMIT = MAXSUB* (2*NDIM+2*NUMFUN+2) + 17*MDIV*NUMFUN + 1
X      IF (NW.LT.LIMIT) THEN
X          IFAIL = 11
X          GO TO 999
X      END IF
C
C    Check on legal RESTAR.
C
X      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
X          IFAIL = 12
X          GO TO 999
X      END IF
999   RETURN
C
C***END SCHHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'schhre.f' &&
  chmod 0644 'schhre.f' ||
  echo 'restore of schhre.f failed'
  shar_count="`wc -c < 'schhre.f'`"
  test 8316 -eq "$shar_count" ||
    echo "schhre.f: original size 8316, current size $shar_count"
fi
# ============= sadhre.f ==============
if test -f 'sadhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping sadhre.f (file already exists)'
else
  echo 'x - extracting sadhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'sadhre.f' &&
X      SUBROUTINE SADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,
X     +                  EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,
X     +                  RESULT,ABSERR,NEVAL,NSUB,IFAIL,VALUES,
X     +                  ERRORS,CENTRS,HWIDTS,GREATE,DIR,OLDRES,WORK,G,W,
X     +                  RULPTS,CENTER,HWIDTH,X,SCALES,NORMS)
C***BEGIN PROLOGUE SADHRE
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals, I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            SADHRE repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with  greatest
C            estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The regions are divided in two equally sized parts along
C            the direction with greatest absolute fourth divided
C            difference.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            Defines the number of new subregions that are divided
C            in each subdivision step.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINSUB Integer.
C            The computations proceed until there are at least
C            MINSUB subregions in the data structure.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from SADHRE
C            must not be changed since the last
C            exit from SADHRE.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*MDIV*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to SADHRE.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by SADHRE.
C     NSUB   Integer.
C            Number of stored subregions.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for SADHRE
C              to obtain the required accuracy. In this case SADHRE
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB.
C            DIR is used to store the directions for
C            further subdivision.
C     OLDRES Real array of dimension (NUMFUN,MDIV).
C            Used to store old estimates of the integrals over the
C            subregions.
C     WORK   Real array of dimension LENW.
C            Used  in SRLHRE and STRHRE.
C     G      Real array of dimension (NDIM,WTLENG,2*MDIV).
C            The fully symmetric sum generators for the rules.
C            G(1,J,1),...,G(NDIM,J,1) are the generators for the
C            points associated with the Jth weights.
C            When MDIV subregions are divided in 2*MDIV
C            subregions, the subregions may be processed on different
C            processors and we must make a copy of the generators
C            for each processor.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in SINHRE.
C     CENTER Real array of dimension NDIM.
C            Work array used in STRHRE.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in STRHRE.
C     X      Real array of dimension (NDIM,2*MDIV).
C            Work array used in SRLHRE.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by SINHRE and SRLHRE.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by SINHRE and SRLHRE.
C
C***REFERENCES
C
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
C
C***ROUTINES CALLED STRHRE,SINHRE,SRLHRE
C***END PROLOGUE SADHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN,MDIV,MINSUB,MAXSUB,KEY,LENW,RESTAR
X      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG
X      REAL A(NDIM),B(NDIM),EPSABS,EPSREL
X      REAL RESULT(NUMFUN),ABSERR(NUMFUN)
X      REAL VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB)
X      REAL CENTRS(NDIM,MAXSUB)
X      REAL HWIDTS(NDIM,MAXSUB)
X      REAL GREATE(MAXSUB),DIR(MAXSUB)
X      REAL OLDRES(NUMFUN,MDIV)
X      REAL WORK(LENW),RULPTS(WTLENG)
X      REAL G(NDIM,WTLENG,2*MDIV),W(5,WTLENG)
X      REAL CENTER(NDIM),HWIDTH(NDIM),X(NDIM,2*MDIV)
X      REAL SCALES(3,WTLENG),NORMS(3,WTLENG)
C
C   Local variables.
C
C   INTSGN is used to get correct sign on the integral.
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coeff. defined in SINHRE and used by SRLHRE
C          and SADHRE.
C
X      INTEGER I,J,K
X      INTEGER INTSGN,SBRGNS
X      INTEGER L1
X      INTEGER NDIV,POINTR,DIRECT,INDEX
X      REAL OLDCEN,EST1,EST2,ERRCOF(6)
C
C***FIRST EXECUTABLE STATEMENT SADHRE
C
C   Get the correct sign on the integral.
C
X      INTSGN = 1
X      DO 10 J = 1,NDIM
X          IF (B(J).LT.A(J)) THEN
X              INTSGN = - INTSGN
X          END IF
10    CONTINUE
C
C   Call SINHRE to compute the weights and abscissas of
C   the function evaluation points.
C
X      CALL SINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C
C   If RESTAR = 1, then this is a restart run.
C
X      IF (RESTAR.EQ.1) THEN
X          SBRGNS = NSUB
X          GO TO 110
X      END IF
C
C   Initialize the SBRGNS, CENTRS and HWIDTS.
C
X      SBRGNS = 1
X      DO 15 J = 1,NDIM
X          CENTRS(J,1) = (A(J)+B(J))/2
X          HWIDTS(J,1) = ABS(B(J)-A(J))/2
15    CONTINUE
C
C   Initialize RESULT, ABSERR and NEVAL.
C
X      DO 20 J = 1,NUMFUN
X          RESULT(J) = 0
X          ABSERR(J) = 0
20    CONTINUE
X      NEVAL = 0
C
C   Apply SRLHRE over the whole region.
C
X      CALL SRLHRE(NDIM,CENTRS(1,1),HWIDTS(1,1),WTLENG,G,W,ERRCOF,NUMFUN,
X     +            FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,1),ERRORS(1,1),
X     +            DIR(1))
X      NEVAL = NEVAL + NUM
C
C   Add the computed values to RESULT and ABSERR.
C
X      DO 55 J = 1,NUMFUN
X          RESULT(J) = RESULT(J) + VALUES(J,1)
55    CONTINUE
X      DO 65 J = 1,NUMFUN
X          ABSERR(J) = ABSERR(J) + ERRORS(J,1)
65    CONTINUE
C
C   Store results in heap.
C
X      INDEX = 1
X      CALL STRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,HWIDTS,
X     +            GREATE,WORK(1),WORK(NUMFUN+1),CENTER,HWIDTH,DIR)
C
C***End initialisation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+1 is less than MAXSUB.
C
110   IF (SBRGNS+1.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over each half of the
C   NDIV subregions with greatest errors.
C   If MAXSUB is great enough, NDIV = MDIV
C
X          IF (MDIV.GT.1) THEN
X              NDIV = MAXSUB - SBRGNS
X              NDIV = MIN(NDIV,MDIV,SBRGNS)
X          ELSE
X              NDIV = 1
X          END IF
C
C   Divide the NDIV subregions in two halves, and compute
C   integral and error over each half.
C
X          DO 150 I = 1,NDIV
X              POINTR = SBRGNS + NDIV + 1 - I
C
C   Adjust RESULT and ABSERR.
C
X              DO 115 J = 1,NUMFUN
X                  RESULT(J) = RESULT(J) - VALUES(J,1)
X                  ABSERR(J) = ABSERR(J) - ERRORS(J,1)
115           CONTINUE
C
C   Compute first half region.
C
X              DO 120 J = 1,NDIM
X                  CENTRS(J,POINTR) = CENTRS(J,1)
X                  HWIDTS(J,POINTR) = HWIDTS(J,1)
120           CONTINUE
X              DIRECT = DIR(1)
X              DIR(POINTR) = DIRECT
X              HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,1)/2
X              OLDCEN = CENTRS(DIRECT,1)
X              CENTRS(DIRECT,POINTR) = OLDCEN - HWIDTS(DIRECT,POINTR)
C
C   Save the computed values of the integrals.
C
X              DO 125 J = 1,NUMFUN
X                  OLDRES(J,NDIV-I+1) = VALUES(J,1)
125           CONTINUE
C
C   Adjust the heap.
C
X              CALL STRHRE(1,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
X     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
X     +                    HWIDTH,DIR)
C
C   Compute second half region.
C
X              DO 130 J = 1,NDIM
X                  CENTRS(J,POINTR-1) = CENTRS(J,POINTR)
X                  HWIDTS(J,POINTR-1) = HWIDTS(J,POINTR)
130           CONTINUE
X              CENTRS(DIRECT,POINTR-1) = OLDCEN + HWIDTS(DIRECT,POINTR)
X              HWIDTS(DIRECT,POINTR-1) = HWIDTS(DIRECT,POINTR)
X              DIR(POINTR-1) = DIRECT
150       CONTINUE
C
C   Make copies of the generators for each processor.
C
X          DO 190 I = 2,2*NDIV
X              DO 190 J = 1,NDIM
X                  DO 190 K = 1,WTLENG
X                      G(J,K,I) = G(J,K,1)
190       CONTINUE
C
C   Apply basic rule.
C
Cvd$l cncall
X          DO 200 I = 1,2*NDIV
X              INDEX = SBRGNS + I
X              L1 = 1 + (I-1)*8*NUMFUN
X              CALL SRLHRE(NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,
X     +                    G(1,1,I),W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,
X     +                    X(1,I),WORK(L1),VALUES(1,INDEX),
X     +                    ERRORS(1,INDEX),DIR(INDEX))
200       CONTINUE
X          NEVAL = NEVAL + 2*NDIV*NUM
C
C   Add new contributions to RESULT.
C
X          DO 220 I = 1,2*NDIV
X              DO 210 J = 1,NUMFUN
X                  RESULT(J) = RESULT(J) + VALUES(J,SBRGNS+I)
210           CONTINUE
220       CONTINUE
C
C   Check consistency of results and if necessary adjust
C   the estimated errors.
C
X          DO 240 I = 1,NDIV
X              GREATE(SBRGNS+2*I-1) = 0
X              GREATE(SBRGNS+2*I) = 0
X              DO 230 J = 1,NUMFUN
X                  EST1 = ABS(OLDRES(J,I)- (VALUES(J,
X     +                   SBRGNS+2*I-1)+VALUES(J,SBRGNS+2*I)))
X                  EST2 = ERRORS(J,SBRGNS+2*I-1) + ERRORS(J,SBRGNS+2*I)
X                  IF (EST2.GT.0) THEN
X                      ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)*
X     +                  (1+ERRCOF(5)*EST1/EST2)
X                      ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)*
X     +                                       (1+ERRCOF(5)*EST1/EST2)
X                  END IF
X                  ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1) +
X     +                                     ERRCOF(6)*EST1
X                  ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I) +
X     +                                   ERRCOF(6)*EST1
X                  IF (ERRORS(J,SBRGNS+2*I-1).GT.
X     +                GREATE(SBRGNS+2*I-1)) THEN
X                      GREATE(SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)
X                  END IF
X                  IF (ERRORS(J,SBRGNS+2*I).GT.GREATE(SBRGNS+2*I)) THEN
X                      GREATE(SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)
X                  END IF
X                  ABSERR(J) = ABSERR(J) + ERRORS(J,SBRGNS+2*I-1) +
X     +                        ERRORS(J,SBRGNS+2*I)
230           CONTINUE
240       CONTINUE
C
C   Store results in heap.
C
X          DO 250 I = 1,2*NDIV
X              INDEX = SBRGNS + I
X              CALL STRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,
X     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
X     +                    HWIDTH,DIR)
250       CONTINUE
X          SBRGNS = SBRGNS + 2*NDIV
C
C   Check for termination.
C
X          IF (SBRGNS.LT.MINSUB) THEN
X              GO TO 110
X          END IF
X          DO 255 J = 1,NUMFUN
X              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
X     +            ABSERR(J).GT.EPSABS) THEN
X                  GO TO 110
X              END IF
255       CONTINUE
X          IFAIL = 0
X          GO TO 499
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
X      ELSE
X          IFAIL = 1
X      END IF
C
C   Compute more accurate values of RESULT and ABSERR.
C
499   CONTINUE
X      DO 500 J = 1,NUMFUN
X          RESULT(J) = 0
X          ABSERR(J) = 0
500   CONTINUE
X      DO 510 I = 1,SBRGNS
X          DO 505 J = 1,NUMFUN
X              RESULT(J) = RESULT(J) + VALUES(J,I)
X              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
505       CONTINUE
510   CONTINUE
C
C   Compute correct sign on the integral.
C
X      DO 600 J = 1,NUMFUN
X          RESULT(J) = RESULT(J)*INTSGN
600   CONTINUE
X      NSUB = SBRGNS
X      RETURN
C
C***END SADHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'sadhre.f' &&
  chmod 0644 'sadhre.f' ||
  echo 'restore of sadhre.f failed'
  shar_count="`wc -c < 'sadhre.f'`"
  test 16487 -eq "$shar_count" ||
    echo "sadhre.f: original size 16487, current size $shar_count"
fi
# ============= strhre.f ==============
if test -f 'strhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping strhre.f (file already exists)'
else
  echo 'x - extracting strhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'strhre.f' &&
X      SUBROUTINE STRHRE(DVFLAG,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
X     +                  HWIDTS,GREATE,ERROR,VALUE,CENTER,HWIDTH,DIR)
C***BEGIN PROLOGUE STRHRE
C***PURPOSE STRHRE maintains a heap of subregions.
C***DESCRIPTION STRHRE maintains a heap of subregions.
C            The subregions are ordered according to the size
C            of the greatest error estimates of each subregion(GREATE).
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     NDIM   Integer.
C            Number of variables.
C     NUMFUN Integer.
C            Number of components of the integral.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     VALUES Real array of dimension (NUMFUN,SBRGNS).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,SBRGNS).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,SBRGNS).
C            Used to store the center limits of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,SBRGNS).
C            Used to store the hwidth limits of the stored subregions.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     ERROR  Real array of dimension NUMFUN.
C            Used as intermediate storage for the error of a subregion.
C     VALUE  Real array of dimension NUMFUN.
C            Used as intermediate storage for the estimate
C            of the integral over a subregion.
C     CENTER Real array of dimension NDIM.
C            Used as intermediate storage for the center of
C            the subregion.
C     HWIDTH Real array of dimension NDIM.
C            Used as intermediate storage for the half width of
C            the subregion.
C     DIR    Integer array of dimension SBRGNS.
C            DIR is used to store the directions for
C            further subdivision.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE STRHRE
C
C   Global variables.
C
X      INTEGER DVFLAG,NDIM,NUMFUN,SBRGNS
X      REAL VALUES(NUMFUN,*),ERRORS(NUMFUN,*)
X      REAL CENTRS(NDIM,*)
X      REAL HWIDTS(NDIM,*)
X      REAL GREATE(*)
X      REAL ERROR(NUMFUN),VALUE(NUMFUN)
X      REAL CENTER(NDIM),HWIDTH(NDIM)
X      REAL DIR(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   DIRECT is used as intermediate storage for the direction of further
C          subdivision.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
C
X      INTEGER J,SUBRGN,SUBTMP
X      REAL GREAT,DIRECT
C
C***FIRST EXECUTABLE STATEMENT STRHRE
C
C   Save values to be stored in their correct place in the heap.
C
X      GREAT = GREATE(SBRGNS)
X      DIRECT = DIR(SBRGNS)
X      DO 5 J = 1,NUMFUN
X          ERROR(J) = ERRORS(J,SBRGNS)
X          VALUE(J) = VALUES(J,SBRGNS)
5     CONTINUE
X      DO 10 J = 1,NDIM
X          CENTER(J) = CENTRS(J,SBRGNS)
X          HWIDTH(J) = HWIDTS(J,SBRGNS)
10    CONTINUE
C
C    If DVFLAG = 1, we will remove the region
C    with greatest estimated error from the heap.
C
X      IF (DVFLAG.EQ.1) THEN
X          SBRGNS = SBRGNS - 1
X          SUBRGN = 1
20        SUBTMP = 2*SUBRGN
X          IF (SUBTMP.LE.SBRGNS) THEN
X              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
X                  IF (GREATE(SUBTMP).LT.GREATE(SUBTMP+1)) THEN
X                      SUBTMP = SUBTMP + 1
X                  END IF
X              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
X              IF (GREAT.LT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp up the heap.
C
X                  GREATE(SUBRGN) = GREATE(SUBTMP)
X                  DO 25 J = 1,NUMFUN
X                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
X                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
25                CONTINUE
X                  DIR(SUBRGN) = DIR(SUBTMP)
X                  DO 30 J = 1,NDIM
X                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
X                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
30                CONTINUE
X                  SUBRGN = SUBTMP
X                  GO TO 20
X              END IF
X          END IF
X      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, then insert new region in the heap.
C
X          SUBRGN = SBRGNS
40        SUBTMP = SUBRGN/2
X          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
X              IF (GREAT.GT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp down the heap.
C
X                  GREATE(SUBRGN) = GREATE(SUBTMP)
X                  DO 45 J = 1,NUMFUN
X                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
X                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
45                CONTINUE
X                  DIR(SUBRGN) = DIR(SUBTMP)
X                  DO 50 J = 1,NDIM
X                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
X                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
50                CONTINUE
X                  SUBRGN = SUBTMP
X                  GO TO 40
X              END IF
X          END IF
X      END IF
C
C    Insert the saved values in their correct places.
C
X      IF (SBRGNS.GT.0) THEN
X          GREATE(SUBRGN) = GREAT
X          DO 55 J = 1,NUMFUN
X              ERRORS(J,SUBRGN) = ERROR(J)
X              VALUES(J,SUBRGN) = VALUE(J)
55        CONTINUE
X          DIR(SUBRGN) = DIRECT
X          DO 60 J = 1,NDIM
X              CENTRS(J,SUBRGN) = CENTER(J)
X              HWIDTS(J,SUBRGN) = HWIDTH(J)
60        CONTINUE
X      END IF
C
C***END STRHRE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'strhre.f' &&
  chmod 0644 'strhre.f' ||
  echo 'restore of strhre.f failed'
  shar_count="`wc -c < 'strhre.f'`"
  test 5837 -eq "$shar_count" ||
    echo "strhre.f: original size 5837, current size $shar_count"
fi
# ============= sinhre.f ==============
if test -f 'sinhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping sinhre.f (file already exists)'
else
  echo 'x - extracting sinhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'sinhre.f' &&
X      SUBROUTINE SINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C***BEGIN PROLOGUE SINHRE
C***PURPOSE SINHRE computes abscissas and weights of the integration
C            rule and the null rules to be used in error estimation.
C            These are computed as functions of NDIM and KEY.
C***DESCRIPTION SINHRE will for given values of NDIM and KEY compute or
C            select the correct values of the abscissas and
C            corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W.
C            The heuristic error coefficients ERRCOF
C            will be computed as a function of KEY.
C            Scaling factors SCALES and normalization factors NORMS
C            used in the error estimation are computed.
C
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables.
C     KEY    Integer.
C            Key to selected local integration rule.
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C            It is assumed that the error is computed using:
C             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C               THEN ERROR = ERRCOF(3)*N1
C               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C            where N1-N3 are the null rules, EP is the error for
C            the parent
C            subregion and ES is the error for the sibling subregion.
C     RULPTS Real array of dimension WTLENG.
C            A work array containing the number of points produced by
C            each generator of the selected rule.
C     SCALES Real array of dimension (3,WTLENG).
C            Scaling factors used to construct new null rules,
C            N1, N2 and N3,
C            based on a linear combination of two successive null rules
C            in the sequence of null rules.
C     NORMS  Real array of dimension (3,WTLENG).
C            2**NDIM/(1-norm of the null rule constructed by each of
C            the scaling factors.)
C
C***ROUTINES CALLED  S132RE,S113RE,S07HRE,S09HRE
C***END PROLOGUE SINHRE
C
C   Global variables.
C
X      INTEGER NDIM,KEY,WTLENG
X      REAL G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
X      REAL RULPTS(WTLENG),SCALES(3,WTLENG)
X      REAL NORMS(3,WTLENG)
C
C   Local variables.
C
X      INTEGER I,J,K
X      REAL WE(14)
C
C***FIRST EXECUTABLE STATEMENT SINHRE
C
C   Compute W, G and ERRCOF.
C
X      IF (KEY.EQ.1) THEN
X          CALL S132RE(WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.2) THEN
X          CALL S113RE(WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.3) THEN
X          CALL S09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
X      ELSE IF (KEY.EQ.4) THEN
X          CALL S07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
X      END IF
C
C   Compute SCALES and NORMS.
C
X      DO 100 K = 1,3
X          DO 50 I = 1,WTLENG
X              IF (W(K+1,I).NE.0) THEN
X                  SCALES(K,I) = - W(K+2,I)/W(K+1,I)
X              ELSE
X                  SCALES(K,I) = 100
X              END IF
X              DO 30 J = 1,WTLENG
X                  WE(J) = W(K+2,J) + SCALES(K,I)*W(K+1,J)
30            CONTINUE
X              NORMS(K,I) = 0
X              DO 40 J = 1,WTLENG
X                  NORMS(K,I) = NORMS(K,I) + RULPTS(J)*ABS(WE(J))
40            CONTINUE
X              NORMS(K,I) = 2**NDIM/NORMS(K,I)
50        CONTINUE
100   CONTINUE
X      RETURN
C
C***END SINHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'sinhre.f' &&
  chmod 0644 'sinhre.f' ||
  echo 'restore of sinhre.f failed'
  shar_count="`wc -c < 'sinhre.f'`"
  test 4036 -eq "$shar_count" ||
    echo "sinhre.f: original size 4036, current size $shar_count"
fi
# ============= s132re.f ==============
if test -f 's132re.f' && test X"$1" != X"-c"; then
  echo 'x - skipping s132re.f (file already exists)'
else
  echo 'x - extracting s132re.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 's132re.f' &&
X      SUBROUTINE S132RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE S132RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE S132RE computes abscissas and weights of a 2 dimensional
C            integration rule of degree 13.
C            Two null rules of degree 11, one null rule of degree 9
C            and one null rule of degree 7 to be used in error
C            estimation are also computed.
C ***DESCRIPTION S132RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W. The heuristic error coefficients ERRCOF
C            will also be assigned.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points produced by each generator.
C***REFERENCES S.Eriksen,
C              Thesis of the degree cand.scient, Dept. of Informatics,
C              Univ. of Bergen,Norway, 1984.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE S132RE
C
C   Global variables
C
X      INTEGER WTLENG
X      REAL W(5,WTLENG),G(2,WTLENG),ERRCOF(6)
X      REAL RULPTS(WTLENG)
C
C   Local variables.
C
X      INTEGER I,J
X      REAL DIM2G(16)
X      REAL DIM2W(14,5)
C
X      DATA (DIM2G(I),I=1,16)/0.2517129343453109E+00,
X     +     0.7013933644534266E+00,0.9590960631619962E+00,
X     +     0.9956010478552127E+00,0.5000000000000000E+00,
X     +     0.1594544658297559E+00,0.3808991135940188E+00,
X     +     0.6582769255267192E+00,0.8761473165029315E+00,
X     +     0.9982431840531980E+00,0.9790222658168462E+00,
X     +     0.6492284325645389E+00,0.8727421201131239E+00,
X     +     0.3582614645881228E+00,0.5666666666666666E+00,
X     +     0.2077777777777778E+00/
C
X      DATA (DIM2W(I,1),I=1,14)/0.3379692360134460E-01,
X     +     0.9508589607597761E-01,0.1176006468056962E+00,
X     +     0.2657774586326950E-01,0.1701441770200640E-01,
X     +     0.0000000000000000E+00,0.1626593098637410E-01,
X     +     0.1344892658526199E+00,0.1328032165460149E+00,
X     +     0.5637474769991870E-01,0.3908279081310500E-02,
X     +     0.3012798777432150E-01,0.1030873234689166E+00,
X     +     0.6250000000000000E-01/
C
X      DATA (DIM2W(I,2),I=1,14)/0.3213775489050763E+00,
X     +     - .1767341636743844E+00,0.7347600537466072E-01,
X     +     - .3638022004364754E-01,0.2125297922098712E-01,
X     +     0.1460984204026913E+00,0.1747613286152099E-01,
X     +     0.1444954045641582E+00,0.1307687976001325E-03,
X     +     0.5380992313941161E-03,0.1042259576889814E-03,
X     +     - .1401152865045733E-02,0.8041788181514763E-02,
X     +     - .1420416552759383E+00/
C
X      DATA (DIM2W(I,3),I=1,14)/0.3372900883288987E+00,
X     +     - .1644903060344491E+00,0.7707849911634622E-01,
X     +     - .3804478358506310E-01,0.2223559940380806E-01,
X     +     0.1480693879765931E+00,0.4467143702185814E-05,
X     +     0.1508944767074130E+00,0.3647200107516215E-04,
X     +     0.5777198999013880E-03,0.1041757313688177E-03,
X     +     - .1452822267047819E-02,0.8338339968783705E-02,
X     +     - .1472796329231960E+00/
C
X      DATA (DIM2W(I,4),I=1,14)/ - .8264123822525677E+00,
X     +     0.3065838614094360E+00,0.2389292538329435E-02,
X     +     - .1343024157997222E+00,0.8833366840533900E-01,
X     +     0.0000000000000000E+00,0.9786283074168292E-03,
X     +     - .1319227889147519E+00,0.7990012200150630E-02,
X     +     0.3391747079760626E-02,0.2294915718283264E-02,
X     +     - .1358584986119197E-01,0.4025866859057809E-01,
X     +     0.3760268580063992E-02/
C
X      DATA (DIM2W(I,5),I=1,14)/0.6539094339575232E+00,
X     +     - .2041614154424632E+00, - .1746981515794990E+00,
X     +     0.3937939671417803E-01,0.6974520545933992E-02,
X     +     0.0000000000000000E+00,0.6667702171778258E-02,
X     +     0.5512960621544304E-01,0.5443846381278607E-01,
X     +     0.2310903863953934E-01,0.1506937747477189E-01,
X     +     - .6057021648901890E-01,0.4225737654686337E-01,
X     +     0.2561989142123099E-01/
C
C***FIRST EXECUTABLE STATEMENT S132RE
C
C   Assign values to W.
C
X      DO 10 I = 1,14
X          DO 10 J = 1,5
X              W(J,I) = DIM2W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
X      DO 20 I = 1,2
X          DO 20 J = 1,14
X              G(I,J) = 0
20    CONTINUE
X      G(1,2) = DIM2G(1)
X      G(1,3) = DIM2G(2)
X      G(1,4) = DIM2G(3)
X      G(1,5) = DIM2G(4)
X      G(1,6) = DIM2G(5)
X      G(1,7) = DIM2G(6)
X      G(2,7) = G(1,7)
X      G(1,8) = DIM2G(7)
X      G(2,8) = G(1,8)
X      G(1,9) = DIM2G(8)
X      G(2,9) = G(1,9)
X      G(1,10) = DIM2G(9)
X      G(2,10) = G(1,10)
X      G(1,11) = DIM2G(10)
X      G(2,11) = G(1,11)
X      G(1,12) = DIM2G(11)
X      G(2,12) = DIM2G(12)
X      G(1,13) = DIM2G(13)
X      G(2,13) = DIM2G(14)
X      G(1,14) = DIM2G(15)
X      G(2,14) = DIM2G(16)
C
C   Assign values to RULPTS.
C
X      RULPTS(1) = 1
X      DO 30 I = 2,11
X          RULPTS(I) = 4
30    CONTINUE
X      RULPTS(12) = 8
X      RULPTS(13) = 8
X      RULPTS(14) = 8
C
C   Assign values to ERRCOF.
C
X      ERRCOF(1) = 10
X      ERRCOF(2) = 10
X      ERRCOF(3) = 1.
X      ERRCOF(4) = 5.
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END S132RE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 's132re.f' &&
  chmod 0644 's132re.f' ||
  echo 'restore of s132re.f failed'
  shar_count="`wc -c < 's132re.f'`"
  test 5904 -eq "$shar_count" ||
    echo "s132re.f: original size 5904, current size $shar_count"
fi
# ============= s113re.f ==============
if test -f 's113re.f' && test X"$1" != X"-c"; then
  echo 'x - skipping s113re.f (file already exists)'
else
  echo 'x - extracting s113re.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 's113re.f' &&
X      SUBROUTINE S113RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE S113RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE S113RE computes abscissas and weights of a 3 dimensional
C            integration rule of degree 11.
C            Two null rules of degree 9, one null rule of degree 7
C            and one null rule of degree 5 to be used in error
C            estimation are also computed.
C***DESCRIPTION S113RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to G
C            and W.
C            The heuristic error coefficients ERRCOF
C            will also be computed.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points used by each generator.
C
C***REFERENCES  J.Berntsen, Cautious adaptive numerical integration
C               over the 3-cube, Reports in Informatics 17, Dept. of
C               Inf.,Univ. of Bergen, Norway, 1985.
C               J.Berntsen and T.O.Espelid, On the construction of
C               higher degree three-dimensional embedded integration
C               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234,
C               1988.
C***ROUTINES CALLED-NONE
C***END PROLOGUE S113RE
C
C   Global variables.
C
X      INTEGER WTLENG
X      REAL W(5,WTLENG),G(3,WTLENG),ERRCOF(6)
X      REAL RULPTS(WTLENG)
C
C   Local variables.
C
X      INTEGER I,J
X      REAL DIM3G(14)
X      REAL DIM3W(13,5)
C
X      DATA (DIM3G(I),I=1,14)/0.1900000000000000E+00,
X     +     0.5000000000000000E+00,0.7500000000000000E+00,
X     +     0.8000000000000000E+00,0.9949999999999999E+00,
X     +     0.9987344998351400E+00,0.7793703685672423E+00,
X     +     0.9999698993088767E+00,0.7902637224771788E+00,
X     +     0.4403396687650737E+00,0.4378478459006862E+00,
X     +     0.9549373822794593E+00,0.9661093133630748E+00,
X     +     0.4577105877763134E+00/
C
X      DATA (DIM3W(I,1),I=1,13)/0.7923078151105734E-02,
X     +     0.6797177392788080E-01,0.1086986538805825E-02,
X     +     0.1838633662212829E+00,0.3362119777829031E-01,
X     +     0.1013751123334062E-01,0.1687648683985235E-02,
X     +     0.1346468564512807E+00,0.1750145884600386E-02,
X     +     0.7752336383837454E-01,0.2461864902770251E+00,
X     +     0.6797944868483039E-01,0.1419962823300713E-01/
C
X      DATA (DIM3W(I,2),I=1,13)/0.1715006248224684E+01,
X     +     - .3755893815889209E+00,0.1488632145140549E+00,
X     +     - .2497046640620823E+00,0.1792501419135204E+00,
X     +     0.3446126758973890E-02, - .5140483185555825E-02,
X     +     0.6536017839876425E-02, - .6513454939229700E-03,
X     +     - .6304672433547204E-02,0.1266959399788263E-01,
X     +     - .5454241018647931E-02,0.4826995274768427E-02/
C
X      DATA (DIM3W(I,3),I=1,13)/0.1936014978949526E+01,
X     +     - .3673449403754268E+00,0.2929778657898176E-01,
X     +     - .1151883520260315E+00,0.5086658220872218E-01,
X     +     0.4453911087786469E-01, - .2287828257125900E-01,
X     +     0.2908926216345833E-01, - .2898884350669207E-02,
X     +     - .2805963413307495E-01,0.5638741361145884E-01,
X     +     - .2427469611942451E-01,0.2148307034182882E-01/
C
X      DATA (DIM3W(I,4),I=1,13)/0.5170828195605760E+00,
X     +     0.1445269144914044E-01, - .3601489663995932E+00,
X     +     0.3628307003418485E+00,0.7148802650872729E-02,
X     +     - .9222852896022966E-01,0.1719339732471725E-01,
X     +     - .1021416537460350E+00, - .7504397861080493E-02,
X     +     0.1648362537726711E-01,0.5234610158469334E-01,
X     +     0.1445432331613066E-01,0.3019236275367777E-02/
C
X      DATA (DIM3W(I,5),I=1,13)/0.2054404503818520E+01,
X     +     0.1377759988490120E-01, - .5768062917904410E+00,
X     +     0.3726835047700328E-01,0.6814878939777219E-02,
X     +     0.5723169733851849E-01, - .4493018743811285E-01,
X     +     0.2729236573866348E-01,0.3547473950556990E-03,
X     +     0.1571366799739551E-01,0.4990099219278567E-01,
X     +     0.1377915552666770E-01,0.2878206423099872E-02/
C
C***FIRST EXECUTABLE STATEMENT S113RE
C
C   Assign values to W.
C
X      DO 10 I = 1,13
X          DO 10 J = 1,5
X              W(J,I) = DIM3W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
X      DO 20 I = 1,3
X          DO 20 J = 1,13
X              G(I,J) = 0
20    CONTINUE
X      G(1,2) = DIM3G(1)
X      G(1,3) = DIM3G(2)
X      G(1,4) = DIM3G(3)
X      G(1,5) = DIM3G(4)
X      G(1,6) = DIM3G(5)
X      G(1,7) = DIM3G(6)
X      G(2,7) = G(1,7)
X      G(1,8) = DIM3G(7)
X      G(2,8) = G(1,8)
X      G(1,9) = DIM3G(8)
X      G(2,9) = G(1,9)
X      G(3,9) = G(1,9)
X      G(1,10) = DIM3G(9)
X      G(2,10) = G(1,10)
X      G(3,10) = G(1,10)
X      G(1,11) = DIM3G(10)
X      G(2,11) = G(1,11)
X      G(3,11) = G(1,11)
X      G(1,12) = DIM3G(12)
X      G(2,12) = DIM3G(11)
X      G(3,12) = G(2,12)
X      G(1,13) = DIM3G(13)
X      G(2,13) = G(1,13)
X      G(3,13) = DIM3G(14)
C
C   Assign values to RULPTS.
C
X      RULPTS(1) = 1
X      RULPTS(2) = 6
X      RULPTS(3) = 6
X      RULPTS(4) = 6
X      RULPTS(5) = 6
X      RULPTS(6) = 6
X      RULPTS(7) = 12
X      RULPTS(8) = 12
X      RULPTS(9) = 8
X      RULPTS(10) = 8
X      RULPTS(11) = 8
X      RULPTS(12) = 24
X      RULPTS(13) = 24
C
C   Assign values to ERRCOF.
C
X      ERRCOF(1) = 4
X      ERRCOF(2) = 4.
X      ERRCOF(3) = 0.5
X      ERRCOF(4) = 3.
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END S113RE
C
X      RETURN
X      END
SHAR_EOF
  $shar_touch -am 0531091995 's113re.f' &&
  chmod 0644 's113re.f' ||
  echo 'restore of s113re.f failed'
  shar_count="`wc -c < 's113re.f'`"
  test 6155 -eq "$shar_count" ||
    echo "s113re.f: original size 6155, current size $shar_count"
fi
# ============= s09hre.f ==============
if test -f 's09hre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping s09hre.f (file already exists)'
else
  echo 'x - extracting s09hre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 's09hre.f' &&
X      SUBROUTINE S09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE S09HRE
C***KEYWORDS basic integration rule, degree 9
C***PURPOSE  To initialize a degree 9 basic rule and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-20
C***DESCRIPTION  S09HRE initializes a degree 9 integration rule,
C            two degree 7 null rules, one degree 5 null rule and one
C            degree 3 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE S09HRE
C
C   Global variables
C
X      INTEGER NDIM,WTLENG
X      REAL W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
X      REAL RULPTS(WTLENG)
C
C   Local Variables
C
X      REAL RATIO,LAM0,LAM1,LAM2,LAM3,LAMP,TWONDM
X      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT S09HRE
C
C
C     Initialize generators, weights and RULPTS
C
X      DO 30 J = 1,WTLENG
X          DO 10 I = 1,NDIM
X              G(I,J) = 0
10        CONTINUE
X          DO 20 I = 1,5
X              W(I,J) = 0
20        CONTINUE
X          RULPTS(J) = 2*NDIM
30    CONTINUE
X      TWONDM = 2**NDIM
X      RULPTS(WTLENG) = TWONDM
X      IF (NDIM.GT.2) RULPTS(8) = (4*NDIM* (NDIM-1)* (NDIM-2))/3
X      RULPTS(7) = 4*NDIM* (NDIM-1)
X      RULPTS(6) = 2*NDIM* (NDIM-1)
X      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
X      LAM0 = 0.4707
X      LAM1 = 4/ (15-5/LAM0)
X      RATIO = (1-LAM1/LAM0)/27
X      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
X      RATIO = RATIO* (1-LAM2/LAM0)/3
X      LAM3 = (7-9* (LAM2+LAM1)+63*LAM2*LAM1/5-63*RATIO)/
X     +       (9-63* (LAM2+LAM1)/5+21*LAM2*LAM1-63*RATIO/LAM0)
X      LAMP = 0.0625
C
C     Compute degree 9 rule weights
C
X      W(1,WTLENG) = 1/ (3*LAM0)**4/TWONDM
X      IF (NDIM.GT.2) W(1,8) = (1-1/ (3*LAM0))/ (6*LAM1)**3
X      W(1,7) = (1-7* (LAM0+LAM1)/5+7*LAM0*LAM1/3)/
X     +         (84*LAM1*LAM2* (LAM2-LAM0)* (LAM2-LAM1))
X      W(1,6) = (1-7* (LAM0+LAM2)/5+7*LAM0*LAM2/3)/
X     +         (84*LAM1*LAM1* (LAM1-LAM0)* (LAM1-LAM2)) -
X     +         W(1,7)*LAM2/LAM1 - 2* (NDIM-2)*W(1,8)
X      W(1,4) = (1-9* ((LAM0+LAM1+LAM2)/7- (LAM0*LAM1+LAM0*LAM2+
X     +         LAM1*LAM2)/5)-3*LAM0*LAM1*LAM2)/
X     +         (18*LAM3* (LAM3-LAM0)* (LAM3-LAM1)* (LAM3-LAM2))
X      W(1,3) = (1-9* ((LAM0+LAM1+LAM3)/7- (LAM0*LAM1+LAM0*LAM3+
X     +         LAM1*LAM3)/5)-3*LAM0*LAM1*LAM3)/
X     +         (18*LAM2* (LAM2-LAM0)* (LAM2-LAM1)* (LAM2-LAM3)) -
X     +         2* (NDIM-1)*W(1,7)
X      W(1,2) = (1-9* ((LAM0+LAM2+LAM3)/7- (LAM0*LAM2+LAM0*LAM3+
X     +         LAM2*LAM3)/5)-3*LAM0*LAM2*LAM3)/
X     +         (18*LAM1* (LAM1-LAM0)* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(1,7)+W(1,6)+ (NDIM-2)*W(1,8))
C
C     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
C
X      W(2,WTLENG) = 1/ (108*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(2,8) = (1-27*TWONDM*W(2,9)*LAM0**3)/ (6*LAM1)**3
X      W(2,7) = (1-5*LAM1/3-15*TWONDM*W(2,WTLENG)*LAM0**2* (LAM0-LAM1))/
X     +          (60*LAM1*LAM2* (LAM2-LAM1))
X      W(2,6) = (1-9* (8*LAM1*LAM2*W(2,7)+TWONDM*W(2,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(2,8)* (NDIM-2)
X      W(2,4) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
X      W(2,3) = (1-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(2,7)
X      W(2,2) = (1-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(2,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(2,7)+W(2,6)+ (NDIM-2)*W(2,8))
X      W(3,WTLENG) = 5/ (324*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(3,8) = (1-27*TWONDM*W(3,9)*LAM0**3)/ (6*LAM1)**3
X      W(3,7) = (1-5*LAM1/3-15*TWONDM*W(3,WTLENG)*LAM0**2* (LAM0-LAM1))/
X     +          (60*LAM1*LAM2* (LAM2-LAM1))
X      W(3,6) = (1-9* (8*LAM1*LAM2*W(3,7)+TWONDM*W(3,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(3,8)* (NDIM-2)
X      W(3,5) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAMP* (LAMP-LAM1)* (LAMP-LAM2))
X      W(3,3) = (1-7* ((LAM1+LAMP)/5-LAM1*LAMP/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAMP)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAMP)) - 2* (NDIM-1)*W(3,7)
X      W(3,2) = (1-7* ((LAM2+LAMP)/5-LAM2*LAMP/3+TWONDM*W(3,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAMP)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAMP)) -
X     +         2* (NDIM-1)* (W(3,7)+W(3,6)+ (NDIM-2)*W(3,8))
X      W(4,WTLENG) = 2/ (81*LAM0**4)/TWONDM
X      IF (NDIM.GT.2) W(4,8) = (2-27*TWONDM*W(4,9)*LAM0**3)/ (6*LAM1)**3
X      W(4,7) = (2-15*LAM1/9-15*TWONDM*W(4,WTLENG)*LAM0* (LAM0-LAM1))/
X     +         (60*LAM1*LAM2* (LAM2-LAM1))
X      W(4,6) = (1-9* (8*LAM1*LAM2*W(4,7)+TWONDM*W(4,WTLENG)*LAM0**2))/
X     +         (36*LAM1*LAM1) - 2*W(4,8)* (NDIM-2)
X      W(4,4) = (2-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
X     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
X      W(4,3) = (2-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
X     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(4,7)
X      W(4,2) = (2-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(4,
X     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
X     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
X     +         2* (NDIM-1)* (W(4,7)+W(4,6)+ (NDIM-2)*W(4,8))
X      W(5,2) = 1/ (6*LAM1)
C
C     Set generator values
C
X      LAM0 = SQRT(LAM0)
X      LAM1 = SQRT(LAM1)
X      LAM2 = SQRT(LAM2)
X      LAM3 = SQRT(LAM3)
X      LAMP = SQRT(LAMP)
X      DO 40 I = 1,NDIM
X          G(I,WTLENG) = LAM0
40    CONTINUE
X      IF (NDIM.GT.2) THEN
X          G(1,8) = LAM1
X          G(2,8) = LAM1
X          G(3,8) = LAM1
X      END IF
X      G(1,7) = LAM1
X      G(2,7) = LAM2
X      G(1,6) = LAM1
X      G(2,6) = LAM1
X      G(1,5) = LAMP
X      G(1,4) = LAM3
X      G(1,3) = LAM2
X      G(1,2) = LAM1
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 9 rule weights and lower degree rule weights.
C
X      W(1,1) = TWONDM
X      DO 70 J = 2,5
X          DO 50 I = 2,WTLENG
X              W(J,I) = W(J,I) - W(1,I)
X              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
X      DO 80 I = 2,WTLENG
X          W(1,I) = TWONDM*W(1,I)
X          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
X      ERRCOF(1) = 5
X      ERRCOF(2) = 5
X      ERRCOF(3) = 1.
X      ERRCOF(4) = 5
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END S09HRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 's09hre.f' &&
  chmod 0644 's09hre.f' ||
  echo 'restore of s09hre.f failed'
  shar_count="`wc -c < 's09hre.f'`"
  test 7827 -eq "$shar_count" ||
    echo "s09hre.f: original size 7827, current size $shar_count"
fi
# ============= s07hre.f ==============
if test -f 's07hre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping s07hre.f (file already exists)'
else
  echo 'x - extracting s07hre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 's07hre.f' &&
X      SUBROUTINE S07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE S07HRE
C***KEYWORDS basic integration rule, degree 7
C***PURPOSE  To initialize a degree 7 basic rule, and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-31
C***DESCRIPTION  S07HRE initializes a degree 7 integration rule,
C            two degree 5 null rules, one degree 3 null rule and one
C            degree 1 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C          WTLENG MUST be set equal to 6.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE S07HRE
C
C   Global variables
C
X      INTEGER NDIM,WTLENG
X      REAL W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
X      REAL RULPTS(WTLENG)
C
C   Local Variables
C
X      REAL RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
X      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT S07HRE
C
C
C     Initialize generators, weights and RULPTS
C
X      DO 30 J = 1,WTLENG
X          DO 10 I = 1,NDIM
X              G(I,J) = 0
10        CONTINUE
X          DO 20 I = 1,5
X              W(I,J) = 0
20        CONTINUE
X          RULPTS(J) = 2*NDIM
30    CONTINUE
X      TWONDM = 2**NDIM
X      RULPTS(WTLENG) = TWONDM
X      RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
X      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
X      LAM0 = 0.4707
X      LAMP = 0.5625
X      LAM1 = 4/ (15-5/LAM0)
X      RATIO = (1-LAM1/LAM0)/27
X      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
C
C     Compute degree 7 rule weights
C
X      W(1,6) = 1/ (3*LAM0)**3/TWONDM
X      W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
X      W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/
X     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
X      W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAM2* (LAM2-LAM1))
C
C     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
C
X      W(2,6) = 1/ (36*LAM0**3)/TWONDM
X      W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
X      W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/
X     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
X      W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAM2* (LAM2-LAM1))
X      W(3,6) = 5/ (108*LAM0**3)/TWONDM
X      W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
X      W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/
X     +         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
X      W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/
X     +         (10*LAMP* (LAMP-LAM1))
X      W(4,6) = 1/ (54*LAM0**3)/TWONDM
X      W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
X      W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/
X     +         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
X      W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/
X     +         (20*LAM2* (LAM2-LAM1))
C
C     Set generator values
C
X      LAM0 = SQRT(LAM0)
X      LAM1 = SQRT(LAM1)
X      LAM2 = SQRT(LAM2)
X      LAMP = SQRT(LAMP)
X      DO 40 I = 1,NDIM
X          G(I,WTLENG) = LAM0
40    CONTINUE
X      G(1,WTLENG-1) = LAM1
X      G(2,WTLENG-1) = LAM1
X      G(1,WTLENG-4) = LAM2
X      G(1,WTLENG-3) = LAM1
X      G(1,WTLENG-2) = LAMP
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 7 rule weights and lower degree rule weights.
C
X      W(1,1) = TWONDM
X      DO 70 J = 2,5
X          DO 50 I = 2,WTLENG
X              W(J,I) = W(J,I) - W(1,I)
X              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
X      DO 80 I = 2,WTLENG
X          W(1,I) = TWONDM*W(1,I)
X          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
X      ERRCOF(1) = 5
X      ERRCOF(2) = 5
X      ERRCOF(3) = 1
X      ERRCOF(4) = 5
X      ERRCOF(5) = 0.5
X      ERRCOF(6) = 0.25
C
C***END S07HRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 's07hre.f' &&
  chmod 0644 's07hre.f' ||
  echo 'restore of s07hre.f failed'
  shar_count="`wc -c < 's07hre.f'`"
  test 4885 -eq "$shar_count" ||
    echo "s07hre.f: original size 4885, current size $shar_count"
fi
# ============= srlhre.f ==============
if test -f 'srlhre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping srlhre.f (file already exists)'
else
  echo 'x - extracting srlhre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'srlhre.f' &&
X      SUBROUTINE SRLHRE(NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,
X     +                  FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT)
C***BEGIN PROLOGUE SRLHRE
C***KEYWORDS basic numerical integration rule
C***PURPOSE  To compute basic integration rule values.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 90-02-06
C***DESCRIPTION SRLHRE computes basic integration rule values for a
C            vector of integrands over a hyper-rectangular region.
C            These are estimates for the integrals. SRLHRE also computes
C            estimates for the errors and determines the coordinate axis
C            where the fourth difference for the integrands is largest.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   WTLENG Integer.
C          The number of weights in the basic integration rule.
C   G      Real array of dimension (NDIM,WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1,J), ..., G(NDIM,J) are the are the generators for the
C          points associated with the Jth weights.
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   ERRCOF Real array of dimension 6.
C          The error coefficients for the rules.
C          It is assumed that the error is computed using:
C           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C             THEN ERROR = ERRCOF(3)*N1
C             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C          where N1-N4 are the null rules, EP is the error
C          for the parent
C          subregion and ES is the error for the sibling subregion.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   SCALES Real array of dimension (3,WTLENG).
C          Scaling factors used to construct new null rules based
C          on a linear combination of two successive null rules
C          in the sequence of null rules.
C   NORMS  Real array of dimension (3,WTLENG).
C          2**NDIM/(1-norm of the null rule constructed by each of the
C          scaling factors.)
C   X      Real Array of dimension NDIM.
C          A work array.
C   NULL   Real array of dimension (NUMFUN, 8)
C          A work array.
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   DIRECT Real.
C          The coordinate axis where the fourth difference of the
C          integrand values is largest.
C
C***REFERENCES
C   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
C   integration over an N-dimensional rectangular region,
C   J.Comp.Appl.Math., 6:295-302, 1980.
C
C   T.O.Espelid, Integration Rules, Null Rules and Error
C   Estimation, Reports in Informatics 33, Dept. of Informatics,
C   Univ. of Bergen, 1988.
C
C***ROUTINES CALLED: SFSHRE, FUNSUB
C
C***END PROLOGUE SRLHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER WTLENG,NUMFUN,NDIM
X      REAL CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),
X     +                 RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),
X     +                 G(NDIM,WTLENG),ERRCOF(6),DIRECT,SCALES(3,WTLENG),
X     +                 NORMS(3,WTLENG)
C
C   Local variables.
C
X      REAL RGNVOL,DIFSUM,DIFMAX,FRTHDF
X      INTEGER I,J,K,DIVAXN
X      REAL SEARCH,RATIO
C
C***FIRST EXECUTABLE STATEMENT SRLHRE
C
C
C       Compute volume of subregion, initialize DIVAXN and rule sums;
C       compute fourth differences and new DIVAXN (RGNERR is used
C       for a work array here). The integrand values used for the
C       fourth divided differences are accumulated in rule arrays.
C
X      RGNVOL = 1
X      DIVAXN = 1
X      DO 10 I = 1,NDIM
X          RGNVOL = RGNVOL*HWIDTH(I)
X          X(I) = CENTER(I)
X          IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN = I
10    CONTINUE
X      CALL FUNSUB(NDIM,X,NUMFUN,RGNERR)
X      DO 30 J = 1,NUMFUN
X          BASVAL(J) = W(1,1)*RGNERR(J)
X          DO 20 K = 1,4
X              NULL(J,K) = W(K+1,1)*RGNERR(J)
20        CONTINUE
30    CONTINUE
X      DIFMAX = 0
X      RATIO = (G(1,3)/G(1,2))**2
X      DO 60 I = 1,NDIM
X          X(I) = CENTER(I) - HWIDTH(I)*G(1,2)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,5))
X          X(I) = CENTER(I) + HWIDTH(I)*G(1,2)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,6))
X          X(I) = CENTER(I) - HWIDTH(I)*G(1,3)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,7))
X          X(I) = CENTER(I) + HWIDTH(I)*G(1,3)
X          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,8))
X          X(I) = CENTER(I)
X          DIFSUM = 0
X          DO 50 J = 1,NUMFUN
X              FRTHDF = 2* (1-RATIO)*RGNERR(J) - (NULL(J,7)+NULL(J,8)) +
X     +                 RATIO* (NULL(J,5)+NULL(J,6))
C
C           Ignore differences below roundoff
C
X              IF (RGNERR(J)+FRTHDF/4.NE.RGNERR(J)) DIFSUM = DIFSUM +
X     +            ABS(FRTHDF)
X              DO 40 K = 1,4
X                  NULL(J,K) = NULL(J,K) + W(K+1,2)*
X     +                        (NULL(J,5)+NULL(J,6)) +
X     +                        W(K+1,3)* (NULL(J,7)+NULL(J,8))
40            CONTINUE
X              BASVAL(J) = BASVAL(J) + W(1,2)* (NULL(J,5)+NULL(J,6)) +
X     +                    W(1,3)* (NULL(J,7)+NULL(J,8))
50        CONTINUE
X          IF (DIFSUM.GT.DIFMAX) THEN
X              DIFMAX = DIFSUM
X              DIVAXN = I
X          END IF
60    CONTINUE
X      DIRECT = DIVAXN
C
C    Finish computing the rule values.
C
X      DO 90 I = 4,WTLENG
X          CALL SFSHRE(NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,
X     +                NULL(1,5))
X          DO 80 J = 1,NUMFUN
X              BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
X              DO 70 K = 1,4
X                  NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
70            CONTINUE
80        CONTINUE
90    CONTINUE
C
C    Compute errors.
C
X      DO 130 J = 1,NUMFUN
C
C    We search for the null rule, in the linear space spanned by two
C    successive null rules in our sequence, which gives the greatest
C    error estimate among all normalized (1-norm) null rules in this
C    space.
C
X          DO 110 I = 1,3
X              SEARCH = 0
X              DO 100 K = 1,WTLENG
X                  SEARCH = MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,
X     +                     K)*NULL(J,I))*NORMS(I,K))
100           CONTINUE
X              NULL(J,I) = SEARCH
110       CONTINUE
X          IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2) .AND.
X     +        ERRCOF(2)*NULL(J,2).LE.NULL(J,3)) THEN
X              RGNERR(J) = ERRCOF(3)*NULL(J,1)
X          ELSE
X              RGNERR(J) = ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
X          END IF
X          RGNERR(J) = RGNVOL*RGNERR(J)
X          BASVAL(J) = RGNVOL*BASVAL(J)
130   CONTINUE
C
C***END SRLHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'srlhre.f' &&
  chmod 0644 'srlhre.f' ||
  echo 'restore of srlhre.f failed'
  shar_count="`wc -c < 'srlhre.f'`"
  test 7902 -eq "$shar_count" ||
    echo "srlhre.f: original size 7902, current size $shar_count"
fi
# ============= sfshre.f ==============
if test -f 'sfshre.f' && test X"$1" != X"-c"; then
  echo 'x - skipping sfshre.f (file already exists)'
else
  echo 'x - extracting sfshre.f (text)'
  sed 's/^X//' << 'SHAR_EOF' > 'sfshre.f' &&
X      SUBROUTINE SFSHRE(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
X     +                  FUNVLS)
C***BEGIN PROLOGUE SFSHRE
C***KEYWORDS fully symmetric sum
C***PURPOSE  To compute fully symmetric basic rule sums
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION SFSHRE computes a fully symmetric sum for a vector
C            of integrand values over a hyper-rectangular region.
C            The sum is fully symmetric with respect to the center of
C            the region and is taken over all sign changes and
C            permutations of the generators for the sum.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   X      Real Array of dimension NDIM.
C          A work array.
C   G      Real Array of dimension NDIM.
C          The generators for the fully symmetric sum. These MUST BE
C          non-negative and non-increasing.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   ON RETURN
C
C   FULSMS Real array of dimension NUMFUN.
C          The values for the fully symmetric sums for each component
C          of the integrand.
C   FUNVLS Real array of dimension NUMFUN.
C          A work array.
C
C***ROUTINES CALLED: FUNSUB
C
C***END PROLOGUE SFSHRE
C
C   Global variables.
C
X      EXTERNAL FUNSUB
X      INTEGER NDIM,NUMFUN
X      REAL CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
X     +                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
C
C   Local variables.
C
X      INTEGER IXCHNG,LXCHNG,I,J,L
X      REAL GL,GI
C
C***FIRST EXECUTABLE STATEMENT SFSHRE
C
X      DO 10 J = 1,NUMFUN
X          FULSMS(J) = 0
10    CONTINUE
C
C     Compute centrally symmetric sum for permutation of G
C
20    DO 30 I = 1,NDIM
X          X(I) = CENTER(I) + G(I)*HWIDTH(I)
30    CONTINUE
40    CALL FUNSUB(NDIM,X,NUMFUN,FUNVLS)
X      DO 50 J = 1,NUMFUN
X          FULSMS(J) = FULSMS(J) + FUNVLS(J)
50    CONTINUE
X      DO 60 I = 1,NDIM
X          G(I) = - G(I)
X          X(I) = CENTER(I) + G(I)*HWIDTH(I)
X          IF (G(I).LT.0) GO TO 40
60    CONTINUE
C
C       Find next distinct permutation of G and loop back for next sum.
C       Permutations are generated in reverse lexicographic order.
C
X      DO 80 I = 2,NDIM
X          IF (G(I-1).GT.G(I)) THEN
X              GI = G(I)
X              IXCHNG = I - 1
X              DO 70 L = 1, (I-1)/2
X                  GL = G(L)
X                  G(L) = G(I-L)
X                  G(I-L) = GL
X                  IF (GL.LE.GI) IXCHNG = IXCHNG - 1
X                  IF (G(L).GT.GI) LXCHNG = L
70            CONTINUE
X              IF (G(IXCHNG).LE.GI) IXCHNG = LXCHNG
X              G(I) = G(IXCHNG)
X              G(IXCHNG) = GI
X              GO TO 20
X          END IF
80    CONTINUE
C
C     Restore original order to generators
C
X      DO 90 I = 1,NDIM/2
X          GI = G(I)
X          G(I) = G(NDIM-I+1)
X          G(NDIM-I+1) = GI
90    CONTINUE
C
C***END SFSHRE
C
X      END
SHAR_EOF
  $shar_touch -am 0531091995 'sfshre.f' &&
  chmod 0644 'sfshre.f' ||
  echo 'restore of sfshre.f failed'
  shar_count="`wc -c < 'sfshre.f'`"
  test 3822 -eq "$shar_count" ||
    echo "sfshre.f: original size 3822, current size $shar_count"
fi
exit 0

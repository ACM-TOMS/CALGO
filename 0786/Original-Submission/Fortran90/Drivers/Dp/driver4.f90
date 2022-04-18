    PROGRAM sample

!                             David M. Smith              6-17-96

!  This is a test program for FMLIB 1.1, a multiple-precision real
!  arithmetic package.  A few example FM calculations are carried
!  out using 60 significant digit precision.

!  The output is saved in file FMSAMPLE.LOG.  A comparison file,
!  FMSAMPLE.CHK, is provided showing the expected output from 32-bit
!  (IEEE arithmetic) machines.  When run on other computers, all the
!  numerical results should still be the same, but the number of terms
!  needed for some of the results might be slightly different.  The
!  program checks all the results and the last line of the log file
!  should be "All results were ok."

!-----------------------------------------------------------------------

!  These four common blocks contain information that must be saved
!  between calls, so they should be declared in the main program.
!  The parameter statement defines various array sizes.

!-----------------------------------------------------------------------

! .. Intrinsic Functions ..
      INTRINSIC mod
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Local Scalars ..
      INTEGER :: iter, j, k, klog, nerror

!             Character string used for input and output.

      CHARACTER (80) :: st1

!             Declare arrays for FM variables.  All are in
!             unpacked format.

! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp, imcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmaddi, fmdiv, fmdivi, fmeq, fmform, fmi2m, &
        fmmpy, fmmpyi, fmset, fmsqr, fmst2m, fmsub, imadd, imdivi, imform, &
        imi2m, immpyi, impmod, impwr, imst2m, imsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Set precision to give at least 60 significant digits
!             and initialize the FMLIB package.

!             Note that any program using the FM package MUST call
!             FMSET before using the package.

      CALL fmset(60)

      nerror = 0

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file FMSAMPLE.LOG.

      klog = 18
      OPEN (klog,file='FMSAMPLE.LOG')

!             1.  Find a root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Use Newton's method with initial guess x = 3.12.
!                 This version is not tuned for speed.  See the FMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function:
!                 f(x) = ((((x-3)*x+1)*x-4)*x+1)*x-6.

!                 MA is the previous iterate.
!                 MB is the current iterate.

      CALL fmst2m('3.12',ma)

!                 Print the first iteration.

      WRITE (kw,90000)
      WRITE (klog,90000)
      CALL fmform('F65.60',ma,st1)
      WRITE (kw,90010) 0, st1(1:65)
      WRITE (klog,90010) 0, st1(1:65)

      DO 10 iter = 1, 10

!                 MC is f(MA).

        CALL fmeq(ma,mc)
        CALL fmaddi(mc,-3)
        CALL fmmpy(mc,ma,mc)
        CALL fmaddi(mc,1)
        CALL fmmpy(mc,ma,mc)
        CALL fmaddi(mc,-4)
        CALL fmmpy(mc,ma,mc)
        CALL fmaddi(mc,1)
        CALL fmmpy(mc,ma,mc)
        CALL fmaddi(mc,-6)

!                 MD is f'(MA).

        CALL fmmpyi(ma,5,md)
        CALL fmaddi(md,-12)
        CALL fmmpy(md,ma,md)
        CALL fmaddi(md,3)
        CALL fmmpy(md,ma,md)
        CALL fmaddi(md,-8)
        CALL fmmpy(md,ma,md)
        CALL fmaddi(md,1)

        CALL fmdiv(mc,md,mb)
        CALL fmsub(ma,mb,mb)

!                 Print each iteration.

        CALL fmform('F65.60',mb,st1)
        WRITE (kw,90010) iter, st1(1:65)
        WRITE (klog,90010) iter, st1(1:65)

!                 Stop iterating if MA and MB agree to over
!                 60 places.

        CALL fmsub(ma,mb,md)
        CALL fmabs(md,md)
        CALL fmst2m('1.0E-61',mc)
        IF (fmcomp(md,'LT',mc)) GO TO 20

!                 Set MA = MB for the next iteration.

        CALL fmeq(mb,ma)
10    CONTINUE

!                 Check the answer.

20    st1 = '3.120656215326726500470956013523797484654623935599066014' // &
        '9888284358'
      CALL fmst2m(st1,mc)
      CALL fmsub(mc,mb,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-61',mc)
      IF (fmcomp(md,'GT',mc)) THEN
        nerror = nerror + 1
        WRITE (kw,90020)
        WRITE (klog,90020)
      END IF

!             2.  Compute the Riemann Zeta function for s=3.

!                 Use Gosper's formula  Zeta(3) =
!                 (5/4)*Sum[ (-1)**k * (k!)**2 / ((k+1)**2 * (2k+1)!) ]
!                 while k = 0, 1, ....

!                 MA is the current partial sum.
!                 MB is the current term.
!                 MC is k!
!                 MD is (2k+1)!

      CALL fmi2m(1,ma)
      CALL fmeq(ma,mc)
      CALL fmeq(ma,md)
      DO 30 k = 1, 200
        CALL fmmpyi(mc,k,mc)
        j = 2*k*(2*k+1)
        CALL fmmpyi(md,j,md)
        CALL fmsqr(mc,mb)
        j = (k+1)*(k+1)
        CALL fmdivi(mb,j,mb)
        CALL fmdiv(mb,md,mb)
        IF (mod(k,2)==0) THEN
          CALL fmadd(ma,mb,ma)
        ELSE
          CALL fmsub(ma,mb,ma)
        END IF

!                 Test for convergence.  KFLAG will be 1 if the result
!                 of the last add or subtract is the same as one of the
!                 input arguments.

        IF (kflag==1) THEN
          WRITE (kw,90030) k
          WRITE (klog,90030) k
          GO TO 40
        END IF
30    CONTINUE

!                 Print the result.

40    CALL fmmpyi(ma,5,ma)
      CALL fmdivi(ma,4,ma)
      CALL fmform('F65.60',ma,st1)
      WRITE (kw,90040) st1(1:65)
      WRITE (klog,90040) st1(1:65)

!                 Check the answer.

      st1 = '1.20205690315959428539973816151144999076498629234049888' // &
        '1792271555'
      CALL fmst2m(st1,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-61',mc)
      IF (fmcomp(md,'GT',mc)) THEN
        nerror = nerror + 1
        WRITE (kw,90050)
        WRITE (klog,90050)
      END IF

!             3.  Integer multiple precision calculations.

!                 Fermat's theorem says  x**(p-1) mod p = 1
!                 when p is prime and x is not a multiple of p.
!                 If  x**(p-1) mod p  gives 1 for some p with
!                 several different x's, then it is very likely
!                 that p is prime (but it is not certain until
!                 further tests are done).

!                 Find a 70-digit number p that is "probably" prime.

!                 MA is the value p being tested.

      CALL imi2m(10,ma)
      CALL imi2m(69,mb)
      CALL impwr(ma,mb,ma)

!                 To speed up the search, test only values that are
!                 not multiples of 2, 3, 5, 7, 11, 13.

      k = 2*3*5*7*11*13
      CALL imdivi(ma,k,ma)
      CALL immpyi(ma,k,ma)
      CALL imi2m(k,mb)
      CALL imadd(ma,mb,ma)
      CALL imi2m(1,md)
      CALL imadd(ma,md,ma)
      CALL imi2m(3,mc)

      DO 50 j = 1, 100

!                 Compute 3**(p-1) mod p

        CALL imsub(ma,md,mb)
        CALL impmod(mc,mb,ma,mc)
        IF (imcomp(mc,'EQ',md)) THEN

!                 Check that 7**(p-1) mod p is also 1.

          CALL imi2m(7,mc)
          CALL impmod(mc,mb,ma,mc)
          IF (imcomp(mc,'EQ',md)) THEN
            WRITE (kw,90060) j
            WRITE (klog,90060) j
            GO TO 60
          END IF
        END IF

        CALL imi2m(3,mc)
        CALL imi2m(k,mb)
        CALL imadd(ma,mb,ma)
50    CONTINUE

!                 Print the result.

60    CALL imform('I72',ma,st1)
      WRITE (kw,90070) st1(1:72)
      WRITE (klog,90070) st1(1:72)

!                 Check the answer.

      st1 = '1000000000000000000000000000000000000000000000000000' // &
        '000000000000659661'
      CALL imst2m(st1,mc)
      IF (imcomp(ma,'NE',mc)) THEN
        nerror = nerror + 1
        WRITE (kw,90080)
        WRITE (klog,90080)
      END IF

      IF (nerror==0) THEN
        WRITE (kw,90090) ' All results were ok.'
        WRITE (klog,90090) ' All results were ok.'
      END IF
      STOP
90000 FORMAT (//' Sample 1.  Find a root of f(x) = x**5 - 3x**4 + ', &
        'x**3 - 4x**2 + x - 6 = 0.'///' Iteration       Newton Approximation')
90010 FORMAT (/I10,4X,A)
90020 FORMAT (/' Error in sample case number 1.'/)
90030 FORMAT (///' Sample 2.',8X,I5,' terms were added'/)
90040 FORMAT (' Zeta(3) = ',A)
90050 FORMAT (/' Error in sample case number 2.'/)
90060 FORMAT (///' Sample 3.',8X,I5,' values were tested'/)
90070 FORMAT (' p = ',A)
90080 FORMAT (/' Error in sample case number 3.'/)
90090 FORMAT (//A/)
    END PROGRAM sample

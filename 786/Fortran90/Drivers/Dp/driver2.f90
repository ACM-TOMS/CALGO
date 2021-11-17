    PROGRAM sample

!                             David M. Smith              9-17-96

!  This is a test program for ZMLIB 1.1, a multiple-precision real
!  arithmetic package.  A few example ZM calculations are carried
!  out using 30 significant digit precision.

!  This program uses both ZMLIB.f90 and FMLIB.f90.

!  The output is saved in file ZMSAMPLE.LOG.  A comparison file,
!  ZMSAMPLE.CHK, is provided showing the expected output from 32-bit
!  (IEEE arithmetic) machines.  When run on other computers, all the
!  numerical results should still be the same, but the number of terms
!  needed for some of the results might be slightly different.  The
!  program checks all the results and the last line of the log file
!  should be "All results were ok."

!-----------------------------------------------------------------------

!  These five common blocks contain information that must be saved
!  between calls, so they should be declared in the main program.
!  The parameter statement defines array sizes and pointers, and
!  contains the FMLIB parameters, followed by ZMLIB parameters.

!-----------------------------------------------------------------------

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Local Scalars ..
      INTEGER :: iter, k, klog, nerror

!             Character string used for input and output.

      CHARACTER (80) :: st1

!             Declare arrays for ZM variables.  All are in
!             unpacked format.
! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmst2m, zmabs, zmadd, zmaddi, zmdiv, zmdivi, zmeq, zmform, &
        zmi2m, zmmpy, zmmpyi, zmset, zmst2m, zmsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, dlogtw, dpeps, &
        dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli, &
         mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, &
        ngrd52, ntrace
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
      COMMON /zmuser/jformz, jprntz
! ..
!             Set precision to give at least 30 significant digits
!             and initialize both the ZMLIB and FMLIB packages.

!             Note that any program using the ZM package MUST call
!             ZMSET before using the package.

      CALL zmset(30)

      nerror = 0

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file ZMSAMPLE.LOG.

      klog = 18
      OPEN (klog,file='ZMSAMPLE.LOG')

!             1.  Find a complex root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Newton's method with initial guess x = .56 + 1.06 i.
!                 This version is not tuned for speed.  See the ZMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function:
!                 f(x) = ((((x-3)*x+1)*x-4)*x+1)*x-6.

!                 MA is the previous iterate.
!                 MB is the current iterate.

      CALL zmst2m('.56 + 1.06 i',ma)

!                 Print the first iteration.

      WRITE (kw,90000)
      WRITE (klog,90000)
      CALL zmform('F32.30','F32.30',ma,st1)
      WRITE (kw,90010) 0, st1(1:69)
      WRITE (klog,90010) 0, st1(1:69)

      DO 10 iter = 1, 10

!                 MC is f(MA).

        CALL zmeq(ma,mc)
        CALL zmaddi(mc,-3)
        CALL zmmpy(mc,ma,mc)
        CALL zmaddi(mc,1)
        CALL zmmpy(mc,ma,mc)
        CALL zmaddi(mc,-4)
        CALL zmmpy(mc,ma,mc)
        CALL zmaddi(mc,1)
        CALL zmmpy(mc,ma,mc)
        CALL zmaddi(mc,-6)

!                 MD is f'(MA).

        CALL zmmpyi(ma,5,md)
        CALL zmaddi(md,-12)
        CALL zmmpy(md,ma,md)
        CALL zmaddi(md,3)
        CALL zmmpy(md,ma,md)
        CALL zmaddi(md,-8)
        CALL zmmpy(md,ma,md)
        CALL zmaddi(md,1)

        CALL zmdiv(mc,md,mb)
        CALL zmsub(ma,mb,mb)

!                 Print each iteration.

        CALL zmform('F32.30','F32.30',mb,st1)
        WRITE (kw,90010) iter, st1(1:69)
        WRITE (klog,90010) iter, st1(1:69)

!                 Stop iterating if MA and MB agree to over
!                 30 places.

        CALL zmsub(ma,mb,md)
        CALL zmabs(md,mafm)

!                 The ABS result is real -- do a real (FM) compare.

        CALL fmst2m('1.0E-31',mbfm)
        IF (fmcomp(mafm,'LT',mbfm)) GO TO 20

!                 Set MA = MB for the next iteration.

        CALL zmeq(mb,ma)
10    CONTINUE

!                 Check the answer.

20    st1 = '0.561958308335403235498111195347453 +' // &
        '1.061134679604332556983391239058885 i'
      CALL zmst2m(st1,mc)
      CALL zmsub(mc,mb,md)
      CALL zmabs(md,mafm)
      CALL fmst2m('1.0E-31',mbfm)
      IF (fmcomp(mafm,'GT',mbfm)) THEN
        nerror = nerror + 1
        WRITE (kw,90020)
        WRITE (klog,90020)
      END IF

!             2.  Compute Exp(1.23-2.34i).

!                 Use the direct Taylor series.  See the ZMEXP routine
!                 for a faster way to get Exp(x).

!                 MA is x.
!                 MB is the current term, x**n/n!.
!                 MC is the current partial sum.

      CALL zmst2m('1.23-2.34i',ma)
      CALL zmi2m(1,mb)
      CALL zmeq(mb,mc)
      DO 30 k = 1, 100
        CALL zmmpy(mb,ma,mb)
        CALL zmdivi(mb,k,mb)
        CALL zmadd(mc,mb,mc)

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

40    CALL zmform('F33.30','F32.30',mc,st1)
      WRITE (kw,90040) st1(1:70)
      WRITE (klog,90040) st1(1:70)

!                 Check the answer.

      st1 = '-2.379681796854777515745457977696745 -' // &
        '2.458032970832342652397461908326042 i'
      CALL zmst2m(st1,md)
      CALL zmsub(md,mc,md)
      CALL zmabs(md,mafm)
      CALL fmst2m('1.0E-31',mbfm)
      IF (fmcomp(mafm,'GT',mbfm)) THEN
        nerror = nerror + 1
        WRITE (kw,90050)
        WRITE (klog,90050)
      END IF

      IF (nerror==0) THEN
        WRITE (kw,90060) ' All results were ok.'
        WRITE (klog,90060) ' All results were ok.'
      END IF
      STOP
90000 FORMAT (//' Sample 1.  Find a root of f(x) = x**5 - 3x**4 + ', &
        'x**3 - 4x**2 + x - 6 = 0.'///' Iteration       Newton Approximation')
90010 FORMAT (/I6,4X,A)
90020 FORMAT (/' Error in sample case number 1.'/)
90030 FORMAT (///' Sample 2.',8X,I5,' terms were added to get ', &
        'Exp(1.23-2.34i)'/)
90040 FORMAT (' Result= ',A)
90050 FORMAT (/' Error in sample case number 2.'/)
90060 FORMAT (//A/)
    END PROGRAM sample

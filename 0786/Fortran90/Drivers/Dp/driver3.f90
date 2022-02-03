    PROGRAM test

!                             David M. Smith               6-14-96

!  This is a test program for FMLIB 1.1, a multiple-precision real
!  arithmetic package.  Most of the FM (floating-point) routines
!  are tested, and the results are checked to 50 significant digits.
!  Most of the IM (integer) routines are tested, with exact results
!  required to pass the tests.

!  This program uses FMLIB.f90.

!  The four common blocks contain information that must be saved
!  between calls, so they should be declared in the main program.
!  The parameter statement defines various array sizes.

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Local Scalars ..
      INTEGER :: klog, ncase, nerror

!             Character strings used for input and output.

      CHARACTER (80) :: st1, st2

!             Declare arrays for FM variables.  All are in
!             unpacked format.

! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmset, test1, test10, test11, test12, test13, test14, test15, &
        test2, test3, test4, test5, test6, test7, test8, test9
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
!             Set precision to give at least 50 significant digits
!             and initialize the FMLIB package.

      CALL fmset(50)

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file TESTFM.LOG.

      klog = 18
      OPEN (klog,file='TESTFM.LOG')

!             NERROR is the number of errors found.
!             NCASE is the number of cases tested.

      nerror = 0

!             Test input and output conversion.

      CALL test1(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test add and subtract.

      CALL test2(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test multiply, divide and square root.

      CALL test3(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test stored constants.

      CALL test4(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test exponentials.

      CALL test5(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test logarithms.

      CALL test6(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test trigonometric functions.

      CALL test7(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test inverse trigonometric functions.

      CALL test8(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test hyperbolic functions.

      CALL test9(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test integer input and output conversion.

      CALL test10(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test integer add and subtract.

      CALL test11(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test integer multiply and divide.

      CALL test12(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test conversions between FM and IM format.

      CALL test13(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test integer power and GCD functions.

      CALL test14(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             Test integer modular functions.

      CALL test15(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!             End of tests.

      IF (nerror==0) THEN
        WRITE (kw,90000) ncase
        WRITE (klog,90000) ncase
      ELSE

!             Write some of the initialized values in common.

        WRITE (klog,*) ' NDIG,MBASE,JFORM1,JFORM2,KRAD = '
        WRITE (klog,*) ndig, mbase, jform1, jform2, krad
        WRITE (klog,*) ' KW,NTRACE,LVLTRC,KFLAG,KWARN,KROUND = '
        WRITE (klog,*) kw, ntrace, lvltrc, kflag, kwarn, kround
        WRITE (klog,*) ' NCALL,MXEXP,MXEXP2,KACCSW,MEXPUN,MEXPOV'
        WRITE (klog,*) ncall, mxexp, mxexp2, kaccsw, mexpun, mexpov
        WRITE (klog,*) ' MUNKNO,IUNKNO,RUNKNO,MXBASE,NDG2MX = '
        WRITE (klog,*) munkno, iunkno, runkno, mxbase, ndg2mx
        WRITE (klog,*) ' MAXINT,INTMAX,SPMAX,DPMAX = '
        WRITE (klog,*) maxint, intmax, spmax, dpmax
        WRITE (klog,*) ' ALOGMB,ALOGM2,ALOGMX,ALOGMT,DLOGMB,DLOGTN ='
        WRITE (klog,*) alogmb, alogm2, alogmx, alogmt, dlogmb, dlogtn
        WRITE (klog,*) ' DLOGTW,DLOGTP,DLOGPI,DPPI ='
        WRITE (klog,*) dlogtw, dlogtp, dlogpi, dppi
        WRITE (klog,*) ' DPEPS,DLOGEB ='
        WRITE (klog,*) dpeps, dlogeb

        WRITE (kw,90010) ncase, nerror
        WRITE (klog,90010) ncase, nerror
      END IF
      WRITE (kw,*) ' End of run.'

      STOP
90000 FORMAT (///1X,I5,' cases tested.  No errors were found.'/)
90010 FORMAT (///1X,I5,' cases tested.',I4,' error(s) found.'/)
    END PROGRAM test
    SUBROUTINE test1(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Input and output testing.

!             Logical function for comparing FM numbers.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmdiv, fmform, fmi2m, fmipwr, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 1
      CALL fmst2m('123',ma)
      CALL fmi2m(123,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmi2m(10,mb)
      CALL fmipwr(mb,-48,mb)

!             Use the .NOT. because FMCOMP returns FALSE for special
!             cases like MD = UNKNOWN, and these should be treated
!             as errors for these tests.

      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 2
      st1 = '1.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmi2m(131,mb)
      CALL fmi2m(97,mc)
      CALL fmdiv(mb,mc,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 3
      st1 = '1.3505154639175257731958762886597938144329896907216495E-2'
      CALL fmst2m(st1,ma)
      CALL fmi2m(131,mb)
      CALL fmi2m(9700,mc)
      CALL fmdiv(mb,mc,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-52',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 4
      st1 = '1.3505154639175257731958762886597938144329896907216495E-2'
      CALL fmst2m(st1,ma)
      CALL fmform('F40.30',ma,st2)
      CALL fmst2m(st2,ma)
      st1 = '         .013505154639175257731958762887'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF (( .NOT. fmcomp(md,'LE',mb)) .OR. st1/=st2) THEN
        CALL errprt('FMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 5
      st1 = '1.3505154639175257731958762886597938144329896907216495E+16'
      CALL fmst2m(st1,ma)
      CALL fmform('F53.33',ma,st2)
      CALL fmst2m(st2,ma)
      st1 = '13505154639175257.731958762886597938144329896907216'
      CALL fmst2m(st1,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 6
      st1 = '1.3505154639175257731958762886597938144329896907216495E+16'
      CALL fmst2m(st1,ma)
      CALL fmform('I24',ma,st2)
      CALL fmst2m(st2,ma)
      st1 = '13505154639175258'
      CALL fmst2m(st1,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 7
      st1 = '-1.3505154639175257731958762886597938144329896907216495E+16'
      CALL fmst2m(st1,ma)
      CALL fmform('E55.49',ma,st2)
      CALL fmst2m(st2,ma)
      st1 = '-1.350515463917525773195876288659793814432989690722D16'
      CALL fmst2m(st1,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 8
      st1 = '-1.3505154639175257731958762886597938144329896907216495E+16'
      CALL fmst2m(st1,ma)
      CALL fmform('1PE54.46',ma,st2)
      CALL fmst2m(st2,ma)
      st1 = '-1.350515463917525773195876288659793814432989691M+16'
      CALL fmst2m(st1,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing input and output routines.')
    END SUBROUTINE test1
    SUBROUTINE test2(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test add and subtract.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmadd, fmaddi, fmi2m, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 9
      CALL fmst2m('123',ma)
      CALL fmst2m('789',mb)
      CALL fmadd(ma,mb,ma)
      CALL fmi2m(912,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMADD ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 10
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmadd(ma,mb,ma)
      st2 = '1.0824742268041237113402061855670103092783505154639175'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMADD ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 11
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmsub(ma,mb,ma)
      st2 = '-.3814432989690721649484536082474226804123711340206185'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSUB ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 12
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.3505154639175257731443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmsub(ma,mb,ma)
      st2 = '5.15463917525773195876288659793815M-20'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSUB ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 13
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmaddi(ma,1)
      st2 = '1.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMADDI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 14
      st1 = '4.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmaddi(ma,5)
      st2 = '9.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMADDI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing add and subtract routines.')
    END SUBROUTINE test2
    SUBROUTINE test3(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test multiply, divide and square root.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmdiv, fmdivi, fmi2m, fmmpy, fmmpyi, fmsqr, &
        fmsqrt, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 15
      CALL fmst2m('123',ma)
      CALL fmst2m('789',mb)
      CALL fmmpy(ma,mb,ma)
      CALL fmi2m(97047,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMMPY ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 16
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmmpy(ma,mb,ma)
      st2 = '0.2565628653416941226485280051014985652035285365075991'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMMPY ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 17
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmdiv(ma,mb,ma)
      st2 = '0.4788732394366197183098591549295774647887323943661972'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMDIV ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 18
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,ma)
      CALL fmmpyi(ma,14,ma)
      st2 = '10.2474226804123711340206185567010309278350515463917526'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMMPYI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 19
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,ma)
      CALL fmdivi(ma,24,ma)
      st2 = '0.0304982817869415807560137457044673539518900343642612'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMDIVI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 20
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsqr(ma,ma)
      st2 = '0.1228610904453183122542246784993091720692953555106813'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSQR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 21
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsqrt(ma,ma)
      st2 = '0.5920434645509785316136003710368759268547372945659987'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSQRT',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing multiply, divide and square root routines.')
    END SUBROUTINE test3
    SUBROUTINE test4(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test stored constants.

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mbsave
      INTEGER :: j, jexp, ndgsav
! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: mlnsv2(0:lunpck), mlnsv3(0:lunpck), mlnsv5(0:lunpck), &
        mlnsv7(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmcons, fmeq, fmexp, fmi2m, fmipwr, fmln, fmpi, &
        fmsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, dlogtw,  &
        dpeps, dppi
      REAL (KIND(0.0D0)) :: mbase, mblogs, mbse, mbslb, mbsli, mbspi, mexpab
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, &
        ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

!             Switch to base 10 and check the stored digits.

      mbsave = mbase
      ndgsav = ndig
      ncase = 22
      mbase = 10
      ndig = 200
      CALL fmcons
      CALL fmi2m(1,mb)
      CALL fmexp(mb,mc)
      DO 10 j = 142, 144
        ndig = j
        ndige = 0
        CALL fmi2m(1,mb)
        CALL fmexp(mb,ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j + 1
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' e    ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 20
        END IF
10    CONTINUE

20    ncase = 23
      mbase = 10
      ndig = 200
      CALL fmi2m(2,mb)
      CALL fmln(mb,mc)
      CALL fmeq(mln1,mlnsv2)
      CALL fmeq(mln2,mlnsv3)
      CALL fmeq(mln3,mlnsv5)
      CALL fmeq(mln4,mlnsv7)
      WRITE (kw,90010)
      DO 30 j = 142, 144
        ndig = j
        ndigli = 0
        CALL fmi2m(2,mb)
        CALL fmln(mb,ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' ln(2)',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 40
        END IF
30    CONTINUE

40    ncase = 24
      mbase = 10
      ndig = 200
      WRITE (kw,90020)
      CALL fmeq(mlnsv3,mc)
      DO 50 j = 142, 144
        ndig = j
        ndigli = 0
        CALL fmi2m(3,mb)
        CALL fmln(mb,ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j + 1
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' ln(3)',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 60
        END IF
50    CONTINUE

60    ncase = 25
      mbase = 10
      ndig = 200
      WRITE (kw,90030)
      CALL fmeq(mlnsv5,mc)
      DO 70 j = 142, 144
        ndig = j
        ndigli = 0
        CALL fmi2m(5,mb)
        CALL fmln(mb,ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j + 1
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' ln(5)',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 80
        END IF
70    CONTINUE

80    ncase = 26
      mbase = 10
      ndig = 200
      WRITE (kw,90040)
      CALL fmeq(mlnsv7,mc)
      DO 90 j = 142, 144
        ndig = j
        ndigli = 0
        CALL fmi2m(7,mb)
        CALL fmln(mb,ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j + 1
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' ln(7)',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 100
        END IF
90    CONTINUE

100   ncase = 27
      mbase = 10
      ndig = 200
      WRITE (kw,90050)
      CALL fmpi(mc)
      DO 110 j = 142, 144
        ndig = j
        ndigpi = 0
        CALL fmpi(ma)
        CALL fmsub(ma,mc,md)
        CALL fmabs(md,md)
        CALL fmi2m(10,mb)
        jexp = -j + 1
        CALL fmipwr(mb,jexp,mb)
        IF ( .NOT. fmcomp(md,'LE',mb)) THEN
          CALL errprt(' pi   ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
          GO TO 120
        END IF
110   CONTINUE

!             Restore base and precision.

120   mbase = mbsave
      ndig = ndgsav
      CALL fmcons
      RETURN
90000 FORMAT (/' Testing stored constants.'//'     Check e.'/)
90010 FORMAT ('     Check ln(2).'/)
90020 FORMAT ('     Check ln(3).'/)
90030 FORMAT ('     Check ln(5).'/)
90040 FORMAT ('     Check ln(7).'/)
90050 FORMAT ('     Check pi.')
    END SUBROUTINE test4
    SUBROUTINE test5(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test exponentials.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmexp, fmipwr, fmpwr, fmrpwr, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 28
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmexp(ma,ma)
      st2 = '0.7043249420381570899426746185150096342459216636010743'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMEXP ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 29
      st1 = '5.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmexp(ma,ma)
      st2 = '210.7168868293979289717186453717687341395104929999527672'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-48',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMEXP ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 30
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmipwr(ma,13,ma)
      st2 = '1.205572620050170403854527299272882946980306577287581E-6'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-56',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMIPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 31
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,ma)
      CALL fmipwr(ma,-1234,ma)
      st2 = '1.673084074011006302103793189789209370839697748745938E167'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E+120',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMIPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 32
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmpwr(ma,mb,ma)
      st2 = '0.4642420045002127676457665673753493595170650613692580'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMPWR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 33
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      st1 = '-34.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,mb)
      CALL fmpwr(ma,mb,ma)
      st2 = '6.504461581246879800523526109766882955934341922848773E15'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-34',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMPWR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 34
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmrpwr(ma,1,3,ma)
      st2 = '0.7050756680967220302067310420367584779561732592049823'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMRPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 35
      st1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL fmst2m(st1,ma)
      CALL fmrpwr(ma,-17,5,ma)
      st2 = '2.8889864895853344043562747681699203201333872009477318'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMRPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing exponential routines.')
    END SUBROUTINE test5
    SUBROUTINE test6(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test logarithms.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmlg10, fmln, fmlni, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 36
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmln(ma,ma)
      st2 = '-1.0483504538872214324499548823726586101452117557127813'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMLN  ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 37
      st1 = '0.3505154639175257731958762886597938144329896907216495E123'
      CALL fmst2m(st1,ma)
      CALL fmln(ma,ma)
      st2 = '282.1696159843803977017629940438041389247902713456262947'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-47',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMLN  ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 38
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmlg10(ma,ma)
      st2 = '-0.4552928172239897280304530226127473926500843247517120'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMLG10',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 39
      CALL fmlni(210,ma)
      st2 = '5.3471075307174686805185894350500696418856767760333836'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMIPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 40
      CALL fmlni(211,ma)
      st2 = '5.3518581334760664957419562654542801180411581735816684'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMPWR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing logarithm routines.')
    END SUBROUTINE test6
    SUBROUTINE test7(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test trigonometric functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmcos, fmcssn, fmsin, fmst2m, fmsub, fmtan
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 41
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcos(ma,ma)
      st2 = '0.9391958366109693586000906984500978377093121163061328'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCOS ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 42
      st1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcos(ma,ma)
      st2 = '0.8069765551968063243992244125871029909816207609700968'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCOS ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 43
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsin(ma,ma)
      st2 = '-0.3433819746180939949443652360333010581867042625893927'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSIN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 44
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsin(ma,ma)
      st2 = '-0.5905834736620182429243173169772978155668602154136946'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSIN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 45
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmtan(ma,ma)
      st2 = '0.3656127521360899712035823015565426347554405301360773'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMTAN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 46
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmtan(ma,ma)
      st2 = '-0.7318471272291003544610122296764031536071117330470298'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMTAN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 47
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcssn(ma,ma,mc)
      st2 = '0.9391958366109693586000906984500978377093121163061328'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 48
      st1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcssn(ma,ma,mc)
      st2 = '0.8069765551968063243992244125871029909816207609700968'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 49
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcssn(ma,mc,ma)
      st2 = '-0.3433819746180939949443652360333010581867042625893927'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 50
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcssn(ma,mc,ma)
      st2 = '-0.5905834736620182429243173169772978155668602154136946'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing trigonometric routines.')
    END SUBROUTINE test7
    SUBROUTINE test8(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test inverse trigonometric functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmacos, fmasin, fmatan, fmst2m, fmsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 51
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmacos(ma,ma)
      st2 = '1.2126748979730954046873545995574544481988102502510807'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMACOS',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 52
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmacos(ma,ma)
      st2 = '1.9289177556166978337752887837220484359983591491240252'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMACOS',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 53
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmasin(ma,ma)
      st2 = '0.3581214288218012145439670920822969938997744494364723'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMASIN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 54
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmasin(ma,ma)
      st2 = '-0.3581214288218012145439670920822969938997744494364723'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMASIN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 55
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmatan(ma,ma)
      st2 = '0.3371339561772373443347761845672381725353758541616570'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMATAN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 56
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmatan(ma,ma)
      st2 = '1.5477326406586162039457549832092678908202994134569781'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMATAN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing inverse trigonometric routines.')
    END SUBROUTINE test8
    SUBROUTINE test9(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test hyperbolic functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmabs, fmchsh, fmcosh, fmsinh, fmst2m, fmsub, fmtanh
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 57
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcosh(ma,ma)
      st2 = '1.0620620786534654254819884264931372964608741056397718'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCOSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 58
      st1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmcosh(ma,ma)
      st2 = '3.356291383454381441662669560464886179346554730604556E+18'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-31',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCOSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 59
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsinh(ma,ma)
      st2 = '-0.3577371366153083355393138079781276622149524420386975'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSINH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 60
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmsinh(ma,ma)
      st2 = '3.356291383454381441662669560464886179197580776059111E+18'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-31',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMSINH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 61
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmtanh(ma,ma)
      st2 = '0.3368326049912874057089491946232983472275659538703038'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMTANH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 62
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmtanh(ma,ma)
      st2 = '0.9999999999999999999999999999999999999556135217341837'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMTANH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 63
      st1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmchsh(ma,ma,mc)
      st2 = '1.0620620786534654254819884264931372964608741056397718'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-49',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 64
      st1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmchsh(ma,ma,mc)
      st2 = '3.356291383454381441662669560464886179346554730604556E+18'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-31',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 65
      st1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmchsh(ma,mc,ma)
      st2 = '-0.3577371366153083355393138079781276622149524420386975'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-50',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 66
      st1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL fmst2m(st1,ma)
      CALL fmchsh(ma,mc,ma)
      st2 = '3.356291383454381441662669560464886179197580776059111E+18'
      CALL fmst2m(st2,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('1.0E-31',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('FMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing hyperbolic routines.')
    END SUBROUTINE test9
    SUBROUTINE test10(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Input and output testing for IM routines.

!             Logical function for comparing IM numbers.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, imform, imi2m, impwr, imst2m
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 67
      CALL imst2m('123',ma)
      CALL imi2m(123,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMST2M',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 68
      st1 = '-350515'
      CALL imst2m(st1,ma)
      CALL imi2m(-350515,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMST2M',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 69
      st1 = '19895113660064588580108197261066338165074766609'
      CALL imst2m(st1,ma)
      CALL imi2m(23,mb)
      CALL imi2m(34,mc)
      CALL impwr(mb,mc,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMST2M',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 70
      st1 = '-20800708073664542533904165663516279809808659679033703'
      CALL imst2m(st1,ma)
      CALL imi2m(-567,mb)
      CALL imi2m(19,mc)
      CALL impwr(mb,mc,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMST2M',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 71
      st1 = '19895113660064588580108197261066338165074766609'
      CALL imst2m(st1,ma)
      CALL imform('I53',ma,st2)
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMFORM',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 72
      st1 = '-20800708073664542533904165663516279809808659679033703'
      CALL imst2m(st1,ma)
      CALL imform('I73',ma,st2)
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMFORM',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing integer input and output routines.')
    END SUBROUTINE test10
    SUBROUTINE test11(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test add and subtract for IM routines.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, imadd, imi2m, imst2m, imsub
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 73
      CALL imst2m('123',ma)
      CALL imst2m('789',mb)
      CALL imadd(ma,mb,ma)
      CALL imi2m(912,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMADD ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 74
      st1 = '3505154639175257731958762886597938144329896907216495'
      CALL imst2m(st1,ma)
      st1 = '7319587628865979381443298969072164948453608247422680'
      CALL imst2m(st1,mb)
      CALL imadd(ma,mb,ma)
      st2 = '10824742268041237113402061855670103092783505154639175'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMADD ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 75
      st1 = '3505154639175257731958762886597938144329896907216495'
      CALL imst2m(st1,ma)
      st1 = '7319587628865979381443298969072164948453608247422680'
      CALL imst2m(st1,mb)
      CALL imsub(ma,mb,ma)
      st2 = '-3814432989690721649484536082474226804123711340206185'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMSUB ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 76
      st1 = '3505154639175257731958762886597938144329896907216495'
      CALL imst2m(st1,ma)
      st1 = '3505154639175257731443298969072164948453608247422680'
      CALL imst2m(st1,mb)
      CALL imsub(ma,mb,ma)
      st2 = '515463917525773195876288659793815'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMSUB ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing integer add and subtract routines.')
    END SUBROUTINE test11
    SUBROUTINE test12(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test integer multiply and divide.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: irem
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, imdiv, imdivi, imdivr, imdvir, imi2m, immod, immpy, &
        immpyi, imsqr, imst2m
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 77
      CALL imst2m('123',ma)
      CALL imst2m('789',mb)
      CALL immpy(ma,mb,ma)
      CALL imi2m(97047,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPY ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 78
      st1 = '10430738374625018354698'
      CALL imst2m(st1,ma)
      st1 = '2879494424799214514791045985'
      CALL imst2m(st1,mb)
      CALL immpy(ma,mb,ma)
      st2 = '30035252996271960952238822892375588336807158787530'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPY ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 79
      CALL imst2m('12347',ma)
      CALL imst2m('47',mb)
      CALL imdiv(ma,mb,ma)
      CALL imst2m('262',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIV ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 80
      st1 = '2701314697583086005158008013691015597308949443159762'
      CALL imst2m(st1,ma)
      st1 = '-978132616472842669976589722394'
      CALL imst2m(st1,mb)
      CALL imdiv(ma,mb,ma)
      CALL imst2m('-2761705981469115610382',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIV ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 81
      CALL imst2m('12368',ma)
      CALL imst2m('67',mb)
      CALL immod(ma,mb,mb)
      CALL imst2m('40',mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMMOD ',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 82
      st1 = '2701314697583086005158008013691015597308949443159762'
      CALL imst2m(st1,ma)
      st1 = '-978132616472842669976589722394'
      CALL imst2m(st1,mb)
      CALL immod(ma,mb,mb)
      CALL imst2m('450750319653685523300198865254',mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMMOD ',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 83
      CALL imst2m('1234',ma)
      CALL imst2m('17',mb)
      CALL imdivr(ma,mb,ma,mb)
      CALL imst2m('72',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIVR',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF
      CALL imst2m('10',mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMDIVR',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 84
      st1 = '34274652243817531418235301715935108945364446765801943'
      CALL imst2m(st1,ma)
      st1 = '-54708769795848731641842224621693'
      CALL imst2m(st1,mb)
      CALL imdivr(ma,mb,ma,mb)
      CALL imst2m('-626492834178447772323',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIVR',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF
      CALL imst2m('31059777254296217822749494999104',mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMDIVR',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 85
      CALL imst2m('4866',ma)
      CALL immpyi(ma,14,ma)
      CALL imst2m('68124',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPYI',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 86
      CALL imst2m('270131469758308600515800801369101559730894',ma)
      CALL immpyi(ma,-2895,ma)
      CALL imst2m('-782030604950303398493243319963549015420938130',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPYI ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 87
      CALL imst2m('-37179',ma)
      CALL imdivi(ma,129,ma)
      CALL imst2m('-288',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIVI',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 88
      st1 = '8267538919383255454483790743961990401918726073065738'
      CALL imst2m(st1,ma)
      CALL imdivi(ma,1729,ma)
      st2 = '4781688212483085861471249707323302719444028960708'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDIVI',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 89
      CALL imst2m('-71792',ma)
      CALL imdvir(ma,65,ma,irem)
      CALL imst2m('-1104',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDVIR',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF
      CALL imi2m(irem,mb)
      CALL imi2m(-32,mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMDVIR',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 90
      st1 = '97813261647284266997658972239417958580120170263408655'
      CALL imst2m(st1,ma)
      CALL imdvir(ma,826,ma,irem)
      st2 = '118417992309060855929369215786220288837917881674828'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMDVIR',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF
      CALL imi2m(irem,mb)
      CALL imi2m(727,mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMDVIR',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 91
      CALL imst2m('538',ma)
      CALL imsqr(ma,ma)
      CALL imst2m('289444',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMSQR ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 92
      CALL imst2m('-47818191879814587168242632',ma)
      CALL imsqr(ma,ma)
      st2 = '2286579474654765721668058416662636606051551222287424'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMSQR ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing integer multiply, divide and square routines.')
    END SUBROUTINE test12
    SUBROUTINE test13(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test conversions between FM and IM format.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp, imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, errprt, fmabs, fmi2m, fmst2m, fmsub, imfm2i, imi2fm, &
        imi2m, imst2m
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 93
      CALL imst2m('123',ma)
      CALL imi2fm(ma,mb)
      CALL fmi2m(123,mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('IMI2FM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 94
      CALL imst2m('979282999076598337488362000995916',ma)
      CALL imi2fm(ma,mb)
      CALL fmst2m('979282999076598337488362000995916',mc)
      CALL fmsub(ma,mc,md)
      CALL fmabs(md,md)
      CALL fmst2m('0',mb)
      IF ( .NOT. fmcomp(md,'LE',mb)) THEN
        CALL errprt('IMI2FM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 95
      CALL fmst2m('123.4',ma)
      CALL imfm2i(ma,mb)
      CALL imi2m(123,mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMFM2I',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 96
      CALL fmst2m('979282999076598337488362000995916',ma)
      CALL imfm2i(ma,mb)
      CALL imst2m('979282999076598337488362000995916',mc)
      IF ( .NOT. imcomp(mb,'EQ',mc)) THEN
        CALL errpr2('IMFM2I',mb,'MB',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing conversions between FM and IM format.')
    END SUBROUTINE test13
    SUBROUTINE test14(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test integer power and GCD functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, imgcd, imi2m, impwr, imst2m
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 97
      CALL imst2m('123',ma)
      CALL imst2m('789',mb)
      CALL imgcd(ma,mb,ma)
      CALL imi2m(3,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMGCD ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 98
      st1 = '431134020618556701030927835051546391752577319587628885'
      CALL imst2m(st1,ma)
      st1 = '900309278350515463917525773195876288659793814432989640'
      CALL imst2m(st1,mb)
      CALL imgcd(ma,mb,ma)
      CALL imst2m('615',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMGCD ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 99
      st1 = '5877631675869176172956662762822298812326084745145447940'
      CALL imst2m(st1,ma)
      st1 = '10379997509886032090765062511740075746391432253007667'
      CALL imst2m(st1,mb)
      CALL imgcd(ma,mb,ma)
      CALL imst2m('1',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMGCD ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 100
      CALL imst2m('47',ma)
      CALL imst2m('34',mb)
      CALL impwr(ma,mb,ma)
      st2 = '710112520079088427392020925014421733344154169313556279969'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPWR ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 101
      CALL imst2m('2',ma)
      CALL imst2m('187',mb)
      CALL impwr(ma,mb,ma)
      st2 = '196159429230833773869868419475239575503198607639501078528'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPWR ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 102
      CALL imst2m('-3',ma)
      CALL imst2m('101',mb)
      CALL impwr(ma,mb,ma)
      st2 = '-1546132562196033993109383389296863818106322566003'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPWR ',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing integer GCD and power routines.')
    END SUBROUTINE test14
    SUBROUTINE test15(ma,mb,mc,md,st1,st2,ncase,nerror,klog)

!  Test integer modular functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (80) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL errpr2, imi2m, immpym, impmod, imst2m
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      WRITE (kw,90000)

      ncase = 103
      CALL imst2m('123',ma)
      CALL imst2m('789',mb)
      CALL imst2m('997',mc)
      CALL immpym(ma,mb,mc,ma)
      CALL imi2m(338,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPYM',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 104
      st1 = '431134020618556701030927835051546391752577319587628885'
      CALL imst2m(st1,ma)
      st1 = '36346366019557973241042306587666640486264616086971724'
      CALL imst2m(st1,mb)
      st1 = '900309278350515463917525773195876288659793814432989640'
      CALL imst2m(st1,mc)
      CALL immpym(ma,mb,mc,ma)
      st2 = '458279704440780378752997531208983184411293504187816380'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPYM',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 105
      st1 = '914726194238000125985765939883182'
      CALL imst2m(st1,ma)
      st1 = '-75505764717193044779376979508186553225192'
      CALL imst2m(st1,mb)
      st1 = '18678872625055834600521936'
      CALL imst2m(st1,mc)
      CALL immpym(ma,mb,mc,ma)
      st2 = '-7769745969769966093344960'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMMPYM',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 106
      CALL imst2m('123',ma)
      CALL imst2m('789',mb)
      CALL imst2m('997',mc)
      CALL impmod(ma,mb,mc,ma)
      CALL imi2m(240,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPMOD',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 107
      st1 = '431134020618556701030927835051546391752577319587628885'
      CALL imst2m(st1,ma)
      st1 = '36346366019557973241042306587666640486264616086971724'
      CALL imst2m(st1,mb)
      st1 = '900309278350515463917525773195876288659793814432989640'
      CALL imst2m(st1,mc)
      CALL impmod(ma,mb,mc,ma)
      st2 = '755107893576299697276281907390144058060594744720442385'
      CALL imst2m(st2,mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPMOD',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      ncase = 108
      CALL imst2m('314159',ma)
      CALL imst2m('1411695892374393248272691827763664225585897550',mb)
      CALL imst2m('1411695892374393248272691827763664225585897551',mc)
      CALL impmod(ma,mb,mc,ma)
      CALL imst2m('1',mc)
      IF ( .NOT. imcomp(ma,'EQ',mc)) THEN
        CALL errpr2('IMPMOD',ma,'MA',mc,'MC',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing integer modular routines.')
    END SUBROUTINE test15
    SUBROUTINE errprt(nrout,m1,name1,m2,name2,m3,name3,ncase,nerror,klog)

!  Print error messages.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using FMST2M.
!  M3 is ABS(M1-M2), and ERRPRT is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in main
!                    correspond to M1,M2,M3.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (2) :: name1, name2, name3
      CHARACTER (6) :: nrout
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: m1(0:lunpck), m2(0:lunpck), m3(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwsave
! ..
! .. External Subroutines ..
      EXTERNAL fmprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      nerror = nerror + 1
      WRITE (kw,90000) ncase, nrout
      WRITE (klog,90000) ncase, nrout

!              Temporarily change KW to KLOG so FMPRNT
!              will write to the log file.

      kwsave = kw
      kw = klog
      WRITE (klog,90010) name1
      CALL fmprnt(m1)
      WRITE (klog,90010) name2
      CALL fmprnt(m2)
      WRITE (klog,90010) name3
      CALL fmprnt(m3)
      kw = kwsave
      RETURN
90000 FORMAT (//' Error in case',I3,'.  The routine',' being tested was ',A6)
90010 FORMAT (1X,A2,' =')
    END SUBROUTINE errprt
    SUBROUTINE errpr2(nrout,m1,name1,m2,name2,ncase,nerror,klog)

!  Print error messages for testing of integer (IM) routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using IMST2M.
!  NAME1,NAME2 are strings identifying which variables in main
!              correspond to M1,M2.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (2) :: name1, name2
      CHARACTER (6) :: nrout
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: m1(0:lunpck), m2(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwsave
! ..
! .. External Subroutines ..
      EXTERNAL imprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      nerror = nerror + 1
      WRITE (kw,90000) ncase, nrout
      WRITE (klog,90000) ncase, nrout

!              Temporarily change KW to KLOG so IMPRNT
!              will write to the log file.

      kwsave = kw
      kw = klog
      WRITE (klog,90010) name1
      CALL imprnt(m1)
      WRITE (klog,90010) name2
      CALL imprnt(m2)
      kw = kwsave
      RETURN
90000 FORMAT (//' Error in case',I3,'.  The routine',' being tested was ',A6)
90010 FORMAT (1X,A2,' =')
    END SUBROUTINE errpr2

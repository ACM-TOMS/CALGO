    PROGRAM test

!                             David M. Smith               6-14-96

!  This is a test program for ZMLIB 1.1, a multiple-precision complex
!  arithmetic package.  Most of the ZM routines are tested, and the
!  results are checked to 50 significant digits.

!  This program uses both ZMLIB.f90 and FMLIB.f90.

!  These five common blocks contain information that must be saved
!  between calls, so they should be declared in the main program.
!  The parameter statement defines array sizes and pointers, and
!  contains the FMLIB parameters, followed by ZMLIB parameters.

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
      INTEGER :: klog, ncase, nerror

!             Character strings used for input and output.

      CHARACTER (160) :: st1, st2

!             Declare arrays for ZM complex variables (MA, MB, MC, MD)
!             and for FM real variables (MAFM, MBFM).  All are in
!             unpacked format.
! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL test1, test2, test3, test4, test5, test6, test7, test8, zmset
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
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
!             Set precision to give at least 50 significant digits
!             and initialize the FMLIB package.

      CALL zmset(50)

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file TESTZM.LOG.

      klog = 18
      OPEN (klog,file='TESTZM.LOG')

!             NERROR is the number of errors found.
!             NCASE is the number of cases tested.

      nerror = 0

!             Test input and output conversion.

      CALL test1(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test add and subtract.

      CALL test2(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test multiply, divide and square root.

      CALL test3(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test exponentials.

      CALL test4(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test logarithms.

      CALL test5(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test trigonometric functions.

      CALL test6(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test inverse trigonometric functions.

      CALL test7(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!             Test hyperbolic functions.

      CALL test8(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

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
    SUBROUTINE test1(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Input and output testing.

!             Logical function for comparing FM numbers.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zm2i2m, zmabs, zmdivi, zmform, zmmpyi, &
        zmst2m, zmsub
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
      CALL zmst2m('123 + 456 i',ma)
      CALL zm2i2m(123,456,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)

!             Use the .NOT. because FMCOMP returns FALSE for special
!             cases like MD = UNKNOWN, and these should be treated
!             as errors for these tests.

      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 2
      st1 = '0.3505154639175257731958762886597938144329896907216495 + ' // &
        '0.7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      CALL zm2i2m(34,71,mc)
      CALL zmdivi(mc,97,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 3
      st1 = '0.3505154639175257731958762886597938144329896907216495E-5 ' // &
        '+ 0.7319587628865979381443298969072164948453608247422680D-5 i'
      CALL zmst2m(st1,ma)
      CALL zm2i2m(34,71,mc)
      CALL zmdivi(mc,9700000,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-55,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 4
      st1 = '7.699115044247787610619469026548672566371681415929204e 03 ' // &
        '- 5.221238938053097345132743362831858407079646017699115M 03 I'
      CALL zmst2m(st1,ma)
      CALL zm2i2m(87,-59,mc)
      CALL zmdivi(mc,113,mc)
      CALL zmmpyi(mc,10000,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-47,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMST2M',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 5
      st1 = '7.699115044247787610619469026548672566371681415929204e+3 ' // &
        '- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL zmst2m(st1,ma)
      CALL zmform('F53.33','F50.30',ma,st2)
      CALL zmst2m(st2,ma)
      st1 = '7699.115044247787610619469026548673 ' // &
        '-5221.238938053097345132743362831858 i'
      CALL zmst2m(st1,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-30,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 6
      st1 = '7.699115044247787610619469026548672566371681415929204e+3 ' // &
        '- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL zmst2m(st1,ma)
      CALL zmform('I9','I7',ma,st2)
      CALL zmst2m(st2,ma)
      st1 = '7699 -5221 i'
      CALL zmst2m(st1,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(0,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 7
      st1 = '7.699115044247787610619469026548672566371681415929204e+3 ' // &
        '- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL zmst2m(st1,ma)
      CALL zmform('E59.50','E58.49',ma,st2)
      CALL zmst2m(st2,ma)
      st1 = '7.6991150442477876106194690265486725663716814159292E3' // &
        '- 5.221238938053097345132743362831858407079646017699E3 i'
      CALL zmst2m(st1,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 8
      st1 = '7.699115044247787610619469026548672566371681415929204e+3 ' // &
        '- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL zmst2m(st1,ma)
      CALL zmform('1PE59.50','1PE58.49',ma,st2)
      CALL zmst2m(st2,ma)
      CALL zmst2m(st1,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-44,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMFORM',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing input and output routines.')
    END SUBROUTINE test1
    SUBROUTINE test2(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test add and subtract.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zm2i2m, zmabs, zmadd, zmst2m, zmsub
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
      CALL zmst2m('123 + 456 i',ma)
      CALL zmst2m('789 - 543 i',mb)
      CALL zmadd(ma,mb,ma)
      CALL zm2i2m(912,-87,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(0,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMADD ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 10
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmadd(ma,mb,ma)
      st2 = '1.1204269683423045342578231913146610710701578323145698 ' // &
        '+ 0.2098348690812882036310555606240306541373962229723565 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMADD ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 11
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmsub(ma,mb,ma)
      st2 = '0.4193960405072529878660706139950734422041784508712709 ' // &
        '- 1.2540826566919076726576042331904023355533254265121795 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSUB ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 12
      st1 = '.7699115044247787610619469026548672566371681415929204E3 ' // &
        '- .5221238938053097345132743362831858407079646017699115E3 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmsub(ma,mb,ma)
      st2 = '769.5609889608612352887510263662074628227351519021987045 ' // &
        '- 522.8558525681963324514186661800930572028099625946537725 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-47,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSUB ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing add and subtract routines.')
    END SUBROUTINE test2
    SUBROUTINE test3(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test multiply, divide and square root.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zm2i2m, zmabs, zmdiv, zmdivi, zmmpy, &
        zmmpyi, zmsqr, zmsqrt, zmst2m, zmsub
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

      ncase = 13
      CALL zmst2m('123 + 456 i',ma)
      CALL zmst2m('789 - 543 i',mb)
      CALL zmmpy(ma,mb,ma)
      CALL zm2i2m(344655,292995,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(0,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMMPY ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 14
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmmpy(ma,mb,ma)
      st2 = '0.6520390475321594745005017790347596022260742632971444 ' // &
        '+ 0.3805309734513274336283185840707964601769911504424779 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMMPY ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 15
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmdiv(ma,mb,ma)
      st2 = '-.1705178497731560089737969128653459210208765017614861 ' // &
        '- 1.1335073636829696356072949942949842987114804337239972 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMDIV ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 16
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmmpyi(ma,36,ma)
      st2 = '27.7168141592920353982300884955752212389380530973451327 ' // &
        '- 18.7964601769911504424778761061946902654867256637168142 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMMPYI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 17
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmdivi(ma,37,ma)
      st2 = '2.080841903850753408275532169337479071992346328629514E-2 ' // &
        '- 1.411145658933269552738579287251853623535039464243004E-2 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-52,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMDIVI',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 18
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsqr(ma,ma)
      st2 = '0.3201503641632077688150990680554467851828647505677813 ' // &
        '- 0.8039783851515388832328295089670295246299631921058814 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSQR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 19
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsqrt(ma,ma)
      st2 = '0.9219999909012323458336720551458583330580388434229845 ' // &
        '- 0.2831474506279259570386845864488094697732718981999941 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSQRT',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing multiply, divide and square root routines.')
    END SUBROUTINE test3
    SUBROUTINE test4(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test exponentials.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zmabs, zmexp, zmipwr, zmpwr, zmrpwr, &
        zmst2m, zmsub
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

      ncase = 20
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmexp(ma,ma)
      st2 = '1.8718374504057787925867989348073888855260008469310002 ' // &
        '- 1.0770279996847678711699041910427261417963102075889234 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMEXP ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 21
      st1 = '5.7699115044247787610619469026548672566371681415929204 ' // &
        '- 4.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmexp(ma,ma)
      st2 = '-60.6144766542152809520229386164396710991242264070603612 ' // &
        '+ 314.7254994809539691403004121118801578835669635535466592 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-47,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMEXP ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 22
      st1 = '1.7699115044247787610619469026548672566371681415929204 ' // &
        '- 1.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmipwr(ma,45,ma)
      st2 = '31595668743300099.70429472191424818167262151605608585179 ' // &
        '- 19209634448276799.67717448173630165852744930837930753788 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-33,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMIPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 23
      st1 = '1.7699115044247787610619469026548672566371681415929204 ' // &
        '- 1.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmipwr(ma,-122,ma)
      st2 = '3.1000215641022021714480000129414241564868699479432E-46 ' // &
        '- 1.1687846789859477815450163510927243367234863123667E-45 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-93,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMIPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 24
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,mb)
      CALL zmpwr(ma,mb,ma)
      st2 = '1.4567089343012352449621841355636496276866203747888724 ' // &
        '- 0.3903177712261966292764255714390622205129978923650749 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMPWR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 25
      st1 = '.3505154639175257731958762886597938144329896907216495 ' // &
        '+ .7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      st1 = '2.7699115044247787610619469026548672566371681415929204 ' // &
        '- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,mb)
      CALL zmpwr(ma,mb,ma)
      st2 = '-1.0053105716678380336247948739245187868180079734997482 ' // &
        '- 0.0819537653234704467729051473979237153087038930127116 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMPWR ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 26
      st1 = '0.7699115044247787610619469026548672566371681415929204 ' // &
        '- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmrpwr(ma,2,7,ma)
      st2 = '0.9653921326136512316639621651337975772631340364271270 ' // &
        '- 0.1659768285667051396562270035411852432430188906482848 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMRPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 27
      st1 = '0.7699115044247787610619469026548672566371681415929204 ' // &
        '- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmrpwr(ma,-19,7,ma)
      st2 = '-0.0567985880053556315170006800325686036902111276420647 ' // &
        '+ 1.2154793972711356706410882510363594270389067962568571 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMRPWR',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing exponential routines.')
    END SUBROUTINE test4
    SUBROUTINE test5(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test logarithms.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zmabs, zmlg10, zmln, zmst2m, zmsub
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
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmln(ma,ma)
      st2 = '-0.0722949652393911311212450699415231782692434885813725 ' // &
        '-  0.5959180055163009910007765127008371205749515965219804 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMLN  ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 29
      st1 = '.7699115044247787610619469026548672566371681415929204E28 ' // &
        '- .5221238938053097345132743362831858407079646017699115E28 i'
      CALL zmst2m(st1,ma)
      CALL zmln(ma,ma)
      st2 = '64.4000876385938880213825156612206746345615981930242708 ' // &
        '-  0.5959180055163009910007765127008371205749515965219804 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMLN  ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 30
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmlg10(ma,ma)
      st2 = '-0.0313973044728549715287589498363619677438302809470943 ' // &
        '-  0.2588039014625211035392823012785304771809982053965284 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMLG10',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 31
      st1 = '.7699115044247787610619469026548672566371681415929204E82 ' // &
        '- .5221238938053097345132743362831858407079646017699115E82 i'
      CALL zmst2m(st1,ma)
      CALL zmlg10(ma,ma)
      st2 = '81.9686026955271450284712410501636380322561697190529057 ' // &
        '-  0.2588039014625211035392823012785304771809982053965284 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMLG10',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing logarithm routines.')
    END SUBROUTINE test5
    SUBROUTINE test6(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test trigonometric functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zmabs, zmcos, zmcssn, zmsin, zmst2m, &
        zmsub, zmtan
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

      ncase = 32
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmcos(ma,ma)
      st2 = '0.8180802525254482451348613286211514555816444253416895 ' // &
        '+  0.3801751200076938035500853542125525088505055292851393 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCOS ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 33
      st1 = '34.7699115044247787610619469026548672566371681415929204 ' // &
        '- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmcos(ma,ma)
      st2 = '-1432925478410268113.5816466154230974355002592549420099 ' // &
        '-  309002816679456015.00151246245263842483282458519462258 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-31,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCOS ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 34
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsin(ma,ma)
      st2 = '0.7931260548991613428648822413402447097755865697557818 ' // &
        '-  0.3921366045897070762848927655743167937790944353110710 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSIN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 35
      st1 = '34.7699115044247787610619469026548672566371681415929204 ' // &
        '- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsin(ma,ma)
      st2 = '-3.090028166794560150015124624526384249047272360765358E17 ' // &
        '+  1.432925478410268113581646615423097435166828182950161E18 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-31,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSIN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 36
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmtan(ma,ma)
      st2 = '0.6141156219447569167198437040270236055089243090199979 ' // &
        '-  0.7647270337230070156308196055474639461102792169274526 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMTAN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 37
      st1 = '35.7699115044247787610619469026548672566371681415929204 ' // &
        '- 43.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmtan(ma,ma)
      st2 = '2.068934241218867332441292427642153175237611151321340E-38 ' // &
        '-  1.000000000000000000000000000000000000023741659169354 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMTAN ',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 38
      st1 = '0.3505154639175257731958762886597938144329896907216495 ' // &
        '+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      CALL zmcssn(ma,ma,mc)
      st2 = '1.2022247452809115256533054407001508718694617802593324 ' // &
        '-  0.2743936538120352873902095801531325075994392065668943 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 39
      st1 = '0.3505154639175257731958762886597938144329896907216495 ' // &
        '+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      CALL zmcssn(ma,mc,ma)
      st2 = '0.4395486978082638069281369170831952476351663772871008 ' // &
        '+  0.7505035100906417134864779281080728222900154610025883 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCSSN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing trigonometric routines.')
    END SUBROUTINE test6
    SUBROUTINE test7(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test inverse trigonometric functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zmabs, zmacos, zmasin, zmatan, zmst2m, &
        zmsub
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

      ncase = 40
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmacos(ma,ma)
      st2 = '0.8797127900868121872960714368309657795959216549012347 ' // &
        '+  0.6342141347945396859119941874681961111936156338608130 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMACOS',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 41
      st1 = '.7699115044247787610619469026548672566371681415929204E12 ' // &
        '- .5221238938053097345132743362831858407079646017699115E12 i'
      CALL zmst2m(st1,ma)
      CALL zmacos(ma,ma)
      st2 = '0.5959180055163009910007767810953294528367807973983794 ' // &
        '+28.2518733312491023865118844008522768856672089946951468 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMACOS',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 42
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmasin(ma,ma)
      st2 = '0.6910835367080844319352502548087856625026630447863182 ' // &
        '-  0.6342141347945396859119941874681961111936156338608130 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMASIN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 43
      st1 = '.7699115044247787610619469026548672566371681415929204E13 ' // &
        '- .5221238938053097345132743362831858407079646017699115E13 i'
      CALL zmst2m(st1,ma)
      CALL zmasin(ma,ma)
      st2 = '0.9748783212785956282305451762549693982010148111568094 ' // &
        '-30.5544584242431480705298759613446206186670533428066404 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-48,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMASIN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 44
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmatan(ma,ma)
      st2 = '0.7417952692265900376512911713942700568648670953521258 ' // &
        '- 0.3162747143126729004878357203292329539837025170484857 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMATAN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 45
      st1 = '.7699115044247787610619469026548672566371681415929204E13 ' // &
        '- .5221238938053097345132743362831858407079646017699115E13 i'
      CALL zmst2m(st1,ma)
      CALL zmatan(ma,ma)
      st2 = '   1.570796326794807650905529836436131532596233124329403 ' // &
        '-6.033484162895927601809954710695221401671437742867605E-14 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMATAN',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing inverse trigonometric routines.')
    END SUBROUTINE test7
    SUBROUTINE test8(ma,mb,mc,md,mafm,mbfm,st1,st2,ncase,nerror,klog)

!  Test hyperbolic functions.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (160) :: st1, st2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mafm(0:lunpck), mb(0:lunpkz), mbfm(0:lunpck), &
        mc(0:lunpkz), md(0:lunpkz)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL errprt, fmi2m, fmipwr, zmabs, zmchsh, zmcosh, zmsinh, zmst2m, &
        zmsub, zmtanh
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

      ncase = 46
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmcosh(ma,ma)
      st2 = '1.1365975275870879962259716562608779977957563621412079 ' // &
        '-  0.4230463404769118342540441830446134405410543954181579 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-49,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCOSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 47
      st1 = '34.7699115044247787610619469026548672566371681415929204 ' // &
        '- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmcosh(ma,ma)
      st2 = '69552104658681.7558589320148420094288419217262200765435 ' // &
        '+ 626163773308016.884007302915197616300902876551542156676 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-35,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCOSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 48
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsinh(ma,ma)
      st2 = '0.7352399228186907963608272785465108877302444847897922 ' // &
        '-  0.6539816592078560369158600079981127012552558121707655 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSINH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 49
      st1 = '34.7699115044247787610619469026548672566371681415929204 ' // &
        '- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmsinh(ma,ma)
      st2 = '6.955210465868175585893201484192181376093291191637290E 13 ' // &
        '+ 6.261637733080168840073029151984050820616907795167046E 14 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-35,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMSINH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 50
      st1 = '.7699115044247787610619469026548672566371681415929204 ' // &
        '- .5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmtanh(ma,ma)
      st2 = '0.7562684782933185240709480231996041186654551038993505 ' // &
        '-  0.2938991498221693198532255749292372853685311106820169 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMTANH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 51
      st1 = '35.7699115044247787610619469026548672566371681415929204 ' // &
        '- 43.5221238938053097345132743362831858407079646017699115 i'
      CALL zmst2m(st1,ma)
      CALL zmtanh(ma,ma)
      st2 = '9.999999999999999999999999999998967653135180689424497E-01 ' // &
        '+ 1.356718776492102400812550018433337461876455254467192E-31 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMTANH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 52
      st1 = '0.3505154639175257731958762886597938144329896907216495 ' // &
        '+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      CALL zmchsh(ma,ma,mc)
      st2 = '0.7900326499280864816444807620997665088044412803737969 ' // &
        '+ 0.2390857359988804105051429301542214823277594407302781 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      ncase = 53
      st1 = '0.3505154639175257731958762886597938144329896907216495 ' // &
        '+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL zmst2m(st1,ma)
      CALL zmchsh(ma,mc,ma)
      st2 = '0.2661087555034471983220879532235334422670297141428191 ' // &
        '+  0.7098057980612199357870532628105009808447460332437714 i'
      CALL zmst2m(st2,mc)
      CALL zmsub(ma,mc,md)
      CALL zmabs(md,mafm)
      CALL fmi2m(10,mbfm)
      CALL fmipwr(mbfm,-50,mbfm)
      IF ( .NOT. fmcomp(mafm,'LE',mbfm)) THEN
        CALL errprt('ZMCHSH',ma,'MA',mc,'MC',md,'MD',ncase,nerror,klog)
      END IF

      RETURN
90000 FORMAT (/' Testing hyperbolic routines.')
    END SUBROUTINE test8
    SUBROUTINE errprt(nrout,m1,name1,m2,name2,m3,name3,ncase,nerror,klog)

!  Print error messages.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using ZMST2M.
!  M3 is ABS(M1-M2), and ERRPRT is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in main
!                    correspond to M1,M2,M3.

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
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
! .. Scalar Arguments ..
      INTEGER :: klog, ncase, nerror
      CHARACTER (2) :: name1, name2, name3
      CHARACTER (6) :: nrout
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: m1(0:lunpkz), m2(0:lunpkz), m3(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: kwsave
! ..
! .. External Subroutines ..
      EXTERNAL zmprnt
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

!              Temporarily change KW to KLOG so ZMPRNT
!              will write to the log file.

      kwsave = kw
      kw = klog
      WRITE (klog,90010) name1
      CALL zmprnt(m1)
      WRITE (klog,90010) name2
      CALL zmprnt(m2)
      WRITE (klog,90010) name3
      CALL zmprnt(m3)
      kw = kwsave
      RETURN
90000 FORMAT (//' Error in case',I3,'.  The routine',' being tested was ',A6)
90010 FORMAT (1X,A2,' =')
    END SUBROUTINE errprt

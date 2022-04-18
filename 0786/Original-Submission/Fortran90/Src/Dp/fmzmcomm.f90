    MODULE fmzmcommon

!             Define the array sizes:

!  Here are the common blocks used in FM and ZM.

!  /FMUSER/, /FM/, /FMSAVE/, /FMBUFF/, and /ZMUSER/ should be declared
!  in the main program.

!             FMUSER contains values that may need to be
!                    changed by the calling program.

!             FM contains the work array used by the low-level
!                arithmetic routines, definitions for overflow
!                and underflow thresholds, and other
!                machine-dependent values.

!             FMSAVE contains information about saved constants.

!             MJSUMS is an array that can contain several FM numbers
!             being used to accumulate concurrent sums in exponential
!             and trigonometric functions.  When NDIGMX = 256, eight is
!             about the maximum number of sums needed (but this depends
!             on MBASE).  For larger NDIGMX, dimensioning MJSUMS to hold
!             more than eight FM numbers could increase the speed of the
!             functions.

!             FMWA contains two work arrays similar to MWA.  They are
!             used in routines FMDIVD, FMMPYD, and FMMPYE.

!             CMBUFF is a character array used by FMPRNT for printing
!                    output from FMOUT.  This array may also be used
!                    for calls to FMOUT from outside the FM package.
!             CMCHAR is the letter used before the exponent field
!                    in FMOUT.  It is defined in FMSET.
!             NAMEST is a stack for names of the routines.  It is
!                    used for trace printing and error messages.

!             FM1 contains scratch arrays for temporary storage of FM
!             numbers while computing various functions.

!             FMPCK contains scratch arrays used to hold input arguments
!             in unpacked format when the packed versions of functions
!             are used.

!             Common blocks used by ZMLIB for complex operations.

!             ZMUSER contains values that may need to be
!                    changed by the calling program.



!             ZM1 contains scratch arrays for temporary storage of ZM
!             numbers while computing various functions.



!             ZMPCK contains scratch arrays used to hold input arguments
!             in unpacked format when the packed versions of functions
!             are used.



!             ZMBUFF contains the complex i/o buffer.

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
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (kind(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, dlogtw, dpeps, &
        dpmax, dppi
      REAL (kind(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli, &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, &
        ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (kind(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa), mwd(lmwa), mwe(lmwa), mz01(0:lunpkz), &
        mz02(0:lunpkz), mz03(0:lunpkz), mz04(0:lunpkz), mzx(0:lunpkz), &
        mzy(0:lunpkz)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /fmwa/mwd, mwe
      COMMON /zm1/mz01, mz02, mz03, mz04
      COMMON /zmbuff/cmbufz
      COMMON /zmpck/mzx, mzy
      COMMON /zmuser/jformz, jprntz
! ..
    END MODULE fmzmcommon

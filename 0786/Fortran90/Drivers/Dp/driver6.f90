    PROGRAM test90

!                             David M. Smith              9-17-96

!    Program using the FM Fortran-90 modules for doing
!    arithmetic using the FM, IM, and ZM derived types.

!    This program does the same calculations as FMSAMPLE and ZMSAMPLE.

!    The output is saved in file SAMPLE90.LOG.  A comparison file,
!    SAMPLE90.CHK, is provided showing the expected output from 32-bit
!    (IEEE arithmetic) machines.  When run on other computers, all the
!    numerical results should still be the same, but the number of terms
!    needed for some of the results might be slightly different.  The
!    program checks all the results and the last line of the log file
!    should be "All results were ok."

!    In a few places, an explicit call is made to an FM or ZM routine.
!    For a call like CALL FM_FORM('F65.60',MAFM,ST1), note that the
!    "FM_" form is used since MAFM is a TYPE (FM) variable and not just
!    an array.  See the discussion in FMZM90.f90.

! .. Use Statements ..
      USE fmzm
! ..
! .. Local Structures ..
      TYPE (fm) :: mafm, mbfm, mcfm, mdfm
      TYPE (im) :: maim, mbim, mcim
      TYPE (zm) :: mazm, mbzm, mczm, mdzm
! ..
! .. Local Scalars ..
      INTEGER :: iter, j, k, klog, nerror

!             Character string used to format multiple-precision output.

      CHARACTER (80) :: st1
! ..
!             Note that any program using the FM package MUST call
!             FM_SET before using the package.

!             Since the argument to FM_SET is not an FM number,
!             calling FMSET would have the same effect.  The "FM_"
!             version is available so that all calls in a program
!             using the derived types can have the "FM_" form.

!             Later in this program complex arithmetic is be used,
!             and ZM_SET is called there to initialize the ZM package.

!             Set precision to give at least 60 significant digits
!             and initialize the FMLIB package.

      CALL fm_set(60)

      nerror = 0

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file SAMPLE90.LOG.

      klog = 18
      OPEN (klog,file='SAMPLE90.LOG')

!             1.  Find a root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Use Newton's method with initial guess x = 3.12.
!                 This version is not tuned for speed.  See the FMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function.

!                 MAFM is the previous iterate.
!                 MBFM is the current iterate.

!             TO_FM is a function for converting other types of numbers
!             to type FM.  Note that TO_FM(3.12) converts the REAL
!             constant to FM, but it is accurate only to single
!             precision.  TO_FM(3.12D0) agrees with 3.12 to double
!             precision accuracy, and TO_FM('3.12') or
!             TO_FM(312)/TO_FM(100) agrees to full FM accuracy.

      mafm = to_fm('3.12')

!                 Print the first iteration.

      WRITE (kw,90000)
      WRITE (klog,90000)
      CALL fm_form('F65.60',mafm,st1)
      WRITE (kw,90010) 0, st1(1:65)
      WRITE (klog,90010) 0, st1(1:65)

      DO iter = 1, 10

!                 MCFM is f(MAFM).

        mcfm = ((((mafm-3)*mafm+1)*mafm-4)*mafm+1)*mafm - 6

!                 MDFM is f'(MAFM).

        mdfm = (((5*mafm-12)*mafm+3)*mafm-8)*mafm + 1

        mbfm = mafm - mcfm/mdfm

!                 Print each iteration.

        CALL fm_form('F65.60',mbfm,st1)
        WRITE (kw,90010) iter, st1(1:65)
        WRITE (klog,90010) iter, st1(1:65)

!                 Stop iterating if MAFM and MBFM agree to over
!                 60 places.

        mdfm = abs(mafm-mbfm)
        IF (mdfm<1.0D-61) EXIT

!                 Set MAFM = MBFM for the next iteration.

        mafm = mbfm
      END DO

!                 Check the answer.

      mcfm = to_fm('3.120656215326726500470956013523797484654623'// &
        '9355990660149888284358')
      IF (abs(mcfm-mbfm)>1.0D-61) THEN
        nerror = nerror + 1
        WRITE (kw,90020)
        WRITE (klog,90020)
      END IF

!             2.  Compute the Riemann Zeta function for s=3.

!                 Use Gosper's formula:  Zeta(3) =
!                 (5/4)*Sum[ (-1)**k * (k!)**2 / ((k+1)**2 * (2k+1)!) ]
!                 while k = 0, 1, ....

!                 MAFM is the current partial sum.
!                 MBFM is the current term.
!                 MCFM is k!
!                 MDFM is (2k+1)!

      mafm = 1
      mcfm = 1
      mdfm = 1
      DO k = 1, 200
        mcfm = k*mcfm
        j = 2*k*(2*k+1)
        mdfm = j*mdfm
        mbfm = mcfm**2
        j = (k+1)*(k+1)
        mbfm = (mbfm/j)/mdfm
        IF (mod(k,2)==0) THEN
          mafm = mafm + mbfm
        ELSE
          mafm = mafm - mbfm
        END IF

!                 Test for convergence.

        IF (mafm-mbfm==mafm) THEN
          WRITE (kw,90030) k
          WRITE (klog,90030) k
          EXIT
        END IF
      END DO

!                 Print the result.

      mafm = (5*mafm)/4
      CALL fm_form('F65.60',mafm,st1)
      WRITE (kw,90040) st1(1:65)
      WRITE (klog,90040) st1(1:65)

!                 Check the answer.

      mcfm = to_fm('1.20205690315959428539973816151144999076498'// &
        '6292340498881792271555')
      IF (abs(mafm-mcfm)>1.0D-61) THEN
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

!                 MAIM is the value p being tested.

      maim = to_im(10)**69

!                 To speed up the search, test only values that are
!                 not multiples of 2, 3, 5, 7, 11, 13.

      k = 2*3*5*7*11*13
      maim = (maim/k)*k + k + 1
      mcim = 3

      DO j = 1, 100

!                 Compute 3**(p-1) mod p

        mbim = maim - 1
        CALL im_pmod(mcim,mbim,maim,mcim)
        IF (mcim==1) THEN

!                 Check that 7**(p-1) mod p is also 1.

          mcim = 7
          CALL im_pmod(mcim,mbim,maim,mcim)
          IF (mcim==1) THEN
            WRITE (kw,90060) j
            WRITE (klog,90060) j
            EXIT
          END IF
        END IF

        mcim = 3
        maim = maim + k
      END DO

!                 Print the result.

      CALL im_form('I72',maim,st1)
      WRITE (kw,90070) st1(1:72)
      WRITE (klog,90070) st1(1:72)

!                 Check the answer.

      mcim = to_im('1000000000000000000000000000000000000000000'// &
        '000000000000000000000659661')
      IF (maim/=mcim) THEN
        nerror = nerror + 1
        WRITE (kw,90080)
        WRITE (klog,90080)
      END IF


!             Complex arithmetic.

!             Set precision to give at least 30 significant digits
!             and initialize the ZMLIB package.  Both FM and ZM
!             operations will now have this precision.

!             Note that any program using the ZM package MUST call
!             ZM_SET before using the package.

      CALL zm_set(30)

!             4.  Find a complex root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Newton's method with initial guess x = .56 + 1.06 i.
!                 This version is not tuned for speed.  See the ZMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function.

!                 MAZM is the previous iterate.
!                 MBZM is the current iterate.

      mazm = to_zm('.56 + 1.06 i')

!                 Print the first iteration.

      WRITE (kw,90090)
      WRITE (klog,90090)
      CALL zm_form('F32.30','F32.30',mazm,st1)
      WRITE (kw,90100) 0, st1(1:69)
      WRITE (klog,90100) 0, st1(1:69)

      DO iter = 1, 10

!                 MCZM is f(MAZM).

        mczm = ((((mazm-3)*mazm+1)*mazm-4)*mazm+1)*mazm - 6

!                 MDZM is f'(MAZM).

        mdzm = (((5*mazm-12)*mazm+3)*mazm-8)*mazm + 1

        mbzm = mazm - mczm/mdzm

!                 Print each iteration.

        CALL zm_form('F32.30','F32.30',mbzm,st1)
        WRITE (kw,90100) iter, st1(1:69)
        WRITE (klog,90100) iter, st1(1:69)

!                 Stop iterating if MAZM and MBZM agree to over
!                 30 places.

        IF (abs(mazm-mbzm)<1.0D-31) EXIT

!                 Set MAZM = MBZM for the next iteration.

        mazm = mbzm
      END DO

!                 Check the answer.

      mczm = to_zm('0.561958308335403235498111195347453 +'// &
        '1.061134679604332556983391239058885 i')
      IF (abs(mczm-mbzm)>1.0D-31) THEN
        nerror = nerror + 1
        WRITE (kw,90110)
        WRITE (klog,90110)
      END IF

!             5.  Compute Exp(1.23-2.34i).

!                 Use the direct Taylor series.  See the ZMEXP routine
!                 for a faster way to get Exp(x).

!                 MAZM is x.
!                 MBZM is the current term, x**n/n!.
!                 MCZM is the current partial sum.

      mazm = to_zm('1.23-2.34i')
      mbzm = 1
      mczm = 1
      DO k = 1, 100
        mbzm = mbzm*mazm/k
        mdzm = mczm + mbzm

!                 Test for convergence.

        IF (mdzm==mczm) THEN
          WRITE (kw,90120) k
          WRITE (klog,90120) k
          EXIT
        END IF
        mczm = mdzm
      END DO

!                 Print the result.

      CALL zm_form('F33.30','F32.30',mczm,st1)
      WRITE (kw,90130) st1(1:70)
      WRITE (klog,90130) st1(1:70)

!                 Check the answer.

      mdzm = to_zm('-2.379681796854777515745457977696745 -'// &
        ' 2.458032970832342652397461908326042 i')
      IF (abs(mdzm-mczm)>1.0D-31) THEN
        nerror = nerror + 1
        WRITE (kw,90140)
        WRITE (klog,90140)
      END IF

      IF (nerror==0) THEN
        WRITE (kw,90150) ' All results were ok.'
        WRITE (klog,90150) ' All results were ok.'
      END IF
90000 FORMAT (//' Sample 1.  Real root of f(x) = x**5 - 3x**4 + ', &
        'x**3 - 4x**2 + x - 6 = 0.'///' Iteration       Newton Approximation')
90010 FORMAT (/I10,4X,A)
90020 FORMAT (/' Error in sample case number 1.'/)
90030 FORMAT (///' Sample 2.',8X,I5,' terms were added'/)
90040 FORMAT (' Zeta(3) = ',A)
90050 FORMAT (/' Error in sample case number 2.'/)
90060 FORMAT (///' Sample 3.',8X,I5,' values were tested'/)
90070 FORMAT (' p = ',A)
90080 FORMAT (/' Error in sample case number 3.'/)
90090 FORMAT (//' Sample 4.  Complex root of f(x) = x**5 - 3x**4 + ', &
        'x**3 - 4x**2 + x - 6 = 0.'///' Iteration       Newton Approximation')
90100 FORMAT (/I6,4X,A)
90110 FORMAT (/' Error in sample case number 4.'/)
90120 FORMAT (///' Sample 5.',8X,I5,' terms were added ', &
        'to get Exp(1.23-2.34i)'/)
90130 FORMAT (' Result= ',A)
90140 FORMAT (/' Error in sample case number 5.'/)
90150 FORMAT (//A/)
    END PROGRAM test90

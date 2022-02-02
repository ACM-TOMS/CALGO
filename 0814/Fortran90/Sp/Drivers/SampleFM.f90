      PROGRAM SAMPLE

!                             David M. Smith

!    This is a sample program using the FM Fortran-90 modules for
!    doing arithmetic using the FM, IM, and ZM derived types.

!    The output is saved in file SAMPLE.LOG.  A comparison file,
!    SAMPLE.CHK, is provided showing the expected output from 32-bit
!    (IEEE arithmetic) machines.  When run on other computers, all the
!    numerical results should still be the same, but the number of terms
!    needed for some of the results might be slightly different.  The
!    program checks all the results and the last line of the log file
!    should be "All results were ok."

!    In a few places, an explicit call is made to an FM or ZM routine.
!    For a call like CALL FM_FORM('F65.60',MAFM,ST1), note that the
!    "FM_" form is used since MAFM is a TYPE (FM) variable and not just
!    an array.  See the discussion in FMZM90.f.


      USE FMZM

      IMPLICIT NONE

      TYPE ( FM ) MAFM,MBFM,MCFM,MDFM
      TYPE ( IM ) MAIM,MBIM,MCIM
      TYPE ( ZM ) MAZM,MBZM,MCZM,MDZM

      CHARACTER(80)  :: ST1
      CHARACTER(175) :: FMT
      INTEGER ITER,J,K,KLOG,LFLAG,NERROR
      INTEGER SEED(7)
      DOUBLE PRECISION VALUE

!             Write output to the screen (unit *), and also to the
!             file SAMPLE.LOG.

      KLOG = 18
      OPEN (KLOG,FILE='SAMPLE.LOG')

      NERROR = 0


!             1.  Find a root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Set precision to give at least 60 significant digits.

      CALL FM_SET(60)

!                 Use Newton's method with initial guess x = 3.12.
!                 This version is not tuned for speed.  See the FMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function.

!                 MAFM is the previous iterate.
!                 MBFM is the current iterate.

!                 TO_FM is a function for converting other types of numbers
!                 to type FM.  Note that TO_FM(3.12) converts the REAL
!                 constant to FM, but it is accurate only to single
!                 precision.  TO_FM(3.12D0) agrees with 3.12 to double
!                 precision accuracy, and TO_FM('3.12') or
!                 TO_FM(312)/TO_FM(100) agrees to full FM accuracy.
!                 Here, TO_FM(3.12) would be ok, since Newton iteration
!                 will correct the error coming from single precision,
!                 but it is a good habit to use the more accurate form.

      MAFM = TO_FM('3.12')

!                 Print the first iteration.

      FMT = "(//' Sample 1.  Real root of f(x) = x**5 - 3x**4 + ',"// &
            "'x**3 - 4x**2 + x - 6 = 0.'///"                       // &
            "' Iteration       Newton Approximation')"
      WRITE (*,FMT)
      WRITE (KLOG,FMT)

!                 FM_FORMAT is a formatting function that returns a
!                 character string (of length 200).
!                 Avoid using FM_FORMAT in the write list, since this
!                 function itself does internal WRITE operations, and
!                 some compilers object to recursive WRITE references.

      ST1 = FM_FORMAT('F65.60',MAFM)
      WRITE (*   ,"(/I10,4X,A)") 0,TRIM(ST1)
      WRITE (KLOG,"(/I10,4X,A)") 0,TRIM(ST1)

      DO ITER = 1, 10

!                 MCFM is f(MAFM).

         MCFM = ((((MAFM-3)*MAFM+1)*MAFM-4)*MAFM+1)*MAFM-6

!                 MDFM is f'(MAFM).

         MDFM = (((5*MAFM-12)*MAFM+3)*MAFM-8)*MAFM+1

         MBFM = MAFM - MCFM/MDFM

!                 Print each iteration.

!                 FM_FORM is a formatting subroutine.  FM_FORM can
!                 handle output strings longer that 200 characters.

         CALL FM_FORM('F65.60',MBFM,ST1)
         WRITE (*   ,"(/I10,4X,A)") ITER,TRIM(ST1)
         WRITE (KLOG,"(/I10,4X,A)") ITER,TRIM(ST1)

!                 Stop iterating if MAFM and MBFM agree to over 60 places.

         MDFM = ABS(MAFM-MBFM)
         IF (MDFM < 1.0D-61) EXIT

!                 Set MAFM = MBFM for the next iteration.

         MAFM = MBFM
      ENDDO

!                 Check the answer.

      MCFM = TO_FM('3.120656215326726500470956013523797484654623'// &
                   '9355990660149888284358')

!                 It is slightly safer to do this test with .NOT. instead of
!                       IF (ABS(MCFM-MBFM) >= 1.0D-61) THEN
!                 because if the result of ABS(MCFM-MBFM) is FM's UNKNOWN value,
!                 the comparison returns false for all comparisons.

      IF (.NOT.(ABS(MCFM-MBFM) < 1.0D-61)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 1.'/)")
          WRITE (KLOG,"(/' Error in sample case number 1.'/)")
      ENDIF


!             2.  Higher Precision.  Compute the root above to 300 decimal places.

      CALL FM_SET(300)

!                 It is tempting to just say MAFM = MCFM here to initialize the
!                 start of the higher precision iterations to be the check value
!                 defined above.  That will not work, because precision has
!                 changed.  Most of the digits of MCFM may be undefined at the
!                 new precision.
!                 The usual way to pad a lower precision value with zeros when
!                 raising precision is to use subroutine FM_EQU, but here it is
!                 easier to define MAFM from scratch at the new precision.

      MAFM = TO_FM('3.120656215326726500470956013523797484654623'// &
                   '9355990660149888284358')

      DO ITER = 1, 10

!                 MCFM is f(MAFM).

         MCFM = ((((MAFM-3)*MAFM+1)*MAFM-4)*MAFM+1)*MAFM-6

!                 MDFM is f'(MAFM).

         MDFM = (((5*MAFM-12)*MAFM+3)*MAFM-8)*MAFM+1

         MBFM = MAFM - MCFM/MDFM

!                 Stop iterating if MAFM and MBFM agree to over 300 places.

         MDFM = ABS(MAFM-MBFM)
         IF (MDFM < TO_FM('1.0E-301')) EXIT

!                 Set MAFM = MBFM for the next iteration.

         MAFM = MBFM
      ENDDO

!                 For very high precision output, it is sometimes more
!                 convenient to use FM_PRNT to format and print the numbers,
!                 since the line breaks are handled automatically.
!                 The unit number for the output, KW, and the format codes
!                 to be used, JFORM1 and JFORM2, are internal FM variables.
!                 Subroutine FMSETVAR is used to re-define these, and the
!                 new values will remain in effect for any further calls
!                 to FM_PRNT.

!                 Other variables that can be changed and the  options they
!                 control are listed in the documentation at the top of file
!                 FM.f.

!                 Set the format to F305.300

      CALL FMSETVAR(' JFORM1 = 2 ')
      CALL FMSETVAR(' JFORM2 = 300 ')

!                 Set the output screen width to 90 columns.

      CALL FMSETVAR(' KSWIDE = 90 ')

      WRITE (*   ,"(///' Sample 2.  Find the root above to 300 decimal places.'/)")
      WRITE (KLOG,"(///' Sample 2.  Find the root above to 300 decimal places.'/)")

!                 Write to the log file.

      CALL FMSETVAR(' KW = 18 ')
      CALL FM_PRNT(MBFM)

!                 Write to the screen (unit 6).

      CALL FMSETVAR(' KW = 6 ')
      CALL FM_PRNT(MBFM)

!                 Check the answer.

      MCFM = TO_FM('3.12065621532672650047095601352379748465462393559906601'// &
                   '4988828435819026499951795468978325745001715109581192343'// &
                   '1332682839420040840535954560118152245371792881305271951'// &
                   '0171189388982124036620583073039835473769132820001100582'// &
                   '7350420283867070989561927541348452154928259189115694520'// &
                   '0789415818387529512010999602155131321076797099026664236')

      IF (.NOT.(ABS(MCFM-MBFM) < TO_FM('1.0E-301'))) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 2.'/)")
          WRITE (KLOG,"(/' Error in sample case number 2.'/)")
      ENDIF


!             3.  Compute the Riemann Zeta function for s=3.

!                 Use Gosper's formula:  Zeta(3) =
!                 (5/4)*Sum[ (-1)**k * (k!)**2 / ((k+1)**2 * (2k+1)!) ]
!                 while k = 0, 1, ....

!                 MAFM is the current partial sum.
!                 MBFM is the current term.
!                 MCFM is k!
!                 MDFM is (2k+1)!

      CALL FM_SET(60)
      MAFM = 1
      MCFM = 1
      MDFM = 1
      DO K = 1, 200
         MCFM = K*MCFM
         J = 2*K*(2*K+1)
         MDFM = J*MDFM
         MBFM = MCFM**2
         J = (K+1)*(K+1)
         MBFM = (MBFM/J)/MDFM
         IF (MOD(K,2) == 0) THEN
             MAFM = MAFM + MBFM
         ELSE
             MAFM = MAFM - MBFM
         ENDIF

!                 Test for convergence.

         IF (MAFM-MBFM == MAFM) THEN
             WRITE (*   ,  &
             "(///' Sample 3.',8X,I5,' terms were added in the Zeta sum'/)") K
             WRITE (KLOG,  &
             "(///' Sample 3.',8X,I5,' terms were added in the Zeta sum'/)") K
             EXIT
         ENDIF
      ENDDO

!                 Print the result.

      MAFM = (5*MAFM)/4
      CALL FM_FORM('F62.60',MAFM,ST1)
      WRITE (*   ,"(' Zeta(3) = ',A)") TRIM(ST1)
      WRITE (KLOG,"(' Zeta(3) = ',A)") TRIM(ST1)

!                 Check the answer.

      MCFM = TO_FM('1.20205690315959428539973816151144999076498'// &
                   '6292340498881792271555')
      IF (.NOT.(ABS(MAFM-MCFM) < 1.0D-61)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 3.'/)")
          WRITE (KLOG,"(/' Error in sample case number 3.'/)")
      ENDIF


!             4.  Integer multiple precision calculations.

!                 Fermat's theorem says  x**(p-1) mod p = 1
!                 when p is prime and x is not a multiple of p.
!                 If  x**(p-1) mod p  gives 1 for some p with
!                 several different x's, then it is very likely
!                 that p is prime (but it is not certain until
!                 further tests are done).

!                 Find a 70-digit number p that is "probably" prime.

!                 Use FM_RANDOM_NUMBER to generate a random 70-digit
!                 starting value and search for a prime from that point.
!                 Initialize the generator.
!                 Note that VALUE is double precision, unlike the similar
!                 Fortran intrinsic random number routine, which returns
!                 a single-precision result.

      CALL FM_SET(80)
      SEED = (/ 2718281,8284590,4523536,0287471,3526624,9775724,7093699 /)
      CALL FM_RANDOM_SEED(PUT=SEED)

!                 MAIM is the value p being tested.

      MAIM = 0
      MCIM = TO_IM(10)**13
      DO J = 1, 6
         CALL FM_RANDOM_NUMBER(VALUE)
         MBIM = 1.0D+13*VALUE
         MAIM = MAIM*MCIM + MBIM
      ENDDO
      MCIM = TO_IM(10)**70
      MAIM = MOD(MAIM,MCIM)

!                 To speed up the search, test only values that are
!                 not multiples of 2, 3, 5, 7, 11, 13.

      K = 2*3*5*7*11*13
      MAIM = (MAIM/K)*K + K + 1
      MCIM = 3

      DO J = 1, 100
         MBIM = MAIM - 1

!                 Compute 3**(p-1) mod p

         MCIM = POWER_MOD(MCIM,MBIM,MAIM)
         IF (MCIM == 1) THEN

!                 Check that 7**(p-1) mod p is also 1.

             MCIM = 7
             MCIM = POWER_MOD(MCIM,MBIM,MAIM)
             IF (MCIM == 1) THEN
                 FMT = "(///' Sample 4.',8X,I5,' values were"//  &
                       " checked before finding a prime p.'/)"
                 WRITE (*   ,FMT) J
                 WRITE (KLOG,FMT) J
                 EXIT
             ENDIF
         ENDIF

         MCIM = 3
         MAIM = MAIM + K
      ENDDO

!                 Print the result.

      CALL IM_FORM('I72',MAIM,ST1)
      WRITE (*   ,"(' p =',A)") TRIM(ST1)
      WRITE (KLOG,"(' p =',A)") TRIM(ST1)

!                 Check the answer.

      MCIM = TO_IM('546831788457201910369201221205379315384'// &
                   '5065543480825746529998049913561')
      IF (.NOT.(MAIM == MCIM)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 4.'/)")
          WRITE (KLOG,"(/' Error in sample case number 4.'/)")
      ENDIF


!             5.  Gamma function.

!                 Check that Gamma(1/2) is Sqrt(pi)

      CALL FM_SET(60)
      WRITE (*   ,"(///' Sample 5.  Check that Gamma(1/2) = Sqrt(pi)'/)")
      WRITE (KLOG,"(///' Sample 5.  Check that Gamma(1/2) = Sqrt(pi)'/)")

      MBFM = GAMMA(TO_FM('0.5'))

!                 Print the result.

      CALL FM_FORM('F62.60',MBFM,ST1)
      WRITE (*   ,"(' Gamma(1/2) = ',A)") TRIM(ST1)
      WRITE (KLOG,"(' Gamma(1/2) = ',A)") TRIM(ST1)

!                 Check the answer.

      MCFM = SQRT(4*ATAN(TO_FM(1)))
      IF (.NOT.(ABS(MCFM-MBFM) < 1.0D-61)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 5.'/)")
          WRITE (KLOG,"(/' Error in sample case number 5.'/)")
      ENDIF


!             6.  Psi and Polygamma functions.

!                 Rational series can often be summed using these functions.
!                 Sum (n=1 to infinity) 1/(n**2 * (8n+1)**2) =
!                 16*(Psi(1) - Psi(9/8)) + Polygamma(1,1) + Polygamma(1,9/8)
!                 Ref: Abramowitz & Stegun, Handbook of Mathematical Functions,
!                      chapter 6, Example 10.

      WRITE (*   ,"(///' Sample 6.  Psi and Polygamma functions.'/)")
      WRITE (KLOG,"(///' Sample 6.  Psi and Polygamma functions.'/)")

      MBFM = 16*(PSI(TO_FM(1)) - PSI(TO_FM(9)/8)) +  &
             POLYGAMMA(1,TO_FM(1)) + POLYGAMMA(1,TO_FM(9)/8)

!                 Print the result.

      CALL FM_FORM('F65.60',MBFM,ST1)
      WRITE (*   ,"(' Sum (n=1 to infinity) 1/(n**2 * (8n+1)**2) = '/9X,A)") TRIM(ST1)
      WRITE (KLOG,"(' Sum (n=1 to infinity) 1/(n**2 * (8n+1)**2) = '/9X,A)") TRIM(ST1)

!                 Check the answer.

      MCFM = TO_FM('1.34994861454130247551078291050351479506449786'//  &
                   '35837270816327396M-2')
      IF (.NOT.(ABS(MCFM-MBFM) < 1.0D-61)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 6.'/)")
          WRITE (KLOG,"(/' Error in sample case number 6.'/)")
      ENDIF


!             7.  Incomplete gamma and Gamma functions.

!                 Find the probability that an observed chi-square for a correct
!                 model should be less that 2.3 when the number of degrees of
!                 freedom is 5.
!                 Ref: Knuth, Volume 2, 3rd ed., Page 56, and Press, Flannery,
!                      Teukolsky, Vetterling, Numerical Recipes, 1st ed., Page 165.

      WRITE (*   ,"(///' Sample 7.  Incomplete gamma and Gamma functions.'/)")
      WRITE (KLOG,"(///' Sample 7.  Incomplete gamma and Gamma functions.'/)")

      MAFM = TO_FM(5)/2
      MBFM = INCOMPLETE_GAMMA1(MAFM,TO_FM('2.3')/2) / GAMMA(MAFM)

!                 Print the result.

      CALL FM_FORM('F61.60',MBFM,ST1)
      WRITE (*   ,"(' Probability = ',A)") TRIM(ST1)
      WRITE (KLOG,"(' Probability = ',A)") TRIM(ST1)

!                 Check the answer.

      MCFM = TO_FM('0.193733130114871446327510259182505999534723186'//  &
                   '07121386973066283739')
      IF (.NOT.(ABS(MCFM-MBFM) < 1.0D-61)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 7.'/)")
          WRITE (KLOG,"(/' Error in sample case number 7.'/)")
      ENDIF


!             Complex arithmetic.

!             Set precision to give at least 30 significant digits.

      CALL FM_SET(30)


!             8.  Find a complex root of the equation
!                 f(x) = x**5 - 3x**4 + x**3 - 4x**2 + x - 6 = 0.

!                 Newton's method with initial guess x = .56 + 1.06 i.
!                 This version is not tuned for speed.  See the ZMSQRT
!                 routine for possible ways to increase speed.
!                 Horner's rule is used to evaluate the function.

!                 MAZM is the previous iterate.
!                 MBZM is the current iterate.

      MAZM = TO_ZM('.56 + 1.06 i')

!                 Print the first iteration.

      FMT = "(///' Sample 8.  Complex root of f(x) = x**5 - 3x**4 + ',"// &
            "'x**3 - 4x**2 + x - 6 = 0.'///"                           // &
            "' Iteration       Newton Approximation')"
      WRITE (*,FMT)
      WRITE (KLOG,FMT)
      CALL ZM_FORM('F32.30','F32.30',MAZM,ST1)
      WRITE (*   ,"(/I6,4X,A)") 0,TRIM(ST1)
      WRITE (KLOG,"(/I6,4X,A)") 0,TRIM(ST1)

      DO ITER = 1, 10

!                 MCZM is f(MAZM).

         MCZM = ((((MAZM-3)*MAZM+1)*MAZM-4)*MAZM+1)*MAZM-6

!                 MDZM is f'(MAZM).

         MDZM = (((5*MAZM-12)*MAZM+3)*MAZM-8)*MAZM+1

         MBZM = MAZM - MCZM/MDZM

!                 Print each iteration.

         CALL ZM_FORM('F32.30','F32.30',MBZM,ST1)
         WRITE (*   ,"(/I6,4X,A)") ITER,TRIM(ST1)
         WRITE (KLOG,"(/I6,4X,A)") ITER,TRIM(ST1)

!                 Stop iterating if MAZM and MBZM agree to over 30 places.

         IF (ABS(MAZM-MBZM) < 1.0D-31) EXIT

!                 Set MAZM = MBZM for the next iteration.

         MAZM = MBZM
      ENDDO

!                 Check the answer.

      MCZM = TO_ZM('0.561958308335403235498111195347453 +'// &
                   '1.061134679604332556983391239058885 i')
      IF (.NOT.(ABS(MCZM-MBZM) < 1.0D-31)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 8.'/)")
          WRITE (KLOG,"(/' Error in sample case number 8.'/)")
      ENDIF


!             9.  Compute Exp(1.23-2.34i).

!                 Use the direct Taylor series.  See the ZMEXP routine
!                 for a faster way to get Exp(x).

!                 MAZM is x.
!                 MBZM is the current term, x**n/n!.
!                 MCZM is the current partial sum.

      MAZM = TO_ZM('1.23-2.34i')
      MBZM = 1
      MCZM = 1
      DO K = 1, 100
         MBZM = MBZM*MAZM/K
         MDZM = MCZM + MBZM

!                 Test for convergence.

         IF (MDZM == MCZM) THEN
             FMT = "(///' Sample 9.',8X,I5,' terms were added ',"// &
                   "'to get Exp(1.23-2.34i)'/)"
             WRITE (*   ,FMT) K
             WRITE (KLOG,FMT) K
             EXIT
         ENDIF
         MCZM = MDZM
      ENDDO

!                 Print the result.

      CALL ZM_FORM('F33.30','F32.30',MCZM,ST1)
      WRITE (*   ,"(' Result= ',A)") TRIM(ST1)
      WRITE (KLOG,"(' Result= ',A)") TRIM(ST1)

!                 Check the answer.

      MDZM = TO_ZM('-2.379681796854777515745457977696745 -'// &
                   ' 2.458032970832342652397461908326042 i')
      IF (.NOT.(ABS(MDZM-MCZM) < 1.0D-31)) THEN
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 9.'/)")
          WRITE (KLOG,"(/' Error in sample case number 9.'/)")
      ENDIF


!            10.  Exception handling.
!                 Iterate (real) Exp(x) starting at 1.0 until overflow occurs.
!
!                 Testing type FM numbers directly using an IF can
!                 be tricky.  When MAFM is +overflow, the statement
!                       IF (MAFM == TO_FM(' +OVERFLOW ')) THEN
!                 will return false, since the comparison routine cannot be
!                 sure that two different overflowed results would have been
!                 equal if the overflow threshold had been higher.
!
!                 In this case, calling subroutine FMFLAG will tell when
!                 an exception has happened.
!
!                 However, for a complicated expression that generates several
!                 FM calls using the derived type numbers, note that the FM
!                 result flag may be zero at the end of the expression even if
!                 an exception occurred.  For example, if EXP(A) overflows in
!                       X = (3 + 1/EXP(A))*2
!                 then the result is 6 with a flag of zero, since the exception
!                 caused no loss of accuracy in the final result.  A warning
!                 message will still appear because of the overflow.
!
!                 The FM warning message is written on unit KW, so in this test
!                 it appears on the screen and not in the log file.
!
!                 The final result is checked by formatting the result and finding
!                 that the output string is '+ OVERFLOW'.

      CALL FM_SET(60)

      MAFM = TO_FM(1)

      FMT = "(///' Sample 10.  Exception handling.'//"                         // &
            "12X,' Iterate Exp(x) starting at 1.0 until overflow occurs.'//"  // &
            "12X,' An FM warning message will be printed before the last iteration.'/)"
      WRITE (*,FMT)
      FMT = "(///' Sample 10.  Exception handling.'//"                         // &
            "12X,' Iterate Exp(x) starting at 1.0 until overflow occurs.'/)"
      WRITE (KLOG,FMT)

      DO J = 1, 10
         MAFM = EXP(MAFM)
         CALL FMFLAG(LFLAG)
         CALL FM_FORM('1PE60.40',MAFM,ST1)
         WRITE (*   ,"(/' Iteration',I3,5X,A)") J,TRIM(ST1)
         WRITE (KLOG,"(/' Iteration',I3,5X,A)") J,TRIM(ST1)
         IF (LFLAG < 0) EXIT
      ENDDO

!             Check that the last result was +overflow.

      IF (FM_FORMAT('E60.40',MAFM) == FM_FORMAT('E60.40',TO_FM('+OVERFLOW'))) THEN
          WRITE (*   ,"(/' Overflow was correctly detected.')")
          WRITE (KLOG,"(/' Overflow was correctly detected.')")
      ELSE
          NERROR = NERROR + 1
          WRITE (*   ,"(/' Error in sample case number 10.'/)")
          WRITE (*   ,"(/' Overflow was not correctly detected.')")
          WRITE (KLOG   ,"(/' Error in sample case number 10.'/)")
          WRITE (KLOG   ,"(/' Overflow was not correctly detected.')")
      ENDIF

      IF (NERROR == 0) THEN
          WRITE (*   ,"(//A/)") ' All results were ok.'
          WRITE (KLOG,"(//A/)") ' All results were ok.'
      ELSE
          WRITE (*   ,"(//I3,A/)") NERROR,' error(s) found.'
          WRITE (KLOG,"(//I3,A/)") NERROR,' error(s) found.'
      ENDIF

      END PROGRAM SAMPLE


!                        David M. Smith

!  This is a test program for FMLIB 1.2, a multiple-precision
!  arithmetic package.  Most of the FM (floating-point real)
!  and ZM (floating-point complex) routines are tested.
!  Precision is set to 50 significant digits and the results
!  are checked to that accuracy.
!  Most of the IM (integer) routines are tested, with exact
!  results required to pass the tests.
!  Most of the USE FMZM derived type interface routines are
!  tested in the same manner as those described above.

!  If all tests are completed successfully, this line is printed:

!  935 cases tested.  No errors were found.

      MODULE TEST_VARS

      USE FMVALS
      USE FMZM

!             Declare arrays for FM variables.

      REAL (KIND(1.0D0)) :: MA(-1:LUNPCK),MB(-1:LUNPCK),MC(-1:LUNPCK),  &
                            MD(-1:LUNPCK),ME(-1:LUNPCK),MP1(-1:LPACK),  &
                            MP2(-1:LPACK),MP3(-1:LPACK)
      REAL (KIND(1.0D0)) :: ZA(-1:LUNPKZ),ZB(-1:LUNPKZ),ZC(-1:LUNPKZ),  &
                            ZD(-1:LUNPKZ),ZE(-1:LUNPKZ)
      REAL (KIND(1.0D0)) :: MLNSV2(-1:LUNPCK),MLNSV3(-1:LUNPCK),  &
                            MLNSV5(-1:LUNPCK),MLNSV7(-1:LUNPCK)

!             Declare derived type variables.

      TYPE (FM), SAVE :: M_A,M_B,M_C,M_D,MFM1,MFM2,MFM3,MFM4,MFM5,MFM6,  &
                         MSMALL,MFMV1(3),MFMV2(3),MFMA(3,3),MFMB(3,3),MFMC(3,3)
      TYPE (IM), SAVE :: M_J,M_K,M_L,MIM1,MIM2,MIM3,MIM4,MIM5,MIMV1(3),  &
                         MIMV2(3),MIMA(2,2),MIMB(2,2),MIMC(2,2)
      TYPE (ZM), SAVE :: M_X,M_Y,M_Z,MZM1,MZM2,MZM3,MZM4,MZM5,MZMV1(3),  &
                         MZMV2(3),MZMA(2,3),MZMB(3,4),MZMC(2,4)

      INTEGER, SAVE :: J1,J2,J3,J4,J5
      REAL, SAVE :: R1,R2,R3,R4,R5,RSMALL
      DOUBLE PRECISION, SAVE :: D1,D2,D3,D4,D5,DSMALL
      COMPLEX, SAVE :: C1,C2,C3,C4,C5
      COMPLEX (KIND(0.0D0)), SAVE :: CD1,CD2,CD3,CD4

      END MODULE TEST_VARS

      PROGRAM TEST

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

!             Character strings used for input and output.

      CHARACTER(80) :: ST1,ST2
      CHARACTER(160) :: STZ1,STZ2

      INTEGER KLOG,KWSAVE,NCASE,NERROR
      REAL TIME1,TIME2

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file TESTFM.LOG.

      KLOG = 18
      OPEN (KLOG,FILE='TESTFM.LOG')
      KWSAVE = KW
      KW = KLOG

!             Set precision to give at least 50 significant digits
!             and initialize the FM package.
!             This call also checks many of the initialization values
!             used in module FMVALS (file FMSAVE.f90).  Set KW = KLOG for
!             this call so that any messages concerning these values will
!             appear in file TESTFM.LOG.

      CALL FMSET(50)
      KW = KWSAVE

      CALL TIMEIT(TIME1)

      J2 = 131
      R2 = 241.21
      D2 = 391.61D0
      C2 = ( 411.11D0 , 421.21D0 )
      CD2 = ( 431.11D0 , 441.21D0 )
      CALL FM_ST2M('581.21',MFM1)
      CALL FM_ST2M('-572.42',MFM2)
      CALL IM_ST2M('661',MIM1)
      CALL IM_ST2M('-602',MIM2)
      CALL ZM_ST2M('731.51 + 711.41 i',MZM1)
      CALL ZM_ST2M('-762.12 - 792.42 i',MZM2)

!             NERROR is the number of errors found.
!             NCASE is the number of cases tested.

      NERROR = 0

!             Test input and output conversion.

      CALL TEST1(ST1,ST2,NCASE,NERROR,KLOG)

!             Test add and subtract.

      CALL TEST2(ST1,ST2,NCASE,NERROR,KLOG)

!             Test multiply, divide and square root.

      CALL TEST3(ST1,ST2,NCASE,NERROR,KLOG)

!             Test stored constants.

      CALL TEST4(NCASE,NERROR,KLOG)

!             Test exponentials.

      CALL TEST5(ST1,ST2,NCASE,NERROR,KLOG)

!             Test logarithms.

      CALL TEST6(ST1,ST2,NCASE,NERROR,KLOG)

!             Test trigonometric functions.

      CALL TEST7(ST1,ST2,NCASE,NERROR,KLOG)

!             Test inverse trigonometric functions.

      CALL TEST8(ST1,ST2,NCASE,NERROR,KLOG)

!             Test hyperbolic functions.

      CALL TEST9(ST1,ST2,NCASE,NERROR,KLOG)

!             Test integer input and output conversion.

      CALL TEST10(ST1,ST2,NCASE,NERROR,KLOG)

!             Test integer add and subtract.

      CALL TEST11(ST1,ST2,NCASE,NERROR,KLOG)

!             Test integer multiply and divide.

      CALL TEST12(ST1,ST2,NCASE,NERROR,KLOG)

!             Test conversions between FM and IM format.

      CALL TEST13(NCASE,NERROR,KLOG)

!             Test integer power and GCD functions.

      CALL TEST14(ST1,ST2,NCASE,NERROR,KLOG)

!             Test integer modular functions.

      CALL TEST15(ST1,ST2,NCASE,NERROR,KLOG)

!             Test complex input and output conversion.

      CALL TEST16(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex add and subtract.

      CALL TEST17(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex multiply, divide and square root.

      CALL TEST18(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex exponentials.

      CALL TEST19(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex logarithms.

      CALL TEST20(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex trigonometric functions.

      CALL TEST21(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex inverse trigonometric functions.

      CALL TEST22(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test complex hyperbolic functions.

      CALL TEST23(STZ1,STZ2,NCASE,NERROR,KLOG)

!             Test the derived type = interface.

      CALL TEST24(NCASE,NERROR,KLOG)

!             Test the derived type == interface.

      CALL TEST25(NCASE,NERROR,KLOG)

!             Test the derived type /= interface.

      CALL TEST26(NCASE,NERROR,KLOG)

!             Test the derived type > interface.

      CALL TEST27(NCASE,NERROR,KLOG)

!             Test the derived type >= interface.

      CALL TEST28(NCASE,NERROR,KLOG)

!             Test the derived type < interface.

      CALL TEST29(NCASE,NERROR,KLOG)

!             Test the derived type <= interface.

      CALL TEST30(NCASE,NERROR,KLOG)

!             Test the derived type + interface.

      CALL TEST31(NCASE,NERROR,KLOG)

!             Test the derived type - interface.

      CALL TEST32(NCASE,NERROR,KLOG)

!             Test the derived type * interface.

      CALL TEST33(NCASE,NERROR,KLOG)

!             Test the derived type / interface.

      CALL TEST34(NCASE,NERROR,KLOG)

!             Test the derived type ** interface.

      CALL TEST35(NCASE,NERROR,KLOG)

!             Test the derived type functions ABS, ..., CEILING interface.

      CALL TEST36(NCASE,NERROR,KLOG)

!             Test the derived type functions CMPLX, ..., EXPONENT interface.

      CALL TEST37(NCASE,NERROR,KLOG)

!             Test the derived type functions FLOOR, ..., MIN interface.

      CALL TEST38(NCASE,NERROR,KLOG)

!             Test the derived type functions MINEXPONENT, ..., RRSPACING interface.

      CALL TEST39(NCASE,NERROR,KLOG)

!             Test the derived type functions SCALE, ..., TINY interface.

      CALL TEST40(NCASE,NERROR,KLOG)

!             Test the derived type functions TO_FM, TO_IM, TO_ZM, ..., TO_DPZ interface.

      CALL TEST41(NCASE,NERROR,KLOG)

!             Test the derived type functions ADDI, ..., Z2M interface.

      CALL TEST42(NCASE,NERROR,KLOG)

!             Test Bernoulli numbers, Pochhammer's function, Euler's constant.

      CALL TEST43(NCASE,NERROR,KLOG)

!             Test Gamma, Factorial, Log(Gamma), Beta, Binomial.

      CALL TEST44(NCASE,NERROR,KLOG)

!             Test Incomplete Gamma, Incomplete Beta.

      CALL TEST45(NCASE,NERROR,KLOG)

!             Test Polygamma, Psi.

      CALL TEST46(NCASE,NERROR,KLOG)

!             Test the different rounding modes.

      CALL TEST47(NCASE,NERROR,KLOG)

!             End of tests.

      CALL TIMEIT(TIME2)

      IF (NERROR == 0) THEN
          WRITE (KW,  &
                 "(///1X,I5,' cases tested.  No errors were found.'/)"  &
          ) NCASE
          WRITE (KLOG,  &
                 "(///1X,I5,' cases tested.  No errors were found.'/)"  &
          ) NCASE
      ELSE IF (NERROR == 1) THEN
          WRITE (KW,  &
                 "(///1X,I5,' cases tested.  1 error was found.'/)"  &
          ) NCASE
          WRITE (KLOG,  &
                 "(///1X,I5,' cases tested.  1 error was found.'/)"  &
          ) NCASE
      ELSE
          WRITE (KW,  &
                 "(///1X,I5,' cases tested.',I4,' errors were found.'/)"  &
          ) NCASE,NERROR
          WRITE (KLOG,  &
                 "(///1X,I5,' cases tested.',I4,' errors were found.'/)"  &
          ) NCASE,NERROR
      ENDIF

      IF (NERROR >= 1) THEN
          KWSAVE = KW
          KW = KLOG

!             Write some of the initialized values in common.

          CALL FMVARS
          KW = KWSAVE
      ENDIF

      WRITE (KW,*) ' '
      WRITE (KW,"(F10.2,A)") TIME2-TIME1,' Seconds for TestFM.'
      WRITE (KW,*) ' '
      WRITE (KLOG,*) ' '
      WRITE (KLOG,"(F10.2,A)") TIME2-TIME1,' Seconds for TestFM.'
      WRITE (KLOG,*) ' '

      WRITE (KW,*)' End of run.'

      STOP
      END PROGRAM TEST

      SUBROUTINE TEST1(ST1,ST2,NCASE,NERROR,KLOG)

!  Input and output testing.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

!             Logical function for comparing FM numbers.

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing input and output routines.')")

      NCASE = 1
      CALL FMST2M('123',MA)
      CALL FMI2M(123,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)

!             Use the .NOT. because FMCOMP returns FALSE for special
!             cases like MD = UNKNOWN, and these should be treated
!             as errors for these tests.

      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMST2M',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 2
      ST1 = '1.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMI2M(131,MB)
      CALL FMI2M(97,MC)
      CALL FMDIV(MB,MC,ME)
      CALL FMEQ(ME,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMST2M',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 3
      ST1 = '1.3505154639175257731958762886597938144329896907216495E-2'
      CALL FMST2M(ST1,MA)
      CALL FMI2M(131,MB)
      CALL FMI2M(9700,MC)
      CALL FMDIV(MB,MC,ME)
      CALL FMEQ(ME,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-52',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMST2M',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 4
      ST1 = '1.3505154639175257731958762886597938144329896907216495E-2'
      CALL FMST2M(ST1,MA)
      CALL FMFORM('F40.30',MA,ST2)
      CALL FMST2M(ST2,MA)
      ST1 = '         .013505154639175257731958762887'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF ((.NOT.FMCOMP(MD,'LE',MB)) .OR. ST1 /= ST2) THEN
          CALL ERRPRTFM('FMFORM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 5
      ST1 = '1.3505154639175257731958762886597938144329896907216495E+16'
      CALL FMST2M(ST1,MA)
      CALL FMFORM('F53.33',MA,ST2)
      CALL FMST2M(ST2,MA)
      ST1 = '13505154639175257.731958762886597938144329896907216'
      CALL FMST2M(ST1,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMFORM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 6
      ST1 = '1.3505154639175257731958762886597938144329896907216495E+16'
      CALL FMST2M(ST1,MA)
      CALL FMFORM('I24',MA,ST2)
      CALL FMST2M(ST2,MA)
      ST1 = '13505154639175258'
      CALL FMST2M(ST1,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMFORM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 7
      ST1 ='-1.3505154639175257731958762886597938144329896907216495E+16'
      CALL FMST2M(ST1,MA)
      CALL FMFORM('E55.49',MA,ST2)
      CALL FMST2M(ST2,MA)
      ST1 = '-1.350515463917525773195876288659793814432989690722D16'
      CALL FMST2M(ST1,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMFORM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 8
      ST1 ='-1.3505154639175257731958762886597938144329896907216495E+16'
      CALL FMST2M(ST1,MA)
      CALL FMFORM('1PE54.46',MA,ST2)
      CALL FMST2M(ST2,MA)
      ST1 = '-1.350515463917525773195876288659793814432989691M+16'
      CALL FMST2M(ST1,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMFORM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST1

      SUBROUTINE TEST2(ST1,ST2,NCASE,NERROR,KLOG)

!  Test add and subtract.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing add and subtract routines.')")

      NCASE = 9
      CALL FMST2M('123',MA)
      CALL FMST2M('789',MB)
      CALL FMADD(MA,MB,ME)
      CALL FMEQ(ME,MA)
      CALL FMI2M(912,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMADD ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 10
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMADD(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.0824742268041237113402061855670103092783505154639175'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMADD ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 11
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMSUB(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-.3814432989690721649484536082474226804123711340206185'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSUB ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 12
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.3505154639175257731443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMSUB(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '5.15463917525773195876288659793815M-20'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSUB ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 13
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMADDI(MA,1)
      ST2 = '1.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMADDI',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 14
      ST1 = '4.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMADDI(MA,5)
      ST2 = '9.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMADDI',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST2

      SUBROUTINE TEST3(ST1,ST2,NCASE,NERROR,KLOG)

!  Test multiply, divide and square root.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing multiply, divide and square root routines.')")

      NCASE = 15
      CALL FMST2M('123',MA)
      CALL FMST2M('789',MB)
      CALL FMMPY(MA,MB,ME)
      CALL FMEQ(ME,MA)
      CALL FMI2M(97047,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMMPY ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 16
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMMPY(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.2565628653416941226485280051014985652035285365075991'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMMPY ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 17
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMDIV(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.4788732394366197183098591549295774647887323943661972'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMDIV ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 18
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MA)
      CALL FMMPYI(MA,14,ME)
      CALL FMEQ(ME,MA)
      ST2 = '10.2474226804123711340206185567010309278350515463917526'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMMPYI',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 19
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MA)
      CALL FMDIVI(MA,24,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.0304982817869415807560137457044673539518900343642612'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMDIVI',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 20
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSQR(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.1228610904453183122542246784993091720692953555106813'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSQR ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 21
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSQRT(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.5920434645509785316136003710368759268547372945659987'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSQRT',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST3

      SUBROUTINE TEST4(NCASE,NERROR,KLOG)

!  Test stored constants.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      REAL (KIND(1.0D0)) :: MBSAVE
      INTEGER J,JEXP,KLOG,NCASE,NDGSAV,NERROR

      WRITE (KW,"(/' Testing stored constants.'//'     Check e.'/)")

!             Switch to base 10 and check the stored digits.

      IF (NDIGMX < 55) THEN
          WRITE (KLOG,*) ' '
          WRITE (KLOG,*)  &
                  ' To test these constants at their stored precision requires'
          WRITE (KLOG,*)  &
                  ' setting NDIG=55 (number of digits).  The current maximum'
          WRITE (KLOG,*) ' for NDIG is NDIGMX = ',NDIGMX
          WRITE (KLOG,*) ' Skip the tests for stored constants.'
          RETURN
      ENDIF

      MBSAVE = MBASE
      NDGSAV = NDIG
      NCASE = 22
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      CALL FMCONS
      CALL FMI2M(1,MB)
      CALL FMEXP(MB,MC)
      DO J = 49, 51
         NDIG = J
         NDIGE = 0
         CALL FMI2M(1,MB)
         CALL FMEXP(MB,MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J + 1
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' e    ',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

      NCASE = 23
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      CALL FMI2M(2,MB)
      CALL FMLN(MB,MC)
      CALL FMEQ(MLN1,MLNSV2)
      CALL FMEQ(MLN2,MLNSV3)
      CALL FMEQ(MLN3,MLNSV5)
      CALL FMEQ(MLN4,MLNSV7)
      WRITE (KW,"('     Check ln(2).'/)")
      DO J = 49, 51
         NDIG = J
         NDIGLI = 0
         CALL FMI2M(2,MB)
         CALL FMLN(MB,MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' ln(2)',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

      NCASE = 24
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      WRITE (KW,"('     Check ln(3).'/)")
      CALL FMEQ(MLNSV3,MC)
      DO J = 49, 51
         NDIG = J
         NDIGLI = 0
         CALL FMI2M(3,MB)
         CALL FMLN(MB,MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J + 1
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' ln(3)',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

      NCASE = 25
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      WRITE (KW,"('     Check ln(5).'/)")
      CALL FMEQ(MLNSV5,MC)
      DO J = 49, 51
         NDIG = J
         NDIGLI = 0
         CALL FMI2M(5,MB)
         CALL FMLN(MB,MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J + 1
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' ln(5)',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

      NCASE = 26
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      WRITE (KW,"('     Check ln(7).'/)")
      CALL FMEQ(MLNSV7,MC)
      DO J = 49, 51
         NDIG = J
         NDIGLI = 0
         CALL FMI2M(7,MB)
         CALL FMLN(MB,MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J + 1
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' ln(7)',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

      NCASE = 27
      CALL FMSETVAR(' MBASE = 1000 ')
      CALL FMSETVAR(' NDIG = 55 ')
      WRITE (KW,"('     Check pi.')")
      CALL FMPI(MC)
      DO J = 49, 51
         NDIG = J
         NDIGPI = 0
         CALL FMPI(MA)
         CALL FMSUB(MA,MC,MD)
         CALL FMABS(MD,ME)
         CALL FMEQ(ME,MD)
         CALL FMI2M(1000,MB)
         JEXP = -J + 1
         CALL FMIPWR(MB,JEXP,ME)
         CALL FMEQ(ME,MB)
         IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
             CALL ERRPRTFM(' pi   ',MA,'MA',MC,'MC',MD,'MD',  &
                            NCASE,NERROR,KLOG)
             EXIT
         ENDIF
      ENDDO

!             Restore base and precision.

      MBASE = MBSAVE
      NDIG = NDGSAV
      CALL FMCONS
      RETURN
      END SUBROUTINE TEST4

      SUBROUTINE TEST5(ST1,ST2,NCASE,NERROR,KLOG)

!  Test exponentials.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing exponential routines.')")

      NCASE = 28
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMEXP(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.7043249420381570899426746185150096342459216636010743'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMEXP ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 29
      ST1 = '5.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMEXP(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '210.7168868293979289717186453717687341395104929999527672'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-48',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMEXP ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 30
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMIPWR(MA,13,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.205572620050170403854527299272882946980306577287581E-6'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-56',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMIPWR',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 31
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MA)
      CALL FMIPWR(MA,-1234,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.673084074011006302103793189789209370839697748745938E167'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E+120',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMIPWR',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 32
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMPWR(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.4642420045002127676457665673753493595170650613692580'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMPWR ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 33
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      ST1 = '-34.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MB)
      CALL FMPWR(MA,MB,ME)
      CALL FMEQ(ME,MA)
      ST2 = '6.504461581246879800523526109766882955934341922848773E15'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-34',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMPWR ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 34
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMRPWR(MA,1,3,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.7050756680967220302067310420367584779561732592049823'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMRPWR',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 35
      ST1 = '0.7319587628865979381443298969072164948453608247422680'
      CALL FMST2M(ST1,MA)
      CALL FMRPWR(MA,-17,5,ME)
      CALL FMEQ(ME,MA)
      ST2 = '2.8889864895853344043562747681699203201333872009477318'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMRPWR',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST5

      SUBROUTINE TEST6(ST1,ST2,NCASE,NERROR,KLOG)

!  Test logarithms.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing logarithm routines.')")

      NCASE = 36
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMLN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-1.0483504538872214324499548823726586101452117557127813'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMLN  ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 37
      ST1 = '0.3505154639175257731958762886597938144329896907216495E123'
      CALL FMST2M(ST1,MA)
      CALL FMLN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '282.1696159843803977017629940438041389247902713456262947'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-47',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMLN  ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 38
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMLG10(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.4552928172239897280304530226127473926500843247517120'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMLG10',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 39
      CALL FMLNI(210,MA)
      ST2 = '5.3471075307174686805185894350500696418856767760333836'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMIPWR',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 40
      CALL FMLNI(211,MA)
      ST2 = '5.3518581334760664957419562654542801180411581735816684'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMPWR ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST6

      SUBROUTINE TEST7(ST1,ST2,NCASE,NERROR,KLOG)

!  Test trigonometric functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing trigonometric routines.')")

      NCASE = 41
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCOS(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.9391958366109693586000906984500978377093121163061328'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCOS ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 42
      ST1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCOS(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.8069765551968063243992244125871029909816207609700968'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCOS ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 43
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSIN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.3433819746180939949443652360333010581867042625893927'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSIN ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 44
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSIN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.5905834736620182429243173169772978155668602154136946'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSIN ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 45
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMTAN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.3656127521360899712035823015565426347554405301360773'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMTAN ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 46
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMTAN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.7318471272291003544610122296764031536071117330470298'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMTAN ',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 47
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCSSN(MA,ME,MC)
      CALL FMEQ(ME,MA)
      ST2 = '0.9391958366109693586000906984500978377093121163061328'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCSSN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 48
      ST1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCSSN(MA,ME,MC)
      CALL FMEQ(ME,MA)
      ST2 = '0.8069765551968063243992244125871029909816207609700968'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCSSN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 49
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCSSN(MA,MC,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.3433819746180939949443652360333010581867042625893927'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCSSN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 50
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCSSN(MA,MC,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.5905834736620182429243173169772978155668602154136946'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCSSN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST7

      SUBROUTINE TEST8(ST1,ST2,NCASE,NERROR,KLOG)

!  Test inverse trigonometric functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing inverse trigonometric routines.')")

      NCASE = 51
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMACOS(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.2126748979730954046873545995574544481988102502510807'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMACOS',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 52
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMACOS(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.9289177556166978337752887837220484359983591491240252'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMACOS',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 53
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMASIN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.3581214288218012145439670920822969938997744494364723'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMASIN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 54
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMASIN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.3581214288218012145439670920822969938997744494364723'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMASIN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 55
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMATAN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.3371339561772373443347761845672381725353758541616570'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMATAN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 56
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMATAN(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.5477326406586162039457549832092678908202994134569781'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMATAN',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST8

      SUBROUTINE TEST9(ST1,ST2,NCASE,NERROR,KLOG)

!  Test hyperbolic functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing hyperbolic routines.')")

      NCASE = 57
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCOSH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '1.0620620786534654254819884264931372964608741056397718'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCOSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 58
      ST1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCOSH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '3.356291383454381441662669560464886179346554730604556E+18'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-31',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCOSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 59
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSINH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.3577371366153083355393138079781276622149524420386975'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSINH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 60
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMSINH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '3.356291383454381441662669560464886179197580776059111E+18'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-31',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMSINH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 61
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMTANH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.3368326049912874057089491946232983472275659538703038'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMTANH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 62
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMTANH(MA,ME)
      CALL FMEQ(ME,MA)
      ST2 = '0.9999999999999999999999999999999999999556135217341837'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMTANH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 63
      ST1 = '0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCHSH(MA,ME,MC)
      CALL FMEQ(ME,MA)
      ST2 = '1.0620620786534654254819884264931372964608741056397718'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-49',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCHSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 64
      ST1 = '-43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCHSH(MA,ME,MC)
      CALL FMEQ(ME,MA)
      ST2 = '3.356291383454381441662669560464886179346554730604556E+18'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-31',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCHSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 65
      ST1 = '-0.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCHSH(MA,MC,ME)
      CALL FMEQ(ME,MA)
      ST2 = '-0.3577371366153083355393138079781276622149524420386975'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-50',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCHSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 66
      ST1 = '43.3505154639175257731958762886597938144329896907216495'
      CALL FMST2M(ST1,MA)
      CALL FMCHSH(MA,MC,ME)
      CALL FMEQ(ME,MA)
      ST2 = '3.356291383454381441662669560464886179197580776059111E+18'
      CALL FMST2M(ST2,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('1.0E-31',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('FMCHSH',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST9

      SUBROUTINE TEST10(ST1,ST2,NCASE,NERROR,KLOG)

!  Input and output testing for IM routines.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

!             Logical function for comparing IM numbers.

      LOGICAL IMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing integer input and output routines.')")

      NCASE = 67
      CALL IMST2M('123',MA)
      CALL IMI2M(123,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMST2M',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 68
      ST1 = '-350515'
      CALL IMST2M(ST1,MA)
      CALL IMI2M(-350515,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMST2M',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 69
      ST1 = '19895113660064588580108197261066338165074766609'
      CALL IMST2M(ST1,MA)
      CALL IMI2M(23,MB)
      CALL IMI2M(34,MC)
      CALL IMPWR(MB,MC,ME)
      CALL IMEQ(ME,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPWR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 70
      ST1 = '-20800708073664542533904165663516279809808659679033703'
      CALL IMST2M(ST1,MA)
      CALL IMI2M(-567,MB)
      CALL IMI2M(19,MC)
      CALL IMPWR(MB,MC,ME)
      CALL IMEQ(ME,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPWR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 71
      ST1 = '19895113660064588580108197261066338165074766609'
      CALL IMST2M(ST1,MA)
      CALL IMFORM('I53',MA,ST2)
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMFORM',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 72
      ST1 = '-20800708073664542533904165663516279809808659679033703'
      CALL IMST2M(ST1,MA)
      CALL IMFORM('I73',MA,ST2)
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMFORM',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST10

      SUBROUTINE TEST11(ST1,ST2,NCASE,NERROR,KLOG)

!  Test add and subtract for IM routines.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL IMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing integer add and subtract routines.')")

      NCASE = 73
      CALL IMST2M('123',MA)
      CALL IMST2M('789',MB)
      CALL IMADD(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMI2M(912,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMADD ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 74
      ST1 = '3505154639175257731958762886597938144329896907216495'
      CALL IMST2M(ST1,MA)
      ST1 = '7319587628865979381443298969072164948453608247422680'
      CALL IMST2M(ST1,MB)
      CALL IMADD(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '10824742268041237113402061855670103092783505154639175'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMADD ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 75
      ST1 = '3505154639175257731958762886597938144329896907216495'
      CALL IMST2M(ST1,MA)
      ST1 = '7319587628865979381443298969072164948453608247422680'
      CALL IMST2M(ST1,MB)
      CALL IMSUB(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '-3814432989690721649484536082474226804123711340206185'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMSUB ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 76
      ST1 = '3505154639175257731958762886597938144329896907216495'
      CALL IMST2M(ST1,MA)
      ST1 = '3505154639175257731443298969072164948453608247422680'
      CALL IMST2M(ST1,MB)
      CALL IMSUB(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '515463917525773195876288659793815'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMSUB ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST11

      SUBROUTINE TEST12(ST1,ST2,NCASE,NERROR,KLOG)

!  Test integer multiply and divide.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL IMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER IREM,KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing integer multiply, divide and square routines.')")

      NCASE = 77
      CALL IMST2M('123',MA)
      CALL IMST2M('789',MB)
      CALL IMMPY(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMI2M(97047,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPY ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 78
      ST1 = '10430738374625018354698'
      CALL IMST2M(ST1,MA)
      ST1 = '2879494424799214514791045985'
      CALL IMST2M(ST1,MB)
      CALL IMMPY(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '30035252996271960952238822892375588336807158787530'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPY ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 79
      CALL IMST2M('12347',MA)
      CALL IMST2M('47',MB)
      CALL IMDIV(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('262',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIV ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 80
      ST1 = '2701314697583086005158008013691015597308949443159762'
      CALL IMST2M(ST1,MA)
      ST1 = '-978132616472842669976589722394'
      CALL IMST2M(ST1,MB)
      CALL IMDIV(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('-2761705981469115610382',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIV ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 81
      CALL IMST2M('12368',MA)
      CALL IMST2M('67',MB)
      CALL IMMOD(MA,MB,ME)
      CALL IMEQ(ME,MB)
      CALL IMST2M('40',MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMOD ',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 82
      ST1 = '2701314697583086005158008013691015597308949443159762'
      CALL IMST2M(ST1,MA)
      ST1 = '-978132616472842669976589722394'
      CALL IMST2M(ST1,MB)
      CALL IMMOD(MA,MB,ME)
      CALL IMEQ(ME,MB)
      CALL IMST2M('450750319653685523300198865254',MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMOD ',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 83
      CALL IMST2M('1234',MA)
      CALL IMST2M('17',MB)
      CALL IMDIVR(MA,MB,MC,MD)
      CALL IMEQ(MC,MA)
      CALL IMEQ(MD,MB)
      CALL IMST2M('72',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVR',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF
      CALL IMST2M('10',MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVR',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 84
      ST1 = '34274652243817531418235301715935108945364446765801943'
      CALL IMST2M(ST1,MA)
      ST1 = '-54708769795848731641842224621693'
      CALL IMST2M(ST1,MB)
      CALL IMDIVR(MA,MB,MC,MD)
      CALL IMEQ(MC,MA)
      CALL IMEQ(MD,MB)
      CALL IMST2M('-626492834178447772323',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVR',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF
      CALL IMST2M('31059777254296217822749494999104',MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVR',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 85
      CALL IMST2M('4866',MA)
      CALL IMMPYI(MA,14,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('68124',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPYI',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 86
      CALL IMST2M('270131469758308600515800801369101559730894',MA)
      CALL IMMPYI(MA,-2895,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('-782030604950303398493243319963549015420938130',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPYI ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 87
      CALL IMST2M('-37179',MA)
      CALL IMDIVI(MA,129,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('-288',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVI',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 88
      ST1 = '8267538919383255454483790743961990401918726073065738'
      CALL IMST2M(ST1,MA)
      CALL IMDIVI(MA,1729,ME)
      CALL IMEQ(ME,MA)
      ST2 = '4781688212483085861471249707323302719444028960708'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDIVI',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 89
      CALL IMST2M('-71792',MA)
      CALL IMDVIR(MA,65,MC,IREM)
      CALL IMEQ(MC,MA)
      CALL IMST2M('-1104',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDVIR',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF
      CALL IMI2M(IREM,MB)
      CALL IMI2M(-32,MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDVIR',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 90
      ST1 = '97813261647284266997658972239417958580120170263408655'
      CALL IMST2M(ST1,MA)
      CALL IMDVIR(MA,826,MC,IREM)
      CALL IMEQ(MC,MA)
      ST2 = '118417992309060855929369215786220288837917881674828'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDVIR',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF
      CALL IMI2M(IREM,MB)
      CALL IMI2M(727,MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMDVIR',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 91
      CALL IMST2M('538',MA)
      CALL IMSQR(MA,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('289444',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMSQR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 92
      CALL IMST2M('-47818191879814587168242632',MA)
      CALL IMSQR(MA,ME)
      CALL IMEQ(ME,MA)
      ST2 = '2286579474654765721668058416662636606051551222287424'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMSQR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST12

      SUBROUTINE TEST13(NCASE,NERROR,KLOG)

!  Test conversions between FM and IM format.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP,IMCOMP

      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing conversions between FM and IM format.')")

      NCASE = 93
      CALL IMST2M('123',MA)
      CALL IMI2FM(MA,MB)
      CALL FMI2M(123,MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('IMI2FM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 94
      CALL IMST2M('979282999076598337488362000995916',MA)
      CALL IMI2FM(MA,MB)
      CALL FMST2M('979282999076598337488362000995916',MC)
      CALL FMSUB(MA,MC,MD)
      CALL FMABS(MD,ME)
      CALL FMEQ(ME,MD)
      CALL FMST2M('0',MB)
      IF (.NOT.FMCOMP(MD,'LE',MB)) THEN
          CALL ERRPRTFM('IMI2FM',MA,'MA',MC,'MC',MD,'MD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 95
      CALL FMST2M('123.4',MA)
      CALL IMFM2I(MA,MB)
      CALL IMI2M(123,MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMFM2I',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 96
      CALL FMST2M('979282999076598337488362000995916',MA)
      CALL IMFM2I(MA,MB)
      CALL IMST2M('979282999076598337488362000995916',MC)
      IF (.NOT.IMCOMP(MB,'EQ',MC)) THEN
          CALL ERRPRTIM('IMFM2I',MB,'MB',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST13

      SUBROUTINE TEST14(ST1,ST2,NCASE,NERROR,KLOG)

!  Test integer power and GCD functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL IMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing integer GCD and power routines.')")

      NCASE = 97
      CALL IMST2M('123',MA)
      CALL IMST2M('789',MB)
      CALL IMGCD(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMI2M(3,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMGCD ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 98
      ST1 = '431134020618556701030927835051546391752577319587628885'
      CALL IMST2M(ST1,MA)
      ST1 = '900309278350515463917525773195876288659793814432989640'
      CALL IMST2M(ST1,MB)
      CALL IMGCD(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('615',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMGCD ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 99
      ST1 = '5877631675869176172956662762822298812326084745145447940'
      CALL IMST2M(ST1,MA)
      ST1 = '10379997509886032090765062511740075746391432253007667'
      CALL IMST2M(ST1,MB)
      CALL IMGCD(MA,MB,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('1',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMGCD ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 100
      CALL IMST2M('47',MA)
      CALL IMST2M('34',MB)
      CALL IMPWR(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '710112520079088427392020925014421733344154169313556279969'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPWR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 101
      CALL IMST2M('2',MA)
      CALL IMST2M('187',MB)
      CALL IMPWR(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '196159429230833773869868419475239575503198607639501078528'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPWR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 102
      CALL IMST2M('-3',MA)
      CALL IMST2M('101',MB)
      CALL IMPWR(MA,MB,ME)
      CALL IMEQ(ME,MA)
      ST2 = '-1546132562196033993109383389296863818106322566003'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPWR ',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST14

      SUBROUTINE TEST15(ST1,ST2,NCASE,NERROR,KLOG)

!  Test integer modular functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL IMCOMP

      CHARACTER(80) :: ST1,ST2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing integer modular routines.')")

      NCASE = 103
      CALL IMST2M('123',MA)
      CALL IMST2M('789',MB)
      CALL IMST2M('997',MC)
      CALL IMMPYM(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      CALL IMI2M(338,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPYM',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 104
      ST1 = '431134020618556701030927835051546391752577319587628885'
      CALL IMST2M(ST1,MA)
      ST1 = '36346366019557973241042306587666640486264616086971724'
      CALL IMST2M(ST1,MB)
      ST1 = '900309278350515463917525773195876288659793814432989640'
      CALL IMST2M(ST1,MC)
      CALL IMMPYM(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      ST2 = '458279704440780378752997531208983184411293504187816380'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPYM',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 105
      ST1 = '914726194238000125985765939883182'
      CALL IMST2M(ST1,MA)
      ST1 = '-75505764717193044779376979508186553225192'
      CALL IMST2M(ST1,MB)
      ST1 = '18678872625055834600521936'
      CALL IMST2M(ST1,MC)
      CALL IMMPYM(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      ST2 = '-7769745969769966093344960'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMMPYM',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 106
      CALL IMST2M('123',MA)
      CALL IMST2M('789',MB)
      CALL IMST2M('997',MC)
      CALL IMPMOD(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      CALL IMI2M(240,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPMOD',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 107
      ST1 = '431134020618556701030927835051546391752577319587628885'
      CALL IMST2M(ST1,MA)
      ST1 = '36346366019557973241042306587666640486264616086971724'
      CALL IMST2M(ST1,MB)
      ST1 = '900309278350515463917525773195876288659793814432989640'
      CALL IMST2M(ST1,MC)
      CALL IMPMOD(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      ST2 = '755107893576299697276281907390144058060594744720442385'
      CALL IMST2M(ST2,MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPMOD',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 108
      CALL IMST2M('314159',MA)
      CALL IMST2M('1411695892374393248272691827763664225585897550',MB)
      CALL IMST2M('1411695892374393248272691827763664225585897551',MC)
      CALL IMPMOD(MA,MB,MC,ME)
      CALL IMEQ(ME,MA)
      CALL IMST2M('1',MC)
      IF (.NOT.IMCOMP(MA,'EQ',MC)) THEN
          CALL ERRPRTIM('IMPMOD',MA,'MA',MC,'MC',NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST15

      SUBROUTINE TEST16(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Complex input and output testing.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

!             Logical function for comparing FM numbers.

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex input and output routines.')")

      NCASE = 109
      CALL ZMST2M('123 + 456 i',ZA)
      CALL ZM2I2M(123,456,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)

!             Use the .NOT. because FMCOMP returns FALSE for special
!             cases like ZD = UNKNOWN, and these should be treated
!             as errors for these tests.

      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMST2M',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 110
      STZ1 = '0.3505154639175257731958762886597938144329896907216495 + '  &
         // '0.7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZM2I2M(34,71,ZC)
      CALL ZMDIVI(ZC,97,ZE)
      CALL ZMEQ(ZE,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMST2M',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 111
      STZ1 = '0.3505154639175257731958762886597938144329896907216495E-5 '  &
       //'+ 0.7319587628865979381443298969072164948453608247422680D-5 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZM2I2M(34,71,ZC)
      CALL ZMDIVI(ZC,9700000,ZE)
      CALL ZMEQ(ZE,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-55,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMST2M',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 112
      STZ1 = '7.699115044247787610619469026548672566371681415929204e 03 '  &
       //'- 5.221238938053097345132743362831858407079646017699115M 03 I'
      CALL ZMST2M(STZ1,ZA)
      CALL ZM2I2M(87,-59,ZC)
      CALL ZMDIVI(ZC,113,ZE)
      CALL ZMEQ(ZE,ZC)
      CALL ZMMPYI(ZC,10000,ZE)
      CALL ZMEQ(ZE,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-47,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMST2M',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 113
      STZ1 = '7.699115044247787610619469026548672566371681415929204e+3 '  &
       //'- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMFORM('F53.33','F50.30',ZA,STZ2)
      CALL ZMST2M(STZ2,ZA)
      STZ1 = '7699.115044247787610619469026548673 '  &
       // '-5221.238938053097345132743362831858 i'
      CALL ZMST2M(STZ1,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-30,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMFORM',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 114
      STZ1 = '7.699115044247787610619469026548672566371681415929204e+3 '  &
       //'- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMFORM('I9','I7',ZA,STZ2)
      CALL ZMST2M(STZ2,ZA)
      STZ1 = '7699 -5221 i'
      CALL ZMST2M(STZ1,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(0,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMFORM',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 115
      STZ1 = '7.699115044247787610619469026548672566371681415929204e+3 '  &
       //'- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMFORM('E59.50','E58.49',ZA,STZ2)
      CALL ZMST2M(STZ2,ZA)
      STZ1 = '7.6991150442477876106194690265486725663716814159292E3'  &
       //'- 5.221238938053097345132743362831858407079646017699E3 i'
      CALL ZMST2M(STZ1,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMFORM',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 116
      STZ1 = '7.699115044247787610619469026548672566371681415929204e+3 '  &
       //'- 5.221238938053097345132743362831858407079646017699115M+3 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMFORM('1PE59.50','1PE58.49',ZA,STZ2)
      CALL ZMST2M(STZ2,ZA)
      CALL ZMST2M(STZ1,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-44,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMFORM',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST16

      SUBROUTINE TEST17(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex add and subtract.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex add and subtract routines.')")

      NCASE = 117
      CALL ZMST2M('123 + 456 i',ZA)
      CALL ZMST2M('789 - 543 i',ZB)
      CALL ZMADD(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      CALL ZM2I2M(912,-87,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(0,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMADD ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 118
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMADD(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '1.1204269683423045342578231913146610710701578323145698 '  &
       //'+ 0.2098348690812882036310555606240306541373962229723565 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMADD ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 119
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMSUB(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.4193960405072529878660706139950734422041784508712709 '  &
       //'- 1.2540826566919076726576042331904023355533254265121795 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSUB ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 120
      STZ1 = '.7699115044247787610619469026548672566371681415929204E3 '  &
       //'- .5221238938053097345132743362831858407079646017699115E3 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMSUB(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '769.5609889608612352887510263662074628227351519021987045 '  &
       //'- 522.8558525681963324514186661800930572028099625946537725 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-47,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSUB ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST17

      SUBROUTINE TEST18(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex multiply, divide and square root.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,  &
            "(/' Testing complex multiply, divide and square root routines.')")

      NCASE = 121
      CALL ZMST2M('123 + 456 i',ZA)
      CALL ZMST2M('789 - 543 i',ZB)
      CALL ZMMPY(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      CALL ZM2I2M(344655,292995,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(0,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMMPY ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 122
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMMPY(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.6520390475321594745005017790347596022260742632971444 '  &
       //'+ 0.3805309734513274336283185840707964601769911504424779 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMMPY ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 123
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMDIV(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-.1705178497731560089737969128653459210208765017614861 '  &
       //'- 1.1335073636829696356072949942949842987114804337239972 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMDIV ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 124
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMMPYI(ZA,36,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '27.7168141592920353982300884955752212389380530973451327 '  &
       //'- 18.7964601769911504424778761061946902654867256637168142 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMMPYI',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 125
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMDIVI(ZA,37,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '2.080841903850753408275532169337479071992346328629514E-2 '  &
       //'- 1.411145658933269552738579287251853623535039464243004E-2 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-52,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMDIVI',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 126
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSQR(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.3201503641632077688150990680554467851828647505677813 '  &
       //'- 0.8039783851515388832328295089670295246299631921058814 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSQR ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 127
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSQRT(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.9219999909012323458336720551458583330580388434229845 '  &
       //'- 0.2831474506279259570386845864488094697732718981999941 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSQRT',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST18

      SUBROUTINE TEST19(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex exponentials.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex exponential routines.')")

      NCASE = 128
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMEXP(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '1.8718374504057787925867989348073888855260008469310002 '  &
       //'- 1.0770279996847678711699041910427261417963102075889234 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMEXP ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 129
      STZ1 = '5.7699115044247787610619469026548672566371681415929204 '  &
       //'- 4.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMEXP(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-60.6144766542152809520229386164396710991242264070603612 '  &
       //'+ 314.7254994809539691403004121118801578835669635535466592 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-47,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMEXP ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 130
      STZ1 = '1.7699115044247787610619469026548672566371681415929204 '  &
       //'- 1.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMIPWR(ZA,45,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '31595668743300099.70429472191424818167262151605608585179 '  &
       //'- 19209634448276799.67717448173630165852744930837930753788 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-33,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMIPWR',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 131
      STZ1 = '1.7699115044247787610619469026548672566371681415929204 '  &
       //'- 1.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMIPWR(ZA,-122,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '3.1000215641022021714480000129414241564868699479432E-46 '  &
       //'- 1.1687846789859477815450163510927243367234863123667E-45 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-93,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMIPWR',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 132
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMPWR(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '1.4567089343012352449621841355636496276866203747888724 '  &
       //'- 0.3903177712261966292764255714390622205129978923650749 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMPWR ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 133
      STZ1 = '.3505154639175257731958762886597938144329896907216495 '  &
       //'+ .7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      STZ1 = '2.7699115044247787610619469026548672566371681415929204 '  &
       //'- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZB)
      CALL ZMPWR(ZA,ZB,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-1.0053105716678380336247948739245187868180079734997482 '  &
       // '- 0.0819537653234704467729051473979237153087038930127116 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMPWR ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 134
      STZ1 = '0.7699115044247787610619469026548672566371681415929204 '  &
       //'- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMRPWR(ZA,2,7,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.9653921326136512316639621651337975772631340364271270 '  &
       //'- 0.1659768285667051396562270035411852432430188906482848 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMRPWR',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 135
      STZ1 = '0.7699115044247787610619469026548672566371681415929204 '  &
       //'- 0.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMRPWR(ZA,-19,7,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-0.0567985880053556315170006800325686036902111276420647 '  &
       // '+ 1.2154793972711356706410882510363594270389067962568571 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMRPWR',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST19

      SUBROUTINE TEST20(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex logarithms.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex logarithm routines.')")

      NCASE = 136
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMLN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-0.0722949652393911311212450699415231782692434885813725 '  &
       //'-  0.5959180055163009910007765127008371205749515965219804 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMLN  ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 137
      STZ1 = '.7699115044247787610619469026548672566371681415929204E28 '  &
       //'- .5221238938053097345132743362831858407079646017699115E28 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMLN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '64.4000876385938880213825156612206746345615981930242708 '  &
       //'-  0.5959180055163009910007765127008371205749515965219804 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMLN  ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 138
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMLG10(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-0.0313973044728549715287589498363619677438302809470943 '  &
       //'-  0.2588039014625211035392823012785304771809982053965284 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMLG10',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 139
      STZ1 = '.7699115044247787610619469026548672566371681415929204E82 '  &
       //'- .5221238938053097345132743362831858407079646017699115E82 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMLG10(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '81.9686026955271450284712410501636380322561697190529057 '  &
       //'-  0.2588039014625211035392823012785304771809982053965284 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMLG10',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST20

      SUBROUTINE TEST21(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex trigonometric functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex trigonometric routines.')")

      NCASE = 140
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCOS(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.8180802525254482451348613286211514555816444253416895 '  &
       //'+  0.3801751200076938035500853542125525088505055292851393 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCOS ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 141
      STZ1 = '34.7699115044247787610619469026548672566371681415929204 '  &
       //'- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCOS(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-1432925478410268113.5816466154230974355002592549420099 '  &
       //'-  309002816679456015.00151246245263842483282458519462258 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-31,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCOS ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 142
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSIN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.7931260548991613428648822413402447097755865697557818 '  &
       //'-  0.3921366045897070762848927655743167937790944353110710 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSIN ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 143
      STZ1 = '34.7699115044247787610619469026548672566371681415929204 '  &
       //'- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSIN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '-3.090028166794560150015124624526384249047272360765358E17 '  &
       //'+  1.432925478410268113581646615423097435166828182950161E18 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-31,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSIN ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 144
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMTAN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.6141156219447569167198437040270236055089243090199979 '  &
       //'-  0.7647270337230070156308196055474639461102792169274526 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMTAN ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 145
      STZ1 = '35.7699115044247787610619469026548672566371681415929204 '  &
       //'- 43.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMTAN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '2.068934241218867332441292427642153175237611151321340E-38 '  &
       //'-  1.000000000000000000000000000000000000023741659169354 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMTAN ',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 146
      STZ1 = '0.3505154639175257731958762886597938144329896907216495 '  &
       //'+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCSSN(ZA,ZE,ZC)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '1.2022247452809115256533054407001508718694617802593324 '  &
       //'-  0.2743936538120352873902095801531325075994392065668943 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCSSN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 147
      STZ1 = '0.3505154639175257731958762886597938144329896907216495 '  &
       //'+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCSSN(ZA,ZC,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.4395486978082638069281369170831952476351663772871008 '  &
       //'+  0.7505035100906417134864779281080728222900154610025883 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCSSN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST21

      SUBROUTINE TEST22(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex inverse trigonometric functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex inverse trigonometric routines.')")

      NCASE = 148
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMACOS(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.8797127900868121872960714368309657795959216549012347 '  &
       //'+  0.6342141347945396859119941874681961111936156338608130 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMACOS',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 149
      STZ1 = '.7699115044247787610619469026548672566371681415929204E12 '  &
       //'- .5221238938053097345132743362831858407079646017699115E12 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMACOS(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.5959180055163009910007767810953294528367807973983794 '  &
       //'+28.2518733312491023865118844008522768856672089946951468 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMACOS',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 150
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMASIN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.6910835367080844319352502548087856625026630447863182 '  &
       //'-  0.6342141347945396859119941874681961111936156338608130 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMASIN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 151
      STZ1 = '.7699115044247787610619469026548672566371681415929204E13 '  &
       //'- .5221238938053097345132743362831858407079646017699115E13 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMASIN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.9748783212785956282305451762549693982010148111568094 '  &
       //'-30.5544584242431480705298759613446206186670533428066404 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-48,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMASIN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 152
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMATAN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.7417952692265900376512911713942700568648670953521258 '  &
       //'- 0.3162747143126729004878357203292329539837025170484857 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMATAN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 153
      STZ1 = '.7699115044247787610619469026548672566371681415929204E13 '  &
       //'- .5221238938053097345132743362831858407079646017699115E13 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMATAN(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '   1.570796326794807650905529836436131532596233124329403 '  &
       //'-6.033484162895927601809954710695221401671437742867605E-14 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMATAN',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST22

      SUBROUTINE TEST23(STZ1,STZ2,NCASE,NERROR,KLOG)

!  Test complex hyperbolic functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      LOGICAL FMCOMP

      CHARACTER(160) :: STZ1,STZ2
      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing complex hyperbolic routines.')")

      NCASE = 154
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCOSH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '1.1365975275870879962259716562608779977957563621412079 '  &
       //'-  0.4230463404769118342540441830446134405410543954181579 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-49,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCOSH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 155
      STZ1 = '34.7699115044247787610619469026548672566371681415929204 '  &
       //'- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCOSH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '69552104658681.7558589320148420094288419217262200765435 '  &
       //'+ 626163773308016.884007302915197616300902876551542156676 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-35,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCOSH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 156
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSINH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.7352399228186907963608272785465108877302444847897922 '  &
       //'-  0.6539816592078560369158600079981127012552558121707655 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSINH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 157
      STZ1 = '34.7699115044247787610619469026548672566371681415929204 '  &
       //'- 42.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMSINH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '6.955210465868175585893201484192181376093291191637290E 13 '  &
       //'+ 6.261637733080168840073029151984050820616907795167046E 14 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-35,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMSINH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 158
      STZ1 = '.7699115044247787610619469026548672566371681415929204 '  &
       //'- .5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMTANH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.7562684782933185240709480231996041186654551038993505 '  &
       //'-  0.2938991498221693198532255749292372853685311106820169 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMTANH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 159
      STZ1 = '35.7699115044247787610619469026548672566371681415929204 '  &
       //'- 43.5221238938053097345132743362831858407079646017699115 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMTANH(ZA,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '9.999999999999999999999999999998967653135180689424497E-01 '  &
       //'+ 1.356718776492102400812550018433337461876455254467192E-31 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMTANH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 160
      STZ1 = '0.3505154639175257731958762886597938144329896907216495 '  &
       //'+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCHSH(ZA,ZE,ZC)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.7900326499280864816444807620997665088044412803737969 '  &
       //'+ 0.2390857359988804105051429301542214823277594407302781 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCHSH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 161
      STZ1 = '0.3505154639175257731958762886597938144329896907216495 '  &
       //'+  0.7319587628865979381443298969072164948453608247422680 i'
      CALL ZMST2M(STZ1,ZA)
      CALL ZMCHSH(ZA,ZC,ZE)
      CALL ZMEQ(ZE,ZA)
      STZ2 = '0.2661087555034471983220879532235334422670297141428191 '  &
       //'+  0.7098057980612199357870532628105009808447460332437714 i'
      CALL ZMST2M(STZ2,ZC)
      CALL ZMSUB(ZA,ZC,ZD)
      CALL ZMABS(ZD,MA)
      CALL FMI2M(10,MB)
      CALL FMIPWR(MB,-50,ME)
      CALL FMEQ(ME,MB)
      IF (.NOT.FMCOMP(MA,'LE',MB)) THEN
          CALL ERRPRTZM('ZMCHSH',ZA,'ZA',ZC,'ZC',ZD,'ZD',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST23

      SUBROUTINE TEST24(NCASE,NERROR,KLOG)

!             Test the = assignment interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER NERROR,NCASE,KLOG
      LOGICAL FM_COMP,IM_COMP

      WRITE (KW,"(/' Testing the derived type = interface.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0
      MSMALL = EPSILON(TO_FM(1))*10000.0
      NCASE = 162
      J4 = MFM1
      IF (J4 /= 581) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 163
      J4 = MIM1
      IF (J4 /= 661) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 164
      J4 = MZM1
      IF (J4 /= 731) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 165
      R4 = MFM1
      IF (ABS((R4-581.21)/581.21) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 166
      R4 = MIM1
      IF (ABS((R4-661.0)/661.0) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 167
      R4 = MZM1
      IF (ABS((R4-731.51)/731.51) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 168
      D4 = MFM1
      IF (ABS((D4-581.21D0)/581.21D0) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 169
      D4 = MIM1
      IF (ABS((D4-661.0D0)/661.0D0) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 170
      D4 = MZM1
      IF (ABS((D4-731.51D0)/731.51D0) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 171
      C4 = MFM1
      IF (ABS((C4-581.21)/581.21) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 172
      C4 = MIM1
      IF (ABS((C4-661.0)/661.0) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 173
      C4 = MZM1
      IF (ABS((C4-(731.51,711.41))/(731.51,711.41)) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 174
      CD4 = MFM1
      IF (ABS((CD4-581.21D0)/581.21D0) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 175
      CD4 = MIM1
      IF (ABS((CD4-661.0D0)/661.0D0) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 176
      CD4 = MZM1
      IF (ABS((CD4-(731.51D0,711.41D0))/(731.51D0,711.41D0)) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 177
      MFM3 = J2
      CALL FM_I2M(131,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ST2M('0',MFM3)
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 178
      MFM3 = R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 179
      MFM3 = D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 180
      MFM3 = C2
      CALL FM_ST2M('411.11',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 181
      MFM3 = CD2
      CALL FM_ST2M('431.11',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 182
      MFM3 = MFM1
      CALL FM_ST2M('581.21',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_EQ(MSMALL,MFM3)
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 183
      MFM3 = MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ST2M('0',MFM3)
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 184
      MFM3 = MZM1
      CALL FM_ST2M('731.51',MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 185
      MIM3 = J2
      CALL IM_I2M(131,MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 186
      MIM3 = R2
      CALL IM_ST2M('241',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 187
      MIM3 = D2
      CALL IM_ST2M('391',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 188
      MIM3 = C2
      CALL IM_ST2M('411',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 189
      MIM3 = CD2
      CALL IM_ST2M('431',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 190
      MIM3 = MFM1
      CALL IM_ST2M('581',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 191
      MIM3 = MIM1
      CALL IM_ST2M('661',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 192
      MIM3 = MZM1
      CALL IM_ST2M('731',MIM4)
      CALL IM_SUB(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      CALL IM_ST2M('0',MIM3)
      IF (IM_COMP(MIM4,'GT',MIM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 193
      MZM3 = J2
      CALL ZM_I2M(131,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      CALL FM_ST2M('0',MFM3)
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 194
      MZM3 = R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 195
      MZM3 = D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 196
      MZM3 = C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 197
      MZM3 = CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 198
      MZM3 = MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = MSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 199
      MZM3 = MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      CALL FM_ST2M('0',MFM3)
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 200
      MZM3 = MZM1
      CALL ZM_ST2M('731.51 + 711.41 i',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = MSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST24

      SUBROUTINE TEST25(NCASE,NERROR,KLOG)

!  Test the derived type == interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type == interface.')")

      NCASE = 201
      M_A = 123
      M_B = M_A
      IF (.NOT.FM_COMP(M_A,'==',M_B)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 202
      M_A = 123
      M_B = M_A
      IF (.NOT.FM_COMP(M_A,'EQ',M_B)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 203
      J1 = 123
      M_A = J1
      IF (.NOT.(M_A == J1)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 204
      J1 = 123
      M_A = J1
      IF (.NOT.(J1 == M_A)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 205
      J1 = 123
      M_J = J1
      IF (.NOT.(M_J == J1)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 206
      J1 = 123
      M_J = J1
      IF (.NOT.(J1 == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 207
      J1 = 123
      M_Z = J1
      IF (.NOT.(M_Z == J1)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 208
      J1 = 123
      M_Z = J1
      IF (.NOT.(J1 == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 209
      R1 = 12.3
      M_A = R1
      IF (.NOT.(M_A == R1)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 210
      R1 = 12.3
      M_A = R1
      IF (.NOT.(R1 == M_A)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 211
      R1 = 123
      M_J = R1
      IF (.NOT.(M_J == R1)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 212
      R1 = 123
      M_J = R1
      IF (.NOT.(R1 == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 213
      R1 = 12.3
      M_Z = R1
      IF (.NOT.(M_Z == R1)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 214
      R1 = 12.3
      M_Z = R1
      IF (.NOT.(R1 == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 215
      D1 = 12.3
      M_A = D1
      IF (.NOT.(M_A == D1)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 216
      D1 = 12.3
      M_A = D1
      IF (.NOT.(D1 == M_A)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 217
      D1 = 123
      M_J = D1
      IF (.NOT.(M_J == D1)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 218
      D1 = 123
      M_J = D1
      IF (.NOT.(D1 == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 219
      D1 = 12.3
      M_Z = D1
      IF (.NOT.(M_Z == D1)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 220
      D1 = 12.3
      M_Z = D1
      IF (.NOT.(D1 == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 221
      C1 = 12.3
      M_A = C1
      IF (.NOT.(M_A == C1)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 222
      C1 = 12.3
      M_A = C1
      IF (.NOT.(C1 == M_A)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 223
      C1 = 123
      M_J = C1
      IF (.NOT.(M_J == C1)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 224
      C1 = 123
      M_J = C1
      IF (.NOT.(C1 == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 225
      C1 = (12.3 , 45.6)
      M_Z = C1
      IF (.NOT.(M_Z == C1)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 226
      C1 = (12.3 , 45.6)
      M_Z = C1
      IF (.NOT.(C1 == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 227
      CD1 = 12.3
      M_A = CD1
      IF (.NOT.(M_A == CD1)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 228
      CD1 = 12.3
      M_A = CD1
      IF (.NOT.(CD1 == M_A)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 229
      CD1 = 123
      M_J = CD1
      IF (.NOT.(M_J == CD1)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 230
      CD1 = 123
      M_J = CD1
      IF (.NOT.(CD1 == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 231
      CD1 = (12.3 , 45.6)
      M_Z = CD1
      IF (.NOT.(M_Z == CD1)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 232
      CD1 = (12.3 , 45.6)
      M_Z = CD1
      IF (.NOT.(CD1 == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 233
      M_B = 12.3
      M_A = M_B
      IF (.NOT.(M_A == M_B)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 234
      M_B = 123
      M_J = M_B
      IF (.NOT.(M_J == M_B)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 235
      M_B = 123
      M_J = M_B
      IF (.NOT.(M_B == M_J)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 236
      M_B = (12.3 , 45.6)
      M_Z = M_B
      IF (.NOT.(M_Z == M_B)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 237
      M_B = (12.3 , 45.6)
      M_Z = M_B
      IF (.NOT.(M_B == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 238
      M_K = 123
      M_J = M_K
      IF (.NOT.(M_J == M_K)) THEN
          CALL ERRPRT_IM('  ==  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 239
      M_K = (12.3 , 45.6)
      M_Z = M_K
      IF (.NOT.(M_Z == M_K)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 240
      M_K = (12.3 , 45.6)
      M_Z = M_K
      IF (.NOT.(M_K == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 241
      M_Y = (12.3 , 45.6)
      M_Z = M_Y
      IF (.NOT.(M_Y == M_Z)) THEN
          CALL ERRPRT_ZM('  ==  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST25

      SUBROUTINE TEST26(NCASE,NERROR,KLOG)

!  Test the derived type /= interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type /= interface.')")

      NCASE = 242
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'/=',M_B)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 243
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'NE',M_B)) THEN
          CALL ERRPRT_FM('  ==  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 244
      J1 = 123
      M_A = 1 + J1
      IF (.NOT.(M_A /= J1)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 245
      J1 = 123
      M_A = 1 + J1
      IF (.NOT.(J1 /= M_A)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 246
      J1 = 123
      M_J = 1 + J1
      IF (.NOT.(M_J /= J1)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 247
      J1 = 123
      M_J = 1 + J1
      IF (.NOT.(J1 /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 248
      J1 = 123
      M_Z = 1 + J1
      IF (.NOT.(M_Z /= J1)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 249
      J1 = 123
      M_Z = 1 + J1
      IF (.NOT.(J1 /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 250
      R1 = 12.3
      M_A = 1 + R1
      IF (.NOT.(M_A /= R1)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 251
      R1 = 12.3
      M_A = 1 + R1
      IF (.NOT.(R1 /= M_A)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 252
      R1 = 123
      M_J = 1 + R1
      IF (.NOT.(M_J /= R1)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 253
      R1 = 123
      M_J = 1 + R1
      IF (.NOT.(R1 /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 254
      R1 = 12.3
      M_Z = 1 + R1
      IF (.NOT.(M_Z /= R1)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 255
      R1 = 12.3
      M_Z = 1 + R1
      IF (.NOT.(R1 /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 256
      D1 = 12.3
      M_A = 1 + D1
      IF (.NOT.(M_A /= D1)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 257
      D1 = 12.3
      M_A = 1 + D1
      IF (.NOT.(D1 /= M_A)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 258
      D1 = 123
      M_J = 1 + D1
      IF (.NOT.(M_J /= D1)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 259
      D1 = 123
      M_J = 1 + D1
      IF (.NOT.(D1 /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 260
      D1 = 12.3
      M_Z = 1 + D1
      IF (.NOT.(M_Z /= D1)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 261
      D1 = 12.3
      M_Z = 1 + D1
      IF (.NOT.(D1 /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 262
      C1 = 12.3
      M_A = 1 + C1
      IF (.NOT.(M_A /= C1)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 263
      C1 = 12.3
      M_A = 1 + C1
      IF (.NOT.(C1 /= M_A)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 264
      C1 = 123
      M_J = 1 + C1
      IF (.NOT.(M_J /= C1)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 265
      C1 = 123
      M_J = 1 + C1
      IF (.NOT.(C1 /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 266
      C1 = (12.3 , 45.6)
      M_Z = 1 + C1
      IF (.NOT.(M_Z /= C1)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 267
      C1 = (12.3 , 45.6)
      M_Z = 1 + C1
      IF (.NOT.(C1 /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 268
      CD1 = 12.3
      M_A = 1 + CD1
      IF (.NOT.(M_A /= CD1)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 269
      CD1 = 12.3
      M_A = 1 + CD1
      IF (.NOT.(CD1 /= M_A)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 270
      CD1 = 123
      M_J = 1 + CD1
      IF (.NOT.(M_J /= CD1)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 271
      CD1 = 123
      M_J = 1 + CD1
      IF (.NOT.(CD1 /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 272
      CD1 = (12.3 , 45.6)
      M_Z = 1 + CD1
      IF (.NOT.(M_Z /= CD1)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 273
      CD1 = (12.3 , 45.6)
      M_Z = 1 + CD1
      IF (.NOT.(CD1 /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 274
      M_B = 12.3
      M_A = 1 + M_B
      IF (.NOT.(M_A /= M_B)) THEN
          CALL ERRPRT_FM('  /=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 275
      M_B = 123
      M_J = 1 + M_B
      IF (.NOT.(M_J /= M_B)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 276
      M_B = 123
      M_J = 1 + M_B
      IF (.NOT.(M_B /= M_J)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 277
      M_B = (12.3 , 45.6)
      M_Z = 1 + M_B
      IF (.NOT.(M_Z /= M_B)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 278
      M_B = (12.3 , 45.6)
      M_Z = 1 + M_B
      IF (.NOT.(M_B /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 279
      M_K = 123
      M_J = 1 + M_K
      IF (.NOT.(M_J /= M_K)) THEN
          CALL ERRPRT_IM('  /=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 280
      M_K = (12.3 , 45.6)
      M_Z = 1 + M_K
      IF (.NOT.(M_Z /= M_K)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 281
      M_K = (12.3 , 45.6)
      M_Z = 1 + M_K
      IF (.NOT.(M_K /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 282
      M_Y = (12.3 , 45.6)
      M_Z = 1 + M_Y
      IF (.NOT.(M_Y /= M_Z)) THEN
          CALL ERRPRT_ZM('  /=  ',M_Z,'M_Z',M_Z,'M_Z',M_Z,'M_Z',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST26

      SUBROUTINE TEST27(NCASE,NERROR,KLOG)

!  Test the derived type > interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type > interface.')")

      NCASE = 283
      M_A = 125
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'>',M_B)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 284
      M_A = 125
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'GT',M_B)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 285
      J1 = 123
      M_A = J1 + 1
      IF (.NOT.(M_A > J1)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 286
      J1 = 123
      M_A = J1 - 1
      IF (.NOT.(J1 > M_A)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 287
      J1 = 123
      M_J = J1 + 1
      IF (.NOT.(M_J > J1)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 288
      J1 = 123
      M_J = J1 - 1
      IF (.NOT.(J1 > M_J)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 289
      R1 = 12.3
      M_A = R1 + 1
      IF (.NOT.(M_A > R1)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 290
      R1 = 12.3
      M_A = R1 - 1
      IF (.NOT.(R1 > M_A)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 291
      R1 = 123
      M_J = R1 + 1
      IF (.NOT.(M_J > R1)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 292
      R1 = 123
      M_J = R1 - 1
      IF (.NOT.(R1 > M_J)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 293
      D1 = 12.3
      M_A = D1 + 1
      IF (.NOT.(M_A > D1)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 294
      D1 = 12.3
      M_A = D1 - 1
      IF (.NOT.(D1 > M_A)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 295
      D1 = 123
      M_J = D1 + 1
      IF (.NOT.(M_J > D1)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 296
      D1 = 123
      M_J = D1 - 1
      IF (.NOT.(D1 > M_J)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 297
      M_B = 12.3
      M_A = M_B + 1
      IF (.NOT.(M_A > M_B)) THEN
          CALL ERRPRT_FM('   >  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 298
      M_B = 123
      M_J = M_B + 1
      IF (.NOT.(M_J > M_B)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 299
      M_B = 123
      M_J = M_B - 1
      IF (.NOT.(M_B > M_J)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 300
      M_K = 123
      M_J = M_K + 1
      IF (.NOT.(M_J > M_K)) THEN
          CALL ERRPRT_IM('   >  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST27

      SUBROUTINE TEST28(NCASE,NERROR,KLOG)

!  Test the derived type >= interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type >= interface.')")

      NCASE = 301
      M_A = 125
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'>=',M_B)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 302
      M_A = 125
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'GE',M_B)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 303
      J1 = 123
      M_A = J1 + 1
      IF (.NOT.(M_A >= J1)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 304
      J1 = 123
      M_A = J1 - 1
      IF (.NOT.(J1 >= M_A)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 305
      J1 = 123
      M_J = J1 + 1
      IF (.NOT.(M_J >= J1)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 306
      J1 = 123
      M_J = J1 - 1
      IF (.NOT.(J1 >= M_J)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 307
      R1 = 12.3
      M_A = R1 + 1
      IF (.NOT.(M_A >= R1)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 308
      R1 = 12.3
      M_A = R1 - 1
      IF (.NOT.(R1 >= M_A)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 309
      R1 = 123
      M_J = R1 + 1
      IF (.NOT.(M_J >= R1)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 310
      R1 = 123
      M_J = R1 - 1
      IF (.NOT.(R1 >= M_J)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 311
      D1 = 12.3
      M_A = D1 + 1
      IF (.NOT.(M_A >= D1)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 312
      D1 = 12.3
      M_A = D1 - 1
      IF (.NOT.(D1 >= M_A)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 313
      D1 = 123
      M_J = D1 + 1
      IF (.NOT.(M_J >= D1)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 314
      D1 = 123
      M_J = D1 - 1
      IF (.NOT.(D1 >= M_J)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 315
      M_B = 12.3
      M_A = M_B + 1
      IF (.NOT.(M_A >= M_B)) THEN
          CALL ERRPRT_FM('  >=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 316
      M_B = 123
      M_J = M_B + 1
      IF (.NOT.(M_J >= M_B)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 317
      M_B = 123
      M_J = M_B - 1
      IF (.NOT.(M_B >= M_J)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 318
      M_K = 123
      M_J = M_K + 1
      IF (.NOT.(M_J >= M_K)) THEN
          CALL ERRPRT_IM('  >=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST28

      SUBROUTINE TEST29(NCASE,NERROR,KLOG)

!  Test the derived type < interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type < interface.')")

      NCASE = 319
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'<',M_B)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 320
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'LT',M_B)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 321
      J1 = 123
      M_A = J1 - 2
      IF (.NOT.(M_A < J1)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 322
      J1 = 123
      M_A = J1 + 2
      IF (.NOT.(J1 < M_A)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 323
      J1 = 123
      M_J = J1 - 2
      IF (.NOT.(M_J < J1)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 324
      J1 = 123
      M_J = J1 + 2
      IF (.NOT.(J1 < M_J)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 325
      R1 = 12.3
      M_A = R1 - 2
      IF (.NOT.(M_A < R1)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 326
      R1 = 12.3
      M_A = R1 + 2
      IF (.NOT.(R1 < M_A)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 327
      R1 = 123
      M_J = R1 - 2
      IF (.NOT.(M_J < R1)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 328
      R1 = 123
      M_J = R1 + 2
      IF (.NOT.(R1 < M_J)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 329
      D1 = 12.3
      M_A = D1 - 2
      IF (.NOT.(M_A < D1)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 330
      D1 = 12.3
      M_A = D1 + 2
      IF (.NOT.(D1 < M_A)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 331
      D1 = 123
      M_J = D1 - 2
      IF (.NOT.(M_J < D1)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 332
      D1 = 123
      M_J = D1 + 2
      IF (.NOT.(D1 < M_J)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 333
      M_B = 12.3
      M_A = M_B - 2
      IF (.NOT.(M_A < M_B)) THEN
          CALL ERRPRT_FM('   <  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 334
      M_B = 123
      M_J = M_B - 2
      IF (.NOT.(M_J < M_B)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 335
      M_B = 123
      M_J = M_B + 2
      IF (.NOT.(M_B < M_J)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 336
      M_K = 123
      M_J = M_K - 2
      IF (.NOT.(M_J < M_K)) THEN
          CALL ERRPRT_IM('   <  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST29

      SUBROUTINE TEST30(NCASE,NERROR,KLOG)

!  Test the derived type <= interface.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      LOGICAL FM_COMP

      WRITE (KW,"(/' Testing the derived type <= interface.')")

      NCASE = 337
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'<=',M_B)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 338
      M_A = 123
      M_B = 124
      IF (.NOT.FM_COMP(M_A,'LE',M_B)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_B,'M_B',M_B,'M_B',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 339
      J1 = 123
      M_A = J1 - 2
      IF (.NOT.(M_A <= J1)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 340
      J1 = 123
      M_A = J1 + 2
      IF (.NOT.(J1 <= M_A)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 341
      J1 = 123
      M_J = J1 - 2
      IF (.NOT.(M_J <= J1)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 342
      J1 = 123
      M_J = J1 + 2
      IF (.NOT.(J1 <= M_J)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 343
      R1 = 12.3
      M_A = R1 - 2
      IF (.NOT.(M_A <= R1)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 344
      R1 = 12.3
      M_A = R1 + 2
      IF (.NOT.(R1 <= M_A)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 345
      R1 = 123
      M_J = R1 - 2
      IF (.NOT.(M_J <= R1)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 346
      R1 = 123
      M_J = R1 + 2
      IF (.NOT.(R1 <= M_J)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 347
      D1 = 12.3
      M_A = D1 - 2
      IF (.NOT.(M_A <= D1)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 348
      D1 = 12.3
      M_A = D1 + 2
      IF (.NOT.(D1 <= M_A)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 349
      D1 = 123
      M_J = D1 - 2
      IF (.NOT.(M_J <= D1)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 350
      D1 = 123
      M_J = D1 + 2
      IF (.NOT.(D1 <= M_J)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 351
      M_B = 12.3
      M_A = M_B - 2
      IF (.NOT.(M_A <= M_B)) THEN
          CALL ERRPRT_FM('  <=  ',M_A,'M_A',M_A,'M_A',M_A,'M_A',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 352
      M_B = 123
      M_J = M_B - 2
      IF (.NOT.(M_J <= M_B)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 353
      M_B = 123
      M_J = M_B + 2
      IF (.NOT.(M_B <= M_J)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 354
      M_K = 123
      M_J = M_K - 2
      IF (.NOT.(M_J <= M_K)) THEN
          CALL ERRPRT_IM('  <=  ',M_J,'M_J',M_J,'M_J',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST30

      SUBROUTINE TEST31(NCASE,NERROR,KLOG)

!             Test the '+' arithmetic operator.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type + interface.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0

      NCASE = 355
      MFM3 = J2 + MFM1
      CALL FM_ST2M('131',MFM4)
      CALL FM_ADD(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 356
      MIM3 = J2 + MIM1
      CALL IM_ST2M('131',MIM4)
      CALL IM_ADD(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 357
      MZM3 = J2 + MZM1
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 358
      MFM3 = R2 + MFM1
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ADD(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 359
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_ADD(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = R2 + MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 360
      MZM3 = R2 + MZM1
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 361
      MFM3 = D2 + MFM1
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ADD(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 362
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_ADD(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = D2 + MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 363
      MZM3 = D2 + MZM1
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 364
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 + MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 365
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 + MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 366
      MZM3 = C2 + MZM1
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 367
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 + MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 368
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 + MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 369
      MZM3 = CD2 + MZM1
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 370
      MFM3 = MFM1 + J2
      CALL FM_ST2M('131',MFM4)
      CALL FM_ADD(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 371
      MFM3 = MFM1 + R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ADD(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 372
      MFM3 = MFM1 + D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ADD(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 373
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ADD(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 + C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 374
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 + CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 375
      MFM3 = MFM1 + MFM2
      CALL FM_ADD(MFM1,MFM2,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 376
      MFM3 = MFM1 + MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_ADD(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 377
      MZM3 = MFM1 + MZM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 378
      MIM3 = MIM1 + J2
      CALL IM_ST2M('131',MIM4)
      CALL IM_ADD(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 379
      CALL FM_ST2M('241.21',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_ADD(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 + R2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 380
      CALL FM_ST2M('391.61',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_ADD(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 + D2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 381
      CALL ZM_ST2M('411.11 + 421.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 + C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 382
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_ADD(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 + CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 383
      MFM3 = MIM1 + MFM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_ADD(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 384
      MIM3 = MIM1 + MIM2
      CALL IM_ADD(MIM1,MIM2,MIM4)
      IF (MIM4 /= MIM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 385
      MZM3 = MIM1 + MZM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_ADD(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 386
      MZM3 = MZM1 + J2
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 387
      MZM3 = MZM1 + R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 388
      MZM3 = MZM1 + D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 389
      MZM3 = MZM1 + C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 390
      MZM3 = MZM1 + CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 391
      MZM3 = MZM1 + MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 392
      MZM3 = MZM1 + MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_ADD(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 393
      MZM3 = MZM1 + MZM2
      CALL ZM_ADD(MZM1,MZM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 394
      MFM3 = +MFM1
      CALL FM_EQ(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 395
      MIM3 = +MIM1
      CALL IM_EQ(MIM1,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 396
      MZM3 = +MZM1
      CALL ZM_EQ(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST31

      SUBROUTINE TEST32(NCASE,NERROR,KLOG)

!             Test the '-' arithmetic operator.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type - interface.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0

      NCASE = 397
      MFM3 = J2 - MFM1
      CALL FM_ST2M('131',MFM4)
      CALL FM_SUB(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 398
      MIM3 = J2 - MIM1
      CALL IM_ST2M('131',MIM4)
      CALL IM_SUB(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 399
      MZM3 = J2 - MZM1
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 400
      MFM3 = R2 - MFM1
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_SUB(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 401
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_SUB(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = R2 - MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 402
      MZM3 = R2 - MZM1
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 403
      MFM3 = D2 - MFM1
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_SUB(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 404
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_SUB(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = D2 - MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 405
      MZM3 = D2 - MZM1
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 406
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 - MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 407
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 - MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 408
      MZM3 = C2 - MZM1
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 409
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 - MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 410
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 - MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 411
      MZM3 = CD2 - MZM1
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 412
      MFM3 = MFM1 - J2
      CALL FM_ST2M('131',MFM4)
      CALL FM_SUB(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 413
      MFM3 = MFM1 - R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_SUB(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 414
      MFM3 = MFM1 - D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_SUB(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 415
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 - C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 416
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 - CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 417
      MFM3 = MFM1 - MFM2
      CALL FM_SUB(MFM1,MFM2,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 418
      MFM3 = MFM1 - MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_SUB(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 419
      MZM3 = MFM1 - MZM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 420
      MIM3 = MIM1 - J2
      CALL IM_ST2M('131',MIM4)
      CALL IM_SUB(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 421
      CALL FM_ST2M('241.21',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_SUB(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 - R2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 422
      CALL FM_ST2M('391.61',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_SUB(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 - D2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 423
      CALL ZM_ST2M('411.11 + 421.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 - C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 424
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_SUB(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 - CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 425
      MFM3 = MIM1 - MFM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_SUB(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 426
      MIM3 = MIM1 - MIM2
      CALL IM_SUB(MIM1,MIM2,MIM4)
      IF (MIM4 /= MIM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 427
      MZM3 = MIM1 - MZM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 428
      MZM3 = MZM1 - J2
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 429
      MZM3 = MZM1 - R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 430
      MZM3 = MZM1 - D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 431
      MZM3 = MZM1 - C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 432
      MZM3 = MZM1 - CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 433
      MZM3 = MZM1 - MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 434
      MZM3 = MZM1 - MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_SUB(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 435
      MZM3 = MZM1 - MZM2
      CALL ZM_SUB(MZM1,MZM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 436
      MFM3 = -MFM1
      CALL FM_I2M(0,MFM4)
      CALL FM_SUB(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 437
      MIM3 = -MIM1
      CALL IM_I2M(0,MIM4)
      CALL IM_SUB(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 438
      MZM3 = -MZM1
      CALL ZM_I2M(0,MZM4)
      CALL ZM_SUB(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST32

      SUBROUTINE TEST33(NCASE,NERROR,KLOG)

!             Test the '*' arithmetic operator.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type * interface.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0

      NCASE = 439
      MFM3 = J2 * MFM1
      CALL FM_ST2M('131',MFM4)
      CALL FM_MPY(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 440
      MIM3 = J2 * MIM1
      CALL IM_ST2M('131',MIM4)
      CALL IM_MPY(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 441
      MZM3 = J2 * MZM1
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 442
      MFM3 = R2 * MFM1
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_MPY(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 443
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_MPY(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = R2 * MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 444
      MZM3 = R2 * MZM1
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 445
      MFM3 = D2 * MFM1
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_MPY(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 446
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_MPY(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = D2 * MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 447
      MZM3 = D2 * MZM1
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 448
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 * MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 449
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 * MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 450
      MZM3 = C2 * MZM1
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 451
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 * MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 452
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 * MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 453
      MZM3 = CD2 * MZM1
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 454
      MFM3 = MFM1 * J2
      CALL FM_ST2M('131',MFM4)
      CALL FM_MPY(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 455
      MFM3 = MFM1 * R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_MPY(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 456
      MFM3 = MFM1 * D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_MPY(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 457
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_MPY(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 * C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 458
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 * CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 459
      MFM3 = MFM1 * MFM2
      CALL FM_MPY(MFM1,MFM2,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 460
      MFM3 = MFM1 * MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_MPY(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 461
      MZM3 = MFM1 * MZM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 462
      MIM3 = MIM1 * J2
      CALL IM_ST2M('131',MIM4)
      CALL IM_MPY(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 463
      CALL FM_ST2M('241.21',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_MPY(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 * R2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 464
      CALL FM_ST2M('391.61',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_MPY(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 * D2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 465
      CALL ZM_ST2M('411.11 + 421.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 * C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 466
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_MPY(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 * CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 467
      MFM3 = MIM1 * MFM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_MPY(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 468
      MIM3 = MIM1 * MIM2
      CALL IM_MPY(MIM1,MIM2,MIM4)
      IF (MIM4 /= MIM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 469
      MZM3 = MIM1 * MZM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_MPY(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 470
      MZM3 = MZM1 * J2
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 471
      MZM3 = MZM1 * R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 472
      MZM3 = MZM1 * D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 473
      MZM3 = MZM1 * C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 474
      MZM3 = MZM1 * CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 475
      MZM3 = MZM1 * MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 476
      MZM3 = MZM1 * MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_MPY(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 477
      MZM3 = MZM1 * MZM2
      CALL ZM_MPY(MZM1,MZM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST33

      SUBROUTINE TEST34(NCASE,NERROR,KLOG)

!             Test the '/' arithmetic operator.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type / interface.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0

      NCASE = 478
      MFM3 = J2 / MFM1
      CALL FM_ST2M('131',MFM4)
      CALL FM_DIV(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 479
      MIM3 = J2 / MIM1
      CALL IM_ST2M('131',MIM4)
      CALL IM_DIV(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 480
      MZM3 = J2 / MZM1
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 481
      MFM3 = R2 / MFM1
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_DIV(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 482
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = R2 / MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 483
      MZM3 = R2 / MZM1
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 484
      MFM3 = D2 / MFM1
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_DIV(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 485
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = D2 / MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 486
      MZM3 = D2 / MZM1
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 487
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 / MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 488
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 / MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 489
      MZM3 = C2 / MZM1
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 490
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 / MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 491
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 / MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 492
      MZM3 = CD2 / MZM1
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 493
      MFM3 = MFM1 / J2
      CALL FM_ST2M('131',MFM4)
      CALL FM_DIV(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 494
      MFM3 = MFM1 / R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_DIV(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 495
      MFM3 = MFM1 / D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_DIV(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 496
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_DIV(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 / C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 497
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 / CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 498
      MFM3 = MFM1 / MFM2
      CALL FM_DIV(MFM1,MFM2,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 499
      MFM3 = MFM1 / MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_DIV(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 500
      MZM3 = MFM1 / MZM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 501
      MIM3 = MIM1 / J2
      CALL IM_ST2M('131',MIM4)
      CALL IM_DIV(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 502
      CALL FM_ST2M('241.21',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 / R2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 503
      CALL FM_ST2M('391.61',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 / D2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 504
      CALL ZM_ST2M('411.11 + 421.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 / C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 505
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 / CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 506
      MFM3 = MIM1 / MFM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_DIV(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 507
      MIM3 = MIM1 / MIM2
      CALL IM_DIV(MIM1,MIM2,MIM4)
      IF (MIM4 /= MIM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 508
      MZM3 = MIM1 / MZM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_DIV(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 509
      MZM3 = MZM1 / J2
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 510
      MZM3 = MZM1 / R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 511
      MZM3 = MZM1 / D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 512
      MZM3 = MZM1 / C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 513
      MZM3 = MZM1 / CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 514
      MZM3 = MZM1 / MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 515
      MZM3 = MZM1 / MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_DIV(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 516
      MZM3 = MZM1 / MZM2
      CALL ZM_DIV(MZM1,MZM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST34

      SUBROUTINE TEST35(NCASE,NERROR,KLOG)

!             Test the '**' arithmetic operator.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type ** interface.')")

!             Use a larger error tolerance for large exponents.

      RSMALL = EPSILON(1.0)*10000.0
      DSMALL = EPSILON(1.0D0)*10000.0

      NCASE = 517
      MFM3 = J2 ** MFM1
      CALL FM_ST2M('131',MFM4)
      CALL FM_PWR(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 518
      J4 = 2
      MIM3 = J4 ** MIM1
      CALL IM_ST2M('2',MIM4)
      CALL IM_PWR(MIM4,MIM1,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 519
      MZM3 = J2 ** MZM1
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 520
      MFM3 = R2 ** MFM1
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_PWR(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)


      NCASE = 521
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_PWR(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = R2 ** MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 522
      MZM3 = R2 ** MZM1
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 523
      MFM3 = D2 ** MFM1
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_PWR(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 524
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_ST2M('661',MFM3)
      CALL FM_PWR(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = D2 ** MIM1
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 525
      MZM3 = D2 ** MZM1
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 526
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 ** MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 527
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = C2 ** MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 528
      MZM3 = C2 ** MZM1
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 529
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 ** MFM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 530
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_ST2M('661',MZM3)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = CD2 ** MIM1
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 531
      MZM3 = CD2 ** MZM1
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 532
      MFM3 = MFM1 ** J2
      CALL FM_ST2M('131',MFM4)
      CALL FM_PWR(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 533
      MFM3 = MFM1 ** R2
      CALL FM_ST2M('241.21',MFM4)
      CALL FM_PWR(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 534
      MFM3 = MFM1 ** D2
      CALL FM_ST2M('391.61',MFM4)
      CALL FM_PWR(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 535
      CALL ZM_ST2M('581.21',MZM3)
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_PWR(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 ** C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 536
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MFM1 ** CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 537
      MFM3 = MFM1 ** MFM2
      CALL FM_PWR(MFM1,MFM2,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 538
      MFM3 = MFM1 ** MIM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_PWR(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 539
      MZM3 = MFM1 ** MZM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 540
      J4 = 17
      MIM3 = MIM1 ** J4
      CALL IM_ST2M('17',MIM4)
      CALL IM_PWR(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 541
      CALL FM_ST2M('241.21',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_PWR(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 ** R2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 542
      CALL FM_ST2M('391.61',MFM3)
      CALL FM_ST2M('661',MFM4)
      CALL FM_PWR(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = MIM1 ** D2
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 543
      CALL ZM_ST2M('411.11 + 421.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 ** C2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 544
      CALL ZM_ST2M('431.11 + 441.21 i',MZM3)
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_PWR(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      MZM3 = MIM1 ** CD2
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 545
      MFM3 = MIM1 ** MFM1
      CALL FM_ST2M('661',MFM4)
      CALL FM_PWR(MFM4,MFM1,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM4 /= MFM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 546
      MIM4 = 19
      MIM3 = MIM1 ** MIM4
      CALL IM_PWR(MIM1,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM4 /= MIM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 547
      MZM3 = MIM1 ** MZM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_PWR(MZM4,MZM1,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 548
      MZM3 = MZM1 ** J2
      CALL ZM_ST2M('131',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 549
      MZM3 = MZM1 ** R2
      CALL ZM_ST2M('241.21',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 550
      MZM3 = MZM1 ** D2
      CALL ZM_ST2M('391.61',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 551
      MZM3 = MZM1 ** C2
      CALL ZM_ST2M('411.11 + 421.21 i',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > RSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 552
      MZM3 = MZM1 ** CD2
      CALL ZM_ST2M('431.11 + 441.21 i',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      IF (MFM4 > DSMALL) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 553
      MZM3 = MZM1 ** MFM1
      CALL ZM_ST2M('581.21',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 554
      MZM3 = MZM1 ** MIM1
      CALL ZM_ST2M('661',MZM4)
      CALL ZM_PWR(MZM1,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 555
      MZM3 = MZM1 ** MZM2
      CALL ZM_PWR(MZM1,MZM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM4 /= MZM3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST35

      SUBROUTINE TEST36(NCASE,NERROR,KLOG)

!             Test functions ABS, ..., CEILING.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER J,JERR,KLOG,NERROR,NCASE

      WRITE (KW,"(/' Testing the derived type ABS, ..., CEILING interfaces.')")

      NCASE = 556
      MFM3 = ABS(MFM1)
      CALL FM_ABS(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 557
      MIM3 = ABS(MIM1)
      CALL IM_ABS(MIM1,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 558
      MFM3 = ABS(MZM1)
      CALL ZM_ABS(MZM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 559
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = ACOS(MFM4)
      CALL FM_ACOS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 560
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = ACOS(MZM4)
      CALL ZM_ACOS(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 561
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MFM3 = AIMAG(MZM4)
      CALL ZM_IMAG(MZM4,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 562
      MFM3 = AINT(MFM1)
      CALL FM_INT(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 563
      MZM3 = AINT(MZM1)
      CALL ZM_INT(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 564
      MFM3 = ANINT(MFM1)
      CALL FM_NINT(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 565
      MZM3 = ANINT(MZM1)
      CALL ZM_NINT(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 566
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = ASIN(MFM4)
      CALL FM_ASIN(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 567
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = ASIN(MZM4)
      CALL ZM_ASIN(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 568
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = ATAN(MFM4)
      CALL FM_ATAN(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 569
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = ATAN(MZM4)
      CALL ZM_ATAN(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 570
      MFM3 = ATAN2(MFM1,MFM2)
      CALL FM_ATN2(MFM1,MFM2,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 571
      JERR = -1
      DO J = 0, 10
         IF (BTEST(661,J)) THEN
             IF (.NOT.BTEST(MIM1,J)) JERR = J
         ELSE
             IF (BTEST(MIM1,J)) JERR = J
         ENDIF
      ENDDO
      IF (JERR >= 0) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 572
      CALL FM_ST2M('12.37654',MFM4)
      MFM3 = CEILING(MFM4)
      CALL FM_ST2M('13',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 573
      CALL FM_ST2M('-12.7654',MFM4)
      MFM3 = CEILING(MFM4)
      CALL FM_ST2M('-12',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 574
      CALL ZM_ST2M('12.37654 - 22.54 i',MZM4)
      MZM3 = CEILING(MZM4)
      CALL ZM_ST2M('13 - 22 i',MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 575
      CALL ZM_ST2M('-12.7654 + 22.31 i',MZM4)
      MZM3 = CEILING(MZM4)
      CALL ZM_ST2M('-12 + 23 i',MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST36

      SUBROUTINE TEST37(NCASE,NERROR,KLOG)

!             Test functions CMPLX, ..., EXPONENT.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER J,KLOG,NERROR,NCASE

      WRITE (KW,"(/"//  &
      "' Testing the derived type CMPLX, ..., EXPONENT interfaces.')")

      NCASE = 576
      MZM3 = CMPLX(MFM1,MFM2)
      CALL ZM_CMPX(MFM1,MFM2,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 577
      MZM3 = CMPLX(MIM1,MIM2)
      CALL IM_I2FM(MIM1,MFM3)
      CALL IM_I2FM(MIM2,MFM4)
      CALL ZM_CMPX(MFM3,MFM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 578
      MZM3 = CMPLX(MFM1)
      CALL FM_I2M(0,MFM4)
      CALL ZM_CMPX(MFM1,MFM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 579
      MZM3 = CMPLX(MIM1)
      CALL IM_I2FM(MIM1,MFM3)
      CALL FM_I2M(0,MFM4)
      CALL ZM_CMPX(MFM3,MFM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 580
      MZM3 = CONJG(MZM1)
      CALL ZM_CONJ(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 581
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = COS(MFM4)
      CALL FM_COS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 582
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = COS(MZM4)
      CALL ZM_COS(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 583
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = COSH(MFM4)
      CALL FM_COSH(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 584
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = COSH(MZM4)
      CALL ZM_COSH(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 585
      MFM3 = DBLE(MFM1)
      CALL FM_EQ(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 586
      MFM3 = DBLE(MIM1)
      CALL IM_I2FM(MIM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 587
      MFM3 = DBLE(MZM1)
      CALL ZM_REAL(MZM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 588
      J = DIGITS(MFM1)
      IF (J /= NDIG) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 589
      J = DIGITS(MIM1)
      IF (J /= NDIGMX) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 590
      J = DIGITS(MZM1)
      IF (J /= NDIG) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 591
      MFM3 = DIM(MFM1,MFM2)
      CALL FM_DIM(MFM1,MFM2,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 592
      MIM3 = DIM(MIM1,MIM2)
      CALL IM_DIM(MIM1,MIM2,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 593
      MFM3 = DINT (MFM1)
      CALL FM_INT(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 594
      MZM3 = DINT (MZM1)
      CALL ZM_INT(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 595
      CALL FM_ST2M('1.23',MFMV1(1))
      CALL FM_ST2M('2.23',MFMV1(2))
      CALL FM_ST2M('3.23',MFMV1(3))
      CALL FM_ST2M('4.23',MFMV2(1))
      CALL FM_ST2M('5.23',MFMV2(2))
      CALL FM_ST2M('6.23',MFMV2(3))
      MFM3 = DOTPRODUCT(MFMV1,MFMV2)
      MFM4 = 0
      DO J = 1, 3
         MFM4 = MFM4 + MFMV1(J)*MFMV2(J)
      ENDDO
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 596
      CALL IM_ST2M('12',MIMV1(1))
      CALL IM_ST2M('23',MIMV1(2))
      CALL IM_ST2M('34',MIMV1(3))
      CALL IM_ST2M('-14',MIMV2(1))
      CALL IM_ST2M('-5',MIMV2(2))
      CALL IM_ST2M('16',MIMV2(3))
      MIM3 = DOTPRODUCT(MIMV1,MIMV2)
      MIM4 = 0
      DO J = 1, 3
         MIM4 = MIM4 + MIMV1(J)*MIMV2(J)
      ENDDO
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 597
      CALL ZM_ST2M('1.23 + 1.67 i',MZMV1(1))
      CALL ZM_ST2M('2.23 - 2.56 i',MZMV1(2))
      CALL ZM_ST2M('3.23 + 3.45 i',MZMV1(3))
      CALL ZM_ST2M('4.23 - 4.34 i',MZMV2(1))
      CALL ZM_ST2M('5.23 + 5.23 i',MZMV2(2))
      CALL ZM_ST2M('6.23 - 6.12 i',MZMV2(3))
      MZM3 = DOTPRODUCT(MZMV1,MZMV2)
      MZM4 = 0
      DO J = 1, 3
         MZM4 = MZM4 + MZMV1(J)*MZMV2(J)
      ENDDO
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 598
      MFM3 = EPSILON(MFM1)
      CALL FM_I2M(1,MFM4)
      CALL FM_ULP(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 599
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = EXP(MFM4)
      CALL FM_EXP(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 600
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = EXP(MZM4)
      CALL ZM_EXP(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 601
      J = EXPONENT(MFM1)
      IF (J /= INT(MFM1%MFM(1))) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST37

      SUBROUTINE TEST38(NCASE,NERROR,KLOG)

!             Test functions FLOOR, ..., MIN.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER I,J,KLOG,NERROR,NCASE

      WRITE (KW,"(/"//  &
      "' Testing the derived type FLOOR, ..., MIN interfaces.')")

      NCASE = 602
      CALL FM_ST2M('12.37654',MFM4)
      MFM3 = FLOOR(MFM4)
      CALL FM_ST2M('12',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 603
      CALL FM_ST2M('-12.7654',MFM4)
      MFM3 = FLOOR(MFM4)
      CALL FM_ST2M('-13',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 604
      CALL IM_ST2M('12',MIM4)
      MIM3 = FLOOR(MIM4)
      CALL IM_ST2M('12',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 605
      CALL IM_ST2M('-123',MIM4)
      MIM3 = FLOOR(MIM4)
      CALL IM_ST2M('-123',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 606
      CALL ZM_ST2M('12.37654 - 22.54 i',MZM4)
      MZM3 = FLOOR(MZM4)
      CALL ZM_ST2M('12 - 23 i',MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 607
      CALL ZM_ST2M('-12.7654 + 22.31 i',MZM4)
      MZM3 = FLOOR(MZM4)
      CALL ZM_ST2M('-13 + 22 i',MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 608
      CALL FM_ST2M('12.37654',MFM4)
      MFM3 = FRACTION(MFM4)
      MFM4%MFM(1) = 0
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 609
      CALL ZM_ST2M('12.37654 - 22.54',MZM4)
      MZM3 = FRACTION(MZM4)
      MZM4%MZM(1) = 0
      MZM4%MZM(KPTIMU+01) = 0
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 610
      MFM3 = HUGE(MFM1)
      CALL FM_BIG(MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 611
      MIM3 = HUGE(MIM1)
      CALL IM_BIG(MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 612
      MZM3 = HUGE(MZM1)
      CALL FM_BIG(MFM4)
      CALL ZM_CMPX(MFM4,MFM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 613
      MIM3 = INT(MFM1)
      CALL FM_INT(MFM1,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 614
      MIM3 = INT(MIM1)
      CALL IM_EQ(MIM1,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 615
      MIM3 = INT(MZM1)
      CALL ZM_INT(MZM1,MZM4)
      CALL ZM_REAL(MZM4,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 616
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = LOG(MFM4)
      CALL FM_LN(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 617
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = LOG(MZM4)
      CALL ZM_LN(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 618
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = LOG10(MFM4)
      CALL FM_LG10(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 619
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = LOG10(MZM4)
      CALL ZM_LG10(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 620
      DO I = 1, 3
         DO J = 1, 3
            MFMA(I,J) = 3*(J-1) + I
            MFMB(I,J) = 3*(I-1) + J + 10
         ENDDO
      ENDDO
      MFMC = MATMUL(MFMA,MFMB)
      MFM3 = ABS(MFMC(1,1)-186)+ABS(MFMC(1,2)-198)+ABS(MFMC(1,3)-210)+ &
             ABS(MFMC(2,1)-228)+ABS(MFMC(2,2)-243)+ABS(MFMC(2,3)-258)+ &
             ABS(MFMC(3,1)-270)+ABS(MFMC(3,2)-288)+ABS(MFMC(3,3)-306)
      IF (MFM3 /= 0) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 621
      DO I = 1, 2
         DO J = 1, 2
            MIMA(I,J) = 2*(J-1) + I + 20
            MIMB(I,J) = 2*(I-1) + J + 30
         ENDDO
      ENDDO
      MIMC = MATMUL(MIMA,MIMB)
      MIM3 = ABS(MIMC(1,1)-1410) + ABS(MIMC(1,2)-1454) + &
             ABS(MIMC(2,1)-1474) + ABS(MIMC(2,2)-1520)
      IF (MIM3 /= 0) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 622
      DO I = 1, 2
         DO J = 1, 3
            MZMA(I,J) = CMPLX(TO_FM(2*(J-1)+I+10),TO_FM(2*(J-1)+I+20))
         ENDDO
      ENDDO
      DO I = 1, 3
         DO J = 1, 4
            MZMB(I,J) = CMPLX(TO_FM(4*(I-1)+J+50),TO_FM(4*(I-1)+J+30))
         ENDDO
      ENDDO
      MZMC = MATMUL(MZMA,MZMB)
      MFM3 = ABS(MZMC(1,1)-TO_ZM('-270 + 5192 i')) + &
             ABS(MZMC(1,2)-TO_ZM('-300 + 5300 i')) + &
             ABS(MZMC(1,3)-TO_ZM('-330 + 5408 i')) + &
             ABS(MZMC(1,4)-TO_ZM('-360 + 5516 i')) + &
             ABS(MZMC(2,1)-TO_ZM('-210 + 5462 i')) + &
             ABS(MZMC(2,2)-TO_ZM('-240 + 5576 i')) + &
             ABS(MZMC(2,3)-TO_ZM('-270 + 5690 i')) + &
             ABS(MZMC(2,4)-TO_ZM('-300 + 5804 i'))
      IF (MFM3 /= 0) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 623
      MFM3 = MAX(MFM1,MFM2)
      CALL FM_MAX(MFM1,MFM2,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 624
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = MAX(MFM2,MFM1,MFM4)
      CALL FM_MAX(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_MAX(MFM2,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 625
      MIM3 = MAX(MIM1,MIM2)
      CALL IM_MAX(MIM1,MIM2,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 626
      CALL IM_ST2M('7654',MIM4)
      CALL IM_ST2M('-1654',MIM3)
      MIM3 = MAX(MIM2,MIM1,MIM3,MIM4)
      CALL IM_ST2M('7654',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 627
      J = MAXEXPONENT(MFM1)
      IF (J /= INT(MXEXP)+1) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 628
      MFM3 = MIN(MFM1,MFM2)
      CALL FM_MIN(MFM1,MFM2,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 629
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = MIN(MFM2,MFM1,MFM4)
      CALL FM_MIN(MFM1,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_MIN(MFM2,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 630
      MIM3 = MIN(MIM1,MIM2)
      CALL IM_MIN(MIM1,MIM2,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 631
      CALL IM_ST2M('7654',MIM4)
      CALL IM_ST2M('-1654',MIM3)
      MIM3 = MIN(MIM2,MIM1,MIM3,MIM4)
      CALL IM_ST2M('-1654',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST38

      SUBROUTINE TEST39(NCASE,NERROR,KLOG)

!             Test functions MINEXPONENT, ..., RRSPACING.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER J,KLOG,NERROR,NCASE

      WRITE (KW,"(/"//  &
      "' Testing the derived type MINEXPONENT, ..., RRSPACING interfaces.')")

      NCASE = 632
      J = MINEXPONENT(MFM1)
      IF (J /= -INT(MXEXP)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 633
      CALL FM_ST2M('8',MFM3)
      CALL FM_ST2M('5',MFM4)
      MFM3 = MOD(MFM3,MFM4)
      CALL FM_ST2M('3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 634
      CALL FM_ST2M('-8',MFM3)
      CALL FM_ST2M('5',MFM4)
      MFM3 = MOD(MFM3,MFM4)
      CALL FM_ST2M('-3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 635
      CALL FM_ST2M('8',MFM3)
      CALL FM_ST2M('-5',MFM4)
      MFM3 = MOD(MFM3,MFM4)
      CALL FM_ST2M('3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 636
      CALL FM_ST2M('-8',MFM3)
      CALL FM_ST2M('-5',MFM4)
      MFM3 = MOD(MFM3,MFM4)
      CALL FM_ST2M('-3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 637
      CALL IM_ST2M('8',MIM3)
      CALL IM_ST2M('5',MIM4)
      MIM3 = MOD(MIM3,MIM4)
      CALL IM_ST2M('3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 638
      CALL IM_ST2M('-8',MIM3)
      CALL IM_ST2M('5',MIM4)
      MIM3 = MOD(MIM3,MIM4)
      CALL IM_ST2M('-3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 639
      CALL IM_ST2M('8',MIM3)
      CALL IM_ST2M('-5',MIM4)
      MIM3 = MOD(MIM3,MIM4)
      CALL IM_ST2M('3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 640
      CALL IM_ST2M('-8',MIM3)
      CALL IM_ST2M('-5',MIM4)
      MIM3 = MOD(MIM3,MIM4)
      CALL IM_ST2M('-3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 641
      CALL FM_ST2M('8',MFM3)
      CALL FM_ST2M('5',MFM4)
      MFM3 = MODULO(MFM3,MFM4)
      CALL FM_ST2M('3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 642
      CALL FM_ST2M('-8',MFM3)
      CALL FM_ST2M('5',MFM4)
      MFM3 = MODULO(MFM3,MFM4)
      CALL FM_ST2M('2',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 643
      CALL FM_ST2M('8',MFM3)
      CALL FM_ST2M('-5',MFM4)
      MFM3 = MODULO(MFM3,MFM4)
      CALL FM_ST2M('-2',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 644
      CALL FM_ST2M('-8',MFM3)
      CALL FM_ST2M('-5',MFM4)
      MFM3 = MODULO(MFM3,MFM4)
      CALL FM_ST2M('-3',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 645
      CALL IM_ST2M('8',MIM3)
      CALL IM_ST2M('5',MIM4)
      MIM3 = MODULO(MIM3,MIM4)
      CALL IM_ST2M('3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 646
      CALL IM_ST2M('-8',MIM3)
      CALL IM_ST2M('5',MIM4)
      MIM3 = MODULO(MIM3,MIM4)
      CALL IM_ST2M('2',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 647
      CALL IM_ST2M('8',MIM3)
      CALL IM_ST2M('-5',MIM4)
      MIM3 = MODULO(MIM3,MIM4)
      CALL IM_ST2M('-2',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 648
      CALL IM_ST2M('-8',MIM3)
      CALL IM_ST2M('-5',MIM4)
      MIM3 = MODULO(MIM3,MIM4)
      CALL IM_ST2M('-3',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 649
      CALL FM_ST2M('0',MFM4)
      CALL FM_ST2M('1',MFM3)
      CALL FM_BIG(MFM5)
      CALL FM_DIV(MFM3,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM5)
      MFM3 = NEAREST(MFM4,MFM3)
      IF (MFM3 /= MFM5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 650
      CALL FM_ST2M('0',MFM4)
      CALL FM_ST2M('-1',MFM3)
      CALL FM_BIG(MFM5)
      CALL FM_DIV(MFM3,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM5)
      MFM3 = NEAREST(MFM4,MFM3)
      IF (MFM3 /= MFM5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 651
      CALL FM_ST2M('2.345',MFM4)
      CALL FM_ST2M('1',MFM3)
      MFM3 = NEAREST(MFM4,MFM3)
      CALL FM_ULP(MFM4,MFM5)
      CALL FM_ADD(MFM4,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 652
      CALL FM_ST2M('2.345',MFM4)
      CALL FM_ST2M('-1',MFM3)
      MFM3 = NEAREST(MFM4,MFM3)
      CALL FM_ULP(MFM4,MFM5)
      CALL FM_SUB(MFM4,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 653
      CALL FM_ST2M('1',MFM4)
      CALL FM_ST2M('-1',MFM3)
      MFM3 = NEAREST(MFM4,MFM3)
      CALL FM_ST2M('0.99',MFM5)
      CALL FM_ULP(MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM5)
      CALL FM_SUB(MFM4,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 654
      CALL FM_ST2M('-1',MFM4)
      CALL FM_ST2M('12',MFM3)
      MFM3 = NEAREST(MFM4,MFM3)
      CALL FM_ST2M('-0.99',MFM5)
      CALL FM_ULP(MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM5)
      CALL FM_SUB(MFM4,MFM5,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 655
      MIM3 = NINT(MFM1)
      CALL FM_NINT(MFM1,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 656
      MIM3 = NINT(MIM1)
      CALL IM_EQ(MIM1,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 657
      MIM3 = NINT(MZM1)
      CALL ZM_NINT(MZM1,MZM4)
      CALL ZM_REAL(MZM4,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 658
      J = PRECISION(MFM1)
      IF (J /= INT(LOG10(REAL(MBASE))*(NDIG-1) + 1)) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 659
      J = PRECISION(MZM1)
      IF (J /= INT(LOG10(REAL(MBASE))*(NDIG-1) + 1)) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 660
      J = RADIX(MFM1)
      IF (J /= INT(MBASE)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 661
      J = RADIX(MIM1)
      IF (J /= INT(MBASE)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 662
      J = RADIX(MZM1)
      IF (J /= INT(MBASE)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 663
      J = RANGE(MFM1)
      IF (J /= INT(MXEXP*LOG10(REAL(MBASE)))) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 664
      J = RANGE(MIM1)
      IF (J /= INT(NDIGMX*LOG10(REAL(MBASE)))) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 665
      J = RANGE(MZM1)
      IF (J /= INT(MXEXP*LOG10(REAL(MBASE)))) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 666
      MFM3 = REAL(MFM1)
      CALL FM_EQ(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 667
      MFM3 = REAL(MIM1)
      CALL IM_I2FM(MIM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 668
      MFM3 = REAL(MZM1)
      CALL ZM_REAL(MZM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 669
      MFM3 = RRSPACING(MFM1)
      CALL FM_ABS(MFM1,MFM4)
      MFM4%MFM(1) = NDIG
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST39

      SUBROUTINE TEST40(NCASE,NERROR,KLOG)

!             Test functions SCALE, ..., TINY.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/"//  &
      "' Testing the derived type SCALE, ..., TINY interfaces.')")

      NCASE = 670
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SCALE(MFM4,1)
      CALL FM_MPYI(MFM4,INT(MBASE),MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 671
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = SCALE(MZM4,-2)
      CALL ZM_DIVI(MZM4,INT(MBASE),MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIVI(MZM4,INT(MBASE),MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 672
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SETEXPONENT(MFM4,1)
      CALL FM_MPYI(MFM4,INT(MBASE),MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 673
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SIGN(MFM4,MFM2)
      CALL FM_SIGN(MFM4,MFM2,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 674
      CALL IM_ST2M('231',MIM4)
      MIM3 = SIGN(MIM4,MIM2)
      CALL IM_SIGN(MIM4,MIM2,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 675
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SIN(MFM4)
      CALL FM_SIN(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 676
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = SIN(MZM4)
      CALL ZM_SIN(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 677
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SINH(MFM4)
      CALL FM_SINH(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 678
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = SINH(MZM4)
      CALL ZM_SINH(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 679
      CALL FM_ST2M('-0.7654',MFM4)
      MFM3 = SPACING(MFM4)
      CALL FM_ULP(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 680
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = SQRT(MFM4)
      CALL FM_SQRT(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 681
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = SQRT(MZM4)
      CALL ZM_SQRT(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 682
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = TAN(MFM4)
      CALL FM_TAN(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 683
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = TAN(MZM4)
      CALL ZM_TAN(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 684
      CALL FM_ST2M('0.7654',MFM4)
      MFM3 = TANH(MFM4)
      CALL FM_TANH(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 685
      CALL ZM_ST2M('0.7654 - 0.3456 i',MZM4)
      MZM3 = TANH(MZM4)
      CALL ZM_TANH(MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 686
      CALL FM_BIG(MFM4)
      CALL FM_I2M(1,MFM3)
      CALL FM_DIV(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = TINY(MFM1)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 687
      MIM3 = TINY(MIM1)
      CALL IM_I2M(1,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 688
      CALL FM_BIG(MFM4)
      CALL FM_I2M(1,MFM3)
      CALL FM_DIV(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL ZM_CMPX(MFM4,MFM4,MZM4)
      MZM3 = TINY(MZM1)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST40

      SUBROUTINE TEST41(NCASE,NERROR,KLOG)

!             Test functions TO_FM, TO_IM, TO_ZM, ..., TO_DPZ.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      INTEGER KLOG,NERROR,NCASE
      LOGICAL FM_COMP

      WRITE (KW,"(/"//  &
      "' Testing the derived type TO_FM,  ..., TO_DPZ interfaces.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0

      NCASE = 689
      MFM3 = TO_FM(123)
      CALL FM_I2M(123,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 690
      MFM3 = TO_FM(123.4)
      CALL FM_SP2M(123.4,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 691
      MFM3 = TO_FM(123.45D0)
      CALL FM_DP2M(123.45D0,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)


      NCASE = 692
      MFM3 = TO_FM(CMPLX(123.4,567.8))
      CALL FM_SP2M(123.4,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 693
      MFM3 = TO_FM(CMPLX(123.4D0,567.8D0,KIND(1.0D0)))
      CALL FM_DP2M(123.4D0,MFM4)
      CALL FM_SUB(MFM3,MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_DIV(MFM4,MFM3,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      CALL FM_ABS(MFM4,MFM6)
      CALL FM_EQ(MFM6,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 694
      MFM3 = TO_FM(MFM1)
      CALL FM_EQ(MFM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 695
      MFM3 = TO_FM(MIM1)
      CALL IM_I2FM(MIM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 696
      MFM3 = TO_FM(MZM1)
      CALL ZM_REAL(MZM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 697
      MFM3 = TO_FM('-123.654')
      CALL FM_ST2M('-123.654',MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 698
      MIM3 = TO_IM(123)
      CALL IM_I2M(123,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 699
      MIM3 = TO_IM(123.4)
      CALL FM_SP2M(123.4,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 700
      MIM3 = TO_IM(123.45D0)
      CALL FM_DP2M(123.45D0,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 701
      MIM3 = TO_IM(CMPLX(123.4,567.8))
      CALL FM_SP2M(123.4,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 702
      MIM3 = TO_IM(CMPLX(123.4D0,567.8D0,KIND(1.0D0)))
      CALL FM_DP2M(123.4D0,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 703
      MIM3 = TO_IM(MFM1)
      CALL FM_EQ(MFM1,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 704
      MIM3 = TO_IM(MIM1)
      CALL IM_I2FM(MIM1,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 705
      MIM3 = TO_IM(MZM1)
      CALL ZM_REAL(MZM1,MFM4)
      CALL IM_FM2I(MFM4,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 706
      MIM3 = TO_IM('-123654')
      CALL IM_ST2M('-123654',MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 707
      MZM3 = TO_ZM(123)
      CALL ZM_I2M(123,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 708
      MZM3 = TO_ZM(123.4)
      CALL FM_SP2M(123.4,MFM4)
      CALL FM_I2M(0,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 709
      MZM3 = TO_ZM(123.45D0)
      CALL FM_DP2M(123.45D0,MFM4)
      CALL FM_I2M(0,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 710
      MZM3 = TO_ZM(CMPLX(123.4,567.8))
      CALL FM_SP2M(123.4,MFM4)
      CALL FM_SP2M(567.8,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = RSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 711
      MZM3 = TO_ZM(CMPLX(123.4D0,567.8D0,KIND(1.0D0)))
      CALL FM_DP2M(123.4D0,MFM4)
      CALL FM_DP2M(567.8D0,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      CALL ZM_SUB(MZM3,MZM4,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_DIV(MZM4,MZM3,MZM5)
      CALL ZM_EQ(MZM5,MZM4)
      CALL ZM_ABS(MZM4,MFM4)
      MFM3 = DSMALL
      IF (FM_COMP(MFM4,'GT',MFM3)) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 712
      MZM3 = TO_ZM(MFM1)
      CALL FM_EQ(MFM1,MFM4)
      CALL FM_I2M(0,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 713
      MZM3 = TO_ZM(MIM1)
      CALL IM_I2FM(MIM1,MFM4)
      CALL FM_I2M(0,MFM5)
      CALL ZM_CMPX(MFM4,MFM5,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 714
      MZM3 = TO_ZM(MZM1)
      CALL ZM_EQ(MZM1,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 715
      MZM3 = TO_ZM('-123.654 + 98.7 i')
      CALL ZM_ST2M('-123.654 + 98.7 i',MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 716
      CALL FM_M2I(MFM1,J3)
      IF (TO_INT(MFM1) /= J3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 717
      CALL IM_M2I(MIM1,J3)
      IF (TO_INT(MIM1) /= J3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 718
      CALL ZM_M2I(MZM1,J3)
      IF (TO_INT(MZM1) /= J3) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 719
      CALL FM_M2SP(MFM1,R3)
      IF (ABS((TO_SP(MFM1)-R3)/R3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 720
      CALL IM_M2DP(MIM1,D3)
      R3 = D3
      IF (ABS((TO_SP(MIM1)-R3)/R3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 721
      CALL ZM_REAL(MZM1,MFM4)
      CALL FM_M2SP(MFM4,R3)
      IF (ABS((TO_SP(MZM1)-R3)/R3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 722
      CALL FM_M2DP(MFM1,D3)
      IF (ABS((TO_DP(MFM1)-D3)/D3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 723
      CALL IM_M2DP(MIM1,D3)
      IF (ABS((TO_DP(MIM1)-D3)/D3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 724
      CALL ZM_REAL(MZM1,MFM4)
      CALL FM_M2DP(MFM4,D3)
      IF (ABS((TO_DP(MZM1)-D3)/D3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 725
      CALL FM_M2SP(MFM1,R3)
      C3 = R3
      IF (ABS((TO_SPZ(MFM1)-C3)/C3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 726
      CALL IM_M2DP(MIM1,D3)
      C3 = D3
      IF (ABS((TO_SPZ(MIM1)-C3)/C3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 727
      CALL ZM_M2Z(MZM1,C3)
      IF (ABS((TO_SPZ(MZM1)-C3)/C3) > RSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 728
      CALL FM_M2DP(MFM1,D3)
      CD3 = D3
      IF (ABS((TO_DPZ(MFM1)-CD3)/CD3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 729
      CALL IM_M2DP(MIM1,D3)
      CD3 = D3
      IF (ABS((TO_DPZ(MIM1)-CD3)/CD3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      NCASE = 730
      CALL ZM_REAL(MZM1,MFM4)
      CALL FM_M2DP(MFM4,D3)
      CALL ZM_IMAG(MZM1,MFM4)
      CALL FM_M2DP(MFM4,D4)
      CD3 = CMPLX( D3 , D4 , KIND(0.0D0) )
      IF (ABS((TO_DPZ(MZM1)-CD3)/CD3) > DSMALL) THEN
          CALL PRTERR(KW,KLOG,NCASE,NERROR)
      ENDIF

      END SUBROUTINE TEST41

      SUBROUTINE TEST42(NCASE,NERROR,KLOG)

!             Test the derived-type interface routines that are not
!             used elsewhere in this program.

      USE FMVALS
      USE FMZM
      USE TEST_VARS
      IMPLICIT NONE

      CHARACTER(80) :: STRING
      INTEGER KLOG,NERROR,NCASE

      WRITE (KW,"(/"//  &
      "' Testing the derived type ADDI, ..., Z2M interfaces.')")

      RSMALL = EPSILON(1.0)*100.0
      DSMALL = EPSILON(1.0D0)*100.0
      MSMALL = EPSILON(TO_FM(1))*10000.0

      NCASE = 731
      MFM3 = MFM1 + 123
      MFM4 = MFM1
      CALL FM_ADDI(MFM4,123)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 732
      CALL FM_CHSH(MFM1,MFM4,MFM3)
      MFM3 = COSH(MFM1)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 733
      CALL FM_CHSH(MFM1,MFM3,MFM4)
      MFM3 = SINH(MFM1)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 734
      CALL FM_CSSN(MFM1,MFM4,MFM3)
      MFM3 = COS(MFM1)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 735
      CALL FM_CSSN(MFM1,MFM3,MFM4)
      MFM3 = SIN(MFM1)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 736
      MFM3 = MFM1 / 123
      CALL FM_DIVI(MFM1,123,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 737
      MFM3 = 123.45D0
      CALL FM_DPM(123.45D0,MFM4)
      IF (ABS((MFM3-MFM4)/MFM4) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 738
      CALL FM_FORM('F70.56',MFM1,STRING)
      CALL FM_ST2M(STRING(1:70),MFM4)
      IF (ABS((MFM1-MFM4)/MFM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 739
      STRING = FM_FORMAT('F70.56',MFM1)
      CALL FM_ST2M(STRING(1:70),MFM4)
      IF (ABS((MFM1-MFM4)/MFM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 740
      MFM3 = MFM1 ** 123
      CALL FM_IPWR(MFM1,123,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 741
      MFM3 = LOG(TO_FM(123))
      CALL FM_LNI(123,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 742
      D4 = MFM1
      CALL FM_M2DP(MFM1,D5)
      IF (ABS((D4-D5)/D4) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 743
      J4 = MFM1
      CALL FM_M2I(MFM1,J5)
      IF (J4 /= J5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 744
      R4 = MFM1
      CALL FM_M2SP(MFM1,R5)
      IF (ABS((R4-R5)/R4) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 745
      MFM3 = 2.67
      CALL FM_MOD(MFM1,MFM3,MFM4)
      MFM3 = MOD(MFM1,MFM3)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 746
      CALL FM_PI(MFM4)
      MFM3 = 4*ATAN(TO_FM(1))
      IF (ABS((MFM3-MFM4)/MFM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 747
      MFM3 = MFM1 ** (TO_FM(1)/TO_FM(3))
      CALL FM_RPWR(MFM1,1,3,MFM4)
      IF (ABS((MFM3-MFM4)/MFM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 748
      CALL FM_SQR(MFM1,MFM4)
      MFM3 = MFM1*MFM1
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 749
      MIM3 = MIM1 / 13
      CALL IM_DIVI(MIM1,13,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 750
      MIM3 = 13
      CALL IM_DIVR(MIM1,MIM3,MIM5,MIM4)
      MIM3 = MOD(MIM1,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 751
      MIM3 = 13
      CALL IM_DIVR(MIM1,MIM3,MIM5,MIM4)
      CALL IM_EQ(MIM5,MIM3)
      MIM4 = MIM1 / 13
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 752
      MIM3 = MIM1 / 13
      CALL IM_DVIR(MIM1,13,MIM4,J5)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 753
      J4 = MOD(MIM1,TO_IM(13))
      CALL IM_DVIR(MIM1,13,MIM4,J5)
      IF (J4 /= J5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 754
      CALL IM_FORM('I70',MIM1,STRING)
      CALL IM_ST2M(STRING(1:70),MIM4)
      IF (MIM1 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 755
      STRING = IM_FORMAT('I70',MIM1)
      CALL IM_ST2M(STRING(1:70),MIM4)
      IF (MIM1 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 756
      MIM3 = 40833
      MIM4 = 16042
      CALL IM_GCD(MIM3,MIM4,MIM5)
      CALL IM_EQ(MIM5,MIM4)
      IF (MIM4 /= 13) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 757
      MIM3 = 40833
      MIM4 = 16042
      MIM4 = GCD(MIM3,MIM4)
      IF (MIM4 /= 13) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 758
      D4 = MIM1
      CALL IM_M2DP(MIM1,D5)
      IF (ABS((D4-D5)/D4) > DSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 759
      J4 = MIM1
      CALL IM_M2I(MIM1,J5)
      IF (J4 /= J5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 760
      MIM3 = 6
      CALL IM_MOD(MIM1,MIM3,MIM4)
      MIM3 = MOD(MIM1,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 761
      MIM3 = MIM1 * 123
      CALL IM_MPYI(MIM1,123,MIM4)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 762
      MIM2 = 3141
      MIM3 = 133
      CALL IM_MPYM(MIM1,MIM2,MIM3,MIM4)
      MIM3 = MOD(MIM1*MIM2,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 763
      MIM2 = 3141
      MIM3 = 133
      MIM4 = MULTIPLY_MOD(MIM1,MIM2,MIM3)
      MIM3 = MOD(MIM1*MIM2,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 764
      MIM2 = 31
      MIM3 = 147
      CALL IM_PMOD(MIM1,MIM2,MIM3,MIM4)
      MIM3 = MOD(MIM1**MIM2,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 765
      MIM2 = 31
      MIM3 = 147
      MIM4 = POWER_MOD(MIM1,MIM2,MIM3)
      MIM3 = MOD(MIM1**MIM2,MIM3)
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 766
      CALL IM_SQR(MIM1,MIM4)
      MIM3 = MIM1*MIM1
      IF (MIM3 /= MIM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 767
      MZM3 = MZM1 + 123
      MZM4 = MZM1
      CALL ZM_ADDI(MZM4,123)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 768
      MFM3 = ATAN2(AIMAG(MZM1),REAL(MZM1))
      CALL ZM_ARG(MZM1,MFM4)
      IF (MFM3 /= MFM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 769
      CALL ZM_CHSH(MZM1,MZM4,MZM3)
      MZM3 = COSH(MZM1)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 770
      CALL ZM_CHSH(MZM1,MZM3,MZM4)
      MZM3 = SINH(MZM1)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 771
      CALL ZM_CSSN(MZM1,MZM4,MZM3)
      MZM3 = COS(MZM1)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 772
      CALL ZM_CSSN(MZM1,MZM3,MZM4)
      MZM3 = SIN(MZM1)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 773
      CALL ZM_FORM('F35.26','F35.26',MZM1,STRING)
      CALL ZM_ST2M(STRING(1:75),MZM4)
      IF (ABS((MZM1-MZM4)/MZM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 774
      STRING = ZM_FORMAT('F35.26','F35.26',MZM1)
      CALL ZM_ST2M(STRING(1:75),MZM4)
      IF (ABS((MZM1-MZM4)/MZM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 775
      MZM3 = TO_ZM('123-456i')
      CALL ZM_2I2M(123,-456,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 776
      MZM3 = MZM1 ** 123
      CALL ZM_IPWR(MZM1,123,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 777
      J4 = MZM1
      CALL ZM_M2I(MZM1,J5)
      IF (J4 /= J5) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 778
      C4 = MZM1
      CALL ZM_M2Z(MZM1,C5)
      IF (ABS((C4-C5)/C4) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 779
      MZM3 = MZM1 * 123
      CALL ZM_MPYI(MZM1,123,MZM4)
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 780
      MZM3 = MZM1 ** (TO_ZM(1)/TO_ZM(3))
      CALL ZM_RPWR(MZM1,1,3,MZM4)
      IF (ABS((MZM3-MZM4)/MZM4) > MSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 781
      CALL ZM_SQR(MZM1,MZM4)
      MZM3 = MZM1*MZM1
      IF (MZM3 /= MZM4) CALL PRTERR(KW,KLOG,NCASE,NERROR)

      NCASE = 782
      MZM3 = C2
      CALL ZM_Z2M(C2,MZM4)
      IF (ABS((MZM3-MZM4)/MZM3) > RSMALL) &
          CALL PRTERR(KW,KLOG,NCASE,NERROR)

      END SUBROUTINE TEST42

      SUBROUTINE TEST43(NCASE,NERROR,KLOG)

!  Test Bernoulli numbers, Pochhammer's function, Euler's constant.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR,NDGSAV

      WRITE (KW,"(/' Testing Bernoulli, Pochhammer, Euler.')")

      NCASE = 783
      M_A = 1
      CALL FM_BERN(10,M_A,M_C)
      M_D = TO_FM('7.5757575757575757575757575757575757575757575757575758M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 784
      M_A = 1
      CALL FM_BERN(0,M_A,M_C)
      M_D = TO_FM('1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 785
      M_A = 1
      CALL FM_BERN(1,M_A,M_C)
      M_D = TO_FM('-0.5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 786
      M_A = 1
      CALL FM_BERN(41,M_A,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 787
      M_A = 0
      CALL FM_BERN(52,M_A,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 788
      M_A = TO_FM('.7699115044247787610619469026548672566371681415929204')
      CALL FM_BERN(102,M_A,M_C)
      M_D = TO_FM('5.7022917356035929245914353639470138260075545712953255M+80')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 789
      M_A = TO_FM('.7699115044247787610619469026548672566371681415929204')
      CALL FM_BERN(76,M_A,M_C)
      M_D = TO_FM('-6.3274121765674850311763600458139008604123253720098077M+50')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 790
      M_A = TO_FM('.7699115044247787610619469026548672566371681415929204')
      M_C = BERNOULLI(76)*M_A
      M_D = TO_FM('-6.3274121765674850311763600458139008604123253720098077M+50')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BERN ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 791
      M_A = TO_FM('769.9115044247787610619469026548672566371681415929204')
      CALL FM_POCH(M_A,10,M_C)
      M_D = TO_FM('7.7568981408767238023000514593534249181767332686451635M+28')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 792
      M_A = TO_FM('7699.115044247787610619469026548672566371681415929204')
      CALL FM_POCH(M_A,2222,M_C)
      M_D = TO_FM('1.3306321985792900130409652455318897459921360351317942M+8763')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 793
      M_A = TO_FM('-7')
      CALL FM_POCH(M_A,12,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 794
      M_A = TO_FM('7.7568981408767238023000514593534249181767332686451635M+281')
      CALL FM_POCH(M_A,6,M_C)
      M_D = TO_FM('2.1783543710019819738631136312604490177244818356538937M+1691')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 795
      M_A = TO_FM('7.7568981408767238023000514593534249181767332686451635M-281')
      CALL FM_POCH(M_A,8,M_C)
      M_D = TO_FM('3.9094766630018687963592259355141261587610735673971624M-277')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 796
      M_A = TO_FM('7.7568981408767238023000514593534249181767332686451635M-281')
      CALL FM_POCH(M_A,1,M_C)
      M_D = TO_FM('7.7568981408767238023000514593534249181767332686451635M-281')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 797
      M_A = TO_FM('7.7568981408767238023000514593534249181767332686451635M-281')
      CALL FM_POCH(M_A,0,M_C)
      M_D = TO_FM('1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 798
      M_A = TO_FM('0')
      CALL FM_POCH(M_A,8,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 799
      M_A = TO_FM('769.9115044247787610619469026548672566371681415929204')
      M_C = POCHHAMMER(M_A,10)
      M_D = TO_FM('7.7568981408767238023000514593534249181767332686451635M+28')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' POCH ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 800
      CALL FM_EULR(M_C)
      M_D = TO_FM('.5772156649015328606065120900824024310421593359399236')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' EULR ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 801
      NDGSAV = NDIG
      NDIG = MIN(NDIGMX,INT(1785*DLOGTN/DLOGMB)+2)
      CALL FM_EULR(M_C)
      M_D = TO_FM(  &
      ' .5772156649015328606065120900824024310421593359399235988057672348848677267'// &
      '776646709369470632917467495146314472498070824809605040144865428362241739976'// &
      '449235362535003337429373377376739427925952582470949160087352039481656708532'// &
      '331517766115286211995015079847937450857057400299213547861466940296043254215'// &
      '190587755352673313992540129674205137541395491116851028079842348775872050384'// &
      '310939973613725530608893312676001724795378367592713515772261027349291394079'// &
      '843010341777177808815495706610750101619166334015227893586796549725203621287'// &
      '922655595366962817638879272680132431010476505963703947394957638906572967929'// &
      '601009015125195950922243501409349871228247949747195646976318506676129063811'// &
      '051824197444867836380861749455169892792301877391072945781554316005002182844'// &
      '096053772434203285478367015177394398700302370339518328690001558193988042707'// &
      '411542227819716523011073565833967348717650491941812300040654693142999297779'// &
      '569303100503086303418569803231083691640025892970890985486825777364288253954'// &
      '925873629596133298574739302373438847070370284412920166417850248733379080562'// &
      '754998434590761643167103146710722370021810745044418664759134803669025532458'// &
      '625442225345181387912434573501361297782278288148945909863846006293169471887'// &
      '149587525492366493520473243641097268276160877595088095126208404544477992299'// &
      '157248292516251278427659657083214610298214617951957959095922704208989627971'// &
      '255363217948873764210660607065982561990102880756125199137511678217643619057'// &
      '058440783573501580056077457934213144988500786415171615194565706170432450750'// &
      '081687052307890937046143066848179164968425491504967243121837838753564894950'// &
      '868454102340601622508515583867234944187880440940770106883795111307872023426'// &
      '395226920971608856908382511378712836820491178925944784861991185293910293099'// &
      '059255266917274468920443869711147174571574573203935209122316085086828')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= MAX(TO_FM('1.0E-1785'),10*EPSILON(M_C)))) THEN
          CALL ERRPRT_FM(' EULR ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF
      NDIG = NDGSAV

      RETURN
      END SUBROUTINE TEST43

      SUBROUTINE TEST44(NCASE,NERROR,KLOG)

!  Test Gamma, Factorial, Log(Gamma), Beta, Binomial.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing Gamma, Factorial, Log(Gamma), Beta, Binomial.')")

      NCASE = 802
      M_A = 19
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('6.402373705728M+15')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 803
      M_A = TO_FM('.7699115044247787610619469026548672566371681415929204')
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('1.1998023858495967876496039855917100290498970370440326')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 804
      M_A = TO_FM('.7699115044247787610619469026548672566371681415929204')
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('1.1998023858495967876496039855917100290498970370440326')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 805
      M_A = TO_FM('5.7699115044247787610619469026548672566371681415929204')
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('8.1434691207877806133071511233406796488474685081500979M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 806
      M_A = TO_FM('7.7568981408767238023000514593534249181767332686451635M-281')
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('1.2891751081921193691625844770542239587773115818085396M+280')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 807
      M_A = TO_FM('2')
      CALL FM_GAM(M_A,M_C)
      M_D = TO_FM('1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 808
      M_A = TO_FM('5.7699115044247787610619469026548672566371681415929204')
      M_C = GAMMA(M_A)
      M_D = TO_FM('8.1434691207877806133071511233406796488474685081500979M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM('  GAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 809
      M_A = 33
      CALL FM_FACT(M_A,M_C)
      M_D = TO_FM('8.68331761881188649551819440128M+36')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' FACT ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 810
      M_A = TO_FM('769.9115044247787610619469026548672566371681415929204')
      CALL FM_FACT(M_A,M_C)
      M_D = TO_FM('5.9982590033571347622193071279165294725603013413394492M+1889')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' FACT ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 811
      M_A = TO_FM('769.9115044247787610619469026548672566371681415929204')
      M_C = FACTORIAL(M_A)
      M_D = TO_FM('5.9982590033571347622193071279165294725603013413394492M+1889')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' FACT ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 812
      M_A = TO_FM('1.0M-222')
      CALL FM_LNGM(M_A,M_C)
      M_D = TO_FM('5.1117389064467814185199410293992885408744453047558760M+2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' LNGM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 813
      M_A = TO_FM('2')
      CALL FM_LNGM(M_A,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' LNGM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 814
      M_A = TO_FM('33')
      CALL FM_LNGM(M_A,M_C)
      M_D = TO_FM('8.1557959456115037178502968666011206687099284403417368M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' LNGM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 815
      M_A = TO_FM('2.00000000000000000001')
      CALL FM_LNGM(M_A,M_C)
      M_D = TO_FM('4.2278433509846713939671258025183870114019600466320121M-21')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' LNGM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 816
      M_C = LOG_GAMMA(TO_FM('33'))
      M_D = TO_FM('8.1557959456115037178502968666011206687099284403417368M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' LNGM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 817
      M_A = TO_FM('2.0706137739520290320140007735608464643737932737070189M-223')
      M_B = TO_FM('.78')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('4.8294858876137637017880452468052846823385248996130407M+222')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 818
      M_A = TO_FM('.78')
      M_B = TO_FM('2.0706137739520290320140007735608464643737932737070189M-223')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('4.8294858876137637017880452468052846823385248996130407M+222')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 819
      M_A = TO_FM('-4.5')
      M_B = TO_FM('4.5')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 820
      M_A = TO_FM('-5.5')
      M_B = TO_FM('4.5')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 821
      M_A = TO_FM('10')
      M_B = TO_FM('4')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('3.4965034965034965034965034965034965034965034965034965M-4')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 822
      M_A = TO_FM('1.0M+1234')
      M_B = TO_FM('2.2')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('1.7462392672319547876554292922652110015806932440139209M-2715')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 823
      M_A = TO_FM('10')
      M_B = TO_FM('5.3')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('7.0836036771097107530120640698518155187687458162734679M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 824
      M_A = TO_FM('10.3')
      M_B = TO_FM('5')
      CALL FM_BETA(M_A,M_B,M_C)
      M_D = TO_FM('8.8146035423244390793072072569173028531206477712519934M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 825
      M_A = TO_FM('10.3')
      M_B = TO_FM('5')
      M_C = BETA(M_A,M_B)
      M_D = TO_FM('8.8146035423244390793072072569173028531206477712519934M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' BETA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 826
      M_A = TO_FM('12.5')
      M_B = TO_FM('0')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 827
      M_A = TO_FM('5')
      M_B = TO_FM('-2')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 828
      M_A = TO_FM('12.5')
      M_B = TO_FM('12.5')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 829
      M_A = TO_FM('-4.5')
      M_B = TO_FM('4.5')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 830
      M_A = TO_FM('-4.5')
      M_B = TO_FM('4.5')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 831
      M_A = TO_FM('-10')
      M_B = TO_FM('3')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('-220')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 832
      M_A = TO_FM('52')
      M_B = TO_FM('5')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('2.59896M+6')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 833
      M_A = TO_FM('1.0M+1234')
      M_B = TO_FM('7')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('1.9841269841269841269841269841269841269841269841269841M+8634')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 834
      M_A = TO_FM('1.0M+123')
      M_B = TO_FM('2.2')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('1.6423797032130683531106846289429264567307029528308099M+270')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 835
      M_A = TO_FM('1.0M-100')
      M_B = TO_FM('4')
      CALL FM_COMB(M_A,M_B,M_C)
      M_D = TO_FM('-2.5M-101')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 836
      M_A = TO_FM('1.0M+123')
      M_B = TO_FM('2.2')
      M_C = BINOMIAL(M_A,M_B)
      M_D = TO_FM('1.6423797032130683531106846289429264567307029528308099M+270')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' COMB ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST44

      SUBROUTINE TEST45(NCASE,NERROR,KLOG)

!  Test Incomplete Gamma, Incomplete Beta.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing Incomplete Gamma, Incomplete Beta.')")

      NCASE = 837
      M_A = TO_FM('2.0706137739520290320140007735608464643737932737070189M-145')
      M_B = TO_FM('.34')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('4.8294858876137637017880452468052846823385248996130407M+144')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 838
      M_A = TO_FM('1.0E-50')
      M_B = TO_FM('1.0E+555')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('9.9999999999999999999999999999999999999999999999999423M+49')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 839
      M_A = TO_FM('1.2')
      M_B = TO_FM('2.3')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('7.9163089830797686672658085698101181778608009481363580M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 840
      M_A = TO_FM('23.4')
      M_B = TO_FM('456.7')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('3.9191215305400046110416169991395759293572844563673750M+21')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 841
      M_A = TO_FM('1.2')
      M_B = TO_FM('0')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 842
      M_A = TO_FM('-1234.5')
      M_B = TO_FM('3.4')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('-2.0892439131810030556730824779643382797767198269736235M-661')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 843
      M_A = TO_FM('10.3')
      M_B = TO_FM('230.7')
      CALL FM_IGM1(M_A,M_B,M_C)
      M_D = TO_FM('7.1643068906237524454762965471616445342244699109269471M+5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 844
      M_A = TO_FM('1.2')
      M_B = TO_FM('2.3')
      M_C = INCOMPLETE_GAMMA1(M_A,M_B)
      M_D = TO_FM('7.9163089830797686672658085698101181778608009481363580M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM1 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 845
      M_A = TO_FM('0')
      M_B = TO_FM('4.5')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('2.0734007547146144328855938695797884889319725701443004M-3')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 846
      M_A = TO_FM('4.5')
      M_B = TO_FM('0')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('1.1631728396567448929144224109426265262108918305803166M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 847
      M_A = TO_FM('1.2')
      M_B = TO_FM('2.3')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('1.2653784409178374391437079820481858290074190484504480M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 848
      M_A = TO_FM('3.4')
      M_B = TO_FM('456.7')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('1.1043526800164195407100289367720949121507981651704628M-192')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 849
      M_A = TO_FM('1.0E-30')
      M_B = TO_FM('40.7')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('5.0619447546123889551107110735110897294460083487536391M-20')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 850
      M_A = TO_FM('-8000.3')
      M_B = TO_FM('1.0e-10')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('1.2499531266327356460522174653022492899665091451890036M+79999')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 851
      M_A = TO_FM('1')
      M_B = TO_FM('-10.7')
      CALL FM_IGM2(M_A,M_B,M_C)
      M_D = TO_FM('4.4355855130297866938628363428602120081387560278336788M+4')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 852
      M_A = TO_FM('1.2')
      M_B = TO_FM('2.3')
      M_C = INCOMPLETE_GAMMA2(M_A,M_B)
      M_D = TO_FM('1.2653784409178374391437079820481858290074190484504480M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IGM2 ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 853
      M_A = TO_FM('0.1')
      M_B = TO_FM('23.4')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('5.8731980918960730463350151650813268739874201571164800M-27')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 854
      M_A = TO_FM('8.115640517330775M-1')
      M_B = TO_FM('2.00853601446773')
      M_C = TO_FM('1.59735792202923')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.0112520048150164306467955877563719782378767062440103M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 855
      M_A = TO_FM('9.01737835258975M-1')
      M_B = TO_FM('2.00853601446773')
      M_C = TO_FM('1.59735792202923')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.2512248738228585976753517954889151150428002974819213M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 856
      M_A = TO_FM('9.6097615596216720E-01')
      M_B = TO_FM('1.970425178583792')
      M_C = TO_FM('5.5680052333367')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.8619456987740165364092968281459448023932520843535423M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 857
      M_A = TO_FM('4.764360371097952E-01')
      M_B = TO_FM('1.161514683661584E+01')
      M_C = TO_FM('2.937801562768354E-01')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.3604503996731113868791517339909092506365724801689105M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 858
      M_A = TO_FM('0.9')
      M_B = TO_FM('23.4')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('7.3148127865937299821246829407023943740949130742928268M-18')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 859
      M_A = TO_FM('9.99496253868099M-1')
      M_B = TO_FM('2.47067979368109M+6')
      M_C = TO_FM('6.09475681774953M-100')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('1.7681753021411259894614747665450637683755190050365931M-544')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 860
      M_A = TO_FM('6.213433771653724M-1')
      M_B = TO_FM('8.854622686031200M-1')
      M_C = TO_FM('5.00000854049816M-121')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('1.1281271573737080091147788530326864610276172049831497M+0')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 861
      M_A = TO_FM('5.304391676698501M-15')
      M_B = TO_FM('4.870186358377400M+2')
      M_C = TO_FM('4.999955247889730M-98')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('8.7892314482956847896604128106803662527479433068750459M-6956')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 862
      M_A = TO_FM('1.882803169800314M-7')
      M_B = TO_FM('1.591547060066600M-169')
      M_C = TO_FM('3.521822614438970M+6')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('6.2831946669434576663925763649227277100409122269443137M+168')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 863
      M_A = TO_FM('.9999999999999')
      M_B = TO_FM('8.591098092677430M+2')
      M_C = TO_FM('1.863210949748253M+1')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('3.9062929191651064065641350979581425238442928803700306M-40')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 864
      M_A = TO_FM('2.531772074701081M-99')
      M_B = TO_FM('3.547571261801072M+2')
      M_C = TO_FM('1.974896958876250M+6')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('4.0957237103166196693191012056689839835950377114705018M-34981')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 865
      M_A = TO_FM('.99999999999999')
      M_B = TO_FM('1.0E-123')
      M_C = TO_FM('1.0E-134')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('1.0M+123')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 866
      M_A = TO_FM('1')
      M_B = TO_FM('2.65')
      M_C = TO_FM('4.88')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('1.5020204575152306127604878970920601604169827852591720M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 867
      M_A = TO_FM('0')
      M_B = TO_FM('2.65')
      M_C = TO_FM('4.88')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('0')
      M_D = ABS(M_C - M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 868
      M_A = TO_FM('.998')
      M_B = TO_FM('759.6')
      M_C = TO_FM('4.95e-57')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('9.7133692099062434492386763673434080317019087637060970M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 869
      M_A = TO_FM('4.764360371097952E-01')
      M_B = TO_FM('1.161514683661584E+01')
      M_C = TO_FM('2.937801562768354E-01')
      M_C = INCOMPLETE_BETA(M_A,M_B,M_C)
      M_D = TO_FM('2.3604503996731113868791517339909092506365724801689105M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST45

      SUBROUTINE TEST46(NCASE,NERROR,KLOG)

!  Test the Polygamma, Psi functions.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR

      WRITE (KW,"(/' Testing Polygamma, Psi.')")


      NCASE = 870
      M_A = TO_FM('4.5')
      CALL FM_PGAM(0,M_A,M_C)
      M_D = TO_FM('1.3888709263595289015114046193821968137592213477205183M+0')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 871
      M_A = TO_FM('1.0E-123')
      CALL FM_PGAM(1,M_A,M_C)
      M_D = TO_FM('1.0M+246')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 872
      M_A = TO_FM('1.0E-123')
      CALL FM_PGAM(2,M_A,M_C)
      M_D = TO_FM('-2.0M+369')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 873
      M_A = TO_FM('2.0706137739520290320140007735608464643737932737070189M-1')
      CALL FM_PGAM(1,M_A,M_C)
      M_D = TO_FM('2.4580954480899934124966756607870377560864828849100481M+1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 874
      M_A = TO_FM('2.0706137739520290320140007735608464643737932737070189M-1')
      CALL FM_PGAM(6,M_A,M_C)
      M_D = TO_FM('-4.4120531379423056741117517146346730469682094212273241M+7')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 875
      M_A = TO_FM('2.0706137739520290320140007735608464643737932737070189M-1')
      CALL FM_PGAM(23,M_A,M_C)
      M_D = TO_FM('6.7006365293376930742991440911935017694098601683947073M+38')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 876
      M_A = TO_FM('1.0E+123')
      CALL FM_PGAM(4,M_A,M_C)
      M_D = TO_FM('-6.0M-492')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 877
      M_A = TO_FM('-6.499999840238790109')
      CALL FM_PGAM(4,M_A,M_C)
      M_D = TO_FM('1.0135142464863270830609416082237513111216512170936928M-16')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 878
      M_C = POLYGAMMA(2,TO_FM('1.0E-123'))
      M_D = TO_FM('-2.0M+369')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PGAM ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 879
      M_A = TO_FM('1.0E-135')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('-1.0M+135')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 880
      M_A = TO_FM('1.2')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('-2.8903989659218829554720796244995210482558827420664281M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 881
      M_A = TO_FM('-3.4')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('2.3844508141180140670320531380285019520468887144980679M+0')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 882
      M_A = TO_FM('57')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('4.0342536898816977739559850955847848905386809772893269M+0')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 883
      M_A = TO_FM('1.0E+56')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('1.2894476520766655830500752146232439562566168336321129M+2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 884
      M_A = TO_FM('1.0')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('-5.7721566490153286060651209008240243104215933593992360M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 885
      M_A = TO_FM('1.0E+23456')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('5.4009435941268335564326007561076446853491436517276499M+4')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 886
      M_A = TO_FM('1.46163214496836234126266')
      CALL FM_PSI(M_A,M_C)
      M_D = TO_FM('4.4287869692570149446165609601581442013784186419176534M-25')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 887
      M_C = PSI(TO_FM('1.2'))
      M_D = TO_FM('-2.8903989659218829554720796244995210482558827420664281M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-49'))) THEN
          CALL ERRPRT_FM(' PSI  ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST46

      SUBROUTINE TEST47(NCASE,NERROR,KLOG)

!  Test the different rounding modes.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      INTEGER KLOG,NCASE,NERROR
      INTEGER SEED(7)

      WRITE (KW,"(/' Testing the different rounding modes.')")


      CALL FMSETVAR(' MBASE = 10 ')
      CALL FMSETVAR(' NDIG = 20 ')
      M_A = 0

      NCASE = 888
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('2')/TO_FM('3')
      M_D = TO_FM('.66666666666666666667')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 889
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('2')/TO_FM('3')
      M_D = TO_FM('.66666666666666666666')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 890
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('2')/TO_FM('3')
      M_D = TO_FM('.66666666666666666666')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 891
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('2')/TO_FM('3')
      M_D = TO_FM('.66666666666666666667')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 892
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('1')/TO_FM('3')
      M_D = TO_FM('.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 893
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('1')/TO_FM('3')
      M_D = TO_FM('.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 894
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('1')/TO_FM('3')
      M_D = TO_FM('.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 895
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('1')/TO_FM('3')
      M_D = TO_FM('.33333333333333333334')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 896
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('-1')/TO_FM('3')
      M_D = TO_FM('-.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 897
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('-1')/TO_FM('3')
      M_D = TO_FM('-.33333333333333333334')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 898
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('-1')/TO_FM('3')
      M_D = TO_FM('-.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 899
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('-1')/TO_FM('3')
      M_D = TO_FM('-.33333333333333333333')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 900
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('-2')/TO_FM('3')
      M_D = TO_FM('-.66666666666666666667')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 901
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('-2')/TO_FM('3')
      M_D = TO_FM('-.66666666666666666667')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 902
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('-2')/TO_FM('3')
      M_D = TO_FM('-.66666666666666666666')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 903
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('-2')/TO_FM('3')
      M_D = TO_FM('-.66666666666666666666')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 904
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('1') + TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 905
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('1') + TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 906
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('1') + TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 907
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('1') + TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000001')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 908
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('1') - TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 909
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('1') - TO_FM('3E-555')
      M_D = TO_FM('.99999999999999999999')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 910
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('1') - TO_FM('3E-555')
      M_D = TO_FM('.99999999999999999999')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 911
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('1') - TO_FM('3E-555')
      M_D = TO_FM('1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 912
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('-1') + TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 913
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('-1') + TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 914
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('-1') + TO_FM('3E-555')
      M_D = TO_FM('-.99999999999999999999')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 915
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('-1') + TO_FM('3E-555')
      M_D = TO_FM('-.99999999999999999999')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 916
      CALL FMSETVAR(' KROUND = 1 ')
      M_C = TO_FM('-1') - TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 917
      CALL FMSETVAR(' KROUND = -1 ')
      M_C = TO_FM('-1') - TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000001')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 918
      CALL FMSETVAR(' KROUND = 0 ')
      M_C = TO_FM('-1') - TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 919
      CALL FMSETVAR(' KROUND = 2 ')
      M_C = TO_FM('-1') - TO_FM('3E-555')
      M_D = TO_FM('-1.0000000000000000000')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      CALL FMSETVAR(' MBASE = 2 ')
      CALL FMSETVAR(' NDIG = 53 ')
      NCASE = 920
      M_A = TO_FM('0.125')
      M_B = TO_FM('23.25')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('6.1345805065305141873M-25')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 921
      M_A = TO_FM('0.52')
      M_B = TO_FM('2.01')
      M_C = TO_FM('1.6')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('1.0304844627978347604M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 922
      M_A = TO_FM('9.01737835258975M-1')
      M_B = TO_FM('2.00853601446773')
      M_C = TO_FM('1.59735792202923')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.2512248738228585986M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 923
      M_A = TO_FM('9.6097615596216720E-01')
      M_B = TO_FM('1.970425178583792')
      M_C = TO_FM('5.5680052333367')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.8619456987740165927M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 924
      M_A = TO_FM('4.764360371097952E-01')
      M_B = TO_FM('1.161514683661584E+01')
      M_C = TO_FM('2.937801562768354E-01')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.3604503996731113869M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 925
      M_A = TO_FM('0.9')
      M_B = TO_FM('23.4')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('7.3148127865937395334M-18')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-15'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      CALL FMSETVAR(' MBASE = 3 ')
      CALL FMSETVAR(' NDIG = 55 ')
      NCASE = 926
      M_A = TO_FM('0.1')
      M_B = TO_FM('23.4')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('5.87319809189607304633501593392681M-27')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 927
      M_A = TO_FM('0.52')
      M_B = TO_FM('2.1')
      M_C = TO_FM('1.6')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('9.25745341552810210762563659429375M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 928
      M_A = TO_FM('9.01737835258975M-1')
      M_B = TO_FM('2.00853601446773')
      M_C = TO_FM('1.59735792202923')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.25122487382285859767535178829535M-1')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 929
      M_A = TO_FM('9.6097615596216720E-01')
      M_B = TO_FM('1.970425178583792')
      M_C = TO_FM('5.5680052333367')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.861945698774016536409296855493M-2')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 930
      M_A = TO_FM('4.764360371097952E-01')
      M_B = TO_FM('1.161514683661584E+01')
      M_C = TO_FM('2.937801562768354E-01')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('2.36045039967311138687915158221269M-5')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 931
      M_A = TO_FM('0.9')
      M_B = TO_FM('23.4')
      M_C = TO_FM('34.5')
      CALL FM_IBTA(M_A,M_B,M_C,MFM6)
      CALL FM_EQ(MFM6,M_C)
      M_D = TO_FM('7.31481278659372998212468424608367M-18')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' IBTA ',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 932
      CALL FPST2M('1.67',MP1)
      CALL FPST2M('2.64',MP2)
      CALL FPADD(MP1,MP2,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPST2M('-3.91',MP2)
      CALL FPSUB(MP1,MP2,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPST2M('4.58',MP2)
      CALL FPMPY(MP1,MP2,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPST2M('0.27',MP2)
      CALL FPDIV(MP1,MP2,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPADDI(MP1,2)
      CALL FPMPYI(MP1,13,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPDIVI(MP1,11,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPLN(MP1,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPSIN(MP1,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPCOS(MP1,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPEXP(MP1,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FPGAM(MP1,MP3)
      CALL FPEQ(MP3,MP1)
      CALL FMUNPK(MP1,M_C%MFM)
      M_D = TO_FM('0.941122001974472326543759839200398')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-25'))) THEN
          CALL ERRPRT_FM(' Pack ',M_C,'M_C',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 933
      SEED = (/ 2718281,8284590,4523536,0287471,3526624,9775724,7093699 /)
      CALL FM_RANDOM_SEED(PUT=SEED)
      DO J1 = 1, 10
         CALL FM_RANDOM_NUMBER(D1)
      ENDDO
      M_C = D1
      M_D = TO_FM('0.945608442536777')
      M_D = ABS((M_C - M_D)/M_D)
      IF (.NOT.(M_D <= TO_FM('1.0E-10'))) THEN
          CALL ERRPRT_FM(' Rand ',M_C,'M_C',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      CALL FMSETVAR(' MBASE = 10000 ')
      CALL FMSETVAR(' NDIG = 5 ')

      NCASE = 934
      CALL FMSETVAR(' KROUND = 1 ')
      CALL FMSETVAR(' KRPERF = 1 ')
      M_C = SQRT( TO_FM('.49841718043038996023') )
      M_D = TO_FM('.70598667156709832621')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      NCASE = 935
      CALL FMSETVAR(' KROUND = 1 ')
      CALL FMSETVAR(' KRPERF = 0 ')
      M_C = SQRT( TO_FM('.49841718043038996023') )
      M_D = TO_FM('.70598667156709832622')
      IF (.NOT.(M_D == M_C)) THEN
          CALL ERRPRT_FM(' Round',M_A,'M_A',M_C,'M_C',M_D,'M_D',  &
                         NCASE,NERROR,KLOG)
      ENDIF

      RETURN
      END SUBROUTINE TEST47

      SUBROUTINE ERRPRTFM(NROUT,M1,NAME1,M2,NAME2,M3,NAME3,  &
                          NCASE,NERROR,KLOG)

!  Print error messages for testing of real (FM) routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using FMST2M.
!  M3 is ABS(M1-M2), and ERRPRT is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in the
!                    calling routine correspond to M1,M2,M3.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      REAL (KIND(1.0D0)) :: M1(-1:LUNPCK),M2(-1:LUNPCK),M3(-1:LUNPCK)

      CHARACTER(2) :: NAME1,NAME2,NAME3
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT
      WRITE (KLOG,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT

!              Temporarily change KW to KLOG so FMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL FMPRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL FMPRNT(M2)
      WRITE (KLOG,"(1X,A,' =')") NAME3
      CALL FMPRNT(M3)
      KW = KWSAVE
      RETURN
      END SUBROUTINE ERRPRTFM

      SUBROUTINE ERRPRTIM(NROUT,M1,NAME1,M2,NAME2,  &
                          NCASE,NERROR,KLOG)

!  Print error messages for testing of integer (IM) routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using IMST2M.
!  NAME1,NAME2 are strings identifying which variables in the calling routine
!              correspond to M1,M2.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      REAL (KIND(1.0D0)) :: M1(-1:LUNPCK),M2(-1:LUNPCK)

      CHARACTER(2) :: NAME1,NAME2
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT
      WRITE (KLOG,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT

!              Temporarily change KW to KLOG so IMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL IMPRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL IMPRNT(M2)
      KW = KWSAVE
      END SUBROUTINE ERRPRTIM

      SUBROUTINE ERRPRTZM(NROUT,M1,NAME1,M2,NAME2,M3,NAME3,  &
                          NCASE,NERROR,KLOG)

!  Print error messages.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using ZMST2M.
!  M3 is ABS(M1-M2), and ERRPRTZM is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in the
!                    calling routine correspond to M1,M2,M3.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      REAL (KIND(1.0D0)) :: M1(-1:LUNPKZ),M2(-1:LUNPKZ),M3(-1:LUNPKZ)

      CHARACTER(2) :: NAME1,NAME2,NAME3
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT
      WRITE (KLOG,  &
          "(//' Error in case',I4,'.  The routine',' being tested was ',A6)"  &
          ) NCASE,NROUT

!              Temporarily change KW to KLOG so ZMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL ZMPRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL ZMPRNT(M2)
      WRITE (KLOG,"(1X,A,' =')") NAME3
      CALL ZMPRNT(M3)
      KW = KWSAVE
      END SUBROUTINE ERRPRTZM

      SUBROUTINE ERRPRT_FM(NROUT,M1,NAME1,M2,NAME2,M3,NAME3,  &
                           NCASE,NERROR,KLOG)

!  Print error messages for testing of TYPE (FM) interface routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using FMST2M.
!  M3 is ABS(M1-M2), and ERRPRT_FM is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in the
!                    calling routine correspond to M1,M2,M3.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      TYPE (FM) M1,M2,M3

      CHARACTER(3) :: NAME1,NAME2,NAME3
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
         "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
         ) NCASE,NROUT
      WRITE (KLOG,  &
         "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
         ) NCASE,NROUT

!              Temporarily change KW to KLOG so FMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL FM_PRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL FM_PRNT(M2)
      WRITE (KLOG,"(1X,A,' =')") NAME3
      CALL FM_PRNT(M3)
      KW = KWSAVE
      END SUBROUTINE ERRPRT_FM

      SUBROUTINE ERRPRT_IM(NROUT,M1,NAME1,M2,NAME2,  &
                           NCASE,NERROR,KLOG)

!  Print error messages for testing of TYPE (IM) interface routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using IMST2M.
!  NAME1,NAME2 are strings identifying which variables in the calling routine
!              correspond to M1,M2.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      TYPE (IM) M1,M2

      CHARACTER(3) :: NAME1,NAME2
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
        "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
        ) NCASE,NROUT
      WRITE (KLOG,  &
        "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
        ) NCASE,NROUT

!              Temporarily change KW to KLOG so IMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL IM_PRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL IM_PRNT(M2)
      KW = KWSAVE
      END SUBROUTINE ERRPRT_IM

      SUBROUTINE ERRPRT_ZM(NROUT,M1,NAME1,M2,NAME2,M3,NAME3,  &
                           NCASE,NERROR,KLOG)

!  Print error messages for testing of TYPE (ZM) interface routines.

!  M1 is the value to be tested, as computed by the routine named NROUT.
!  M2 is the reference value, usually converted using ZMST2M.
!  M3 is ABS(M1-M2), and ERRPRT_ZM is called if this is too big.
!  NAME1,NAME2,NAME3 are strings identifying which variables in the calling routine
!                    correspond to M1,M2,M3.

      USE FMVALS
      USE FMZM
      USE TEST_VARS

      IMPLICIT NONE

      TYPE (ZM) M1,M2,M3

      CHARACTER(3) :: NAME1,NAME2,NAME3
      CHARACTER(6) :: NROUT
      INTEGER KLOG,KWSAVE,NCASE,NERROR

      NERROR = NERROR + 1
      WRITE (KW,  &
        "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
        ) NCASE,NROUT
      WRITE (KLOG,  &
        "(//' Error in case',I4,'.  The interface',' being tested was ',A6)"  &
        ) NCASE,NROUT

!              Temporarily change KW to KLOG so ZMPRNT
!              will write to the log file.

      KWSAVE = KW
      KW = KLOG
      WRITE (KLOG,"(1X,A,' =')") NAME1
      CALL ZM_PRNT(M1)
      WRITE (KLOG,"(1X,A,' =')") NAME2
      CALL ZM_PRNT(M2)
      WRITE (KLOG,"(1X,A,' =')") NAME3
      CALL ZM_PRNT(M3)
      KW = KWSAVE
      END SUBROUTINE ERRPRT_ZM

      SUBROUTINE PRTERR(KW,KLOG,NCASE,NERROR)
      IMPLICIT NONE
      INTEGER KW,KLOG,NCASE,NERROR

      WRITE (KW,*) ' Error in case ',NCASE
      WRITE (KLOG,*) ' '
      WRITE (KLOG,*) ' Error in case ',NCASE
      NERROR = NERROR + 1
      END SUBROUTINE PRTERR

      SUBROUTINE TIMEIT(TIME)

      INTEGER JTIME,JRATE
      REAL TIME

!  Return the system time.  f90 version.

      CALL SYSTEM_CLOCK(JTIME,JRATE)
      TIME = REAL(JTIME)/REAL(JRATE)
      RETURN
      END SUBROUTINE TIMEIT

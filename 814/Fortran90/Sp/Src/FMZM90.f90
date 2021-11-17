 MODULE FMZM_1

!  FMZM 1.2                        David M. Smith

!  This module extends the definition of Fortran-90 arithmetic and function
!  operations so they also apply to multiple precision numbers, using version
!  1.2 of FM.
!  There are three multiple precision data types:
!     FM  (multiple precision real)
!     IM  (multiple precision integer)
!     ZM  (multiple precision complex)

!  Some of the interface routines assume that the precision chosen in the
!  calling program (using FMSET) represents more significant digits than does
!  the machine's double precision.

!  Most of the functions defined in this module are multiple precision versions
!  of standard Fortran-90 functions.  In addition, there are functions for
!  direct conversion, formatting, and some mathematical special functions.

!  TO_FM is a function for converting other types of numbers to type FM.  Note
!  that TO_FM(3.12) converts the REAL constant to FM, but it is accurate only
!  to single precision.  TO_FM(3.12D0) agrees with 3.12 to double precision
!  accuracy, and TO_FM('3.12') or TO_FM(312)/TO_FM(100) agrees to full FM
!  accuracy.

!  TO_IM converts to type IM, and TO_ZM converts to type ZM.

!  Functions are also supplied for converting the three multiple precision types
!  to the other numeric data types:
!     TO_INT   converts to machine precision integer
!     TO_SP    converts to single precision
!     TO_DP    converts to double precision
!     TO_SPZ   converts to single precision complex
!     TO_DPZ   converts to double precision complex

!  WARNING:   When multiple precision type declarations are inserted in an
!             existing program, take care in converting functions like DBLE(X),
!             where X has been declared as a multiple precision type.  If X
!             was single precision in the original program, then replacing
!             the DBLE(X) by TO_DP(X) in the new version could lose accuracy.
!             For this reason, the Fortran type-conversion functions defined
!             in this module assume that results should be multiple precision
!             whenever inputs are.  Examples:
!             DBLE(TO_FM('1.23E+123456')) is type FM
!             REAL(TO_FM('1.23E+123456')) is type FM
!             REAL(TO_ZM('3.12+4.56i'))   is type FM   = TO_FM('3.12')
!             INT(TO_FM('1.23'))          is type IM   = TO_IM(1)
!             INT(TO_IM('1E+23'))         is type IM
!             CMPLX(TO_FM('1.23'),TO_FM('4.56')) is type ZM

!  Programs using this module may sometimes need to call FM, IM, or ZM routines
!  directly.  This is normally the case when routines are needed that are not
!  Fortran-90 intrinsics, such as the formatting subroutine FMFORM.  In a
!  program using this module, suppose MAFM has been declared with
!  TYPE ( FM ) MAFM.  To use the routine FMFORM, which expects the second
!  argument to be an array and not a derived type, the call would have to be
!  CALL FMFORM('F65.60',MAFM%MFM,ST1) so that the array contained in MAFM is
!  passed.

!  As an alternative so the user can refer directly to the FM-, IM-, and ZM-type
!  variables and avoid the cumbersome "%MFM" suffixes, the FM package provides a
!  collection of interface routines to supply any needed argument conversions.
!  For each FM, IM, and ZM routine that is designed to be called by the user,
!  there is also a version that assumes any multiple-precision arguments are
!  derived types instead of arrays.  Each interface routine has the same name as
!  the original with an underscore after the first two letters of the routine
!  name.  To convert the number to a character string with F65.60 format, use
!  CALL FM_FORM('F65.60',MAFM,ST1) if MAFM is of TYPE ( FM ), or use
!  CALL FMFORM('F65.60',MA,ST1) if MA is declared as an array.  All the routines
!  shown below may be used this way.

!  For each of the operations =,  == ,  /= ,  < ,  <= ,  > ,  >= , +, -, *, /,
!  and **, the interface module defines all mixed mode variations involving one
!  of the three multiple precision derived types and another argument having one
!  of the types: { integer, real, double, complex, complex double, FM, IM, ZM }.
!  So mixed mode expressions such as
!        MAFM = 12
!        MAFM = MAFM + 1
!        IF (ABS(MAFM) > 1.0D-23) THEN
!  are handled correctly.

!  Not all the named functions are defined for all three multiple precision
!  derived types, so the list below shows which can be used.  The labels "real",
!  "integer", and "complex" refer to types FM, IM, and ZM respectively, "string"
!  means the function accepts character strings (e.g., TO_FM('3.45')), and
!  "other" means the function can accept any of the machine precision data types
!  integer, real, double, complex, or complex double.  For functions that accept
!  two or more arguments, like ATAN2 or MAX, all the arguments must be of the
!  same type.


!  AVAILABLE OPERATIONS:

!     =
!     +
!     -
!     *
!     /
!     **
!     ==
!     /=
!     <
!     <=
!     >
!     >=
!     ABS          real    integer    complex
!     ACOS         real               complex
!     AIMAG                           complex
!     AINT         real               complex
!     ANINT        real               complex
!     ASIN         real               complex
!     ATAN         real               complex
!     ATAN2        real
!     BTEST                integer
!     CEILING      real               complex
!     CMPLX        real    integer
!     CONJ                            complex
!     COS          real               complex
!     COSH         real               complex
!     DBLE         real    integer    complex
!     DIGITS       real    integer    complex
!     DIM          real    integer
!     DINT         real               complex
!     DOTPRODUCT   real    integer    complex
!     EPSILON      real
!     EXP          real               complex
!     EXPONENT     real
!     FLOOR        real    integer    complex
!     FRACTION     real               complex
!     HUGE         real    integer    complex
!     INT          real    integer    complex
!     LOG          real               complex
!     LOG10        real               complex
!     MATMUL       real    integer    complex
!     MAX          real    integer
!     MAXEXPONENT  real
!     MIN          real    integer
!     MINEXPONENT  real
!     MOD          real    integer
!     MODULO       real    integer
!     NEAREST      real
!     NINT         real    integer    complex
!     PRECISION    real               complex
!     RADIX        real    integer    complex
!     RANGE        real    integer    complex
!     REAL         real    integer    complex
!     RRSPACING    real
!     SCALE        real               complex
!     SETEXPONENT  real
!     SIGN         real    integer
!     SIN          real               complex
!     SINH         real               complex
!     SPACING      real
!     SQRT         real               complex
!     TAN          real               complex
!     TANH         real               complex
!     TINY         real    integer    complex
!     TO_FM        real    integer    complex    string    other
!     TO_IM        real    integer    complex    string    other
!     TO_ZM        real    integer    complex    string    other
!     TO_INT       real    integer    complex
!     TO_SP        real    integer    complex
!     TO_DP        real    integer    complex
!     TO_SPZ       real    integer    complex
!     TO_DPZ       real    integer    complex

!  Some other functions are defined that do not correspond to machine precision
!  intrinsic functions.  These include formatting functions, integer modular
!  functions and GCD, and the Gamma function and its related functions.
!  N below is a machine precision integer, J1, J2, J3 are TYPE (IM), FMT, FMTR,
!  FMTI are character strings, A,B,X are TYPE (FM), and Z is TYPE (ZM).
!  The three formatting functions return a character string containing the
!  formatted number, the three TYPE (IM) functions return a TYPE (IM) result,
!  and the 12 special functions return TYPE (FM) results.

!  Formatting functions:

!     FM_FORMAT(FMT,A)        Put A into FMT (real) format
!     IM_FORMAT(FMT,J1)       Put J1 into FMT (integer) format
!     ZM_FORMAT(FMTR,FMTI,Z)  Put Z into (complex) format, FMTR for the real
!                             part and FMTI for the imaginary part

!     Examples:
!        ST = FM_FORMAT('F65.60',A)
!        WRITE (*,*) ' A = ',TRIM(ST)
!        ST = FM_FORMAT('E75.60',B)
!        WRITE (*,*) ' B = ',ST(1:75)
!        ST = IM_FORMAT('I50',J1)
!        WRITE (*,*) ' J1 = ',ST(1:50)
!        ST = ZM_FORMAT('F35.30','F30.25',Z)
!        WRITE (*,*) ' Z = ',ST(1:70)

!     These functions are used for one-line output.  The returned character
!     strings are of length 200.  Avoid using the formatting function in the
!     write list, as in
!        WRITE (*,*) ' J1 = ',IM_FORMAT('I50',J1)(1:50)
!     since the formatting functions may themselves execute an internal WRITE
!     and that would cause a recursive write reference.

!     For higher precision numbers, the output can be broken onto multiple
!     lines automatically by calling subroutines FM_PRNT, IM_PRNT, ZM_PRNT,
!     or the line breaks can be done by hand after calling one of the
!     subroutines FM_FORM, IM_FORM, ZM_FORM.

!     For ZM_FORMAT the length of the output is 5 more than the sum of the
!     two field widths.

!  Integer functions:

!     GCD(J1,J2)              Greatest Common Divisor of J1 and J2.
!     MULTIPLY_MOD(J1,J2,J3)  J1 * J2 mod J3
!     POWER_MOD(J1,J2,J3)     J1 ** J2 mod J3

!  Special functions:

!     BERNOULLI(N)            Nth Bernoulli number
!     BETA(A,B)               Integral (0 to 1)  t**(A-1) * (1-t)**(B-1)  dt
!     BINOMIAL(A,B)           Binomial Coefficient  A! / ( B! (A-B)! )
!     FACTORIAL(A)            A!
!     GAMMA(A)                Integral (0 to infinity)  t**(A-1) * exp(-t)  dt
!     INCOMPLETE_BETA(X,A,B)  Integral (0 to X)  t**(A-1) * (1-t)**(B-1)  dt
!     INCOMPLETE_GAMMA1(A,X)  Integral (0 to X)  t**(A-1) * exp(-t)  dt
!     INCOMPLETE_GAMMA2(A,X)  Integral (X to infinity)  t**(A-1) * exp(-t)  dt
!     LOG_GAMMA(A)            Ln( GAMMA(A) )
!     POLYGAMMA(N,A)          Nth derivative of Psi(x), evaluated at A
!     POCHHAMMER(A,N)         A*(A+1)*(A+2)*...*(A+N-1)
!     PSI(A)                  Derivative of Ln(Gamma(x)), evaluated at A


!        To keep the FM variables hidden from a program that uses this
!        module, these parameters are set to the same values as the
!        corresponding ones in the FM_VARIABLES module.

    USE FMVALS, ONLY : NDIGMX_2 => NDIGMX
    INTEGER, PARAMETER, PRIVATE :: LUNPCK_2 = (11*NDIGMX_2)/5 + 30
    INTEGER, PARAMETER, PRIVATE :: LUNPKZ_2 = 2*LUNPCK_2+2

    TYPE FM
       REAL (KIND(1.0D0)) :: MFM(-1:LUNPCK_2)
    END TYPE

    TYPE IM
       REAL (KIND(1.0D0)) :: MIM(-1:LUNPCK_2)
    END TYPE

    TYPE ZM
       REAL (KIND(1.0D0)) :: MZM(-1:LUNPKZ_2)
    END TYPE

    REAL (KIND(1.0D0)), SAVE, DIMENSION(-1:LUNPCK_2) :: MTFM,MUFM,MVFM
    REAL (KIND(1.0D0)), SAVE, DIMENSION(-1:LUNPCK_2) :: MTIM,MUIM,MVIM
    REAL (KIND(1.0D0)), SAVE, DIMENSION(-1:LUNPKZ_2) :: MTZM,MUZM,MVZM

 END MODULE FMZM_1

 MODULE FMZM_2
    USE FMZM_1

!   These abbreviations are used for operations
!   on the various data types.

!   I    Integer
!   R    Real
!   D    Double Precision
!   Z    Complex
!   C    Complex Double Precision
!   FM   Multiple precision real
!   IM   Multiple precision integer
!   ZM   Multiple precision complex

!   For example, the "=" procedure FMEQ_FMD is for statements like
!   X = A, where X is type FM and A is type Double Precision.

    INTERFACE ASSIGNMENT (=)
       MODULE PROCEDURE FMEQ_IFM
       MODULE PROCEDURE FMEQ_IIM
       MODULE PROCEDURE FMEQ_IZM
       MODULE PROCEDURE FMEQ_RFM
       MODULE PROCEDURE FMEQ_RIM
       MODULE PROCEDURE FMEQ_RZM
       MODULE PROCEDURE FMEQ_DFM
       MODULE PROCEDURE FMEQ_DIM
       MODULE PROCEDURE FMEQ_DZM
       MODULE PROCEDURE FMEQ_ZFM
       MODULE PROCEDURE FMEQ_ZIM
       MODULE PROCEDURE FMEQ_ZZM
       MODULE PROCEDURE FMEQ_CFM
       MODULE PROCEDURE FMEQ_CIM
       MODULE PROCEDURE FMEQ_CZM
       MODULE PROCEDURE FMEQ_FMI
       MODULE PROCEDURE FMEQ_FMR
       MODULE PROCEDURE FMEQ_FMD
       MODULE PROCEDURE FMEQ_FMZ
       MODULE PROCEDURE FMEQ_FMC
       MODULE PROCEDURE FMEQ_FMFM
       MODULE PROCEDURE FMEQ_FMIM
       MODULE PROCEDURE FMEQ_FMZM
       MODULE PROCEDURE FMEQ_IMI
       MODULE PROCEDURE FMEQ_IMR
       MODULE PROCEDURE FMEQ_IMD
       MODULE PROCEDURE FMEQ_IMZ
       MODULE PROCEDURE FMEQ_IMC
       MODULE PROCEDURE FMEQ_IMFM
       MODULE PROCEDURE FMEQ_IMIM
       MODULE PROCEDURE FMEQ_IMZM
       MODULE PROCEDURE FMEQ_ZMI
       MODULE PROCEDURE FMEQ_ZMR
       MODULE PROCEDURE FMEQ_ZMD
       MODULE PROCEDURE FMEQ_ZMZ
       MODULE PROCEDURE FMEQ_ZMC
       MODULE PROCEDURE FMEQ_ZMFM
       MODULE PROCEDURE FMEQ_ZMIM
       MODULE PROCEDURE FMEQ_ZMZM
    END INTERFACE

    INTERFACE OPERATOR ( == )
       MODULE PROCEDURE FMLEQ_IFM
       MODULE PROCEDURE FMLEQ_IIM
       MODULE PROCEDURE FMLEQ_IZM
       MODULE PROCEDURE FMLEQ_RFM
       MODULE PROCEDURE FMLEQ_RIM
       MODULE PROCEDURE FMLEQ_RZM
       MODULE PROCEDURE FMLEQ_DFM
       MODULE PROCEDURE FMLEQ_DIM
       MODULE PROCEDURE FMLEQ_DZM
       MODULE PROCEDURE FMLEQ_ZFM
       MODULE PROCEDURE FMLEQ_ZIM
       MODULE PROCEDURE FMLEQ_ZZM
       MODULE PROCEDURE FMLEQ_CFM
       MODULE PROCEDURE FMLEQ_CIM
       MODULE PROCEDURE FMLEQ_CZM
       MODULE PROCEDURE FMLEQ_FMI
       MODULE PROCEDURE FMLEQ_FMR
       MODULE PROCEDURE FMLEQ_FMD
       MODULE PROCEDURE FMLEQ_FMZ
       MODULE PROCEDURE FMLEQ_FMC
       MODULE PROCEDURE FMLEQ_FMFM
       MODULE PROCEDURE FMLEQ_FMIM
       MODULE PROCEDURE FMLEQ_FMZM
       MODULE PROCEDURE FMLEQ_IMI
       MODULE PROCEDURE FMLEQ_IMR
       MODULE PROCEDURE FMLEQ_IMD
       MODULE PROCEDURE FMLEQ_IMZ
       MODULE PROCEDURE FMLEQ_IMC
       MODULE PROCEDURE FMLEQ_IMFM
       MODULE PROCEDURE FMLEQ_IMIM
       MODULE PROCEDURE FMLEQ_IMZM
       MODULE PROCEDURE FMLEQ_ZMI
       MODULE PROCEDURE FMLEQ_ZMR
       MODULE PROCEDURE FMLEQ_ZMD
       MODULE PROCEDURE FMLEQ_ZMZ
       MODULE PROCEDURE FMLEQ_ZMC
       MODULE PROCEDURE FMLEQ_ZMFM
       MODULE PROCEDURE FMLEQ_ZMIM
       MODULE PROCEDURE FMLEQ_ZMZM
    END INTERFACE


 CONTAINS

!                                                                   =

   SUBROUTINE FMEQ_IFM(IVAL,MA)
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: IVAL
      INTENT (IN) :: MA
      CALL FMM2I(MA%MFM,IVAL)
   END SUBROUTINE

   SUBROUTINE FMEQ_IIM(IVAL,MA)
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: IVAL
      INTENT (IN) :: MA
      CALL IMM2I(MA%MIM,IVAL)
   END SUBROUTINE

   SUBROUTINE FMEQ_IZM(IVAL,MA)
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: IVAL
      INTENT (IN) :: MA
      CALL ZMM2I(MA%MZM,IVAL)
   END SUBROUTINE

   SUBROUTINE FMEQ_RFM(R,MA)
      TYPE ( FM ) MA
      REAL R
      INTENT (INOUT) :: R
      INTENT (IN) :: MA
      CALL FMM2SP(MA%MFM,R)
   END SUBROUTINE

   SUBROUTINE FMEQ_RIM(R,MA)
      TYPE ( IM ) MA
      REAL R
      INTENT (INOUT) :: R
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMM2SP(MTFM,R)
   END SUBROUTINE

   SUBROUTINE FMEQ_RZM(R,MA)
      TYPE ( ZM ) MA
      REAL R
      INTENT (INOUT) :: R
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2SP(MTFM,R)
   END SUBROUTINE

   SUBROUTINE FMEQ_DFM(D,MA)
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (INOUT) :: D
      INTENT (IN) :: MA
      CALL FMM2DP(MA%MFM,D)
   END SUBROUTINE

   SUBROUTINE FMEQ_DIM(D,MA)
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (INOUT) :: D
      INTENT (IN) :: MA
      CALL IMM2DP(MA%MIM,D)
   END SUBROUTINE

   SUBROUTINE FMEQ_DZM(D,MA)
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (INOUT) :: D
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,D)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZFM(Z,MA)
      TYPE ( FM ) MA
      COMPLEX Z
      REAL R
      INTENT (INOUT) :: Z
      INTENT (IN) :: MA
      CALL FMM2SP(MA%MFM,R)
      Z = CMPLX( R , 0.0 )
   END SUBROUTINE

   SUBROUTINE FMEQ_ZIM(Z,MA)
      TYPE ( IM ) MA
      COMPLEX Z
      DOUBLE PRECISION D
      INTENT (INOUT) :: Z
      INTENT (IN) :: MA
      CALL IMM2DP(MA%MIM,D)
      Z = CMPLX( REAL(D) , 0.0 )
   END SUBROUTINE

   SUBROUTINE FMEQ_ZZM(Z,MA)
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (INOUT) :: Z
      INTENT (IN) :: MA
      CALL ZMM2Z(MA%MZM,Z)
   END SUBROUTINE

   SUBROUTINE FMEQ_CFM(C,MA)
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      INTENT (INOUT) :: C
      INTENT (IN) :: MA
      CALL FMM2DP(MA%MFM,D)
      C = CMPLX( D , 0.0D0 , KIND(0.0D0) )
   END SUBROUTINE

   SUBROUTINE FMEQ_CIM(C,MA)
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      INTENT (INOUT) :: C
      INTENT (IN) :: MA
      CALL IMM2DP(MA%MIM,D)
      C = CMPLX( D , 0.0D0 , KIND(0.0D0) )
   END SUBROUTINE

   SUBROUTINE FMEQ_CZM(C,MA)
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D1,D2
      INTENT (INOUT) :: C
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,D1)
      CALL ZMIMAG(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,D2)
      C = CMPLX( D1 , D2 , KIND(0.0D0) )
   END SUBROUTINE

   SUBROUTINE FMEQ_FMI(MA,IVAL)
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: MA
      INTENT (IN) :: IVAL
      CALL FMI2M(IVAL,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMR(MA,R)
      TYPE ( FM ) MA
      REAL R
      INTENT (INOUT) :: MA
      INTENT (IN) :: R
      CALL FMSP2M(R,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMD(MA,D)
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (INOUT) :: MA
      INTENT (IN) :: D
      CALL FMDP2M(D,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMZ(MA,Z)
      TYPE ( FM ) MA
      COMPLEX Z
      REAL R
      INTENT (INOUT) :: MA
      INTENT (IN) :: Z
      R = REAL(Z)
      CALL FMSP2M(R,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMC(MA,C)
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      INTENT (INOUT) :: MA
      INTENT (IN) :: C
      D = REAL(C,KIND(0.0D0))
      CALL FMDP2M(D,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMFM(MA,MB)
      TYPE ( FM ) MA,MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL FMEQ(MB%MFM,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMIM(MA,MB)
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL IMI2FM(MB%MIM,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_FMZM(MA,MB)
      TYPE ( FM ) MA
      TYPE ( ZM ) MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL ZMREAL(MB%MZM,MA%MFM)
   END SUBROUTINE

   SUBROUTINE FMEQ_IMI(MA,IVAL)
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: MA
      INTENT (IN) :: IVAL
      CALL IMI2M(IVAL,MA%MIM)
   END SUBROUTINE

   SUBROUTINE FMEQ_IMR(MA,R)
      TYPE ( IM ) MA
      INTEGER IVAL
      REAL R
      CHARACTER(25) :: ST
      INTENT (INOUT) :: MA
      INTENT (IN) :: R
      IF (ABS(R) < HUGE(1)) THEN
          IVAL = INT(R)
          CALL IMI2M(IVAL,MA%MIM)
      ELSE
          WRITE (ST,'(E25.16)') R
          CALL IMST2M(ST,MA%MIM)
      ENDIF
   END SUBROUTINE

   SUBROUTINE FMEQ_IMD(MA,D)
      TYPE ( IM ) MA
      INTEGER IVAL
      DOUBLE PRECISION D
      CHARACTER(25) :: ST
      INTENT (INOUT) :: MA
      INTENT (IN) :: D
      IF (ABS(D) < HUGE(1)) THEN
          IVAL = INT(D)
          CALL IMI2M(IVAL,MA%MIM)
      ELSE
          WRITE (ST,'(E25.16)') D
          CALL IMST2M(ST,MA%MIM)
      ENDIF
   END SUBROUTINE

   SUBROUTINE FMEQ_IMZ(MA,Z)
      TYPE ( IM ) MA
      COMPLEX Z
      REAL R
      CHARACTER(25) :: ST
      INTENT (INOUT) :: MA
      INTENT (IN) :: Z
      R = REAL(Z)
      IF (ABS(R) < HUGE(1)) THEN
          IVAL = INT(R)
          CALL IMI2M(IVAL,MA%MIM)
      ELSE
          WRITE (ST,'(E25.16)') R
          CALL IMST2M(ST,MA%MIM)
      ENDIF
   END SUBROUTINE

   SUBROUTINE FMEQ_IMC(MA,C)
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      CHARACTER(25) :: ST
      INTENT (INOUT) :: MA
      INTENT (IN) :: C
      D = REAL(C)
      IF (ABS(D) < HUGE(1)) THEN
          IVAL = INT(D)
          CALL IMI2M(IVAL,MA%MIM)
      ELSE
          WRITE (ST,'(E25.16)') D
          CALL IMST2M(ST,MA%MIM)
      ENDIF
   END SUBROUTINE

   SUBROUTINE FMEQ_IMFM(MA,MB)
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL IMFM2I(MB%MFM,MA%MIM)
   END SUBROUTINE

   SUBROUTINE FMEQ_IMIM(MA,MB)
      TYPE ( IM ) MA,MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL IMEQ(MB%MIM,MA%MIM)
   END SUBROUTINE

   SUBROUTINE FMEQ_IMZM(MA,MB)
      TYPE ( IM ) MA
      TYPE ( ZM ) MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL ZMREAL(MB%MZM,MTFM)
      CALL IMFM2I(MTFM,MA%MIM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMI(MA,IVAL)
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (INOUT) :: MA
      INTENT (IN) :: IVAL
      CALL ZMI2M(IVAL,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMR(MA,R)
      TYPE ( ZM ) MA
      REAL R
      COMPLEX Z
      INTENT (INOUT) :: MA
      INTENT (IN) :: R
      Z = CMPLX(R,0.0)
      CALL ZMZ2M(Z,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMD(MA,D)
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (INOUT) :: MA
      INTENT (IN) :: D
      CALL FMDP2M(D,MTFM)
      CALL FMDP2M(0.0D0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMZ(MA,Z)
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (INOUT) :: MA
      INTENT (IN) :: Z
      CALL ZMZ2M(Z,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMC(MA,C)
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      INTENT (INOUT) :: MA
      INTENT (IN) :: C
      D = REAL(C,KIND(0.0D0))
      CALL FMDP2M(D,MTFM)
      D = AIMAG(C)
      CALL FMDP2M(D,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMFM(MA,MB)
      TYPE ( FM ) MB
      TYPE ( ZM ) MA
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMIM(MA,MB)
      TYPE ( IM ) MB
      TYPE ( ZM ) MA
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MA%MZM)
   END SUBROUTINE

   SUBROUTINE FMEQ_ZMZM(MA,MB)
      TYPE ( ZM ) MA,MB
      INTENT (INOUT) :: MA
      INTENT (IN) :: MB
      CALL ZMEQ(MB%MZM,MA%MZM)
   END SUBROUTINE

!                                                                ==

   FUNCTION FMLEQ_IFM(IVAL,MA)
      LOGICAL FMLEQ_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLEQ_IFM = FMCOMP(MTFM,'EQ',MA%MFM)
   END FUNCTION

   FUNCTION FMLEQ_IIM(IVAL,MA)
      LOGICAL FMLEQ_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLEQ_IIM = IMCOMP(MTIM,'EQ',MA%MIM)
   END FUNCTION

   FUNCTION FMLEQ_IZM(IVAL,MA)
      LOGICAL FMLEQ_IZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_IZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_RFM(R,MA)
      LOGICAL FMLEQ_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLEQ_RFM = FMCOMP(MTFM,'EQ',MA%MFM)
   END FUNCTION

   FUNCTION FMLEQ_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLEQ_RIM,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLEQ_RIM = FMCOMP(MTFM,'EQ',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLEQ_RZM(R,MA)
      LOGICAL FMLEQ_RZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_RZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_DFM(D,MA)
      LOGICAL FMLEQ_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLEQ_DFM = FMCOMP(MTFM,'EQ',MA%MFM)
   END FUNCTION

   FUNCTION FMLEQ_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLEQ_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLEQ_DIM = FMCOMP(MTFM,'EQ',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLEQ_DZM(D,MA)
      LOGICAL FMLEQ_DZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_DZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZFM(Z,MA)
      LOGICAL FMLEQ_ZFM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL FMSP2M(REAL(Z),MTFM)
      L1 = FMCOMP(MTFM,'EQ',MA%MFM)
      L2 = .TRUE.
      IF (AIMAG(Z) /= 0.0) L2 = .FALSE.
      FMLEQ_ZFM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZIM(Z,MA)
      USE FMVALS
      LOGICAL FMLEQ_ZIM,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX Z
      INTEGER KA,NDSAVE
      INTENT (IN) :: Z,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(REAL(Z),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      NDIG = NDSAVE
      L2 = .TRUE.
      IF (AIMAG(Z) /= 0.0) L2 = .FALSE.
      FMLEQ_ZIM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZZM(Z,MA)
      LOGICAL FMLEQ_ZZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMREAL(MTZM,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL ZMIMAG(MTZM,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_CFM(C,MA)
      LOGICAL FMLEQ_CFM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      L1 = FMCOMP(MTFM,'EQ',MA%MFM)
      L2 = .TRUE.
      IF (AIMAG(C) /= 0.0) L2 = .FALSE.
      FMLEQ_CFM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_CIM(C,MA)
      USE FMVALS
      LOGICAL FMLEQ_CIM,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTEGER KA,NDSAVE
      INTENT (IN) :: C,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      NDIG = NDSAVE
      L2 = .TRUE.
      IF (AIMAG(C) /= 0.0) L2 = .FALSE.
      FMLEQ_CIM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_CZM(C,MA)
      LOGICAL FMLEQ_CZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMDP2M(AIMAG(C),MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_CZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_FMI(MA,IVAL)
      LOGICAL FMLEQ_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLEQ_FMI = FMCOMP(MA%MFM,'EQ',MTFM)
   END FUNCTION

   FUNCTION FMLEQ_FMR(MA,R)
      LOGICAL FMLEQ_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLEQ_FMR = FMCOMP(MA%MFM,'EQ',MTFM)
   END FUNCTION

   FUNCTION FMLEQ_FMD(MA,D)
      LOGICAL FMLEQ_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLEQ_FMD = FMCOMP(MA%MFM,'EQ',MTFM)
   END FUNCTION

   FUNCTION FMLEQ_FMZ(MA,Z)
      LOGICAL FMLEQ_FMZ,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL FMSP2M(REAL(Z),MTFM)
      L1 = FMCOMP(MA%MFM,'EQ',MTFM)
      L2 = .TRUE.
      IF (AIMAG(Z) /= 0.0) L2 = .FALSE.
      FMLEQ_FMZ = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_FMC(MA,C)
      LOGICAL FMLEQ_FMC,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      L1 = FMCOMP(MA%MFM,'EQ',MTFM)
      L2 = .TRUE.
      IF (AIMAG(C) /= 0.0) L2 = .FALSE.
      FMLEQ_FMC = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_FMFM(MA,MB)
      LOGICAL FMLEQ_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLEQ_FMFM = FMCOMP(MA%MFM,'EQ',MB%MFM)
   END FUNCTION

   FUNCTION FMLEQ_FMIM(MA,MB)
      LOGICAL FMLEQ_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL FMINT(MA%MFM,MTFM)
      IF (FMCOMP(MA%MFM,'EQ',MTFM)) THEN
          CALL IMI2FM(MB%MIM,MTFM)
          FMLEQ_FMIM = FMCOMP(MA%MFM,'EQ',MTFM)
      ELSE
          FMLEQ_FMIM = .FALSE.
      ENDIF
   END FUNCTION

   FUNCTION FMLEQ_FMZM(MA,MB)
      USE FMVALS
      LOGICAL FMLEQ_FMZM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      TYPE ( ZM ) MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MB%MZM,MTFM)
      L1 = FMCOMP(MA%MFM,'EQ',MTFM)
      L2 = .TRUE.
      IF (MB%MZM(KPTIMU+2) /= 0) L2 = .FALSE.
      FMLEQ_FMZM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_IMI(MA,IVAL)
      LOGICAL FMLEQ_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLEQ_IMI = IMCOMP(MA%MIM,'EQ',MTIM)
   END FUNCTION

   FUNCTION FMLEQ_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLEQ_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLEQ_IMR = FMCOMP(MUFM,'EQ',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLEQ_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLEQ_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLEQ_IMD = FMCOMP(MUFM,'EQ',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLEQ_IMZ(MA,Z)
      USE FMVALS
      LOGICAL FMLEQ_IMZ,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX Z
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,Z
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(REAL(Z),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MUFM,'EQ',MTFM)
      NDIG = NDSAVE
      L2 = .TRUE.
      IF (AIMAG(Z) /= 0.0) L2 = .FALSE.
      FMLEQ_IMZ = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_IMC(MA,C)
      USE FMVALS
      LOGICAL FMLEQ_IMC,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,C
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MUFM,'EQ',MTFM)
      NDIG = NDSAVE
      L2 = .TRUE.
      IF (AIMAG(C) /= 0.0) L2 = .FALSE.
      FMLEQ_IMC = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_IMFM(MA,MB)
      LOGICAL FMLEQ_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTENT (IN) :: MA,MB
      CALL FMINT(MB%MFM,MTFM)
      IF (FMCOMP(MB%MFM,'EQ',MTFM)) THEN
          CALL IMI2FM(MA%MIM,MTFM)
          FMLEQ_IMFM = FMCOMP(MB%MFM,'EQ',MTFM)
      ELSE
          FMLEQ_IMFM = .FALSE.
      ENDIF
   END FUNCTION

   FUNCTION FMLEQ_IMIM(MA,MB)
      LOGICAL FMLEQ_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLEQ_IMIM = IMCOMP(MA%MIM,'EQ',MB%MIM)
   END FUNCTION

   FUNCTION FMLEQ_IMZM(MA,MB)
      USE FMVALS
      LOGICAL FMLEQ_IMZM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( ZM ) MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MB%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MB%MZM(KPTIMU+2) == 0) THEN
          CALL IMI2FM(MA%MIM,MUFM)
          FMLEQ_IMZM = FMCOMP(MUFM,'EQ',MTFM)
      ELSE
          FMLEQ_IMZM = .FALSE.
      ENDIF
   END FUNCTION

   FUNCTION FMLEQ_ZMI(MA,IVAL)
      USE FMVALS
      LOGICAL FMLEQ_ZMI,FMCOMP
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MA%MZM(KPTIMU+2) == 0) THEN
          CALL FMI2M(IVAL,MUFM)
          FMLEQ_ZMI = FMCOMP(MTFM,'EQ',MUFM)
      ELSE
          FMLEQ_ZMI = .FALSE.
      ENDIF
   END FUNCTION

   FUNCTION FMLEQ_ZMR(MA,R)
      LOGICAL FMLEQ_ZMR,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZMR = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZMD(MA,D)
      LOGICAL FMLEQ_ZMD,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZMD = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZMZ(MA,Z)
      LOGICAL FMLEQ_ZMZ,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMREAL(MTZM,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL ZMIMAG(MTZM,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZMZ = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZMC(MA,C)
      LOGICAL FMLEQ_ZMC,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL FMDP2M(AIMAG(C),MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZMC = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLEQ_ZMFM,FMCOMP,L1,L2
      TYPE ( FM ) MB
      TYPE ( ZM ) MA
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      L1 = FMCOMP(MB%MFM,'EQ',MTFM)
      L2 = .TRUE.
      IF (MA%MZM(KPTIMU+2) /= 0) L2 = .FALSE.
      FMLEQ_ZMFM = L1.AND.L2
   END FUNCTION

   FUNCTION FMLEQ_ZMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLEQ_ZMIM,FMCOMP
      TYPE ( IM ) MB
      TYPE ( ZM ) MA
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MA%MZM(KPTIMU+2) == 0) THEN
          CALL IMI2FM(MB%MIM,MUFM)
          FMLEQ_ZMIM = FMCOMP(MUFM,'EQ',MTFM)
      ELSE
          FMLEQ_ZMIM = .FALSE.
      ENDIF
   END FUNCTION

   FUNCTION FMLEQ_ZMZM(MA,MB)
      LOGICAL FMLEQ_ZMZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA,MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      CALL ZMREAL(MB%MZM,MUFM)
      L1 = FMCOMP(MTFM,'EQ',MUFM)
      CALL ZMIMAG(MA%MZM,MTFM)
      CALL ZMIMAG(MB%MZM,MUFM)
      L2 = FMCOMP(MTFM,'EQ',MUFM)
      FMLEQ_ZMZM = L1.AND.L2
   END FUNCTION

 END MODULE FMZM_2

 MODULE FMZM_3
    USE FMZM_1

    INTERFACE OPERATOR ( /= )
       MODULE PROCEDURE FMLNE_IFM
       MODULE PROCEDURE FMLNE_IIM
       MODULE PROCEDURE FMLNE_IZM
       MODULE PROCEDURE FMLNE_RFM
       MODULE PROCEDURE FMLNE_RIM
       MODULE PROCEDURE FMLNE_RZM
       MODULE PROCEDURE FMLNE_DFM
       MODULE PROCEDURE FMLNE_DIM
       MODULE PROCEDURE FMLNE_DZM
       MODULE PROCEDURE FMLNE_ZFM
       MODULE PROCEDURE FMLNE_ZIM
       MODULE PROCEDURE FMLNE_ZZM
       MODULE PROCEDURE FMLNE_CFM
       MODULE PROCEDURE FMLNE_CIM
       MODULE PROCEDURE FMLNE_CZM
       MODULE PROCEDURE FMLNE_FMI
       MODULE PROCEDURE FMLNE_FMR
       MODULE PROCEDURE FMLNE_FMD
       MODULE PROCEDURE FMLNE_FMZ
       MODULE PROCEDURE FMLNE_FMC
       MODULE PROCEDURE FMLNE_FMFM
       MODULE PROCEDURE FMLNE_FMIM
       MODULE PROCEDURE FMLNE_FMZM
       MODULE PROCEDURE FMLNE_IMI
       MODULE PROCEDURE FMLNE_IMR
       MODULE PROCEDURE FMLNE_IMD
       MODULE PROCEDURE FMLNE_IMZ
       MODULE PROCEDURE FMLNE_IMC
       MODULE PROCEDURE FMLNE_IMFM
       MODULE PROCEDURE FMLNE_IMIM
       MODULE PROCEDURE FMLNE_IMZM
       MODULE PROCEDURE FMLNE_ZMI
       MODULE PROCEDURE FMLNE_ZMR
       MODULE PROCEDURE FMLNE_ZMD
       MODULE PROCEDURE FMLNE_ZMZ
       MODULE PROCEDURE FMLNE_ZMC
       MODULE PROCEDURE FMLNE_ZMFM
       MODULE PROCEDURE FMLNE_ZMIM
       MODULE PROCEDURE FMLNE_ZMZM
    END INTERFACE

    INTERFACE OPERATOR ( > )
       MODULE PROCEDURE FMLGT_IFM
       MODULE PROCEDURE FMLGT_IIM
       MODULE PROCEDURE FMLGT_RFM
       MODULE PROCEDURE FMLGT_RIM
       MODULE PROCEDURE FMLGT_DFM
       MODULE PROCEDURE FMLGT_DIM
       MODULE PROCEDURE FMLGT_FMI
       MODULE PROCEDURE FMLGT_FMR
       MODULE PROCEDURE FMLGT_FMD
       MODULE PROCEDURE FMLGT_FMFM
       MODULE PROCEDURE FMLGT_FMIM
       MODULE PROCEDURE FMLGT_IMI
       MODULE PROCEDURE FMLGT_IMR
       MODULE PROCEDURE FMLGT_IMD
       MODULE PROCEDURE FMLGT_IMFM
       MODULE PROCEDURE FMLGT_IMIM
    END INTERFACE

 CONTAINS

!                                                                /=

   FUNCTION FMLNE_IFM(IVAL,MA)
      LOGICAL FMLNE_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLNE_IFM = FMCOMP(MTFM,'NE',MA%MFM)
   END FUNCTION

   FUNCTION FMLNE_IIM(IVAL,MA)
      LOGICAL FMLNE_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLNE_IIM = IMCOMP(MTIM,'NE',MA%MIM)
   END FUNCTION

   FUNCTION FMLNE_IZM(IVAL,MA)
      LOGICAL FMLNE_IZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_IZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_RFM(R,MA)
      LOGICAL FMLNE_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLNE_RFM = FMCOMP(MTFM,'NE',MA%MFM)
   END FUNCTION

   FUNCTION FMLNE_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLNE_RIM,FMCOMP
      TYPE ( IM ) MA
      REAL R
      INTEGER KA,NDSAVE
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLNE_RIM = FMCOMP(MTFM,'NE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLNE_RZM(R,MA)
      LOGICAL FMLNE_RZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_RZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_DFM(D,MA)
      LOGICAL FMLNE_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLNE_DFM = FMCOMP(MTFM,'NE',MA%MFM)
   END FUNCTION

   FUNCTION FMLNE_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLNE_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLNE_DIM = FMCOMP(MTFM,'NE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLNE_DZM(D,MA)
      LOGICAL FMLNE_DZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_DZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZFM(Z,MA)
      LOGICAL FMLNE_ZFM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL FMSP2M(REAL(Z),MTFM)
      L1 = FMCOMP(MTFM,'NE',MA%MFM)
      L2 = .FALSE.
      IF (AIMAG(Z) /= 0.0) L2 = .TRUE.
      FMLNE_ZFM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZIM(Z,MA)
      USE FMVALS
      LOGICAL FMLNE_ZIM,FMCOMP,L1,L2
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      COMPLEX Z
      INTENT (IN) :: Z,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(REAL(Z),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      NDIG = NDSAVE
      L2 = .FALSE.
      IF (AIMAG(Z) /= 0.0) L2 = .TRUE.
      FMLNE_ZIM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZZM(Z,MA)
      LOGICAL FMLNE_ZZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMREAL(MTZM,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL ZMIMAG(MTZM,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_CFM(C,MA)
      LOGICAL FMLNE_CFM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      L1 = FMCOMP(MTFM,'NE',MA%MFM)
      L2 = .FALSE.
      IF (AIMAG(C) /= 0.0) L2 = .TRUE.
      FMLNE_CFM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_CIM(C,MA)
      USE FMVALS
      LOGICAL FMLNE_CIM,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTEGER KA,NDSAVE
      INTENT (IN) :: C,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      NDIG = NDSAVE
      L2 = .FALSE.
      IF (AIMAG(C) /= 0.0) L2 = .TRUE.
      FMLNE_CIM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_CZM(C,MA)
      LOGICAL FMLNE_CZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMDP2M(AIMAG(C),MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_CZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_FMI(MA,IVAL)
      LOGICAL FMLNE_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLNE_FMI = FMCOMP(MA%MFM,'NE',MTFM)
   END FUNCTION

   FUNCTION FMLNE_FMR(MA,R)
      LOGICAL FMLNE_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLNE_FMR = FMCOMP(MA%MFM,'NE',MTFM)
   END FUNCTION

   FUNCTION FMLNE_FMD(MA,D)
      LOGICAL FMLNE_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLNE_FMD = FMCOMP(MA%MFM,'NE',MTFM)
   END FUNCTION

   FUNCTION FMLNE_FMZ(MA,Z)
      LOGICAL FMLNE_FMZ,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL FMSP2M(REAL(Z),MTFM)
      L1 = FMCOMP(MA%MFM,'NE',MTFM)
      L2 = .FALSE.
      IF (AIMAG(Z) /= 0.0) L2 = .TRUE.
      FMLNE_FMZ = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_FMC(MA,C)
      LOGICAL FMLNE_FMC,FMCOMP,L1,L2
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      L1 = FMCOMP(MA%MFM,'NE',MTFM)
      L2 = .FALSE.
      IF (AIMAG(C) /= 0.0) L2 = .TRUE.
      FMLNE_FMC = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_FMFM(MA,MB)
      LOGICAL FMLNE_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLNE_FMFM = FMCOMP(MA%MFM,'NE',MB%MFM)
   END FUNCTION

   FUNCTION FMLNE_FMIM(MA,MB)
      LOGICAL FMLNE_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL FMINT(MA%MFM,MTFM)
      IF (FMCOMP(MA%MFM,'EQ',MTFM)) THEN
          CALL IMI2FM(MB%MIM,MTFM)
          FMLNE_FMIM = FMCOMP(MA%MFM,'NE',MTFM)
      ELSE
          FMLNE_FMIM = .TRUE.
      ENDIF
   END FUNCTION

   FUNCTION FMLNE_FMZM(MA,MB)
      USE FMVALS
      LOGICAL FMLNE_FMZM,FMCOMP,L1,L2
      TYPE ( FM ) MA
      TYPE ( ZM ) MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MB%MZM,MTFM)
      L1 = FMCOMP(MA%MFM,'NE',MTFM)
      L2 = .FALSE.
      IF (MB%MZM(KPTIMU+2) /= 0) L2 = .TRUE.
      FMLNE_FMZM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_IMI(MA,IVAL)
      LOGICAL FMLNE_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLNE_IMI = IMCOMP(MA%MIM,'NE',MTIM)
   END FUNCTION

   FUNCTION FMLNE_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLNE_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLNE_IMR = FMCOMP(MUFM,'NE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLNE_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLNE_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLNE_IMD = FMCOMP(MUFM,'NE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLNE_IMZ(MA,Z)
      USE FMVALS
      LOGICAL FMLNE_IMZ,FMCOMP,L1,L2
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      COMPLEX Z
      INTENT (IN) :: MA,Z
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(REAL(Z),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MUFM,'NE',MTFM)
      NDIG = NDSAVE
      L2 = .FALSE.
      IF (AIMAG(Z) /= 0.0) L2 = .TRUE.
      FMLNE_IMZ = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_IMC(MA,C)
      USE FMVALS
      LOGICAL FMLNE_IMC,FMCOMP,L1,L2
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,C
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      L1 = FMCOMP(MUFM,'NE',MTFM)
      NDIG = NDSAVE
      L2 = .FALSE.
      IF (AIMAG(C) /= 0.0) L2 = .TRUE.
      FMLNE_IMC = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_IMFM(MA,MB)
      LOGICAL FMLNE_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTENT (IN) :: MA,MB
      CALL FMINT(MB%MFM,MTFM)
      IF (FMCOMP(MB%MFM,'EQ',MTFM)) THEN
          CALL IMI2FM(MA%MIM,MTFM)
          FMLNE_IMFM = FMCOMP(MB%MFM,'NE',MTFM)
      ELSE
          FMLNE_IMFM = .TRUE.
      ENDIF
   END FUNCTION

   FUNCTION FMLNE_IMIM(MA,MB)
      LOGICAL FMLNE_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLNE_IMIM = IMCOMP(MA%MIM,'NE',MB%MIM)
   END FUNCTION

   FUNCTION FMLNE_IMZM(MA,MB)
      USE FMVALS
      LOGICAL FMLNE_IMZM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( ZM ) MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MB%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MB%MZM(KPTIMU+2) == 0) THEN
          CALL IMI2FM(MA%MIM,MUFM)
          FMLNE_IMZM = FMCOMP(MUFM,'NE',MTFM)
      ELSE
          FMLNE_IMZM = .TRUE.
      ENDIF
   END FUNCTION

   FUNCTION FMLNE_ZMI(MA,IVAL)
      USE FMVALS
      LOGICAL FMLNE_ZMI,FMCOMP
      TYPE ( ZM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MA%MZM(KPTIMU+2) == 0) THEN
          CALL FMI2M(IVAL,MUFM)
          FMLNE_ZMI = FMCOMP(MTFM,'NE',MUFM)
      ELSE
          FMLNE_ZMI = .TRUE.
      ENDIF
   END FUNCTION

   FUNCTION FMLNE_ZMR(MA,R)
      LOGICAL FMLNE_ZMR,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZMR = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZMD(MA,D)
      LOGICAL FMLNE_ZMD,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMI2M(0,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZMD = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZMZ(MA,Z)
      LOGICAL FMLNE_ZMZ,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMREAL(MTZM,MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL ZMIMAG(MTZM,MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZMZ = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZMC(MA,C)
      LOGICAL FMLNE_ZMC,FMCOMP,L1,L2
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL ZMREAL(MA%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL FMDP2M(AIMAG(C),MTFM)
      CALL ZMIMAG(MA%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZMC = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLNE_ZMFM,FMCOMP,L1,L2
      TYPE ( FM ) MB
      TYPE ( ZM ) MA
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      L1 = FMCOMP(MB%MFM,'NE',MTFM)
      L2 = .FALSE.
      IF (MA%MZM(KPTIMU+2) /= 0) L2 = .TRUE.
      FMLNE_ZMFM = L1.OR.L2
   END FUNCTION

   FUNCTION FMLNE_ZMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLNE_ZMIM,FMCOMP
      TYPE ( IM ) MB
      TYPE ( ZM ) MA
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      IF (FMCOMP(MUFM,'EQ',MTFM).AND.MA%MZM(KPTIMU+2) == 0) THEN
          CALL IMI2FM(MB%MIM,MUFM)
          FMLNE_ZMIM = FMCOMP(MUFM,'NE',MTFM)
      ELSE
          FMLNE_ZMIM = .TRUE.
      ENDIF
   END FUNCTION

   FUNCTION FMLNE_ZMZM(MA,MB)
      LOGICAL FMLNE_ZMZM,FMCOMP,L1,L2
      TYPE ( ZM ) MA,MB
      INTENT (IN) :: MA,MB
      CALL ZMREAL(MA%MZM,MTFM)
      CALL ZMREAL(MB%MZM,MUFM)
      L1 = FMCOMP(MTFM,'NE',MUFM)
      CALL ZMIMAG(MA%MZM,MTFM)
      CALL ZMIMAG(MB%MZM,MUFM)
      L2 = FMCOMP(MTFM,'NE',MUFM)
      FMLNE_ZMZM = L1.OR.L2
   END FUNCTION

!                                                                >

   FUNCTION FMLGT_IFM(IVAL,MA)
      LOGICAL FMLGT_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLGT_IFM = FMCOMP(MTFM,'GT',MA%MFM)
   END FUNCTION

   FUNCTION FMLGT_IIM(IVAL,MA)
      LOGICAL FMLGT_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLGT_IIM = IMCOMP(MTIM,'GT',MA%MIM)
   END FUNCTION

   FUNCTION FMLGT_RFM(R,MA)
      LOGICAL FMLGT_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLGT_RFM = FMCOMP(MTFM,'GT',MA%MFM)
   END FUNCTION

   FUNCTION FMLGT_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLGT_RIM,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGT_RIM = FMCOMP(MTFM,'GT',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_DFM(D,MA)
      LOGICAL FMLGT_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLGT_DFM = FMCOMP(MTFM,'GT',MA%MFM)
   END FUNCTION

   FUNCTION FMLGT_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLGT_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGT_DIM = FMCOMP(MTFM,'GT',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_FMI(MA,IVAL)
      LOGICAL FMLGT_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLGT_FMI = FMCOMP(MA%MFM,'GT',MTFM)
   END FUNCTION

   FUNCTION FMLGT_FMR(MA,R)
      LOGICAL FMLGT_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLGT_FMR = FMCOMP(MA%MFM,'GT',MTFM)
   END FUNCTION

   FUNCTION FMLGT_FMD(MA,D)
      LOGICAL FMLGT_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLGT_FMD = FMCOMP(MA%MFM,'GT',MTFM)
   END FUNCTION

   FUNCTION FMLGT_FMFM(MA,MB)
      LOGICAL FMLGT_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLGT_FMFM = FMCOMP(MA%MFM,'GT',MB%MFM)
   END FUNCTION

   FUNCTION FMLGT_FMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLGT_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MB%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MB%MIM,MTFM)
      FMLGT_FMIM = FMCOMP(MA%MFM,'GT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_IMI(MA,IVAL)
      LOGICAL FMLGT_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLGT_IMI = IMCOMP(MA%MIM,'GT',MTIM)
   END FUNCTION

   FUNCTION FMLGT_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLGT_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGT_IMR = FMCOMP(MUFM,'GT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLGT_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGT_IMD = FMCOMP(MUFM,'GT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_IMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLGT_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MA%MIM,MTFM)
      FMLGT_IMFM = FMCOMP(MTFM,'GT',MB%MFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGT_IMIM(MA,MB)
      LOGICAL FMLGT_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLGT_IMIM = IMCOMP(MA%MIM,'GT',MB%MIM)
   END FUNCTION

 END MODULE FMZM_3

 MODULE FMZM_4
    USE FMZM_1

    INTERFACE OPERATOR ( >= )
       MODULE PROCEDURE FMLGE_IFM
       MODULE PROCEDURE FMLGE_IIM
       MODULE PROCEDURE FMLGE_RFM
       MODULE PROCEDURE FMLGE_RIM
       MODULE PROCEDURE FMLGE_DFM
       MODULE PROCEDURE FMLGE_DIM
       MODULE PROCEDURE FMLGE_FMI
       MODULE PROCEDURE FMLGE_FMR
       MODULE PROCEDURE FMLGE_FMD
       MODULE PROCEDURE FMLGE_FMFM
       MODULE PROCEDURE FMLGE_FMIM
       MODULE PROCEDURE FMLGE_IMI
       MODULE PROCEDURE FMLGE_IMR
       MODULE PROCEDURE FMLGE_IMD
       MODULE PROCEDURE FMLGE_IMFM
       MODULE PROCEDURE FMLGE_IMIM
    END INTERFACE

    INTERFACE OPERATOR ( < )
       MODULE PROCEDURE FMLLT_IFM
       MODULE PROCEDURE FMLLT_IIM
       MODULE PROCEDURE FMLLT_RFM
       MODULE PROCEDURE FMLLT_RIM
       MODULE PROCEDURE FMLLT_DFM
       MODULE PROCEDURE FMLLT_DIM
       MODULE PROCEDURE FMLLT_FMI
       MODULE PROCEDURE FMLLT_FMR
       MODULE PROCEDURE FMLLT_FMD
       MODULE PROCEDURE FMLLT_FMFM
       MODULE PROCEDURE FMLLT_FMIM
       MODULE PROCEDURE FMLLT_IMI
       MODULE PROCEDURE FMLLT_IMR
       MODULE PROCEDURE FMLLT_IMD
       MODULE PROCEDURE FMLLT_IMFM
       MODULE PROCEDURE FMLLT_IMIM
    END INTERFACE

    INTERFACE OPERATOR ( <= )
       MODULE PROCEDURE FMLLE_IFM
       MODULE PROCEDURE FMLLE_IIM
       MODULE PROCEDURE FMLLE_RFM
       MODULE PROCEDURE FMLLE_RIM
       MODULE PROCEDURE FMLLE_DFM
       MODULE PROCEDURE FMLLE_DIM
       MODULE PROCEDURE FMLLE_FMI
       MODULE PROCEDURE FMLLE_FMR
       MODULE PROCEDURE FMLLE_FMD
       MODULE PROCEDURE FMLLE_FMFM
       MODULE PROCEDURE FMLLE_FMIM
       MODULE PROCEDURE FMLLE_IMI
       MODULE PROCEDURE FMLLE_IMR
       MODULE PROCEDURE FMLLE_IMD
       MODULE PROCEDURE FMLLE_IMFM
       MODULE PROCEDURE FMLLE_IMIM
    END INTERFACE

 CONTAINS

!                                                                >=

   FUNCTION FMLGE_IFM(IVAL,MA)
      LOGICAL FMLGE_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLGE_IFM = FMCOMP(MTFM,'GE',MA%MFM)
   END FUNCTION

   FUNCTION FMLGE_IIM(IVAL,MA)
      LOGICAL FMLGE_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLGE_IIM = IMCOMP(MTIM,'GE',MA%MIM)
   END FUNCTION

   FUNCTION FMLGE_RFM(R,MA)
      LOGICAL FMLGE_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLGE_RFM = FMCOMP(MTFM,'GE',MA%MFM)
   END FUNCTION

   FUNCTION FMLGE_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLGE_RIM,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGE_RIM = FMCOMP(MTFM,'GE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_DFM(D,MA)
      LOGICAL FMLGE_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLGE_DFM = FMCOMP(MTFM,'GE',MA%MFM)
   END FUNCTION

   FUNCTION FMLGE_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLGE_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGE_DIM = FMCOMP(MTFM,'GE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_FMI(MA,IVAL)
      LOGICAL FMLGE_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLGE_FMI = FMCOMP(MA%MFM,'GE',MTFM)
   END FUNCTION

   FUNCTION FMLGE_FMR(MA,R)
      LOGICAL FMLGE_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLGE_FMR = FMCOMP(MA%MFM,'GE',MTFM)
   END FUNCTION

   FUNCTION FMLGE_FMD(MA,D)
      LOGICAL FMLGE_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLGE_FMD = FMCOMP(MA%MFM,'GE',MTFM)
   END FUNCTION

   FUNCTION FMLGE_FMFM(MA,MB)
      LOGICAL FMLGE_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLGE_FMFM = FMCOMP(MA%MFM,'GE',MB%MFM)
   END FUNCTION

   FUNCTION FMLGE_FMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLGE_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MB%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MB%MIM,MTFM)
      FMLGE_FMIM = FMCOMP(MA%MFM,'GE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_IMI(MA,IVAL)
      LOGICAL FMLGE_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLGE_IMI = IMCOMP(MA%MIM,'GE',MTIM)
   END FUNCTION

   FUNCTION FMLGE_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLGE_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGE_IMR = FMCOMP(MUFM,'GE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLGE_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLGE_IMD = FMCOMP(MUFM,'GE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_IMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLGE_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MA%MIM,MTFM)
      FMLGE_IMFM = FMCOMP(MTFM,'GE',MB%MFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLGE_IMIM(MA,MB)
      LOGICAL FMLGE_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLGE_IMIM = IMCOMP(MA%MIM,'GE',MB%MIM)
   END FUNCTION

!                                                                <

   FUNCTION FMLLT_IFM(IVAL,MA)
      LOGICAL FMLLT_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLLT_IFM = FMCOMP(MTFM,'LT',MA%MFM)
   END FUNCTION

   FUNCTION FMLLT_IIM(IVAL,MA)
      LOGICAL FMLLT_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLLT_IIM = IMCOMP(MTIM,'LT',MA%MIM)
   END FUNCTION

   FUNCTION FMLLT_RFM(R,MA)
      LOGICAL FMLLT_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLLT_RFM = FMCOMP(MTFM,'LT',MA%MFM)
   END FUNCTION

   FUNCTION FMLLT_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLLT_RIM,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLT_RIM = FMCOMP(MTFM,'LT',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_DFM(D,MA)
      LOGICAL FMLLT_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLLT_DFM = FMCOMP(MTFM,'LT',MA%MFM)
   END FUNCTION

   FUNCTION FMLLT_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLLT_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLT_DIM = FMCOMP(MTFM,'LT',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_FMI(MA,IVAL)
      LOGICAL FMLLT_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLLT_FMI = FMCOMP(MA%MFM,'LT',MTFM)
   END FUNCTION

   FUNCTION FMLLT_FMR(MA,R)
      LOGICAL FMLLT_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLLT_FMR = FMCOMP(MA%MFM,'LT',MTFM)
   END FUNCTION

   FUNCTION FMLLT_FMD(MA,D)
      LOGICAL FMLLT_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLLT_FMD = FMCOMP(MA%MFM,'LT',MTFM)
   END FUNCTION

   FUNCTION FMLLT_FMFM(MA,MB)
      LOGICAL FMLLT_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLLT_FMFM = FMCOMP(MA%MFM,'LT',MB%MFM)
   END FUNCTION

   FUNCTION FMLLT_FMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLLT_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MB%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MB%MIM,MTFM)
      FMLLT_FMIM = FMCOMP(MA%MFM,'LT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_IMI(MA,IVAL)
      LOGICAL FMLLT_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLLT_IMI = IMCOMP(MA%MIM,'LT',MTIM)
   END FUNCTION

   FUNCTION FMLLT_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLLT_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLT_IMR = FMCOMP(MUFM,'LT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLLT_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLT_IMD = FMCOMP(MUFM,'LT',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_IMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLLT_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MA%MIM,MTFM)
      FMLLT_IMFM = FMCOMP(MTFM,'LT',MB%MFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLT_IMIM(MA,MB)
      LOGICAL FMLLT_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLLT_IMIM = IMCOMP(MA%MIM,'LT',MB%MIM)
   END FUNCTION

!                                                                <=

   FUNCTION FMLLE_IFM(IVAL,MA)
      LOGICAL FMLLE_IFM,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      FMLLE_IFM = FMCOMP(MTFM,'LE',MA%MFM)
   END FUNCTION

   FUNCTION FMLLE_IIM(IVAL,MA)
      LOGICAL FMLLE_IIM,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      FMLLE_IIM = IMCOMP(MTIM,'LE',MA%MIM)
   END FUNCTION

   FUNCTION FMLLE_RFM(R,MA)
      LOGICAL FMLLE_RFM,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      FMLLE_RFM = FMCOMP(MTFM,'LE',MA%MFM)
   END FUNCTION

   FUNCTION FMLLE_RIM(R,MA)
      USE FMVALS
      LOGICAL FMLLE_RIM,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: R,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLE_RIM = FMCOMP(MTFM,'LE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_DFM(D,MA)
      LOGICAL FMLLE_DFM,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      FMLLE_DFM = FMCOMP(MTFM,'LE',MA%MFM)
   END FUNCTION

   FUNCTION FMLLE_DIM(D,MA)
      USE FMVALS
      LOGICAL FMLLE_DIM,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: D,MA
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLE_DIM = FMCOMP(MTFM,'LE',MUFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_FMI(MA,IVAL)
      LOGICAL FMLLE_FMI,FMCOMP
      TYPE ( FM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      FMLLE_FMI = FMCOMP(MA%MFM,'LE',MTFM)
   END FUNCTION

   FUNCTION FMLLE_FMR(MA,R)
      LOGICAL FMLLE_FMR,FMCOMP
      TYPE ( FM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      FMLLE_FMR = FMCOMP(MA%MFM,'LE',MTFM)
   END FUNCTION

   FUNCTION FMLLE_FMD(MA,D)
      LOGICAL FMLLE_FMD,FMCOMP
      TYPE ( FM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      FMLLE_FMD = FMCOMP(MA%MFM,'LE',MTFM)
   END FUNCTION

   FUNCTION FMLLE_FMFM(MA,MB)
      LOGICAL FMLLE_FMFM,FMCOMP
      TYPE ( FM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLLE_FMFM = FMCOMP(MA%MFM,'LE',MB%MFM)
   END FUNCTION

   FUNCTION FMLLE_FMIM(MA,MB)
      USE FMVALS
      LOGICAL FMLLE_FMIM,FMCOMP
      TYPE ( FM ) MA
      TYPE ( IM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MB%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MB%MIM,MTFM)
      FMLLE_FMIM = FMCOMP(MA%MFM,'LE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_IMI(MA,IVAL)
      LOGICAL FMLLE_IMI,IMCOMP
      TYPE ( IM ) MA
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      FMLLE_IMI = IMCOMP(MA%MIM,'LE',MTIM)
   END FUNCTION

   FUNCTION FMLLE_IMR(MA,R)
      USE FMVALS
      LOGICAL FMLLE_IMR,FMCOMP
      TYPE ( IM ) MA
      INTEGER KA,NDSAVE
      REAL R
      INTENT (IN) :: MA,R
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLE_IMR = FMCOMP(MUFM,'LE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_IMD(MA,D)
      USE FMVALS
      LOGICAL FMLLE_IMD,FMCOMP
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,D
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      FMLLE_IMD = FMCOMP(MUFM,'LE',MTFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_IMFM(MA,MB)
      USE FMVALS
      LOGICAL FMLLE_IMFM,FMCOMP
      TYPE ( IM ) MA
      TYPE ( FM ) MB
      INTEGER KA,NDSAVE
      INTENT (IN) :: MA,MB
      NDSAVE = NDIG
      KA = MA%MIM(1)
      NDIG = MAX(KA+NGRD52,NDIG)
      CALL IMI2FM(MA%MIM,MTFM)
      FMLLE_IMFM = FMCOMP(MTFM,'LE',MB%MFM)
      NDIG = NDSAVE
   END FUNCTION

   FUNCTION FMLLE_IMIM(MA,MB)
      LOGICAL FMLLE_IMIM,IMCOMP
      TYPE ( IM ) MA,MB
      INTENT (IN) :: MA,MB
      FMLLE_IMIM = IMCOMP(MA%MIM,'LE',MB%MIM)
   END FUNCTION

 END MODULE FMZM_4

 MODULE FMZM_5
    USE FMZM_1

    INTERFACE OPERATOR (+)
       MODULE PROCEDURE FMADD_IFM
       MODULE PROCEDURE FMADD_IIM
       MODULE PROCEDURE FMADD_IZM
       MODULE PROCEDURE FMADD_RFM
       MODULE PROCEDURE FMADD_RIM
       MODULE PROCEDURE FMADD_RZM
       MODULE PROCEDURE FMADD_DFM
       MODULE PROCEDURE FMADD_DIM
       MODULE PROCEDURE FMADD_DZM
       MODULE PROCEDURE FMADD_ZFM
       MODULE PROCEDURE FMADD_ZIM
       MODULE PROCEDURE FMADD_ZZM
       MODULE PROCEDURE FMADD_CFM
       MODULE PROCEDURE FMADD_CIM
       MODULE PROCEDURE FMADD_CZM
       MODULE PROCEDURE FMADD_FMI
       MODULE PROCEDURE FMADD_FMR
       MODULE PROCEDURE FMADD_FMD
       MODULE PROCEDURE FMADD_FMZ
       MODULE PROCEDURE FMADD_FMC
       MODULE PROCEDURE FMADD_FMFM
       MODULE PROCEDURE FMADD_FMIM
       MODULE PROCEDURE FMADD_FMZM
       MODULE PROCEDURE FMADD_IMI
       MODULE PROCEDURE FMADD_IMR
       MODULE PROCEDURE FMADD_IMD
       MODULE PROCEDURE FMADD_IMZ
       MODULE PROCEDURE FMADD_IMC
       MODULE PROCEDURE FMADD_IMFM
       MODULE PROCEDURE FMADD_IMIM
       MODULE PROCEDURE FMADD_IMZM
       MODULE PROCEDURE FMADD_ZMI
       MODULE PROCEDURE FMADD_ZMR
       MODULE PROCEDURE FMADD_ZMD
       MODULE PROCEDURE FMADD_ZMZ
       MODULE PROCEDURE FMADD_ZMC
       MODULE PROCEDURE FMADD_ZMFM
       MODULE PROCEDURE FMADD_ZMIM
       MODULE PROCEDURE FMADD_ZMZM
       MODULE PROCEDURE FMADD_FM
       MODULE PROCEDURE FMADD_IM
       MODULE PROCEDURE FMADD_ZM
    END INTERFACE

    INTERFACE OPERATOR (-)
       MODULE PROCEDURE FMSUB_IFM
       MODULE PROCEDURE FMSUB_IIM
       MODULE PROCEDURE FMSUB_IZM
       MODULE PROCEDURE FMSUB_RFM
       MODULE PROCEDURE FMSUB_RIM
       MODULE PROCEDURE FMSUB_RZM
       MODULE PROCEDURE FMSUB_DFM
       MODULE PROCEDURE FMSUB_DIM
       MODULE PROCEDURE FMSUB_DZM
       MODULE PROCEDURE FMSUB_ZFM
       MODULE PROCEDURE FMSUB_ZIM
       MODULE PROCEDURE FMSUB_ZZM
       MODULE PROCEDURE FMSUB_CFM
       MODULE PROCEDURE FMSUB_CIM
       MODULE PROCEDURE FMSUB_CZM
       MODULE PROCEDURE FMSUB_FMI
       MODULE PROCEDURE FMSUB_FMR
       MODULE PROCEDURE FMSUB_FMD
       MODULE PROCEDURE FMSUB_FMZ
       MODULE PROCEDURE FMSUB_FMC
       MODULE PROCEDURE FMSUB_FMFM
       MODULE PROCEDURE FMSUB_FMIM
       MODULE PROCEDURE FMSUB_FMZM
       MODULE PROCEDURE FMSUB_IMI
       MODULE PROCEDURE FMSUB_IMR
       MODULE PROCEDURE FMSUB_IMD
       MODULE PROCEDURE FMSUB_IMZ
       MODULE PROCEDURE FMSUB_IMC
       MODULE PROCEDURE FMSUB_IMFM
       MODULE PROCEDURE FMSUB_IMIM
       MODULE PROCEDURE FMSUB_IMZM
       MODULE PROCEDURE FMSUB_ZMI
       MODULE PROCEDURE FMSUB_ZMR
       MODULE PROCEDURE FMSUB_ZMD
       MODULE PROCEDURE FMSUB_ZMZ
       MODULE PROCEDURE FMSUB_ZMC
       MODULE PROCEDURE FMSUB_ZMFM
       MODULE PROCEDURE FMSUB_ZMIM
       MODULE PROCEDURE FMSUB_ZMZM
       MODULE PROCEDURE FMSUB_FM
       MODULE PROCEDURE FMSUB_IM
       MODULE PROCEDURE FMSUB_ZM
    END INTERFACE

 CONTAINS

!                                                                   +

   FUNCTION FMADD_IFM(IVAL,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_IFM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMADD(MTFM,MA%MFM,FMADD_IFM%MFM)
   END FUNCTION

   FUNCTION FMADD_IIM(IVAL,MA)
      USE FMVALS
      TYPE ( IM ) MA,FMADD_IIM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      CALL IMADD(MTIM,MA%MIM,FMADD_IIM%MIM)
   END FUNCTION

   FUNCTION FMADD_IZM(IVAL,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_IZM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MTZM,MA%MZM,FMADD_IZM%MZM)
   END FUNCTION

   FUNCTION FMADD_RFM(R,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_RFM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMADD(MTFM,MA%MFM,FMADD_RFM%MFM)
   END FUNCTION

   FUNCTION FMADD_RIM(R,MA)
      USE FMVALS
      TYPE ( FM ) FMADD_RIM
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMADD(MTFM,MUFM,FMADD_RIM%MFM)
   END FUNCTION

   FUNCTION FMADD_RZM(R,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_RZM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MTZM,MA%MZM,FMADD_RZM%MZM)
   END FUNCTION

   FUNCTION FMADD_DFM(D,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_DFM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMADD(MTFM,MA%MFM,FMADD_DFM%MFM)
   END FUNCTION

   FUNCTION FMADD_DIM(D,MA)
      USE FMVALS
      TYPE ( FM ) FMADD_DIM
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMADD(MTFM,MUFM,FMADD_DIM%MFM)
   END FUNCTION

   FUNCTION FMADD_DZM(D,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_DZM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MTZM,MA%MZM,FMADD_DZM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZFM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMADD_ZFM
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMADD(MTZM,MUZM,FMADD_ZFM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZIM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMADD_ZIM
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MTZM,MUZM,FMADD_ZIM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZZM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZZM
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMADD(MTZM,MA%MZM,FMADD_ZZM%MZM)
   END FUNCTION

   FUNCTION FMADD_CFM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMADD_CFM
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMADD(MTZM,MUZM,FMADD_CFM%MZM)
   END FUNCTION

   FUNCTION FMADD_CIM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMADD_CIM
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MTZM,MUZM,FMADD_CIM%MZM)
   END FUNCTION

   FUNCTION FMADD_CZM(C,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_CZM
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MTZM,MA%MZM,FMADD_CZM%MZM)
   END FUNCTION

   FUNCTION FMADD_FMI(MA,IVAL)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_FMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      CALL FMADD(MA%MFM,MTFM,FMADD_FMI%MFM)
   END FUNCTION

   FUNCTION FMADD_FMR(MA,R)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_FMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMADD(MA%MFM,MTFM,FMADD_FMR%MFM)
   END FUNCTION

   FUNCTION FMADD_FMD(MA,D)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_FMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMADD(MA%MFM,MTFM,FMADD_FMD%MFM)
   END FUNCTION

   FUNCTION FMADD_FMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMADD_FMZ
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMADD(MUZM,MTZM,FMADD_FMZ%MZM)
   END FUNCTION

   FUNCTION FMADD_FMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMADD_FMC
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMADD(MUZM,MTZM,FMADD_FMC%MZM)
   END FUNCTION

   FUNCTION FMADD_FMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMADD_FMFM
      INTENT (IN) :: MA,MB
      CALL FMADD(MA%MFM,MB%MFM,FMADD_FMFM%MFM)
   END FUNCTION

   FUNCTION FMADD_FMIM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_FMIM
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMADD(MA%MFM,MTFM,FMADD_FMIM%MFM)
   END FUNCTION

   FUNCTION FMADD_FMZM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( ZM ) MB,FMADD_FMZM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MA%MFM,MTFM,MTZM)
      CALL ZMADD(MTZM,MB%MZM,FMADD_FMZM%MZM)
   END FUNCTION

   FUNCTION FMADD_IMI(MA,IVAL)
      USE FMVALS
      TYPE ( IM ) MA,FMADD_IMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      CALL IMADD(MA%MIM,MTIM,FMADD_IMI%MIM)
   END FUNCTION

   FUNCTION FMADD_IMR(MA,R)
      USE FMVALS
      TYPE ( FM ) FMADD_IMR
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMADD(MUFM,MTFM,FMADD_IMR%MFM)
   END FUNCTION

   FUNCTION FMADD_IMD(MA,D)
      USE FMVALS
      TYPE ( FM ) FMADD_IMD
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMADD(MUFM,MTFM,FMADD_IMD%MFM)
   END FUNCTION

   FUNCTION FMADD_IMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMADD_IMZ
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MUZM,MTZM,FMADD_IMZ%MZM)
   END FUNCTION

   FUNCTION FMADD_IMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMADD_IMC
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MUZM,MTZM,FMADD_IMC%MZM)
   END FUNCTION

   FUNCTION FMADD_IMFM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( FM ) MB,FMADD_IMFM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMADD(MTFM,MB%MFM,FMADD_IMFM%MFM)
   END FUNCTION

   FUNCTION FMADD_IMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMADD_IMIM
      INTENT (IN) :: MA,MB
      CALL IMADD(MA%MIM,MB%MIM,FMADD_IMIM%MIM)
   END FUNCTION

   FUNCTION FMADD_IMZM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( ZM ) MB,FMADD_IMZM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MUZM,MB%MZM,FMADD_IMZM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMI(MA,IVAL)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMI%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMR(MA,R)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMR%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMD(MA,D)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMD%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZMZ
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMZ%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZMC
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMC%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MB
      TYPE ( ZM ) MA,FMADD_ZMFM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MTZM)
      CALL ZMADD(MA%MZM,MTZM,FMADD_ZMFM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MB
      TYPE ( ZM ) MA,FMADD_ZMIM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMADD(MA%MZM,MUZM,FMADD_ZMIM%MZM)
   END FUNCTION

   FUNCTION FMADD_ZMZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA,MB,FMADD_ZMZM
      INTENT (IN) :: MA,MB
      CALL ZMADD(MA%MZM,MB%MZM,FMADD_ZMZM%MZM)
   END FUNCTION

   FUNCTION FMADD_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMADD_FM
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,FMADD_FM%MFM)
   END FUNCTION

   FUNCTION FMADD_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMADD_IM
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,FMADD_IM%MIM)
   END FUNCTION

   FUNCTION FMADD_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMADD_ZM
      INTENT (IN) :: MA
      CALL ZMEQ(MA%MZM,FMADD_ZM%MZM)
   END FUNCTION

!                                                                   -

   FUNCTION FMSUB_IFM(IVAL,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_IFM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMSUB(MTFM,MA%MFM,FMSUB_IFM%MFM)
   END FUNCTION

   FUNCTION FMSUB_IIM(IVAL,MA)
      USE FMVALS
      TYPE ( IM ) MA,FMSUB_IIM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      CALL IMSUB(MTIM,MA%MIM,FMSUB_IIM%MIM)
   END FUNCTION

   FUNCTION FMSUB_IZM(IVAL,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_IZM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MTZM,MA%MZM,FMSUB_IZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_RFM(R,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_RFM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMSUB(MTFM,MA%MFM,FMSUB_RFM%MFM)
   END FUNCTION

   FUNCTION FMSUB_RIM(R,MA)
      USE FMVALS
      TYPE ( FM ) FMSUB_RIM
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMSUB(MTFM,MUFM,FMSUB_RIM%MFM)
   END FUNCTION

   FUNCTION FMSUB_RZM(R,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_RZM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MTZM,MA%MZM,FMSUB_RZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_DFM(D,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_DFM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMSUB(MTFM,MA%MFM,FMSUB_DFM%MFM)
   END FUNCTION

   FUNCTION FMSUB_DIM(D,MA)
      USE FMVALS
      TYPE ( FM ) FMSUB_DIM
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMSUB(MTFM,MUFM,FMSUB_DIM%MFM)
   END FUNCTION

   FUNCTION FMSUB_DZM(D,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_DZM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MTZM,MA%MZM,FMSUB_DZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZFM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMSUB_ZFM
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMSUB(MTZM,MUZM,FMSUB_ZFM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZIM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMSUB_ZIM
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MTZM,MUZM,FMSUB_ZIM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZZM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZZM
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMSUB(MTZM,MA%MZM,FMSUB_ZZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_CFM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMSUB_CFM
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMSUB(MTZM,MUZM,FMSUB_CFM%MZM)
   END FUNCTION

   FUNCTION FMSUB_CIM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMSUB_CIM
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MTZM,MUZM,FMSUB_CIM%MZM)
   END FUNCTION

   FUNCTION FMSUB_CZM(C,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_CZM
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MTZM,MA%MZM,FMSUB_CZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_FMI(MA,IVAL)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_FMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      CALL FMSUB(MA%MFM,MTFM,FMSUB_FMI%MFM)
   END FUNCTION

   FUNCTION FMSUB_FMR(MA,R)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_FMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMSUB(MA%MFM,MTFM,FMSUB_FMR%MFM)
   END FUNCTION

   FUNCTION FMSUB_FMD(MA,D)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_FMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMSUB(MA%MFM,MTFM,FMSUB_FMD%MFM)
   END FUNCTION

   FUNCTION FMSUB_FMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMSUB_FMZ
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMSUB(MUZM,MTZM,FMSUB_FMZ%MZM)
   END FUNCTION

   FUNCTION FMSUB_FMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMSUB_FMC
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMSUB(MUZM,MTZM,FMSUB_FMC%MZM)
   END FUNCTION

   FUNCTION FMSUB_FMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMSUB_FMFM
      INTENT (IN) :: MA,MB
      CALL FMSUB(MA%MFM,MB%MFM,FMSUB_FMFM%MFM)
   END FUNCTION

   FUNCTION FMSUB_FMIM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_FMIM
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMSUB(MA%MFM,MTFM,FMSUB_FMIM%MFM)
   END FUNCTION

   FUNCTION FMSUB_FMZM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( ZM ) MB,FMSUB_FMZM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MA%MFM,MTFM,MTZM)
      CALL ZMSUB(MTZM,MB%MZM,FMSUB_FMZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_IMI(MA,IVAL)
      USE FMVALS
      TYPE ( IM ) MA,FMSUB_IMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      CALL IMSUB(MA%MIM,MTIM,FMSUB_IMI%MIM)
   END FUNCTION

   FUNCTION FMSUB_IMR(MA,R)
      USE FMVALS
      TYPE ( FM ) FMSUB_IMR
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMSUB(MUFM,MTFM,FMSUB_IMR%MFM)
   END FUNCTION

   FUNCTION FMSUB_IMD(MA,D)
      USE FMVALS
      TYPE ( FM ) FMSUB_IMD
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMSUB(MUFM,MTFM,FMSUB_IMD%MFM)
   END FUNCTION

   FUNCTION FMSUB_IMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMSUB_IMZ
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MUZM,MTZM,FMSUB_IMZ%MZM)
   END FUNCTION

   FUNCTION FMSUB_IMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMSUB_IMC
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MUZM,MTZM,FMSUB_IMC%MZM)
   END FUNCTION

   FUNCTION FMSUB_IMFM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( FM ) MB,FMSUB_IMFM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMSUB(MTFM,MB%MFM,FMSUB_IMFM%MFM)
   END FUNCTION

   FUNCTION FMSUB_IMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMSUB_IMIM
      INTENT (IN) :: MA,MB
      CALL IMSUB(MA%MIM,MB%MIM,FMSUB_IMIM%MIM)
   END FUNCTION

   FUNCTION FMSUB_IMZM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( ZM ) MB,FMSUB_IMZM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MUZM,MB%MZM,FMSUB_IMZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMI(MA,IVAL)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMI%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMR(MA,R)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMR%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMD(MA,D)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMD%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZMZ
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMZ%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZMC
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMC%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MB
      TYPE ( ZM ) MA,FMSUB_ZMFM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MTZM)
      CALL ZMSUB(MA%MZM,MTZM,FMSUB_ZMFM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MB
      TYPE ( ZM ) MA,FMSUB_ZMIM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMSUB(MA%MZM,MUZM,FMSUB_ZMIM%MZM)
   END FUNCTION

   FUNCTION FMSUB_ZMZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA,MB,FMSUB_ZMZM
      INTENT (IN) :: MA,MB
      CALL ZMSUB(MA%MZM,MB%MZM,FMSUB_ZMZM%MZM)
   END FUNCTION

   FUNCTION FMSUB_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSUB_FM
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,MTFM)
      IF (MTFM(1) /= MUNKNO .AND. MTFM(2) /= 0)  &
          MTFM(-1) = -MTFM(-1)
      CALL FMEQ(MTFM,FMSUB_FM%MFM)
   END FUNCTION

   FUNCTION FMSUB_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMSUB_IM
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,MTIM)
      IF (MTIM(1) /= MUNKNO .AND. MTIM(2) /= 0)  &
          MTIM(-1) = -MTIM(-1)
      CALL IMEQ(MTIM,FMSUB_IM%MIM)
   END FUNCTION

   FUNCTION FMSUB_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSUB_ZM
      INTENT (IN) :: MA
      CALL ZMEQ(MA%MZM,MTZM)
      IF (MTZM(1) /= MUNKNO .AND. MTZM(2) /= 0)  &
          MTZM(-1) = -MTZM(-1)
      IF (MTZM(KPTIMU+1) /= MUNKNO .AND. MTZM(KPTIMU+2) /= 0) THEN
          MTZM(KPTIMU-1) = -MTZM(KPTIMU-1)
      ENDIF
      CALL ZMEQ(MTZM,FMSUB_ZM%MZM)
   END FUNCTION

 END MODULE FMZM_5

 MODULE FMZM_6
    USE FMZM_1

    INTERFACE OPERATOR (*)
       MODULE PROCEDURE FMMPY_IFM
       MODULE PROCEDURE FMMPY_IIM
       MODULE PROCEDURE FMMPY_IZM
       MODULE PROCEDURE FMMPY_RFM
       MODULE PROCEDURE FMMPY_RIM
       MODULE PROCEDURE FMMPY_RZM
       MODULE PROCEDURE FMMPY_DFM
       MODULE PROCEDURE FMMPY_DIM
       MODULE PROCEDURE FMMPY_DZM
       MODULE PROCEDURE FMMPY_ZFM
       MODULE PROCEDURE FMMPY_ZIM
       MODULE PROCEDURE FMMPY_ZZM
       MODULE PROCEDURE FMMPY_CFM
       MODULE PROCEDURE FMMPY_CIM
       MODULE PROCEDURE FMMPY_CZM
       MODULE PROCEDURE FMMPY_FMI
       MODULE PROCEDURE FMMPY_FMR
       MODULE PROCEDURE FMMPY_FMD
       MODULE PROCEDURE FMMPY_FMZ
       MODULE PROCEDURE FMMPY_FMC
       MODULE PROCEDURE FMMPY_FMFM
       MODULE PROCEDURE FMMPY_FMIM
       MODULE PROCEDURE FMMPY_FMZM
       MODULE PROCEDURE FMMPY_IMI
       MODULE PROCEDURE FMMPY_IMR
       MODULE PROCEDURE FMMPY_IMD
       MODULE PROCEDURE FMMPY_IMZ
       MODULE PROCEDURE FMMPY_IMC
       MODULE PROCEDURE FMMPY_IMFM
       MODULE PROCEDURE FMMPY_IMIM
       MODULE PROCEDURE FMMPY_IMZM
       MODULE PROCEDURE FMMPY_ZMI
       MODULE PROCEDURE FMMPY_ZMR
       MODULE PROCEDURE FMMPY_ZMD
       MODULE PROCEDURE FMMPY_ZMZ
       MODULE PROCEDURE FMMPY_ZMC
       MODULE PROCEDURE FMMPY_ZMFM
       MODULE PROCEDURE FMMPY_ZMIM
       MODULE PROCEDURE FMMPY_ZMZM
    END INTERFACE

    INTERFACE OPERATOR (/)
       MODULE PROCEDURE FMDIV_IFM
       MODULE PROCEDURE FMDIV_IIM
       MODULE PROCEDURE FMDIV_IZM
       MODULE PROCEDURE FMDIV_RFM
       MODULE PROCEDURE FMDIV_RIM
       MODULE PROCEDURE FMDIV_RZM
       MODULE PROCEDURE FMDIV_DFM
       MODULE PROCEDURE FMDIV_DIM
       MODULE PROCEDURE FMDIV_DZM
       MODULE PROCEDURE FMDIV_ZFM
       MODULE PROCEDURE FMDIV_ZIM
       MODULE PROCEDURE FMDIV_ZZM
       MODULE PROCEDURE FMDIV_CFM
       MODULE PROCEDURE FMDIV_CIM
       MODULE PROCEDURE FMDIV_CZM
       MODULE PROCEDURE FMDIV_FMI
       MODULE PROCEDURE FMDIV_FMR
       MODULE PROCEDURE FMDIV_FMD
       MODULE PROCEDURE FMDIV_FMZ
       MODULE PROCEDURE FMDIV_FMC
       MODULE PROCEDURE FMDIV_FMFM
       MODULE PROCEDURE FMDIV_FMIM
       MODULE PROCEDURE FMDIV_FMZM
       MODULE PROCEDURE FMDIV_IMI
       MODULE PROCEDURE FMDIV_IMR
       MODULE PROCEDURE FMDIV_IMD
       MODULE PROCEDURE FMDIV_IMZ
       MODULE PROCEDURE FMDIV_IMC
       MODULE PROCEDURE FMDIV_IMFM
       MODULE PROCEDURE FMDIV_IMIM
       MODULE PROCEDURE FMDIV_IMZM
       MODULE PROCEDURE FMDIV_ZMI
       MODULE PROCEDURE FMDIV_ZMR
       MODULE PROCEDURE FMDIV_ZMD
       MODULE PROCEDURE FMDIV_ZMZ
       MODULE PROCEDURE FMDIV_ZMC
       MODULE PROCEDURE FMDIV_ZMFM
       MODULE PROCEDURE FMDIV_ZMIM
       MODULE PROCEDURE FMDIV_ZMZM
    END INTERFACE

 CONTAINS

!                                                                   *

   FUNCTION FMMPY_IFM(IVAL,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_IFM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMMPYI(MA%MFM,IVAL,FMMPY_IFM%MFM)
   END FUNCTION

   FUNCTION FMMPY_IIM(IVAL,MA)
      USE FMVALS
      TYPE ( IM ) MA,FMMPY_IIM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMMPYI(MA%MIM,IVAL,FMMPY_IIM%MIM)
   END FUNCTION

   FUNCTION FMMPY_IZM(IVAL,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_IZM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL ZMMPYI(MA%MZM,IVAL,FMMPY_IZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_RFM(R,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_RFM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMMPY(MTFM,MA%MFM,FMMPY_RFM%MFM)
   END FUNCTION

   FUNCTION FMMPY_RIM(R,MA)
      USE FMVALS
      TYPE ( FM ) FMMPY_RIM
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMMPY(MTFM,MUFM,FMMPY_RIM%MFM)
   END FUNCTION

   FUNCTION FMMPY_RZM(R,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_RZM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MTZM,MA%MZM,FMMPY_RZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_DFM(D,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_DFM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMMPY(MTFM,MA%MFM,FMMPY_DFM%MFM)
   END FUNCTION

   FUNCTION FMMPY_DIM(D,MA)
      USE FMVALS
      TYPE ( FM ) FMMPY_DIM
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMMPY(MTFM,MUFM,FMMPY_DIM%MFM)
   END FUNCTION

   FUNCTION FMMPY_DZM(D,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_DZM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MTZM,MA%MZM,FMMPY_DZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZFM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMMPY_ZFM
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMMPY(MTZM,MUZM,FMMPY_ZFM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZIM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMMPY_ZIM
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MTZM,MUZM,FMMPY_ZIM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZZM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZZM
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMMPY(MTZM,MA%MZM,FMMPY_ZZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_CFM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMMPY_CFM
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMMPY(MTZM,MUZM,FMMPY_CFM%MZM)
   END FUNCTION

   FUNCTION FMMPY_CIM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMMPY_CIM
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MTZM,MUZM,FMMPY_CIM%MZM)
   END FUNCTION

   FUNCTION FMMPY_CZM(C,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_CZM
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MTZM,MA%MZM,FMMPY_CZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_FMI(MA,IVAL)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_FMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMMPYI(MA%MFM,IVAL,FMMPY_FMI%MFM)
   END FUNCTION

   FUNCTION FMMPY_FMR(MA,R)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_FMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMMPY(MA%MFM,MTFM,FMMPY_FMR%MFM)
   END FUNCTION

   FUNCTION FMMPY_FMD(MA,D)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_FMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMMPY(MA%MFM,MTFM,FMMPY_FMD%MFM)
   END FUNCTION

   FUNCTION FMMPY_FMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMMPY_FMZ
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMMPY(MUZM,MTZM,FMMPY_FMZ%MZM)
   END FUNCTION

   FUNCTION FMMPY_FMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMMPY_FMC
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMMPY(MUZM,MTZM,FMMPY_FMC%MZM)
   END FUNCTION

   FUNCTION FMMPY_FMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMMPY_FMFM
      INTENT (IN) :: MA,MB
      CALL FMMPY(MA%MFM,MB%MFM,FMMPY_FMFM%MFM)
   END FUNCTION

   FUNCTION FMMPY_FMIM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,FMMPY_FMIM
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMMPY(MA%MFM,MTFM,FMMPY_FMIM%MFM)
   END FUNCTION

   FUNCTION FMMPY_FMZM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( ZM ) MB,FMMPY_FMZM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MA%MFM,MTFM,MTZM)
      CALL ZMMPY(MTZM,MB%MZM,FMMPY_FMZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_IMI(MA,IVAL)
      USE FMVALS
      TYPE ( IM ) MA,FMMPY_IMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMMPYI(MA%MIM,IVAL,FMMPY_IMI%MIM)
   END FUNCTION

   FUNCTION FMMPY_IMR(MA,R)
      USE FMVALS
      TYPE ( FM ) FMMPY_IMR
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMMPY(MUFM,MTFM,FMMPY_IMR%MFM)
   END FUNCTION

   FUNCTION FMMPY_IMD(MA,D)
      USE FMVALS
      TYPE ( FM ) FMMPY_IMD
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMMPY(MUFM,MTFM,FMMPY_IMD%MFM)
   END FUNCTION

   FUNCTION FMMPY_IMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMMPY_IMZ
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MUZM,MTZM,FMMPY_IMZ%MZM)
   END FUNCTION

   FUNCTION FMMPY_IMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMMPY_IMC
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MUZM,MTZM,FMMPY_IMC%MZM)
   END FUNCTION

   FUNCTION FMMPY_IMFM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( FM ) MB,FMMPY_IMFM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMMPY(MTFM,MB%MFM,FMMPY_IMFM%MFM)
   END FUNCTION

   FUNCTION FMMPY_IMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMMPY_IMIM
      INTENT (IN) :: MA,MB
      CALL IMMPY(MA%MIM,MB%MIM,FMMPY_IMIM%MIM)
   END FUNCTION

   FUNCTION FMMPY_IMZM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( ZM ) MB,FMMPY_IMZM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MUZM,MB%MZM,FMMPY_IMZM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMI(MA,IVAL)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL ZMMPYI(MA%MZM,IVAL,FMMPY_ZMI%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMR(MA,R)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MA%MZM,MTZM,FMMPY_ZMR%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMD(MA,D)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MA%MZM,MTZM,FMMPY_ZMD%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZMZ
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMMPY(MA%MZM,MTZM,FMMPY_ZMZ%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) MA,FMMPY_ZMC
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMMPY(MA%MZM,MTZM,FMMPY_ZMC%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MB
      TYPE ( ZM ) MA,FMMPY_ZMFM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MTZM)
      CALL ZMMPY(MA%MZM,MTZM,FMMPY_ZMFM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MB
      TYPE ( ZM ) MA,FMMPY_ZMIM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMMPY(MA%MZM,MUZM,FMMPY_ZMIM%MZM)
   END FUNCTION

   FUNCTION FMMPY_ZMZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA,MB,FMMPY_ZMZM
      INTENT (IN) :: MA,MB
      CALL ZMMPY(MA%MZM,MB%MZM,FMMPY_ZMZM%MZM)
   END FUNCTION

!                                                                   /

   FUNCTION FMDIV_IFM(IVAL,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_IFM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMDIV(MTFM,MA%MFM,FMDIV_IFM%MFM)
   END FUNCTION

   FUNCTION FMDIV_IIM(IVAL,MA)
      USE FMVALS
      TYPE ( IM ) MA,FMDIV_IIM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      CALL IMDIV(MTIM,MA%MIM,FMDIV_IIM%MIM)
   END FUNCTION

   FUNCTION FMDIV_IZM(IVAL,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_IZM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MTZM,MA%MZM,FMDIV_IZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_RFM(R,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_RFM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMDIV(MTFM,MA%MFM,FMDIV_RFM%MFM)
   END FUNCTION

   FUNCTION FMDIV_RIM(R,MA)
      USE FMVALS
      TYPE ( FM ) FMDIV_RIM
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMDIV(MTFM,MUFM,FMDIV_RIM%MFM)
   END FUNCTION

   FUNCTION FMDIV_RZM(R,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_RZM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MTZM,MA%MZM,FMDIV_RZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_DFM(D,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_DFM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMDIV(MTFM,MA%MFM,FMDIV_DFM%MFM)
   END FUNCTION

   FUNCTION FMDIV_DIM(D,MA)
      USE FMVALS
      TYPE ( FM ) FMDIV_DIM
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMDIV(MTFM,MUFM,FMDIV_DIM%MFM)
   END FUNCTION

   FUNCTION FMDIV_DZM(D,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_DZM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MTZM,MA%MZM,FMDIV_DZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZFM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMDIV_ZFM
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMDIV(MTZM,MUZM,FMDIV_ZFM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZIM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMDIV_ZIM
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MTZM,MUZM,FMDIV_ZIM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZZM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZZM
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMDIV(MTZM,MA%MZM,FMDIV_ZZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_CFM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMDIV_CFM
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMDIV(MTZM,MUZM,FMDIV_CFM%MZM)
   END FUNCTION

   FUNCTION FMDIV_CIM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMDIV_CIM
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MTZM,MUZM,FMDIV_CIM%MZM)
   END FUNCTION

   FUNCTION FMDIV_CZM(C,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_CZM
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MTZM,MA%MZM,FMDIV_CZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_FMI(MA,IVAL)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_FMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMDIVI(MA%MFM,IVAL,FMDIV_FMI%MFM)
   END FUNCTION

   FUNCTION FMDIV_FMR(MA,R)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_FMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMDIV(MA%MFM,MTFM,FMDIV_FMR%MFM)
   END FUNCTION

   FUNCTION FMDIV_FMD(MA,D)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_FMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMDIV(MA%MFM,MTFM,FMDIV_FMD%MFM)
   END FUNCTION

   FUNCTION FMDIV_FMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMDIV_FMZ
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMDIV(MUZM,MTZM,FMDIV_FMZ%MZM)
   END FUNCTION

   FUNCTION FMDIV_FMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMDIV_FMC
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMDIV(MUZM,MTZM,FMDIV_FMC%MZM)
   END FUNCTION

   FUNCTION FMDIV_FMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMDIV_FMFM
      INTENT (IN) :: MA,MB
      CALL FMDIV(MA%MFM,MB%MFM,FMDIV_FMFM%MFM)
   END FUNCTION

   FUNCTION FMDIV_FMIM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,FMDIV_FMIM
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMDIV(MA%MFM,MTFM,FMDIV_FMIM%MFM)
   END FUNCTION

   FUNCTION FMDIV_FMZM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( ZM ) MB,FMDIV_FMZM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MA%MFM,MTFM,MTZM)
      CALL ZMDIV(MTZM,MB%MZM,FMDIV_FMZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_IMI(MA,IVAL)
      USE FMVALS
      TYPE ( IM ) MA,FMDIV_IMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMDIVI(MA%MIM,IVAL,FMDIV_IMI%MIM)
   END FUNCTION

   FUNCTION FMDIV_IMR(MA,R)
      USE FMVALS
      TYPE ( FM ) FMDIV_IMR
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMDIV(MUFM,MTFM,FMDIV_IMR%MFM)
   END FUNCTION

   FUNCTION FMDIV_IMD(MA,D)
      USE FMVALS
      TYPE ( FM ) FMDIV_IMD
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMDIV(MUFM,MTFM,FMDIV_IMD%MFM)
   END FUNCTION

   FUNCTION FMDIV_IMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMDIV_IMZ
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MUZM,MTZM,FMDIV_IMZ%MZM)
   END FUNCTION

   FUNCTION FMDIV_IMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMDIV_IMC
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MUZM,MTZM,FMDIV_IMC%MZM)
   END FUNCTION

   FUNCTION FMDIV_IMFM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( FM ) MB,FMDIV_IMFM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMDIV(MTFM,MB%MFM,FMDIV_IMFM%MFM)
   END FUNCTION

   FUNCTION FMDIV_IMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMDIV_IMIM
      INTENT (IN) :: MA,MB
      CALL IMDIV(MA%MIM,MB%MIM,FMDIV_IMIM%MIM)
   END FUNCTION

   FUNCTION FMDIV_IMZM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( ZM ) MB,FMDIV_IMZM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MUZM,MB%MZM,FMDIV_IMZM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMI(MA,IVAL)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL ZMDIVI(MA%MZM,IVAL,FMDIV_ZMI%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMR(MA,R)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MA%MZM,MTZM,FMDIV_ZMR%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMD(MA,D)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MA%MZM,MTZM,FMDIV_ZMD%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZMZ
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMDIV(MA%MZM,MTZM,FMDIV_ZMZ%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) MA,FMDIV_ZMC
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMDIV(MA%MZM,MTZM,FMDIV_ZMC%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MB
      TYPE ( ZM ) MA,FMDIV_ZMFM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MTZM)
      CALL ZMDIV(MA%MZM,MTZM,FMDIV_ZMFM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MB
      TYPE ( ZM ) MA,FMDIV_ZMIM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMDIV(MA%MZM,MUZM,FMDIV_ZMIM%MZM)
   END FUNCTION

   FUNCTION FMDIV_ZMZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA,MB,FMDIV_ZMZM
      INTENT (IN) :: MA,MB
      CALL ZMDIV(MA%MZM,MB%MZM,FMDIV_ZMZM%MZM)
   END FUNCTION

 END MODULE FMZM_6

 MODULE FMZM_7
    USE FMZM_1

    INTERFACE OPERATOR (**)
       MODULE PROCEDURE FMPWR_IFM
       MODULE PROCEDURE FMPWR_IIM
       MODULE PROCEDURE FMPWR_IZM
       MODULE PROCEDURE FMPWR_RFM
       MODULE PROCEDURE FMPWR_RIM
       MODULE PROCEDURE FMPWR_RZM
       MODULE PROCEDURE FMPWR_DFM
       MODULE PROCEDURE FMPWR_DIM
       MODULE PROCEDURE FMPWR_DZM
       MODULE PROCEDURE FMPWR_ZFM
       MODULE PROCEDURE FMPWR_ZIM
       MODULE PROCEDURE FMPWR_ZZM
       MODULE PROCEDURE FMPWR_CFM
       MODULE PROCEDURE FMPWR_CIM
       MODULE PROCEDURE FMPWR_CZM
       MODULE PROCEDURE FMPWR_FMI
       MODULE PROCEDURE FMPWR_FMR
       MODULE PROCEDURE FMPWR_FMD
       MODULE PROCEDURE FMPWR_FMZ
       MODULE PROCEDURE FMPWR_FMC
       MODULE PROCEDURE FMPWR_FMFM
       MODULE PROCEDURE FMPWR_FMIM
       MODULE PROCEDURE FMPWR_FMZM
       MODULE PROCEDURE FMPWR_IMI
       MODULE PROCEDURE FMPWR_IMR
       MODULE PROCEDURE FMPWR_IMD
       MODULE PROCEDURE FMPWR_IMZ
       MODULE PROCEDURE FMPWR_IMC
       MODULE PROCEDURE FMPWR_IMFM
       MODULE PROCEDURE FMPWR_IMIM
       MODULE PROCEDURE FMPWR_IMZM
       MODULE PROCEDURE FMPWR_ZMI
       MODULE PROCEDURE FMPWR_ZMR
       MODULE PROCEDURE FMPWR_ZMD
       MODULE PROCEDURE FMPWR_ZMZ
       MODULE PROCEDURE FMPWR_ZMC
       MODULE PROCEDURE FMPWR_ZMFM
       MODULE PROCEDURE FMPWR_ZMIM
       MODULE PROCEDURE FMPWR_ZMZM
    END INTERFACE

   INTERFACE ABS
      MODULE PROCEDURE FMABS_FM
      MODULE PROCEDURE FMABS_IM
      MODULE PROCEDURE FMABS_ZM
   END INTERFACE

   INTERFACE ACOS
      MODULE PROCEDURE FMACOS_FM
      MODULE PROCEDURE FMACOS_ZM
   END INTERFACE

   INTERFACE AIMAG
      MODULE PROCEDURE FMAIMAG_ZM
   END INTERFACE

   INTERFACE AINT
      MODULE PROCEDURE FMAINT_FM
      MODULE PROCEDURE FMAINT_ZM
   END INTERFACE

   INTERFACE ANINT
      MODULE PROCEDURE FMANINT_FM
      MODULE PROCEDURE FMANINT_ZM
   END INTERFACE

   INTERFACE ASIN
      MODULE PROCEDURE FMASIN_FM
      MODULE PROCEDURE FMASIN_ZM
   END INTERFACE

   INTERFACE ATAN
      MODULE PROCEDURE FMATAN_FM
      MODULE PROCEDURE FMATAN_ZM
   END INTERFACE

   INTERFACE ATAN2
      MODULE PROCEDURE FMATAN2_FM
   END INTERFACE

   INTERFACE BTEST
      MODULE PROCEDURE FMBTEST_IM
   END INTERFACE

   INTERFACE CEILING
      MODULE PROCEDURE FMCEILING_FM
      MODULE PROCEDURE FMCEILING_ZM
   END INTERFACE

   INTERFACE CMPLX
      MODULE PROCEDURE FMCMPLX_FM
      MODULE PROCEDURE FMCMPLX_IM
   END INTERFACE

   INTERFACE CONJG
      MODULE PROCEDURE FMCONJG_ZM
   END INTERFACE

   INTERFACE COS
      MODULE PROCEDURE FMCOS_FM
      MODULE PROCEDURE FMCOS_ZM
   END INTERFACE

   INTERFACE COSH
      MODULE PROCEDURE FMCOSH_FM
      MODULE PROCEDURE FMCOSH_ZM
   END INTERFACE

   INTERFACE DBLE
      MODULE PROCEDURE FMDBLE_FM
      MODULE PROCEDURE FMDBLE_IM
      MODULE PROCEDURE FMDBLE_ZM
   END INTERFACE

   INTERFACE DIGITS
      MODULE PROCEDURE FMDIGITS_FM
      MODULE PROCEDURE FMDIGITS_IM
      MODULE PROCEDURE FMDIGITS_ZM
   END INTERFACE

   INTERFACE DIM
      MODULE PROCEDURE FMDIM_FM
      MODULE PROCEDURE FMDIM_IM
   END INTERFACE

   INTERFACE DINT
      MODULE PROCEDURE FMDINT_FM
      MODULE PROCEDURE FMDINT_ZM
   END INTERFACE

   INTERFACE DOTPRODUCT
      MODULE PROCEDURE FMDOTPRODUCT_FM
      MODULE PROCEDURE FMDOTPRODUCT_IM
      MODULE PROCEDURE FMDOTPRODUCT_ZM
   END INTERFACE

 CONTAINS

!                                                                  **

   FUNCTION FMPWR_IFM(IVAL,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_IFM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMPWR(MTFM,MA%MFM,FMPWR_IFM%MFM)
   END FUNCTION

   FUNCTION FMPWR_IIM(IVAL,MA)
      USE FMVALS
      TYPE ( IM ) MA,FMPWR_IIM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL IMI2M(IVAL,MTIM)
      CALL IMPWR(MTIM,MA%MIM,FMPWR_IIM%MIM)
   END FUNCTION

   FUNCTION FMPWR_IZM(IVAL,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_IZM
      INTEGER IVAL
      INTENT (IN) :: IVAL,MA
      CALL FMI2M(IVAL,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MTZM,MA%MZM,FMPWR_IZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_RFM(R,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_RFM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMPWR(MTFM,MA%MFM,FMPWR_RFM%MFM)
   END FUNCTION

   FUNCTION FMPWR_RIM(R,MA)
      USE FMVALS
      TYPE ( FM ) FMPWR_RIM
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMPWR(MTFM,MUFM,FMPWR_RIM%MFM)
   END FUNCTION

   FUNCTION FMPWR_RZM(R,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_RZM
      REAL R
      INTENT (IN) :: R,MA
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MTZM,MA%MZM,FMPWR_RZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_DFM(D,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_DFM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMPWR(MTFM,MA%MFM,FMPWR_DFM%MFM)
   END FUNCTION

   FUNCTION FMPWR_DIM(D,MA)
      USE FMVALS
      TYPE ( FM ) FMPWR_DIM
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMPWR(MTFM,MUFM,FMPWR_DIM%MFM)
   END FUNCTION

   FUNCTION FMPWR_DZM(D,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_DZM
      DOUBLE PRECISION D
      INTENT (IN) :: D,MA
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MTZM,MA%MZM,FMPWR_DZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZFM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMPWR_ZFM
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMPWR(MTZM,MUZM,FMPWR_ZFM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZIM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) FMPWR_ZIM
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MTZM,MUZM,FMPWR_ZIM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZZM(Z,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZZM
      COMPLEX Z
      INTENT (IN) :: Z,MA
      CALL ZMZ2M(Z,MTZM)
      CALL ZMPWR(MTZM,MA%MZM,FMPWR_ZZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_CFM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMPWR_CFM
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMPWR(MTZM,MUZM,FMPWR_CFM%MZM)
   END FUNCTION

   FUNCTION FMPWR_CIM(C,MA)
      USE FMVALS
      TYPE ( ZM ) FMPWR_CIM
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MTZM,MUZM,FMPWR_CIM%MZM)
   END FUNCTION

   FUNCTION FMPWR_CZM(C,MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_CZM
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C,MA
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MTZM,MA%MZM,FMPWR_CZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_FMI(MA,IVAL)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_FMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL FMIPWR(MA%MFM,IVAL,FMPWR_FMI%MFM)
   END FUNCTION

   FUNCTION FMPWR_FMR(MA,R)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_FMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMPWR(MA%MFM,MTFM,FMPWR_FMR%MFM)
   END FUNCTION

   FUNCTION FMPWR_FMD(MA,D)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_FMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMPWR(MA%MFM,MTFM,FMPWR_FMD%MFM)
   END FUNCTION

   FUNCTION FMPWR_FMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMPWR_FMZ
      TYPE ( FM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMPWR(MUZM,MTZM,FMPWR_FMZ%MZM)
   END FUNCTION

   FUNCTION FMPWR_FMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMPWR_FMC
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,MUZM)
      CALL ZMPWR(MUZM,MTZM,FMPWR_FMC%MZM)
   END FUNCTION

   FUNCTION FMPWR_FMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMPWR_FMFM
      INTENT (IN) :: MA,MB
      CALL FMPWR(MA%MFM,MB%MFM,FMPWR_FMFM%MFM)
   END FUNCTION

   FUNCTION FMPWR_FMIM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,FMPWR_FMIM
      TYPE ( IM ) MB
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMPWR(MA%MFM,MTFM,FMPWR_FMIM%MFM)
   END FUNCTION

   FUNCTION FMPWR_FMZM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( ZM ) MB,FMPWR_FMZM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MA%MFM,MTFM,MTZM)
      CALL ZMPWR(MTZM,MB%MZM,FMPWR_FMZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_IMI(MA,IVAL)
      USE FMVALS
      TYPE ( IM ) MA,FMPWR_IMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL IMI2M(IVAL,MTIM)
      CALL IMPWR(MA%MIM,MTIM,FMPWR_IMI%MIM)
   END FUNCTION

   FUNCTION FMPWR_IMR(MA,R)
      USE FMVALS
      TYPE ( FM ) FMPWR_IMR
      TYPE ( IM ) MA
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMPWR(MUFM,MTFM,FMPWR_IMR%MFM)
   END FUNCTION

   FUNCTION FMPWR_IMD(MA,D)
      USE FMVALS
      TYPE ( FM ) FMPWR_IMD
      TYPE ( IM ) MA
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL IMI2FM(MA%MIM,MUFM)
      CALL FMPWR(MUFM,MTFM,FMPWR_IMD%MFM)
   END FUNCTION

   FUNCTION FMPWR_IMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) FMPWR_IMZ
      TYPE ( IM ) MA
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MUZM,MTZM,FMPWR_IMZ%MZM)
   END FUNCTION

   FUNCTION FMPWR_IMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) FMPWR_IMC
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MUZM,MTZM,FMPWR_IMC%MZM)
   END FUNCTION

   FUNCTION FMPWR_IMFM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( FM ) MB,FMPWR_IMFM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMPWR(MTFM,MB%MFM,FMPWR_IMFM%MFM)
   END FUNCTION

   FUNCTION FMPWR_IMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMPWR_IMIM
      INTENT (IN) :: MA,MB
      CALL IMPWR(MA%MIM,MB%MIM,FMPWR_IMIM%MIM)
   END FUNCTION

   FUNCTION FMPWR_IMZM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA
      TYPE ( ZM ) MB,FMPWR_IMZM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MUZM,MB%MZM,FMPWR_IMZM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMI(MA,IVAL)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZMI
      INTEGER IVAL
      INTENT (IN) :: MA,IVAL
      CALL ZMIPWR(MA%MZM,IVAL,FMPWR_ZMI%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMR(MA,R)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZMR
      REAL R
      INTENT (IN) :: MA,R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MA%MZM,MTZM,FMPWR_ZMR%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMD(MA,D)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZMD
      DOUBLE PRECISION D
      INTENT (IN) :: MA,D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MA%MZM,MTZM,FMPWR_ZMD%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMZ(MA,Z)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZMZ
      COMPLEX Z
      INTENT (IN) :: MA,Z
      CALL ZMZ2M(Z,MTZM)
      CALL ZMPWR(MA%MZM,MTZM,FMPWR_ZMZ%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMC(MA,C)
      USE FMVALS
      TYPE ( ZM ) MA,FMPWR_ZMC
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: MA,C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,MTZM)
      CALL ZMPWR(MA%MZM,MTZM,FMPWR_ZMC%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMFM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MB
      TYPE ( ZM ) MA,FMPWR_ZMFM
      INTENT (IN) :: MA,MB
      CALL FMI2M(0,MTFM)
      CALL ZMCMPX(MB%MFM,MTFM,MTZM)
      CALL ZMPWR(MA%MZM,MTZM,FMPWR_ZMFM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMIM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MB
      TYPE ( ZM ) MA,FMPWR_ZMIM
      INTENT (IN) :: MA,MB
      CALL IMI2FM(MB%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,MUZM)
      CALL ZMPWR(MA%MZM,MUZM,FMPWR_ZMIM%MZM)
   END FUNCTION

   FUNCTION FMPWR_ZMZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA,MB,FMPWR_ZMZM
      INTENT (IN) :: MA,MB
      CALL ZMPWR(MA%MZM,MB%MZM,FMPWR_ZMZM%MZM)
   END FUNCTION

!                                                                 ABS

   FUNCTION FMABS_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMABS_FM
      INTENT (IN) :: MA
      CALL FMABS(MA%MFM,FMABS_FM%MFM)
   END FUNCTION

   FUNCTION FMABS_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMABS_IM
      INTENT (IN) :: MA
      CALL IMABS(MA%MIM,FMABS_IM%MIM)
   END FUNCTION

   FUNCTION FMABS_ZM(MA)
      USE FMVALS
      TYPE ( FM ) FMABS_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMABS(MA%MZM,FMABS_ZM%MFM)
   END FUNCTION

!                                                                ACOS

   FUNCTION FMACOS_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMACOS_FM
      INTENT (IN) :: MA
      CALL FMACOS(MA%MFM,FMACOS_FM%MFM)
   END FUNCTION

   FUNCTION FMACOS_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMACOS_ZM
      INTENT (IN) :: MA
      CALL ZMACOS(MA%MZM,FMACOS_ZM%MZM)
   END FUNCTION

!                                                               AIMAG

   FUNCTION FMAIMAG_ZM(MA)
      USE FMVALS
      TYPE ( FM ) FMAIMAG_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMIMAG(MA%MZM,FMAIMAG_ZM%MFM)
   END FUNCTION

!                                                                AINT

   FUNCTION FMAINT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMAINT_FM
      INTENT (IN) :: MA
      CALL FMINT(MA%MFM,FMAINT_FM%MFM)
   END FUNCTION

   FUNCTION FMAINT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMAINT_ZM
      INTENT (IN) :: MA
      CALL ZMINT(MA%MZM,FMAINT_ZM%MZM)
   END FUNCTION

!                                                               ANINT

   FUNCTION FMANINT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMANINT_FM
      INTENT (IN) :: MA
      CALL FMNINT(MA%MFM,FMANINT_FM%MFM)
   END FUNCTION

   FUNCTION FMANINT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMANINT_ZM
      INTENT (IN) :: MA
      CALL ZMNINT(MA%MZM,FMANINT_ZM%MZM)
   END FUNCTION

!                                                                ASIN

   FUNCTION FMASIN_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMASIN_FM
      INTENT (IN) :: MA
      CALL FMASIN(MA%MFM,FMASIN_FM%MFM)
   END FUNCTION

   FUNCTION FMASIN_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMASIN_ZM
      INTENT (IN) :: MA
      CALL ZMASIN(MA%MZM,FMASIN_ZM%MZM)
   END FUNCTION

!                                                                ATAN

   FUNCTION FMATAN_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMATAN_FM
      INTENT (IN) :: MA
      CALL FMATAN(MA%MFM,FMATAN_FM%MFM)
   END FUNCTION

   FUNCTION FMATAN_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMATAN_ZM
      INTENT (IN) :: MA
      CALL ZMATAN(MA%MZM,FMATAN_ZM%MZM)
   END FUNCTION

!                                                               ATAN2

   FUNCTION FMATAN2_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMATAN2_FM
      INTENT (IN) :: MA,MB
      CALL FMATN2(MA%MFM,MB%MFM,FMATAN2_FM%MFM)
   END FUNCTION

!                                                               BTEST

   FUNCTION FMBTEST_IM(MA,POS)
      TYPE ( IM ) MA
      INTEGER POS
      LOGICAL FMBTEST_IM
      INTENT (IN) :: MA,POS
      CALL IMI2M(2,MTIM)
      CALL IMI2M(POS,MUIM)
      CALL IMPWR(MTIM,MUIM,MVIM)
      CALL IMDIV(MA%MIM,MVIM,MUIM)
      MUIM(-1) = 1
      CALL IMMOD(MUIM,MTIM,MVIM)
      IF (MVIM(2) == 0) THEN
          FMBTEST_IM = .FALSE.
      ELSE
          FMBTEST_IM = .TRUE.
      ENDIF
   END FUNCTION

!                                                             CEILING

   FUNCTION FMCEILING_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMCEILING_FM
      INTENT (IN) :: MA
      CALL FMINT(MA%MFM,MTFM)
      CALL FMSUB(MA%MFM,MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MFM,FMCEILING_FM%MFM)
      ELSE IF (MA%MFM(-1) > 0) THEN
          CALL FMADDI(MTFM,1)
          CALL FMEQ(MTFM,FMCEILING_FM%MFM)
      ELSE
          CALL FMEQ(MTFM,FMCEILING_FM%MFM)
      ENDIF
   END FUNCTION

   FUNCTION FMCEILING_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMCEILING_ZM
      INTENT (IN) :: MA
      CALL FMINT(MA%MZM,MTFM)
      CALL FMSUB(MA%MZM,MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MZM,MVFM)
      ELSE IF (MA%MZM(-1) > 0) THEN
          CALL FMADDI(MTFM,1)
          CALL FMEQ(MTFM,MVFM)
      ELSE
          CALL FMEQ(MTFM,MVFM)
      ENDIF
      CALL FMINT(MA%MZM(KPTIMU-1),MTFM)
      CALL FMSUB(MA%MZM(KPTIMU-1),MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MZM(KPTIMU-1),MUFM)
      ELSE IF (MA%MZM(KPTIMU-1) > 0) THEN
          CALL FMADDI(MTFM,1)
          CALL FMEQ(MTFM,MUFM)
      ELSE
          CALL FMEQ(MTFM,MUFM)
      ENDIF
      CALL ZMCMPX(MVFM,MUFM,FMCEILING_ZM%MZM)
   END FUNCTION

!                                                               CMPLX

   FUNCTION FMCMPLX_FM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) FMCMPLX_FM
      TYPE ( FM ) MA
      TYPE ( FM ), OPTIONAL :: MB
      INTENT (IN) :: MA,MB
      IF (PRESENT(MB)) THEN
          CALL ZMCMPX(MA%MFM,MB%MFM,FMCMPLX_FM%MZM)
      ELSE
          CALL FMI2M(0,MTFM)
          CALL ZMCMPX(MA%MFM,MTFM,FMCMPLX_FM%MZM)
      ENDIF
   END FUNCTION

   FUNCTION FMCMPLX_IM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) FMCMPLX_IM
      TYPE ( IM ) MA
      TYPE ( IM ), OPTIONAL :: MB
      INTENT (IN) :: MA,MB
      IF (PRESENT(MB)) THEN
          CALL IMI2FM(MA%MIM,MTFM)
          CALL IMI2FM(MB%MIM,MUFM)
          CALL ZMCMPX(MTFM,MUFM,FMCMPLX_IM%MZM)
      ELSE
          CALL IMI2FM(MA%MIM,MTFM)
          CALL FMI2M(0,MUFM)
          CALL ZMCMPX(MTFM,MUFM,FMCMPLX_IM%MZM)
      ENDIF
   END FUNCTION

!                                                               CONJG

   FUNCTION FMCONJG_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) FMCONJG_ZM,MA
      INTENT (IN) :: MA
      CALL ZMCONJ(MA%MZM,FMCONJG_ZM%MZM)
   END FUNCTION

!                                                                 COS

   FUNCTION FMCOS_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMCOS_FM
      INTENT (IN) :: MA
      CALL FMCOS(MA%MFM,FMCOS_FM%MFM)
   END FUNCTION

   FUNCTION FMCOS_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMCOS_ZM
      INTENT (IN) :: MA
      CALL ZMCOS(MA%MZM,FMCOS_ZM%MZM)
   END FUNCTION

!                                                                COSH

   FUNCTION FMCOSH_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMCOSH_FM
      INTENT (IN) :: MA
      CALL FMCOSH(MA%MFM,FMCOSH_FM%MFM)
   END FUNCTION

   FUNCTION FMCOSH_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMCOSH_ZM
      INTENT (IN) :: MA
      CALL ZMCOSH(MA%MZM,FMCOSH_ZM%MZM)
   END FUNCTION

!                                                                DBLE

   FUNCTION FMDBLE_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMDBLE_FM
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,FMDBLE_FM%MFM)
   END FUNCTION

   FUNCTION FMDBLE_IM(MA)
      USE FMVALS
      TYPE ( FM ) FMDBLE_IM
      TYPE ( IM ) MA
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,FMDBLE_IM%MFM)
   END FUNCTION

   FUNCTION FMDBLE_ZM(MA)
      USE FMVALS
      TYPE ( FM ) FMDBLE_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,FMDBLE_ZM%MFM)
   END FUNCTION

!                                                              DIGITS

   FUNCTION FMDIGITS_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMDIGITS_FM
      INTENT (IN) :: MA
      FMDIGITS_FM = NDIG
   END FUNCTION

   FUNCTION FMDIGITS_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA
      INTEGER FMDIGITS_IM
      INTENT (IN) :: MA
      FMDIGITS_IM = NDIGMX
   END FUNCTION

   FUNCTION FMDIGITS_ZM(MA)
      USE FMVALS
      INTEGER FMDIGITS_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      FMDIGITS_ZM = NDIG
   END FUNCTION

!                                                                 DIM

   FUNCTION FMDIM_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMDIM_FM
      INTENT (IN) :: MA,MB
      CALL FMDIM(MA%MFM,MB%MFM,FMDIM_FM%MFM)
   END FUNCTION

   FUNCTION FMDIM_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMDIM_IM
      INTENT (IN) :: MA,MB
      CALL IMDIM(MA%MIM,MB%MIM,FMDIM_IM%MIM)
   END FUNCTION

!                                                                DINT

   FUNCTION FMDINT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMDINT_FM
      INTENT (IN) :: MA
      CALL FMINT(MA%MFM,FMDINT_FM%MFM)
   END FUNCTION

   FUNCTION FMDINT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMDINT_ZM
      INTENT (IN) :: MA
      CALL ZMINT(MA%MZM,FMDINT_ZM%MZM)
   END FUNCTION

!                                                          DOTPRODUCT

   FUNCTION FMDOTPRODUCT_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA(:),MB(:),FMDOTPRODUCT_FM
      INTEGER J,JA,JB,NDSAVE
      INTENT (IN) :: MA,MB
      IF (SIZE(MA) == SIZE(MB)) THEN
          NDSAVE = NDIG
          J = MAX(NGRD52,2)
          NDIG = MIN(MAX(NDIG+J,2),NDG2MX)
          CALL FMI2M(0,FMDOTPRODUCT_FM%MFM)
          DO J = 1, SIZE(MA)
             JA = LBOUND(MA,DIM=1) + J - 1
             JB = LBOUND(MB,DIM=1) + J - 1
             CALL FMEQ2(MA(JA)%MFM,MUFM,NDSAVE,NDIG)
             CALL FMEQ2(MB(JB)%MFM,MVFM,NDSAVE,NDIG)
             CALL FMMPY(MUFM,MVFM,MTFM)
             CALL FMADD_R1(FMDOTPRODUCT_FM%MFM,MTFM)
          ENDDO
          CALL FMEQ2(FMDOTPRODUCT_FM%MFM,MTFM,NDIG,NDSAVE)
          CALL FMEQ(MTFM,FMDOTPRODUCT_FM%MFM)
          NDIG = NDSAVE
      ELSE
          CALL FMI2M(1,MTFM)
          CALL FMI2M(0,MUFM)
          CALL FMDIV(MTFM,MUFM,FMDOTPRODUCT_FM%MFM)
      ENDIF
   END FUNCTION

   FUNCTION FMDOTPRODUCT_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA(:),MB(:),FMDOTPRODUCT_IM
      INTEGER J,JA,JB
      INTENT (IN) :: MA,MB
      IF (SIZE(MA) == SIZE(MB)) THEN
      CALL IMI2M(0,FMDOTPRODUCT_IM%MIM)
          DO J = 1, SIZE(MA)
             JA = LBOUND(MA,DIM=1) + J - 1
             JB = LBOUND(MB,DIM=1) + J - 1
             CALL IMMPY(MA(JA)%MIM,MB(JB)%MIM,MTIM)
             CALL IMADD(FMDOTPRODUCT_IM%MIM,MTIM,MUIM)
             CALL IMEQ(MUIM,FMDOTPRODUCT_IM%MIM)
          ENDDO
      ELSE
          CALL IMI2M(1,MTIM)
          CALL IMI2M(0,MUIM)
          CALL IMDIV(MTIM,MUIM,FMDOTPRODUCT_IM%MIM)
      ENDIF
   END FUNCTION

   FUNCTION FMDOTPRODUCT_ZM(MA,MB)
      USE FMVALS
      TYPE ( ZM ) MA(:),MB(:),FMDOTPRODUCT_ZM
      INTEGER J,JA,JB,NDSAVE
      INTENT (IN) :: MA,MB
      IF (SIZE(MA) == SIZE(MB)) THEN
          NDSAVE = NDIG
          J = MAX(NGRD52,2)
          NDIG = MIN(MAX(NDIG+J,2),NDG2MX)
          CALL ZMI2M(0,FMDOTPRODUCT_ZM%MZM)
          DO J = 1, SIZE(MA)
             JA = LBOUND(MA,DIM=1) + J - 1
             JB = LBOUND(MB,DIM=1) + J - 1
             CALL ZMEQ2(MA(JA)%MZM,MUZM,NDSAVE,NDIG)
             CALL ZMEQ2(MB(JB)%MZM,MVZM,NDSAVE,NDIG)
             CALL ZMMPY(MUZM,MVZM,MTZM)
             CALL ZMADD(FMDOTPRODUCT_ZM%MZM,MTZM,MUZM)
             CALL ZMEQ(MUZM,FMDOTPRODUCT_ZM%MZM)
          ENDDO
          CALL ZMEQ2(FMDOTPRODUCT_ZM%MZM,MTZM,NDIG,NDSAVE)
          CALL ZMEQ(MTZM,FMDOTPRODUCT_ZM%MZM)
          NDIG = NDSAVE
      ELSE
          CALL ZMI2M(1,MTZM)
          CALL ZMI2M(0,MUZM)
          CALL ZMDIV(MTZM,MUZM,FMDOTPRODUCT_ZM%MZM)
      ENDIF
   END FUNCTION

 END MODULE FMZM_7

 MODULE FMZM_8
    USE FMZM_1

   INTERFACE EPSILON
      MODULE PROCEDURE FMEPSILON_FM
   END INTERFACE

   INTERFACE EXP
      MODULE PROCEDURE FMEXP_FM
      MODULE PROCEDURE FMEXP_ZM
   END INTERFACE

   INTERFACE EXPONENT
      MODULE PROCEDURE FMEXPONENT_FM
   END INTERFACE

   INTERFACE FLOOR
      MODULE PROCEDURE FMFLOOR_FM
      MODULE PROCEDURE FMFLOOR_IM
      MODULE PROCEDURE FMFLOOR_ZM
   END INTERFACE

   INTERFACE FRACTION
      MODULE PROCEDURE FMFRACTION_FM
      MODULE PROCEDURE FMFRACTION_ZM
   END INTERFACE

   INTERFACE HUGE
      MODULE PROCEDURE FMHUGE_FM
      MODULE PROCEDURE FMHUGE_IM
      MODULE PROCEDURE FMHUGE_ZM
   END INTERFACE

   INTERFACE INT
      MODULE PROCEDURE FMINT_FM
      MODULE PROCEDURE FMINT_IM
      MODULE PROCEDURE FMINT_ZM
   END INTERFACE

   INTERFACE LOG
      MODULE PROCEDURE FMLOG_FM
      MODULE PROCEDURE FMLOG_ZM
   END INTERFACE

   INTERFACE LOG10
      MODULE PROCEDURE FMLOG10_FM
      MODULE PROCEDURE FMLOG10_ZM
   END INTERFACE

   INTERFACE MATMUL
      MODULE PROCEDURE FMMATMUL_FM
      MODULE PROCEDURE FMMATMUL_IM
      MODULE PROCEDURE FMMATMUL_ZM
   END INTERFACE

   INTERFACE MAX
      MODULE PROCEDURE FMMAX_FM
      MODULE PROCEDURE FMMAX_IM
   END INTERFACE

   INTERFACE MAXEXPONENT
      MODULE PROCEDURE FMMAXEXPONENT_FM
   END INTERFACE

   INTERFACE MIN
      MODULE PROCEDURE FMMIN_FM
      MODULE PROCEDURE FMMIN_IM
   END INTERFACE

   INTERFACE MINEXPONENT
      MODULE PROCEDURE FMMINEXPONENT_FM
   END INTERFACE

   INTERFACE MOD
      MODULE PROCEDURE FMMOD_FM
      MODULE PROCEDURE FMMOD_IM
   END INTERFACE

   INTERFACE MODULO
      MODULE PROCEDURE FMMODULO_FM
      MODULE PROCEDURE FMMODULO_IM
   END INTERFACE

   INTERFACE NEAREST
      MODULE PROCEDURE FMNEAREST_FM
   END INTERFACE

   INTERFACE NINT
      MODULE PROCEDURE FMNINT_FM
      MODULE PROCEDURE FMNINT_IM
      MODULE PROCEDURE FMNINT_ZM
   END INTERFACE

   INTERFACE PRECISION
      MODULE PROCEDURE FMPRECISION_FM
      MODULE PROCEDURE FMPRECISION_ZM
   END INTERFACE

   INTERFACE RADIX
      MODULE PROCEDURE FMRADIX_FM
      MODULE PROCEDURE FMRADIX_IM
      MODULE PROCEDURE FMRADIX_ZM
   END INTERFACE

   INTERFACE RANGE
      MODULE PROCEDURE FMRANGE_FM
      MODULE PROCEDURE FMRANGE_IM
      MODULE PROCEDURE FMRANGE_ZM
   END INTERFACE

   INTERFACE REAL
      MODULE PROCEDURE FMREAL_FM
      MODULE PROCEDURE FMREAL_IM
      MODULE PROCEDURE FMREAL_ZM
   END INTERFACE

   INTERFACE RRSPACING
      MODULE PROCEDURE FMRRSPACING_FM
   END INTERFACE

   INTERFACE SCALE
      MODULE PROCEDURE FMSCALE_FM
      MODULE PROCEDURE FMSCALE_ZM
   END INTERFACE

   INTERFACE SETEXPONENT
      MODULE PROCEDURE FMSETEXPONENT_FM
   END INTERFACE

   INTERFACE SIGN
      MODULE PROCEDURE FMSIGN_FM
      MODULE PROCEDURE FMSIGN_IM
   END INTERFACE

 CONTAINS

!                                                             EPSILON

   FUNCTION FMEPSILON_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMEPSILON_FM
      INTENT (IN) :: MA
      CALL FMI2M(1,MTFM)
      CALL FMULP(MTFM,FMEPSILON_FM%MFM)
   END FUNCTION

!                                                                 EXP

   FUNCTION FMEXP_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMEXP_FM
      INTENT (IN) :: MA
      CALL FMEXP(MA%MFM,FMEXP_FM%MFM)
   END FUNCTION

   FUNCTION FMEXP_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMEXP_ZM
      INTENT (IN) :: MA
      CALL ZMEXP(MA%MZM,FMEXP_ZM%MZM)
   END FUNCTION

!                                                            EXPONENT

   FUNCTION FMEXPONENT_FM(MA)
      TYPE ( FM ) MA
      INTEGER FMEXPONENT_FM
      INTENT (IN) :: MA
      FMEXPONENT_FM = INT(MA%MFM(1))
   END FUNCTION

!                                                               FLOOR

   FUNCTION FMFLOOR_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMFLOOR_FM
      INTENT (IN) :: MA
      CALL FMINT(MA%MFM,MTFM)
      CALL FMSUB(MA%MFM,MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MFM,FMFLOOR_FM%MFM)
      ELSE IF (MA%MFM(-1) < 0) THEN
          CALL FMADDI(MTFM,-1)
          CALL FMEQ(MTFM,FMFLOOR_FM%MFM)
      ELSE
          CALL FMEQ(MTFM,FMFLOOR_FM%MFM)
      ENDIF
   END FUNCTION

   FUNCTION FMFLOOR_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMFLOOR_IM
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,FMFLOOR_IM%MIM)
   END FUNCTION

   FUNCTION FMFLOOR_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMFLOOR_ZM
      INTENT (IN) :: MA
      CALL FMINT(MA%MZM,MTFM)
      CALL FMSUB(MA%MZM,MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MZM,MVFM)
      ELSE IF (MA%MZM(-1) < 0) THEN
          CALL FMADDI(MTFM,-1)
          CALL FMEQ(MTFM,MVFM)
      ELSE
          CALL FMEQ(MTFM,MVFM)
      ENDIF
      CALL FMINT(MA%MZM(KPTIMU-1),MTFM)
      CALL FMSUB(MA%MZM(KPTIMU-1),MTFM,MUFM)
      IF (MUFM(2) == 0) THEN
          CALL FMEQ(MA%MZM(KPTIMU-1),MUFM)
      ELSE IF (MA%MZM(KPTIMU-1) < 0) THEN
          CALL FMADDI(MTFM,-1)
          CALL FMEQ(MTFM,MUFM)
      ELSE
          CALL FMEQ(MTFM,MUFM)
      ENDIF
      CALL ZMCMPX(MVFM,MUFM,FMFLOOR_ZM%MZM)
   END FUNCTION

!                                                            FRACTION

   FUNCTION FMFRACTION_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMFRACTION_FM
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,MTFM)
      MTFM(1) = 0
      CALL FMEQ(MTFM,FMFRACTION_FM%MFM)
   END FUNCTION

   FUNCTION FMFRACTION_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMFRACTION_ZM
      INTENT (IN) :: MA
      CALL ZMEQ(MA%MZM,MTZM)
      MTZM(1) = 0
      MTZM(KPTIMU+1) = 0
      CALL ZMEQ(MTZM,FMFRACTION_ZM%MZM)
   END FUNCTION

!                                                                HUGE

   FUNCTION FMHUGE_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMHUGE_FM
      INTENT (IN) :: MA
      CALL FMBIG(FMHUGE_FM%MFM)
   END FUNCTION

   FUNCTION FMHUGE_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMHUGE_IM
      INTENT (IN) :: MA
      CALL IMBIG(FMHUGE_IM%MIM)
   END FUNCTION

   FUNCTION FMHUGE_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMHUGE_ZM
      INTENT (IN) :: MA
      CALL FMBIG(MTFM)
      CALL ZMCMPX(MTFM,MTFM,FMHUGE_ZM%MZM)
   END FUNCTION

!                                                                 INT

   FUNCTION FMINT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( IM ) FMINT_FM
      INTENT (IN) :: MA
      CALL FMINT(MA%MFM,MTFM)
      CALL IMFM2I(MTFM,FMINT_FM%MIM)
   END FUNCTION

   FUNCTION FMINT_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMINT_IM
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,FMINT_IM%MIM)
   END FUNCTION

   FUNCTION FMINT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA
      TYPE ( IM ) FMINT_ZM
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMINT(MTFM,MUFM)
      CALL IMFM2I(MUFM,FMINT_ZM%MIM)
   END FUNCTION

!                                                                 LOG

   FUNCTION FMLOG_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMLOG_FM
      INTENT (IN) :: MA
      CALL FMLN(MA%MFM,FMLOG_FM%MFM)
   END FUNCTION

   FUNCTION FMLOG_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMLOG_ZM
      INTENT (IN) :: MA
      CALL ZMLN(MA%MZM,FMLOG_ZM%MZM)
   END FUNCTION

!                                                               LOG10

   FUNCTION FMLOG10_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMLOG10_FM
      INTENT (IN) :: MA
      CALL FMLG10(MA%MFM,FMLOG10_FM%MFM)
   END FUNCTION

   FUNCTION FMLOG10_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMLOG10_ZM
      INTENT (IN) :: MA
      CALL ZMLG10(MA%MZM,FMLOG10_ZM%MZM)
   END FUNCTION

!                                                              MATMUL

   FUNCTION FMMATMUL_FM(MA,MB)        RESULT(MC)
      USE FMVALS
      TYPE ( FM ) MA(:,:),MB(:,:)
      TYPE ( FM ), DIMENSION(SIZE(MA,DIM=1),SIZE(MB,DIM=2)) :: MC
      INTEGER I,J,K,NDSAVE
      INTENT (IN) :: MA,MB
      DO J = 1, SIZE(MA,DIM=1)
         DO K = 1, SIZE(MB,DIM=2)
         ENDDO
      ENDDO
      IF (SIZE(MA,DIM=2) == SIZE(MB,DIM=1)) THEN
          NDSAVE = NDIG
          J = MAX(NGRD52,2)
          NDIG = MIN(MAX(NDIG+J,2),NDG2MX)
          DO I = LBOUND(MA,DIM=1), UBOUND(MA,DIM=1)
             DO J = LBOUND(MB,DIM=2), UBOUND(MB,DIM=2)
                CALL FMI2M(0,MTFM)
                DO K = LBOUND(MA,DIM=2), UBOUND(MA,DIM=2)
                   CALL FMEQ2(MA(I,K)%MFM,MUFM,NDSAVE,NDIG)
                   CALL FMEQ2(MB(K,J)%MFM,MVFM,NDSAVE,NDIG)
                   CALL FMMPY(MUFM,MVFM,M01)
                   CALL FMADD_R1(MTFM,M01)
                ENDDO
                CALL FMEQ2_R1(MTFM,NDIG,NDSAVE)
                CALL FMEQ(MTFM,MC(I,J)%MFM)
             ENDDO
          ENDDO
          NDIG = NDSAVE
      ELSE
          CALL FMI2M(1,MTFM)
          CALL FMI2M(0,MUFM)
          CALL FMDIV(MTFM,MUFM,MVFM)
          DO I = 1, SIZE(MA,DIM=1)
             DO J = 1, SIZE(MB,DIM=2)
                CALL FMEQ(MVFM,MC(I,J)%MFM)
             ENDDO
          ENDDO
      ENDIF
   END FUNCTION

   FUNCTION FMMATMUL_IM(MA,MB)        RESULT(MC)
      USE FMVALS
      TYPE ( IM ) MA(:,:),MB(:,:)
      TYPE ( IM ), DIMENSION(SIZE(MA,DIM=1),SIZE(MB,DIM=2)) :: MC
      INTEGER I,J,K
      INTENT (IN) :: MA,MB
      DO J = 1, SIZE(MA,DIM=1)
         DO K = 1, SIZE(MB,DIM=2)
         ENDDO
      ENDDO
      IF (SIZE(MA,DIM=2) == SIZE(MB,DIM=1)) THEN
          DO I = LBOUND(MA,DIM=1), UBOUND(MA,DIM=1)
             DO J = LBOUND(MB,DIM=2), UBOUND(MB,DIM=2)
                CALL IMI2M(0,MTIM)
                DO K = LBOUND(MA,DIM=2), UBOUND(MA,DIM=2)
                   CALL IMMPY(MA(I,K)%MIM,MB(K,J)%MIM,M01)
                   CALL IMADD(MTIM,M01,MUIM)
                   CALL IMEQ(MUIM,MTIM)
                ENDDO
                CALL IMEQ(MTIM,MC(I,J)%MIM)
             ENDDO
          ENDDO
      ELSE
          CALL IMI2M(1,MTIM)
          CALL IMI2M(0,MUIM)
          CALL IMDIV(MTIM,MUIM,MC(1,1)%MIM)
          DO I = 1, SIZE(MA,DIM=1)
             DO J = 1, SIZE(MB,DIM=2)
                IF (I > 1 .OR. J > 1) CALL IMEQ(MC(1,1)%MIM,MC(I,J)%MIM)
             ENDDO
          ENDDO
      ENDIF
   END FUNCTION

   FUNCTION FMMATMUL_ZM(MA,MB)        RESULT(MC)
      USE FMVALS
      TYPE ( ZM ) MA(:,:),MB(:,:)
      TYPE ( ZM ), DIMENSION(SIZE(MA,DIM=1),SIZE(MB,DIM=2)) :: MC
      INTEGER I,J,K,NDSAVE
      INTENT (IN) :: MA,MB
      DO J = 1, SIZE(MA,DIM=1)
         DO K = 1, SIZE(MB,DIM=2)
         ENDDO
      ENDDO
      IF (SIZE(MA,DIM=2) == SIZE(MB,DIM=1)) THEN
          NDSAVE = NDIG
          J = MAX(NGRD52,2)
          NDIG = MIN(MAX(NDIG+J,2),NDG2MX)
          DO I = LBOUND(MA,DIM=1), UBOUND(MA,DIM=1)
             DO J = LBOUND(MB,DIM=2), UBOUND(MB,DIM=2)
                CALL ZMI2M(0,MTZM)
                DO K = LBOUND(MA,DIM=2), UBOUND(MA,DIM=2)
                   CALL ZMEQ2(MA(I,K)%MZM,MUZM,NDSAVE,NDIG)
                   CALL ZMEQ2(MB(K,J)%MZM,MVZM,NDSAVE,NDIG)
                   CALL ZMMPY(MUZM,MVZM,MZ02)
                   CALL ZMADD(MTZM,MZ02,MUZM)
                   CALL ZMEQ(MUZM,MTZM)
                ENDDO
                CALL ZMEQ2_R1(MTZM,NDIG,NDSAVE)
                CALL ZMEQ(MTZM,MC(I,J)%MZM)
             ENDDO
          ENDDO
          NDIG = NDSAVE
      ELSE
          CALL ZMI2M(1,MTZM)
          CALL ZMI2M(0,MUZM)
          CALL ZMDIV(MTZM,MUZM,MC(1,1)%MZM)
          DO I = 1, SIZE(MA,DIM=1)
             DO J = 1, SIZE(MB,DIM=2)
                IF (I > 1 .OR. J > 1) CALL ZMEQ(MC(1,1)%MZM,MC(I,J)%MZM)
             ENDDO
          ENDDO
      ENDIF
   END FUNCTION

!                                                                 MAX

   FUNCTION FMMAX_FM(MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMMAX_FM
      TYPE ( FM ), OPTIONAL :: MC,MD,ME,MF,MG,MH,MI,MJ
      INTENT (IN) :: MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ
      CALL FMMAX(MA%MFM,MB%MFM,MTFM)
      IF (PRESENT(MC)) THEN
          CALL FMMAX(MTFM,MC%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MD)) THEN
          CALL FMMAX(MTFM,MD%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(ME)) THEN
          CALL FMMAX(MTFM,ME%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MF)) THEN
          CALL FMMAX(MTFM,MF%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MG)) THEN
          CALL FMMAX(MTFM,MG%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MH)) THEN
          CALL FMMAX(MTFM,MH%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MI)) THEN
          CALL FMMAX(MTFM,MI%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MJ)) THEN
          CALL FMMAX(MTFM,MJ%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      CALL FMEQ(MTFM,FMMAX_FM%MFM)
   END FUNCTION

   FUNCTION FMMAX_IM(MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMMAX_IM
      TYPE ( IM ), OPTIONAL :: MC,MD,ME,MF,MG,MH,MI,MJ
      INTENT (IN) :: MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ
      CALL IMMAX(MA%MIM,MB%MIM,MTIM)
      IF (PRESENT(MC)) THEN
          CALL IMMAX(MTIM,MC%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MD)) THEN
          CALL IMMAX(MTIM,MD%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(ME)) THEN
          CALL IMMAX(MTIM,ME%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MF)) THEN
          CALL IMMAX(MTIM,MF%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MG)) THEN
          CALL IMMAX(MTIM,MG%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MH)) THEN
          CALL IMMAX(MTIM,MH%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MI)) THEN
          CALL IMMAX(MTIM,MI%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MJ)) THEN
          CALL IMMAX(MTIM,MJ%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      CALL IMEQ(MTIM,FMMAX_IM%MIM)
   END FUNCTION

!                                                         MAXEXPONENT

   FUNCTION FMMAXEXPONENT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMMAXEXPONENT_FM
      INTENT (IN) :: MA
      FMMAXEXPONENT_FM = INT(MXEXP) + 1
   END FUNCTION

!                                                                 MIN

   FUNCTION FMMIN_FM(MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMMIN_FM
      TYPE ( FM ), OPTIONAL :: MC,MD,ME,MF,MG,MH,MI,MJ
      INTENT (IN) :: MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ
      CALL FMMIN(MA%MFM,MB%MFM,MTFM)
      IF (PRESENT(MC)) THEN
          CALL FMMIN(MTFM,MC%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MD)) THEN
          CALL FMMIN(MTFM,MD%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(ME)) THEN
          CALL FMMIN(MTFM,ME%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MF)) THEN
          CALL FMMIN(MTFM,MF%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MG)) THEN
          CALL FMMIN(MTFM,MG%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MH)) THEN
          CALL FMMIN(MTFM,MH%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MI)) THEN
          CALL FMMIN(MTFM,MI%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      IF (PRESENT(MJ)) THEN
          CALL FMMIN(MTFM,MJ%MFM,MUFM)
          CALL FMEQ(MUFM,MTFM)
      ENDIF
      CALL FMEQ(MTFM,FMMIN_FM%MFM)
   END FUNCTION

   FUNCTION FMMIN_IM(MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMMIN_IM
      TYPE ( IM ), OPTIONAL :: MC,MD,ME,MF,MG,MH,MI,MJ
      INTENT (IN) :: MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ
      CALL IMMIN(MA%MIM,MB%MIM,MTIM)
      IF (PRESENT(MC)) THEN
          CALL IMMIN(MTIM,MC%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MD)) THEN
          CALL IMMIN(MTIM,MD%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(ME)) THEN
          CALL IMMIN(MTIM,ME%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MF)) THEN
          CALL IMMIN(MTIM,MF%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MG)) THEN
          CALL IMMIN(MTIM,MG%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MH)) THEN
          CALL IMMIN(MTIM,MH%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MI)) THEN
          CALL IMMIN(MTIM,MI%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      IF (PRESENT(MJ)) THEN
          CALL IMMIN(MTIM,MJ%MIM,MUIM)
          CALL IMEQ(MUIM,MTIM)
      ENDIF
      CALL IMEQ(MTIM,FMMIN_IM%MIM)
   END FUNCTION

!                                                         MINEXPONENT

   FUNCTION FMMINEXPONENT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMMINEXPONENT_FM
      INTENT (IN) :: MA
      FMMINEXPONENT_FM = -INT(MXEXP)
   END FUNCTION

!                                                                 MOD

   FUNCTION FMMOD_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMMOD_FM
      INTENT (IN) :: MA,MB
      CALL FMMOD(MA%MFM,MB%MFM,FMMOD_FM%MFM)
   END FUNCTION

   FUNCTION FMMOD_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMMOD_IM
      INTENT (IN) :: MA,MB
      CALL IMMOD(MA%MIM,MB%MIM,FMMOD_IM%MIM)
   END FUNCTION

!                                                              MODULO

   FUNCTION FMMODULO_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMMODULO_FM
      INTENT (IN) :: MA,MB
      CALL FMMOD(MA%MFM,MB%MFM,MTFM)
      IF (MTFM(2) /= 0) THEN
          IF ((MA%MFM(2) > 0 .AND. MA%MFM(-1) > 0 .AND.  &
               MB%MFM(2) > 0 .AND. MB%MFM(-1) < 0) .OR.  &
              (MA%MFM(2) > 0 .AND. MA%MFM(-1) < 0 .AND.  &
               MB%MFM(2) > 0 .AND. MB%MFM(-1) > 0)) THEN
              CALL FMADD_R1(MTFM,MB%MFM)
          ENDIF
      ENDIF
      CALL FMEQ(MTFM,FMMODULO_FM%MFM)
   END FUNCTION

   FUNCTION FMMODULO_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMMODULO_IM
      INTENT (IN) :: MA,MB
      CALL IMMOD(MA%MIM,MB%MIM,MTIM)
      IF (MTIM(2) /= 0) THEN
          IF ((MA%MIM(2) > 0 .AND. MA%MIM(-1) > 0 .AND.  &
               MB%MIM(2) > 0 .AND. MB%MIM(-1) < 0) .OR.  &
              (MA%MIM(2) > 0 .AND. MA%MIM(-1) < 0 .AND.  &
               MB%MIM(2) > 0 .AND. MB%MIM(-1) > 0)) THEN
              CALL IMADD(MTIM,MB%MIM,MUIM)
              CALL IMEQ(MUIM,MTIM)
          ENDIF
      ENDIF
      CALL IMEQ(MTIM,FMMODULO_IM%MIM)
   END FUNCTION

!                                                             NEAREST

   FUNCTION FMNEAREST_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMNEAREST_FM
      LOGICAL FMCOMP
      INTENT (IN) :: MA,MB
      IF (MA%MFM(2) == 0) THEN
          IF (MB%MFM(-1) > 0) THEN
              CALL FMBIG(MTFM)
              CALL FMI2M(1,MUFM)
              CALL FMDIV(MUFM,MTFM,FMNEAREST_FM%MFM)
          ELSE
              CALL FMBIG(MTFM)
              CALL FMI2M(-1,MUFM)
              CALL FMDIV(MUFM,MTFM,FMNEAREST_FM%MFM)
          ENDIF
      ELSE
          IF (MB%MFM(-1) > 0) THEN
              CALL FMULP(MA%MFM,MTFM)
              MTFM(-1) = 1
              CALL FMADD(MA%MFM,MTFM,MUFM)
              CALL FMULP(MUFM,MVFM)
              CALL FMABS(MVFM,MUFM)
              IF (FMCOMP(MTFM,'LE',MUFM)) THEN
                  CALL FMADD(MA%MFM,MTFM,FMNEAREST_FM%MFM)
              ELSE
                  CALL FMADD(MA%MFM,MUFM,FMNEAREST_FM%MFM)
              ENDIF
          ELSE
              CALL FMULP(MA%MFM,MTFM)
              MTFM(-1) = 1
              CALL FMSUB(MA%MFM,MTFM,MUFM)
              CALL FMULP(MUFM,MVFM)
              CALL FMABS(MVFM,MUFM)
              IF (FMCOMP(MTFM,'LE',MUFM)) THEN
                  CALL FMSUB(MA%MFM,MTFM,FMNEAREST_FM%MFM)
              ELSE
                  CALL FMSUB(MA%MFM,MUFM,FMNEAREST_FM%MFM)
              ENDIF
          ENDIF
      ENDIF
   END FUNCTION

!                                                                NINT

   FUNCTION FMNINT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      TYPE ( IM ) FMNINT_FM
      INTENT (IN) :: MA
      CALL FMNINT(MA%MFM,MTFM)
      CALL IMFM2I(MTFM,FMNINT_FM%MIM)
   END FUNCTION

   FUNCTION FMNINT_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMNINT_IM
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,FMNINT_IM%MIM)
   END FUNCTION

   FUNCTION FMNINT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA
      TYPE ( IM ) FMNINT_ZM
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMNINT(MTFM,MUFM)
      CALL IMFM2I(MUFM,FMNINT_ZM%MIM)
   END FUNCTION

!                                                           PRECISION

   FUNCTION FMPRECISION_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMPRECISION_FM
      INTENT (IN) :: MA
      FMPRECISION_FM = INT(LOG10(REAL(MBASE))*(NDIG-1) + 1)
   END FUNCTION

   FUNCTION FMPRECISION_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA
      INTEGER FMPRECISION_ZM
      INTENT (IN) :: MA
      FMPRECISION_ZM = INT(LOG10(REAL(MBASE))*(NDIG-1) + 1)
   END FUNCTION

!                                                               RADIX

   FUNCTION FMRADIX_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMRADIX_FM
      INTENT (IN) :: MA
      FMRADIX_FM = INT(MBASE)
   END FUNCTION

   FUNCTION FMRADIX_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA
      INTEGER FMRADIX_IM
      INTENT (IN) :: MA
      FMRADIX_IM = INT(MBASE)
   END FUNCTION

   FUNCTION FMRADIX_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA
      INTEGER FMRADIX_ZM
      INTENT (IN) :: MA
      FMRADIX_ZM = INT(MBASE)
   END FUNCTION

!                                                               RANGE

   FUNCTION FMRANGE_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA
      INTEGER FMRANGE_FM
      INTENT (IN) :: MA
      FMRANGE_FM = INT(MXEXP*LOG10(REAL(MBASE)))
   END FUNCTION

   FUNCTION FMRANGE_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA
      INTEGER FMRANGE_IM
      INTENT (IN) :: MA
      FMRANGE_IM = INT(NDIGMX*LOG10(REAL(MBASE)))
   END FUNCTION

   FUNCTION FMRANGE_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA
      INTEGER FMRANGE_ZM
      INTENT (IN) :: MA
      FMRANGE_ZM = INT(MXEXP*LOG10(REAL(MBASE)))
   END FUNCTION

!                                                                REAL

   FUNCTION FMREAL_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMREAL_FM
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,FMREAL_FM%MFM)
   END FUNCTION

   FUNCTION FMREAL_IM(MA)
      USE FMVALS
      TYPE ( FM ) FMREAL_IM
      TYPE ( IM ) MA
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,FMREAL_IM%MFM)
   END FUNCTION

   FUNCTION FMREAL_ZM(MA)
      USE FMVALS
      TYPE ( FM ) FMREAL_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,FMREAL_ZM%MFM)
   END FUNCTION

!                                                           RRSPACING

   FUNCTION FMRRSPACING_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMRRSPACING_FM
      INTENT (IN) :: MA
      CALL FMABS(MA%MFM,MTFM)
      MTFM(1) = NDIG
      CALL FMEQ(MTFM,FMRRSPACING_FM%MFM)
   END FUNCTION

!                                                               SCALE

   FUNCTION FMSCALE_FM(MA,L)
      USE FMVALS
      TYPE ( FM ) MA,FMSCALE_FM
      INTEGER L
      INTENT (IN) :: MA,L
      CALL FMEQ(MA%MFM,MTFM)
      IF (ABS(MTFM(1)+L) < MXEXP) THEN
          MTFM(1) = MTFM(1) + L
          CALL FMEQ(MTFM,FMSCALE_FM%MFM)
      ELSE
          CALL FMI2M(INT(MBASE),MUFM)
          CALL FMIPWR(MUFM,L,MVFM)
          CALL FMMPY(MA%MFM,MVFM,FMSCALE_FM%MFM)
      ENDIF
   END FUNCTION

   FUNCTION FMSCALE_ZM(MA,L)
      USE FMVALS
      INTEGER L
      TYPE ( ZM ) MA,FMSCALE_ZM
      INTENT (IN) :: MA,L
      CALL ZMEQ(MA%MZM,MTZM)
      IF (ABS(MTZM(1)+L) < MXEXP .AND. &
          ABS(MTZM(KPTIMU+1)+L) < MXEXP) THEN
          MTZM(1) = MTZM(1) + L
          MTZM(KPTIMU+1) = MTZM(KPTIMU+1) + L
          CALL ZMEQ(MTZM,FMSCALE_ZM%MZM)
      ELSE
          CALL ZMI2M(INT(MBASE),MUZM)
          CALL ZMIPWR(MUZM,L,MVZM)
          CALL ZMMPY(MA%MZM,MVZM,FMSCALE_ZM%MZM)
      ENDIF
   END FUNCTION

!                                                         SETEXPONENT

   FUNCTION FMSETEXPONENT_FM(MA,L)
      USE FMVALS
      TYPE ( FM ) MA,FMSETEXPONENT_FM
      INTEGER L
      INTENT (IN) :: MA,L
      CALL FMEQ(MA%MFM,MTFM)
      MTFM(1) = L
      CALL FMEQ(MTFM,FMSETEXPONENT_FM%MFM)
   END FUNCTION

!                                                                SIGN

   FUNCTION FMSIGN_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMSIGN_FM
      INTENT (IN) :: MA,MB
      CALL FMSIGN(MA%MFM,MB%MFM,FMSIGN_FM%MFM)
   END FUNCTION

   FUNCTION FMSIGN_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,FMSIGN_IM
      INTENT (IN) :: MA,MB
      CALL IMSIGN(MA%MIM,MB%MIM,FMSIGN_IM%MIM)
   END FUNCTION

 END MODULE FMZM_8

 MODULE FMZM_9
    USE FMZM_1

   INTERFACE SIN
      MODULE PROCEDURE FMSIN_FM
      MODULE PROCEDURE FMSIN_ZM
   END INTERFACE

   INTERFACE SINH
      MODULE PROCEDURE FMSINH_FM
      MODULE PROCEDURE FMSINH_ZM
   END INTERFACE

   INTERFACE SPACING
      MODULE PROCEDURE FMSPACING_FM
   END INTERFACE

   INTERFACE SQRT
      MODULE PROCEDURE FMSQRT_FM
      MODULE PROCEDURE FMSQRT_ZM
   END INTERFACE

   INTERFACE TAN
      MODULE PROCEDURE FMTAN_FM
      MODULE PROCEDURE FMTAN_ZM
   END INTERFACE

   INTERFACE TANH
      MODULE PROCEDURE FMTANH_FM
      MODULE PROCEDURE FMTANH_ZM
   END INTERFACE

   INTERFACE TINY
      MODULE PROCEDURE FMTINY_FM
      MODULE PROCEDURE FMTINY_IM
      MODULE PROCEDURE FMTINY_ZM
   END INTERFACE

   INTERFACE TO_FM
      MODULE PROCEDURE FM_I
      MODULE PROCEDURE FM_R
      MODULE PROCEDURE FM_D
      MODULE PROCEDURE FM_Z
      MODULE PROCEDURE FM_C
      MODULE PROCEDURE FM_FM
      MODULE PROCEDURE FM_IM
      MODULE PROCEDURE FM_ZM
      MODULE PROCEDURE FM_ST
   END INTERFACE

   INTERFACE TO_IM
      MODULE PROCEDURE IM_I
      MODULE PROCEDURE IM_R
      MODULE PROCEDURE IM_D
      MODULE PROCEDURE IM_Z
      MODULE PROCEDURE IM_C
      MODULE PROCEDURE IM_FM
      MODULE PROCEDURE IM_IM
      MODULE PROCEDURE IM_ZM
      MODULE PROCEDURE IM_ST
   END INTERFACE

   INTERFACE TO_ZM
      MODULE PROCEDURE ZM_I
      MODULE PROCEDURE ZM_R
      MODULE PROCEDURE ZM_D
      MODULE PROCEDURE ZM_Z
      MODULE PROCEDURE ZM_C
      MODULE PROCEDURE ZM_FM
      MODULE PROCEDURE ZM_IM
      MODULE PROCEDURE ZM_ZM
      MODULE PROCEDURE ZM_ST
   END INTERFACE

   INTERFACE TO_INT
      MODULE PROCEDURE FM_2INT
      MODULE PROCEDURE IM_2INT
      MODULE PROCEDURE ZM_2INT
   END INTERFACE

   INTERFACE TO_SP
      MODULE PROCEDURE FM_2SP
      MODULE PROCEDURE IM_2SP
      MODULE PROCEDURE ZM_2SP
   END INTERFACE

   INTERFACE TO_DP
      MODULE PROCEDURE FM_2DP
      MODULE PROCEDURE IM_2DP
      MODULE PROCEDURE ZM_2DP
   END INTERFACE

   INTERFACE TO_SPZ
      MODULE PROCEDURE FM_2SPZ
      MODULE PROCEDURE IM_2SPZ
      MODULE PROCEDURE ZM_2SPZ
   END INTERFACE

   INTERFACE TO_DPZ
      MODULE PROCEDURE FM_2DPZ
      MODULE PROCEDURE IM_2DPZ
      MODULE PROCEDURE ZM_2DPZ
   END INTERFACE


   INTERFACE FM_FORMAT
      MODULE PROCEDURE FMFORMAT_FM
   END INTERFACE

   INTERFACE IM_FORMAT
      MODULE PROCEDURE IMFORMAT_IM
   END INTERFACE

   INTERFACE ZM_FORMAT
      MODULE PROCEDURE ZMFORMAT_ZM
   END INTERFACE

   INTERFACE GCD
      MODULE PROCEDURE GCD_IM
   END INTERFACE

   INTERFACE MULTIPLY_MOD
      MODULE PROCEDURE MULTIPLYMOD_IM
   END INTERFACE

   INTERFACE POWER_MOD
      MODULE PROCEDURE POWERMOD_IM
   END INTERFACE

   INTERFACE FM_RANDOM_SEED
      MODULE PROCEDURE FM_SEED
   END INTERFACE

   INTERFACE BERNOULLI
      MODULE PROCEDURE FMBERNOULLI_FM
   END INTERFACE

   INTERFACE BETA
      MODULE PROCEDURE FMBETA_FM
   END INTERFACE

   INTERFACE BINOMIAL
      MODULE PROCEDURE FMBINOMIAL_FM
   END INTERFACE

   INTERFACE FACTORIAL
      MODULE PROCEDURE FMFACTORIAL_FM
   END INTERFACE

   INTERFACE GAMMA
      MODULE PROCEDURE FMGAMMA_FM
   END INTERFACE

   INTERFACE INCOMPLETE_BETA
      MODULE PROCEDURE FMINCOMPLETE_BETA_FM
   END INTERFACE

   INTERFACE INCOMPLETE_GAMMA1
      MODULE PROCEDURE FMINCOMPLETE_GAMMA1_FM
   END INTERFACE

   INTERFACE INCOMPLETE_GAMMA2
      MODULE PROCEDURE FMINCOMPLETE_GAMMA2_FM
   END INTERFACE

   INTERFACE LOG_GAMMA
      MODULE PROCEDURE FMLOG_GAMMA_FM
   END INTERFACE

   INTERFACE POLYGAMMA
      MODULE PROCEDURE FMPOLYGAMMA_FM
   END INTERFACE

   INTERFACE POCHHAMMER
      MODULE PROCEDURE FMPOCHHAMMER_FM
   END INTERFACE

   INTERFACE PSI
      MODULE PROCEDURE FMPSI_FM
   END INTERFACE

 CONTAINS

!                                                                 SIN

   FUNCTION FMSIN_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSIN_FM
      INTENT (IN) :: MA
      CALL FMSIN(MA%MFM,FMSIN_FM%MFM)
   END FUNCTION

   FUNCTION FMSIN_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSIN_ZM
      INTENT (IN) :: MA
      CALL ZMSIN(MA%MZM,FMSIN_ZM%MZM)
   END FUNCTION

!                                                                SINH

   FUNCTION FMSINH_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSINH_FM
      INTENT (IN) :: MA
      CALL FMSINH(MA%MFM,FMSINH_FM%MFM)
   END FUNCTION

   FUNCTION FMSINH_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSINH_ZM
      INTENT (IN) :: MA
      CALL ZMSINH(MA%MZM,FMSINH_ZM%MZM)
   END FUNCTION

!                                                             SPACING

   FUNCTION FMSPACING_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSPACING_FM
      INTENT (IN) :: MA
      CALL FMABS(MA%MFM,MTFM)
      CALL FMULP(MTFM,FMSPACING_FM%MFM)
   END FUNCTION

!                                                                SQRT

   FUNCTION FMSQRT_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMSQRT_FM
      INTENT (IN) :: MA
      CALL FMSQRT(MA%MFM,FMSQRT_FM%MFM)
   END FUNCTION

   FUNCTION FMSQRT_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMSQRT_ZM
      INTENT (IN) :: MA
      CALL ZMSQRT(MA%MZM,FMSQRT_ZM%MZM)
   END FUNCTION

!                                                                 TAN

   FUNCTION FMTAN_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMTAN_FM
      INTENT (IN) :: MA
      CALL FMTAN(MA%MFM,FMTAN_FM%MFM)
   END FUNCTION

   FUNCTION FMTAN_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMTAN_ZM
      INTENT (IN) :: MA
      CALL ZMTAN(MA%MZM,FMTAN_ZM%MZM)
   END FUNCTION

!                                                                TANH

   FUNCTION FMTANH_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMTANH_FM
      INTENT (IN) :: MA
      CALL FMTANH(MA%MFM,FMTANH_FM%MFM)
   END FUNCTION

   FUNCTION FMTANH_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMTANH_ZM
      INTENT (IN) :: MA
      CALL ZMTANH(MA%MZM,FMTANH_ZM%MZM)
   END FUNCTION

!                                                                TINY

   FUNCTION FMTINY_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMTINY_FM
      INTEGER J
      INTENT (IN) :: MA
      FMTINY_FM%MFM(-1) = 1
      FMTINY_FM%MFM(0) = NINT(NDIG*ALOGM2)
      FMTINY_FM%MFM(1) = -MXEXP
      FMTINY_FM%MFM(2) = 1
      DO J = 3, NDIG+1
         FMTINY_FM%MFM(J) = 0
      ENDDO
   END FUNCTION

   FUNCTION FMTINY_IM(MA)
      USE FMVALS
      TYPE ( IM ) MA,FMTINY_IM
      INTENT (IN) :: MA
      CALL IMI2M(1,FMTINY_IM%MIM)
   END FUNCTION

   FUNCTION FMTINY_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) MA,FMTINY_ZM
      INTEGER J
      INTENT (IN) :: MA
      MTFM(-1) = 1
      MTFM(0) = NINT(NDIG*ALOGM2)
      MTFM(1) = -MXEXP
      MTFM(2) = 1
      DO J = 3, NDIG+1
         MTFM(J) = 0
      ENDDO
      CALL ZMCMPX(MTFM,MTFM,FMTINY_ZM%MZM)
   END FUNCTION

!                                                               TO_FM

   FUNCTION FM_I(IVAL)
      USE FMVALS
      TYPE ( FM ) FM_I
      INTEGER IVAL
      INTENT (IN) :: IVAL
      CALL FMI2M(IVAL,FM_I%MFM)
   END FUNCTION

   FUNCTION FM_R(R)
      USE FMVALS
      TYPE ( FM ) FM_R
      REAL R
      INTENT (IN) :: R
      CALL FMSP2M(R,FM_R%MFM)
   END FUNCTION

   FUNCTION FM_D(D)
      USE FMVALS
      TYPE ( FM ) FM_D
      DOUBLE PRECISION D
      INTENT (IN) :: D
      CALL FMDP2M(D,FM_D%MFM)
   END FUNCTION

   FUNCTION FM_Z(Z)
      USE FMVALS
      TYPE ( FM ) FM_Z
      COMPLEX Z
      INTENT (IN) :: Z
      CALL FMSP2M(REAL(Z),FM_Z%MFM)
   END FUNCTION

   FUNCTION FM_C(C)
      USE FMVALS
      TYPE ( FM ) FM_C
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),FM_C%MFM)
   END FUNCTION

   FUNCTION FM_FM(MA)
      USE FMVALS
      TYPE ( FM ) FM_FM,MA
      INTENT (IN) :: MA
      CALL FMEQ(MA%MFM,FM_FM%MFM)
   END FUNCTION

   FUNCTION FM_IM(MA)
      USE FMVALS
      TYPE ( FM ) FM_IM
      TYPE ( IM ) MA
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,FM_IM%MFM)
   END FUNCTION

   FUNCTION FM_ST(ST)
      USE FMVALS
      TYPE ( FM ) FM_ST
      CHARACTER(*) :: ST
      INTENT (IN) :: ST
      CALL FMST2M(ST,FM_ST%MFM)
   END FUNCTION

   FUNCTION FM_ZM(MA)
      USE FMVALS
      TYPE ( FM ) FM_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,FM_ZM%MFM)
   END FUNCTION

!                                                               TO_IM

   FUNCTION IM_I(IVAL)
      USE FMVALS
      TYPE ( IM ) IM_I
      INTEGER IVAL
      INTENT (IN) :: IVAL
      CALL IMI2M(IVAL,IM_I%MIM)
   END FUNCTION

   FUNCTION IM_R(R)
      USE FMVALS
      TYPE ( IM ) IM_R
      REAL R
      CHARACTER(25) :: ST
      INTENT (IN) :: R
      IF (ABS(R) < HUGE(1)) THEN
          IVAL = INT(R)
          CALL IMI2M(IVAL,IM_R%MIM)
      ELSE
          WRITE (ST,'(E25.16)') R
          CALL IMST2M(ST,IM_R%MIM)
      ENDIF
   END FUNCTION

   FUNCTION IM_D(D)
      USE FMVALS
      TYPE ( IM ) IM_D
      DOUBLE PRECISION D
      CHARACTER(25) :: ST
      INTENT (IN) :: D
      IF (ABS(D) < HUGE(1)) THEN
          IVAL = INT(D)
          CALL IMI2M(IVAL,IM_D%MIM)
      ELSE
          WRITE (ST,'(E25.16)') D
          CALL IMST2M(ST,IM_D%MIM)
      ENDIF
   END FUNCTION

   FUNCTION IM_Z(Z)
      USE FMVALS
      TYPE ( IM ) IM_Z
      COMPLEX Z
      REAL R
      CHARACTER(25) :: ST
      INTENT (IN) :: Z
      R = REAL(Z)
      IF (ABS(R) < HUGE(1)) THEN
          IVAL = INT(R)
          CALL IMI2M(IVAL,IM_Z%MIM)
      ELSE
          WRITE (ST,'(E25.16)') R
          CALL IMST2M(ST,IM_Z%MIM)
      ENDIF
   END FUNCTION

   FUNCTION IM_C(C)
      USE FMVALS
      TYPE ( IM ) IM_C
      COMPLEX (KIND(0.0D0)) :: C
      DOUBLE PRECISION D
      CHARACTER(25) :: ST
      INTENT (IN) :: C
      D = REAL(C)
      IF (ABS(D) < HUGE(1)) THEN
          IVAL = INT(D)
          CALL IMI2M(IVAL,IM_C%MIM)
      ELSE
          WRITE (ST,'(E25.16)') D
          CALL IMST2M(ST,IM_C%MIM)
      ENDIF
   END FUNCTION

   FUNCTION IM_FM(MA)
      USE FMVALS
      TYPE ( IM ) IM_FM
      TYPE ( FM ) MA
      INTENT (IN) :: MA
      CALL IMFM2I(MA%MFM,IM_FM%MIM)
   END FUNCTION

   FUNCTION IM_IM(MA)
      USE FMVALS
      TYPE ( IM ) IM_IM,MA
      INTENT (IN) :: MA
      CALL IMEQ(MA%MIM,IM_IM%MIM)
   END FUNCTION

   FUNCTION IM_ST(ST)
      USE FMVALS
      TYPE ( IM ) IM_ST
      CHARACTER(*) :: ST
      INTENT (IN) :: ST
      CALL IMST2M(ST,IM_ST%MIM)
   END FUNCTION

   FUNCTION IM_ZM(MA)
      USE FMVALS
      TYPE ( IM ) IM_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL IMFM2I(MTFM,IM_ZM%MIM)
   END FUNCTION

!                                                               TO_ZM

   FUNCTION ZM_I(IVAL)
      USE FMVALS
      TYPE ( ZM ) ZM_I
      INTEGER IVAL
      INTENT (IN) :: IVAL
      CALL ZMI2M(IVAL,ZM_I%MZM)
   END FUNCTION

   FUNCTION ZM_R(R)
      USE FMVALS
      TYPE ( ZM ) ZM_R
      REAL R
      INTENT (IN) :: R
      CALL FMSP2M(R,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,ZM_R%MZM)
   END FUNCTION

   FUNCTION ZM_D(D)
      USE FMVALS
      TYPE ( ZM ) ZM_D
      DOUBLE PRECISION D
      INTENT (IN) :: D
      CALL FMDP2M(D,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,ZM_D%MZM)
   END FUNCTION

   FUNCTION ZM_Z(Z)
      USE FMVALS
      TYPE ( ZM ) ZM_Z
      COMPLEX Z
      INTENT (IN) :: Z
      CALL ZMZ2M(Z,ZM_Z%MZM)
   END FUNCTION

   FUNCTION ZM_C(C)
      USE FMVALS
      TYPE ( ZM ) ZM_C
      COMPLEX (KIND(0.0D0)) :: C
      INTENT (IN) :: C
      CALL FMDP2M(REAL(C,KIND(0.0D0)),MTFM)
      CALL FMDP2M(AIMAG(C),MUFM)
      CALL ZMCMPX(MTFM,MUFM,ZM_C%MZM)
   END FUNCTION

   FUNCTION ZM_FM(MA)
      USE FMVALS
      TYPE ( ZM ) ZM_FM
      TYPE ( FM ) MA
      INTENT (IN) :: MA
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MA%MFM,MUFM,ZM_FM%MZM)
   END FUNCTION

   FUNCTION ZM_IM(MA)
      USE FMVALS
      TYPE ( ZM ) ZM_IM
      TYPE ( IM ) MA
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMI2M(0,MUFM)
      CALL ZMCMPX(MTFM,MUFM,ZM_IM%MZM)
   END FUNCTION

   FUNCTION ZM_ST(ST)
      USE FMVALS
      TYPE ( ZM ) ZM_ST
      CHARACTER(*) :: ST
      INTENT (IN) :: ST
      CALL ZMST2M(ST,ZM_ST%MZM)
   END FUNCTION

   FUNCTION ZM_ZM(MA)
      USE FMVALS
      TYPE ( ZM ) ZM_ZM,MA
      INTENT (IN) :: MA
      CALL ZMEQ(MA%MZM,ZM_ZM%MZM)
   END FUNCTION

!                                                              TO_INT

   FUNCTION FM_2INT(MA)
      TYPE ( FM ) MA
      INTEGER FM_2INT
      INTENT (IN) :: MA
      CALL FMM2I(MA%MFM,FM_2INT)
   END FUNCTION

   FUNCTION IM_2INT(MA)
      TYPE ( IM ) MA
      INTEGER IM_2INT
      INTENT (IN) :: MA
      CALL IMM2I(MA%MIM,IM_2INT)
   END FUNCTION

   FUNCTION ZM_2INT(MA)
      TYPE ( ZM ) MA
      INTEGER ZM_2INT
      INTENT (IN) :: MA
      CALL ZMM2I(MA%MZM,ZM_2INT)
   END FUNCTION

!                                                               TO_SP

   FUNCTION FM_2SP(MA)
      TYPE ( FM ) MA
      REAL FM_2SP
      INTENT (IN) :: MA
      CALL FMM2SP(MA%MFM,FM_2SP)
   END FUNCTION

   FUNCTION IM_2SP(MA)
      TYPE ( IM ) MA
      REAL IM_2SP
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMM2SP(MTFM,IM_2SP)
   END FUNCTION

   FUNCTION ZM_2SP(MA)
      TYPE ( ZM ) MA
      REAL ZM_2SP
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2SP(MTFM,ZM_2SP)
   END FUNCTION

!                                                               TO_DP

   FUNCTION FM_2DP(MA)
      TYPE ( FM ) MA
      DOUBLE PRECISION FM_2DP
      INTENT (IN) :: MA
      CALL FMM2DP(MA%MFM,FM_2DP)
   END FUNCTION

   FUNCTION IM_2DP(MA)
      TYPE ( IM ) MA
      DOUBLE PRECISION IM_2DP
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMM2DP(MTFM,IM_2DP)
   END FUNCTION

   FUNCTION ZM_2DP(MA)
      TYPE ( ZM ) MA
      DOUBLE PRECISION ZM_2DP
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,ZM_2DP)
   END FUNCTION

!                                                              TO_SPZ

   FUNCTION FM_2SPZ(MA)
      TYPE ( FM ) MA
      COMPLEX FM_2SPZ
      REAL R
      INTENT (IN) :: MA
      CALL FMM2SP(MA%MFM,R)
      FM_2SPZ = CMPLX( R , 0.0 )
   END FUNCTION

   FUNCTION IM_2SPZ(MA)
      TYPE ( IM ) MA
      COMPLEX IM_2SPZ
      REAL R
      INTENT (IN) :: MA
      CALL IMI2FM(MA%MIM,MTFM)
      CALL FMM2SP(MTFM,R)
      IM_2SPZ = CMPLX( R , 0.0 )
   END FUNCTION

   FUNCTION ZM_2SPZ(MA)
      TYPE ( ZM ) MA
      COMPLEX ZM_2SPZ
      INTENT (IN) :: MA
      CALL ZMM2Z(MA%MZM,ZM_2SPZ)
   END FUNCTION

!                                                              TO_DPZ

   FUNCTION FM_2DPZ(MA)
      TYPE ( FM ) MA
      COMPLEX (KIND(0.0D0)) :: FM_2DPZ
      DOUBLE PRECISION D
      INTENT (IN) :: MA
      CALL FMM2DP(MA%MFM,D)
      FM_2DPZ = CMPLX( D , 0.0D0 , KIND(0.0D0) )
   END FUNCTION

   FUNCTION IM_2DPZ(MA)
      TYPE ( IM ) MA
      COMPLEX (KIND(0.0D0)) :: IM_2DPZ
      DOUBLE PRECISION D
      INTENT (IN) :: MA
      CALL IMM2DP(MA%MIM,D)
      IM_2DPZ = CMPLX( D , 0.0D0 , KIND(0.0D0) )
   END FUNCTION

   FUNCTION ZM_2DPZ(MA)
      TYPE ( ZM ) MA
      COMPLEX (KIND(0.0D0)) :: ZM_2DPZ
      DOUBLE PRECISION D1,D2
      INTENT (IN) :: MA
      CALL ZMREAL(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,D1)
      CALL ZMIMAG(MA%MZM,MTFM)
      CALL FMM2DP(MTFM,D2)
      ZM_2DPZ = CMPLX( D1 , D2 , KIND(0.0D0) )
   END FUNCTION

!                                                         FM_RANDOM_SEED

      SUBROUTINE FM_SEED(PUT,GET,SIZE)

!  Interface routine for FM_RANDOM_SEED, used to initialize the random sequence
!  from FM_RANDOM_NUMBER.

!  Like the Fortran intrinsic function RANDOM_SEED, exactly one of the three
!  arguments must be present, and the call should be with an argument keyword.

!  CALL FM_RANDOM_SEED(SIZE=J)   returns J=7 to the calling program, indicating
!       that the seed array has length 7.

!  CALL FM_RANDOM_SEED(GET=SEED)   returns SEED(1) through SEED(7) as the current
!       seed for the generator, but see the comments in routine FM_RANDOM_NUMBER.

!  CALL FM_RANDOM_SEED(PUT=SEED)   initializes the FM_RANDOM_NUMBER generator.

!  The typical usage is to call FM_RANDOM_SEED once with PUT defined as an
!  integer array of length 7 containing seven seed values used to initialize
!  the generator.  This initializes the table used by the mixed congruential
!  generator.  Then each call to FM_RANDOM_NUMBER gets the next random value.

!  This example seeds the generator and then fills the double precision array R
!  with random values between 0 and 1.

!        SEED = (/ 314159,265358,979323,846264,338327,950288,419716 /)
!        CALL FM_RANDOM_SEED(PUT=SEED)
!        DO J = 1, N
!           CALL FM_RANDOM_NUMBER(R(J))
!        ENDDO

      USE FMVALS

      IMPLICIT NONE

      INTEGER, OPTIONAL, INTENT(IN)  :: PUT(7)
      INTEGER, OPTIONAL, INTENT(OUT) :: GET(7)
      INTEGER, OPTIONAL, INTENT(OUT) :: SIZE

      REAL (KIND(1.0D0)) :: MSAVE
      INTEGER J,K


      IF (PRESENT(SIZE)) THEN
          SIZE = 7
          RETURN
      ENDIF
      MSAVE = MBASE
      MBASE = MBRAND
      IF (PRESENT(PUT)) THEN
          K = 10**7
          CALL IMI2M(ABS(PUT(1)),MRNX)
          DO J = 2, 7
             CALL IMMPYI(MRNX,K,MTIM)
             CALL IMI2M(ABS(PUT(J)),M04)
             CALL IMADD(MTIM,M04,MRNX)
          ENDDO
          CALL IMST2M('2070613773952029032014000773560846464373793273739',M04)
          CALL IMMOD(MRNX,M04,MTIM)
          CALL IMEQ(MTIM,MRNX)
          START_RANDOM_SEQUENCE = 1
          MBASE = MSAVE
          RETURN
      ENDIF
      IF (PRESENT(GET)) THEN
          K = 10**7
          CALL IMI2M(K,M05)
          CALL IMEQ(MRNX,M04)
          DO J = 7, 1, -1
             CALL IMMOD(M04,M05,M06)
             CALL IMM2I(M06,GET(J))
             CALL IMDIVI(M04,K,MTIM)
             CALL IMEQ(MTIM,M04)
          ENDDO
          MBASE = MSAVE
          RETURN
      ENDIF
      END SUBROUTINE

!                                                         FM_FORMAT

   FUNCTION FMFORMAT_FM(FMT,MA)
      USE FMVALS
      CHARACTER(*) :: FMT
      TYPE ( FM ) MA
      CHARACTER(200) :: FMFORMAT_FM
      INTENT (IN) :: FMT,MA
      CALL FMFORM(FMT,MA%MFM,FMFORMAT_FM)
   END FUNCTION

!                                                         IM_FORMAT

   FUNCTION IMFORMAT_IM(FMT,MA)
      USE FMVALS
      CHARACTER(*) :: FMT
      CHARACTER(200) :: IMFORMAT_IM
      TYPE ( IM ) MA
      INTENT (IN) :: FMT,MA
      CALL IMFORM(FMT,MA%MIM,IMFORMAT_IM)
   END FUNCTION

!                                                         ZM_FORMAT

   FUNCTION ZMFORMAT_ZM(FMTR,FMTI,MA)
      USE FMVALS
      CHARACTER(*) :: FMTR,FMTI
      CHARACTER(200) :: ZMFORMAT_ZM
      TYPE ( ZM ) MA
      INTENT (IN) :: FMTR,FMTI,MA
      CALL ZMFORM(FMTR,FMTI,MA%MZM,ZMFORMAT_ZM)
   END FUNCTION

!                                                         GCD

   FUNCTION GCD_IM(MA,MB)
      USE FMVALS
      TYPE ( IM ) MA,MB,GCD_IM
      INTENT (IN) :: MA,MB
      CALL IMGCD(MA%MIM,MB%MIM,GCD_IM%MIM)
   END FUNCTION

!                                                         MULTIPLY_MOD

   FUNCTION MULTIPLYMOD_IM(MA,MB,MC)
      USE FMVALS
      TYPE ( IM ) MA,MB,MC,MULTIPLYMOD_IM
      INTENT (IN) :: MA,MB,MC
      CALL IMMPYM(MA%MIM,MB%MIM,MC%MIM,MULTIPLYMOD_IM%MIM)
   END FUNCTION

!                                                         POWER_MOD

   FUNCTION POWERMOD_IM(MA,MB,MC)
      USE FMVALS
      TYPE ( IM ) MA,MB,MC,POWERMOD_IM
      INTENT (IN) :: MA,MB,MC
      CALL IMPMOD(MA%MIM,MB%MIM,MC%MIM,POWERMOD_IM%MIM)
   END FUNCTION

!                                                         BERNOULLI

   FUNCTION FMBERNOULLI_FM(N)
      USE FMVALS
      TYPE ( FM ) FMBERNOULLI_FM
      INTEGER N
      INTENT (IN) :: N
      CALL FMI2M(1,MTFM)
      CALL FMBERN(N,MTFM,FMBERNOULLI_FM%MFM)
   END FUNCTION

!                                                         BETA

   FUNCTION FMBETA_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMBETA_FM
      INTENT (IN) :: MA,MB
      CALL FMBETA(MA%MFM,MB%MFM,FMBETA_FM%MFM)
   END FUNCTION

!                                                         BINOMIAL

   FUNCTION FMBINOMIAL_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMBINOMIAL_FM
      INTENT (IN) :: MA,MB
      CALL FMCOMB(MA%MFM,MB%MFM,FMBINOMIAL_FM%MFM)
   END FUNCTION

!                                                         FACTORIAL

   FUNCTION FMFACTORIAL_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMFACTORIAL_FM
      INTENT (IN) :: MA
      CALL FMFACT(MA%MFM,FMFACTORIAL_FM%MFM)
   END FUNCTION

!                                                         GAMMA

   FUNCTION FMGAMMA_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMGAMMA_FM
      INTENT (IN) :: MA
      CALL FMGAM(MA%MFM,FMGAMMA_FM%MFM)
   END FUNCTION

!                                                         INCOMPLETE_BETA

   FUNCTION FMINCOMPLETE_BETA_FM(MX,MA,MB)
      USE FMVALS
      TYPE ( FM ) MX,MA,MB,FMINCOMPLETE_BETA_FM
      INTENT (IN) :: MX,MA,MB
      CALL FMIBTA(MX%MFM,MA%MFM,MB%MFM,FMINCOMPLETE_BETA_FM%MFM)
   END FUNCTION

!                                                         INCOMPLETE_GAMMA1

   FUNCTION FMINCOMPLETE_GAMMA1_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMINCOMPLETE_GAMMA1_FM
      INTENT (IN) :: MA,MB
      CALL FMIGM1(MA%MFM,MB%MFM,FMINCOMPLETE_GAMMA1_FM%MFM)
   END FUNCTION

!                                                         INCOMPLETE_GAMMA2

   FUNCTION FMINCOMPLETE_GAMMA2_FM(MA,MB)
      USE FMVALS
      TYPE ( FM ) MA,MB,FMINCOMPLETE_GAMMA2_FM
      INTENT (IN) :: MA,MB
      CALL FMIGM2(MA%MFM,MB%MFM,FMINCOMPLETE_GAMMA2_FM%MFM)
   END FUNCTION

!                                                         LOG_GAMMA

   FUNCTION FMLOG_GAMMA_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMLOG_GAMMA_FM
      INTENT (IN) :: MA
      CALL FMLNGM(MA%MFM,FMLOG_GAMMA_FM%MFM)
   END FUNCTION

!                                                         POLYGAMMA

   FUNCTION FMPOLYGAMMA_FM(N,MA)
      USE FMVALS
      TYPE ( FM ) MA,FMPOLYGAMMA_FM
      INTEGER N
      INTENT (IN) :: N,MA
      CALL FMPGAM(N,MA%MFM,FMPOLYGAMMA_FM%MFM)
   END FUNCTION

!                                                         POCHHAMMER

   FUNCTION FMPOCHHAMMER_FM(MA,N)
      USE FMVALS
      TYPE ( FM ) MA,FMPOCHHAMMER_FM
      INTEGER N
      INTENT (IN) :: N,MA
      CALL FMPOCH(MA%MFM,N,FMPOCHHAMMER_FM%MFM)
   END FUNCTION

!                                                         PSI

   FUNCTION FMPSI_FM(MA)
      USE FMVALS
      TYPE ( FM ) MA,FMPSI_FM
      INTENT (IN) :: MA
      CALL FMPSI(MA%MFM,FMPSI_FM%MFM)
   END FUNCTION

 END MODULE FMZM_9

 MODULE FMZM

   USE FMZM_1
   USE FMZM_2
   USE FMZM_3
   USE FMZM_4
   USE FMZM_5
   USE FMZM_6
   USE FMZM_7
   USE FMZM_8
   USE FMZM_9

 END MODULE FMZM

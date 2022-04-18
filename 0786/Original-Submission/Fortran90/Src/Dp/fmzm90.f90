    MODULE fmzm

!   FMZM 1.1                  David M. Smith               3-23-97

!   This module extends the definition of Fortran-90 arithmetic and
!   function operations so they also apply to multiple precision
!   numbers, using version 1.1 of FMLIB and ZMLIB.
!   There are three multiple precision data types:
!      FM  (multiple precision real)
!      IM  (multiple precision integer)
!      ZM  (multiple precision complex)

!   Some the the interface routines assume that the precision chosen
!   in the calling program (using FMSET or ZMSET) represents more
!   significant digits than does the machine's double precision.

!   All the functions defined in this module are standard Fortran-90
!   functions, except for several direct conversion functions:

!   TO_FM is a function for converting other types of numbers to type
!   FM.  Note that TO_FM(3.12) converts the REAL constant to FM, but
!   it is accurate only to single precision.  TO_FM(3.12D0) agrees
!   with 3.12 to double precision accuracy, and TO_FM('3.12') or
!   TO_FM(312)/TO_FM(100) agrees to full FM accuracy.

!   TO_IM converts to type IM, and TO_ZM converts to type ZM.

!   Functions are also supplied for converting the three multiple
!   precision types to the other numeric data types:
!      TO_INT   converts to machine precision integer
!      TO_SP    converts to single precision
!      TO_DP    converts to double precision
!      TO_SPZ   converts to single precision complex
!      TO_DPZ   converts to double precision complex

!   WARNING:   When multiple precision type declarations are inserted
!              in an existing program, take care in converting functions
!              like DBLE(X), where X has been declared as a multiple
!              precision type.  If X was single precision in the
!              original program, then replacing the DBLE(X) by TO_DP(X)
!              in the new version could lose accuracy.
!              For this reason, the Fortran type-conversion functions
!              defined in this module assume that results should be
!              multiple precision whenever inputs are.  Examples:
!              DBLE(TO_FM('1.23E+123456')) is type FM
!              REAL(TO_FM('1.23E+123456')) is type FM
!              REAL(TO_ZM('3.12+4.56i'))   is type FM   = TO_FM('3.12')
!              INT(TO_FM('1.23'))          is type IM   = TO_IM(1)
!              INT(TO_IM('1E+23'))         is type IM
!              CMPLX(TO_FM('1.23'),TO_FM('4.56')) is type ZM

!   Programs using this module may sometimes need to call FM, IM, or
!   ZM routines directly.  This is normally the case when routines are
!   needed that are not Fortran-90 intrinsics, such as the formatting
!   subroutine FMFORM.  In a program using this module, suppose MAFM
!   has been declared with TYPE ( FM ) MAFM.  To use the routine FMFORM,
!   which expects the second argument to be an array and not a derived
!   type, the call would have to be CALL FMFORM('F65.60',MAFM%MFM,ST1)
!   so that the array contained in MAFM is passed.

!   As an alternative so the user can refer directly to the FM-, IM-,
!   and ZM-type variables and avoid the cumbersome "%MFM" suffixes,
!   this module contains a collection of interface routines to supply
!   any needed argument conversions.  For each FM, IM, and ZM routine
!   that is designed to be called by the user, there is also a version
!   that assumes any multiple-precision arguments are derived types
!   instead of arrays.  Each interface routine has the same name as
!   the original with an underscore after the first two letters of the
!   routine name.  To convert the number to a character string with
!   F65.60 format, use CALL FM_FORM('F65.60',MAFM,ST1) if MAFM is of
!   TYPE ( FM ), or use CALL FMFORM('F65.60',MA,ST1) if MA is declared
!   as an array.  All the routines shown below may be used this way.

!   For each of the operations =, .EQ., .NE., .GT., .GE., .LT., .LE.,
!   +, -, *, /, and **, the interface module defines all mixed mode
!   variations involving one of the three multiple precision derived
!   types and another argument having one of the types:
!   { integer, real, double, complex, complex double, FM, IM, ZM }.
!   So mixed mode expressions such as
!         MAFM = 12
!         MAFM = MAFM + 1
!         IF (ABS(MAFM).LT.1.0D-23) THEN
!   are handled correctly.

!   Not all the named functions are defined for all three multiple
!   precision derived types, so the list below shows which can be used.
!   The labels "real", "integer", and "complex" refer to types FM, IM,
!   and ZM respectively, "string" means the function accepts character
!   strings (e.g., TO_FM('3.45')), and "other" means the function can
!   accept any of the machine precision data types integer, real,
!   double, complex, or complex double.  For functions that accept two
!   or more arguments, like ATAN2 or MAX, all the arguments must be of
!   the same type.

!   AVAILABLE OPERATIONS:

!      =
!      +
!      -
!      *
!      /
!      **
!      .EQ.
!      .NE.
!      .GT.
!      .GE.
!      .LT.
!      .LE.
!      ABS          real    integer    complex
!      ACOS         real               complex
!      AIMAG                           complex
!      AINT         real               complex
!      ANINT        real               complex
!      ASIN         real               complex
!      ATAN         real               complex
!      ATAN2        real
!      BTEST                integer
!      CEILING      real               complex
!      CMPLX        real    integer
!      CONJ                            complex
!      COS          real               complex
!      COSH         real               complex
!      DBLE         real    integer    complex
!      DIGITS       real    integer    complex
!      DIM          real    integer
!      DINT         real               complex
!      DOTPRODUCT   real    integer    complex
!      EPSILON      real
!      EXP          real               complex
!      EXPONENT     real
!      FLOOR        real    integer    complex
!      FRACTION     real               complex
!      HUGE         real    integer    complex
!      INT          real    integer    complex
!      LOG          real               complex
!      LOG10        real               complex
!      MATMUL       real    integer    complex
!      MAX          real    integer
!      MAXEXPONENT  real
!      MIN          real    integer
!      MINEXPONENT  real
!      MOD          real    integer
!      MODULO       real    integer
!      NEAREST      real
!      NINT         real    integer    complex
!      PRECISION    real               complex
!      RADIX        real    integer    complex
!      RANGE        real    integer    complex
!      REAL         real    integer    complex
!      RRSPACING    real
!      SCALE        real               complex
!      SETEXPONENT  real
!      SIGN         real    integer
!      SIN          real               complex
!      SINH         real               complex
!      SPACING      real
!      SQRT         real               complex
!      TAN          real               complex
!      TANH         real               complex
!      TINY         real    integer    complex
!      TO_FM        real    integer    complex    string    other
!      TO_IM        real    integer    complex    string    other
!      TO_ZM        real    integer    complex    string    other
!      TO_INT       real    integer    complex
!      TO_SP        real    integer    complex
!      TO_DP        real    integer    complex
!      TO_SPZ       real    integer    complex
!      TO_DPZ       real    integer    complex

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

! .. Use Statements ..
      USE fmzmcommon
! ..
! .. Generic Interface Blocks ..
      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE fmeq_ifm
        MODULE PROCEDURE fmeq_iim
        MODULE PROCEDURE fmeq_izm
        MODULE PROCEDURE fmeq_rfm
        MODULE PROCEDURE fmeq_rim
        MODULE PROCEDURE fmeq_rzm
        MODULE PROCEDURE fmeq_dfm
        MODULE PROCEDURE fmeq_dim
        MODULE PROCEDURE fmeq_dzm
        MODULE PROCEDURE fmeq_zfm
        MODULE PROCEDURE fmeq_zim
        MODULE PROCEDURE fmeq_zzm
        MODULE PROCEDURE fmeq_cfm
        MODULE PROCEDURE fmeq_cim
        MODULE PROCEDURE fmeq_czm
        MODULE PROCEDURE fmeq_fmi
        MODULE PROCEDURE fmeq_fmr
        MODULE PROCEDURE fmeq_fmd
        MODULE PROCEDURE fmeq_fmz
        MODULE PROCEDURE fmeq_fmc
        MODULE PROCEDURE fmeq_fmfm
        MODULE PROCEDURE fmeq_fmim
        MODULE PROCEDURE fmeq_fmzm
        MODULE PROCEDURE fmeq_imi
        MODULE PROCEDURE fmeq_imr
        MODULE PROCEDURE fmeq_imd
        MODULE PROCEDURE fmeq_imz
        MODULE PROCEDURE fmeq_imc
        MODULE PROCEDURE fmeq_imfm
        MODULE PROCEDURE fmeq_imim
        MODULE PROCEDURE fmeq_imzm
        MODULE PROCEDURE fmeq_zmi
        MODULE PROCEDURE fmeq_zmr
        MODULE PROCEDURE fmeq_zmd
        MODULE PROCEDURE fmeq_zmz
        MODULE PROCEDURE fmeq_zmc
        MODULE PROCEDURE fmeq_zmfm
        MODULE PROCEDURE fmeq_zmim
        MODULE PROCEDURE fmeq_zmzm
      END INTERFACE
      INTERFACE OPERATOR (==)
        MODULE PROCEDURE fmleq_ifm
        MODULE PROCEDURE fmleq_iim
        MODULE PROCEDURE fmleq_izm
        MODULE PROCEDURE fmleq_rfm
        MODULE PROCEDURE fmleq_rim
        MODULE PROCEDURE fmleq_rzm
        MODULE PROCEDURE fmleq_dfm
        MODULE PROCEDURE fmleq_dim
        MODULE PROCEDURE fmleq_dzm
        MODULE PROCEDURE fmleq_zfm
        MODULE PROCEDURE fmleq_zim
        MODULE PROCEDURE fmleq_zzm
        MODULE PROCEDURE fmleq_cfm
        MODULE PROCEDURE fmleq_cim
        MODULE PROCEDURE fmleq_czm
        MODULE PROCEDURE fmleq_fmi
        MODULE PROCEDURE fmleq_fmr
        MODULE PROCEDURE fmleq_fmd
        MODULE PROCEDURE fmleq_fmz
        MODULE PROCEDURE fmleq_fmc
        MODULE PROCEDURE fmleq_fmfm
        MODULE PROCEDURE fmleq_fmim
        MODULE PROCEDURE fmleq_fmzm
        MODULE PROCEDURE fmleq_imi
        MODULE PROCEDURE fmleq_imr
        MODULE PROCEDURE fmleq_imd
        MODULE PROCEDURE fmleq_imz
        MODULE PROCEDURE fmleq_imc
        MODULE PROCEDURE fmleq_imfm
        MODULE PROCEDURE fmleq_imim
        MODULE PROCEDURE fmleq_imzm
        MODULE PROCEDURE fmleq_zmi
        MODULE PROCEDURE fmleq_zmr
        MODULE PROCEDURE fmleq_zmd
        MODULE PROCEDURE fmleq_zmz
        MODULE PROCEDURE fmleq_zmc
        MODULE PROCEDURE fmleq_zmfm
        MODULE PROCEDURE fmleq_zmim
        MODULE PROCEDURE fmleq_zmzm
      END INTERFACE
      INTERFACE OPERATOR (/=)
        MODULE PROCEDURE fmlne_ifm
        MODULE PROCEDURE fmlne_iim
        MODULE PROCEDURE fmlne_izm
        MODULE PROCEDURE fmlne_rfm
        MODULE PROCEDURE fmlne_rim
        MODULE PROCEDURE fmlne_rzm
        MODULE PROCEDURE fmlne_dfm
        MODULE PROCEDURE fmlne_dim
        MODULE PROCEDURE fmlne_dzm
        MODULE PROCEDURE fmlne_zfm
        MODULE PROCEDURE fmlne_zim
        MODULE PROCEDURE fmlne_zzm
        MODULE PROCEDURE fmlne_cfm
        MODULE PROCEDURE fmlne_cim
        MODULE PROCEDURE fmlne_czm
        MODULE PROCEDURE fmlne_fmi
        MODULE PROCEDURE fmlne_fmr
        MODULE PROCEDURE fmlne_fmd
        MODULE PROCEDURE fmlne_fmz
        MODULE PROCEDURE fmlne_fmc
        MODULE PROCEDURE fmlne_fmfm
        MODULE PROCEDURE fmlne_fmim
        MODULE PROCEDURE fmlne_fmzm
        MODULE PROCEDURE fmlne_imi
        MODULE PROCEDURE fmlne_imr
        MODULE PROCEDURE fmlne_imd
        MODULE PROCEDURE fmlne_imz
        MODULE PROCEDURE fmlne_imc
        MODULE PROCEDURE fmlne_imfm
        MODULE PROCEDURE fmlne_imim
        MODULE PROCEDURE fmlne_imzm
        MODULE PROCEDURE fmlne_zmi
        MODULE PROCEDURE fmlne_zmr
        MODULE PROCEDURE fmlne_zmd
        MODULE PROCEDURE fmlne_zmz
        MODULE PROCEDURE fmlne_zmc
        MODULE PROCEDURE fmlne_zmfm
        MODULE PROCEDURE fmlne_zmim
        MODULE PROCEDURE fmlne_zmzm
      END INTERFACE
      INTERFACE OPERATOR (>)
        MODULE PROCEDURE fmlgt_ifm
        MODULE PROCEDURE fmlgt_iim
        MODULE PROCEDURE fmlgt_rfm
        MODULE PROCEDURE fmlgt_rim
        MODULE PROCEDURE fmlgt_dfm
        MODULE PROCEDURE fmlgt_dim
        MODULE PROCEDURE fmlgt_fmi
        MODULE PROCEDURE fmlgt_fmr
        MODULE PROCEDURE fmlgt_fmd
        MODULE PROCEDURE fmlgt_fmfm
        MODULE PROCEDURE fmlgt_fmim
        MODULE PROCEDURE fmlgt_imi
        MODULE PROCEDURE fmlgt_imr
        MODULE PROCEDURE fmlgt_imd
        MODULE PROCEDURE fmlgt_imfm
        MODULE PROCEDURE fmlgt_imim
      END INTERFACE
      INTERFACE OPERATOR (>=)
        MODULE PROCEDURE fmlge_ifm
        MODULE PROCEDURE fmlge_iim
        MODULE PROCEDURE fmlge_rfm
        MODULE PROCEDURE fmlge_rim
        MODULE PROCEDURE fmlge_dfm
        MODULE PROCEDURE fmlge_dim
        MODULE PROCEDURE fmlge_fmi
        MODULE PROCEDURE fmlge_fmr
        MODULE PROCEDURE fmlge_fmd
        MODULE PROCEDURE fmlge_fmfm
        MODULE PROCEDURE fmlge_fmim
        MODULE PROCEDURE fmlge_imi
        MODULE PROCEDURE fmlge_imr
        MODULE PROCEDURE fmlge_imd
        MODULE PROCEDURE fmlge_imfm
        MODULE PROCEDURE fmlge_imim
      END INTERFACE
      INTERFACE OPERATOR (<)
        MODULE PROCEDURE fmllt_ifm
        MODULE PROCEDURE fmllt_iim
        MODULE PROCEDURE fmllt_rfm
        MODULE PROCEDURE fmllt_rim
        MODULE PROCEDURE fmllt_dfm
        MODULE PROCEDURE fmllt_dim
        MODULE PROCEDURE fmllt_fmi
        MODULE PROCEDURE fmllt_fmr
        MODULE PROCEDURE fmllt_fmd
        MODULE PROCEDURE fmllt_fmfm
        MODULE PROCEDURE fmllt_fmim
        MODULE PROCEDURE fmllt_imi
        MODULE PROCEDURE fmllt_imr
        MODULE PROCEDURE fmllt_imd
        MODULE PROCEDURE fmllt_imfm
        MODULE PROCEDURE fmllt_imim
      END INTERFACE
      INTERFACE OPERATOR (<=)
        MODULE PROCEDURE fmlle_ifm
        MODULE PROCEDURE fmlle_iim
        MODULE PROCEDURE fmlle_rfm
        MODULE PROCEDURE fmlle_rim
        MODULE PROCEDURE fmlle_dfm
        MODULE PROCEDURE fmlle_dim
        MODULE PROCEDURE fmlle_fmi
        MODULE PROCEDURE fmlle_fmr
        MODULE PROCEDURE fmlle_fmd
        MODULE PROCEDURE fmlle_fmfm
        MODULE PROCEDURE fmlle_fmim
        MODULE PROCEDURE fmlle_imi
        MODULE PROCEDURE fmlle_imr
        MODULE PROCEDURE fmlle_imd
        MODULE PROCEDURE fmlle_imfm
        MODULE PROCEDURE fmlle_imim
      END INTERFACE
      INTERFACE OPERATOR (+)
        MODULE PROCEDURE fmadd_ifm
        MODULE PROCEDURE fmadd_iim
        MODULE PROCEDURE fmadd_izm
        MODULE PROCEDURE fmadd_rfm
        MODULE PROCEDURE fmadd_rim
        MODULE PROCEDURE fmadd_rzm
        MODULE PROCEDURE fmadd_dfm
        MODULE PROCEDURE fmadd_dim
        MODULE PROCEDURE fmadd_dzm
        MODULE PROCEDURE fmadd_zfm
        MODULE PROCEDURE fmadd_zim
        MODULE PROCEDURE fmadd_zzm
        MODULE PROCEDURE fmadd_cfm
        MODULE PROCEDURE fmadd_cim
        MODULE PROCEDURE fmadd_czm
        MODULE PROCEDURE fmadd_fmi
        MODULE PROCEDURE fmadd_fmr
        MODULE PROCEDURE fmadd_fmd
        MODULE PROCEDURE fmadd_fmz
        MODULE PROCEDURE fmadd_fmc
        MODULE PROCEDURE fmadd_fmfm
        MODULE PROCEDURE fmadd_fmim
        MODULE PROCEDURE fmadd_fmzm
        MODULE PROCEDURE fmadd_imi
        MODULE PROCEDURE fmadd_imr
        MODULE PROCEDURE fmadd_imd
        MODULE PROCEDURE fmadd_imz
        MODULE PROCEDURE fmadd_imc
        MODULE PROCEDURE fmadd_imfm
        MODULE PROCEDURE fmadd_imim
        MODULE PROCEDURE fmadd_imzm
        MODULE PROCEDURE fmadd_zmi
        MODULE PROCEDURE fmadd_zmr
        MODULE PROCEDURE fmadd_zmd
        MODULE PROCEDURE fmadd_zmz
        MODULE PROCEDURE fmadd_zmc
        MODULE PROCEDURE fmadd_zmfm
        MODULE PROCEDURE fmadd_zmim
        MODULE PROCEDURE fmadd_zmzm
        MODULE PROCEDURE fmadd_fm
        MODULE PROCEDURE fmadd_im
        MODULE PROCEDURE fmadd_zm
      END INTERFACE
      INTERFACE OPERATOR (-)
        MODULE PROCEDURE fmsub_ifm
        MODULE PROCEDURE fmsub_iim
        MODULE PROCEDURE fmsub_izm
        MODULE PROCEDURE fmsub_rfm
        MODULE PROCEDURE fmsub_rim
        MODULE PROCEDURE fmsub_rzm
        MODULE PROCEDURE fmsub_dfm
        MODULE PROCEDURE fmsub_dim
        MODULE PROCEDURE fmsub_dzm
        MODULE PROCEDURE fmsub_zfm
        MODULE PROCEDURE fmsub_zim
        MODULE PROCEDURE fmsub_zzm
        MODULE PROCEDURE fmsub_cfm
        MODULE PROCEDURE fmsub_cim
        MODULE PROCEDURE fmsub_czm
        MODULE PROCEDURE fmsub_fmi
        MODULE PROCEDURE fmsub_fmr
        MODULE PROCEDURE fmsub_fmd
        MODULE PROCEDURE fmsub_fmz
        MODULE PROCEDURE fmsub_fmc
        MODULE PROCEDURE fmsub_fmfm
        MODULE PROCEDURE fmsub_fmim
        MODULE PROCEDURE fmsub_fmzm
        MODULE PROCEDURE fmsub_imi
        MODULE PROCEDURE fmsub_imr
        MODULE PROCEDURE fmsub_imd
        MODULE PROCEDURE fmsub_imz
        MODULE PROCEDURE fmsub_imc
        MODULE PROCEDURE fmsub_imfm
        MODULE PROCEDURE fmsub_imim
        MODULE PROCEDURE fmsub_imzm
        MODULE PROCEDURE fmsub_zmi
        MODULE PROCEDURE fmsub_zmr
        MODULE PROCEDURE fmsub_zmd
        MODULE PROCEDURE fmsub_zmz
        MODULE PROCEDURE fmsub_zmc
        MODULE PROCEDURE fmsub_zmfm
        MODULE PROCEDURE fmsub_zmim
        MODULE PROCEDURE fmsub_zmzm
        MODULE PROCEDURE fmsub_fm
        MODULE PROCEDURE fmsub_im
        MODULE PROCEDURE fmsub_zm
      END INTERFACE
      INTERFACE OPERATOR (*)
        MODULE PROCEDURE fmmpy_ifm
        MODULE PROCEDURE fmmpy_iim
        MODULE PROCEDURE fmmpy_izm
        MODULE PROCEDURE fmmpy_rfm
        MODULE PROCEDURE fmmpy_rim
        MODULE PROCEDURE fmmpy_rzm
        MODULE PROCEDURE fmmpy_dfm
        MODULE PROCEDURE fmmpy_dim
        MODULE PROCEDURE fmmpy_dzm
        MODULE PROCEDURE fmmpy_zfm
        MODULE PROCEDURE fmmpy_zim
        MODULE PROCEDURE fmmpy_zzm
        MODULE PROCEDURE fmmpy_cfm
        MODULE PROCEDURE fmmpy_cim
        MODULE PROCEDURE fmmpy_czm
        MODULE PROCEDURE fmmpy_fmi
        MODULE PROCEDURE fmmpy_fmr
        MODULE PROCEDURE fmmpy_fmd
        MODULE PROCEDURE fmmpy_fmz
        MODULE PROCEDURE fmmpy_fmc
        MODULE PROCEDURE fmmpy_fmfm
        MODULE PROCEDURE fmmpy_fmim
        MODULE PROCEDURE fmmpy_fmzm
        MODULE PROCEDURE fmmpy_imi
        MODULE PROCEDURE fmmpy_imr
        MODULE PROCEDURE fmmpy_imd
        MODULE PROCEDURE fmmpy_imz
        MODULE PROCEDURE fmmpy_imc
        MODULE PROCEDURE fmmpy_imfm
        MODULE PROCEDURE fmmpy_imim
        MODULE PROCEDURE fmmpy_imzm
        MODULE PROCEDURE fmmpy_zmi
        MODULE PROCEDURE fmmpy_zmr
        MODULE PROCEDURE fmmpy_zmd
        MODULE PROCEDURE fmmpy_zmz
        MODULE PROCEDURE fmmpy_zmc
        MODULE PROCEDURE fmmpy_zmfm
        MODULE PROCEDURE fmmpy_zmim
        MODULE PROCEDURE fmmpy_zmzm
      END INTERFACE
      INTERFACE OPERATOR (/)
        MODULE PROCEDURE fmdiv_ifm
        MODULE PROCEDURE fmdiv_iim
        MODULE PROCEDURE fmdiv_izm
        MODULE PROCEDURE fmdiv_rfm
        MODULE PROCEDURE fmdiv_rim
        MODULE PROCEDURE fmdiv_rzm
        MODULE PROCEDURE fmdiv_dfm
        MODULE PROCEDURE fmdiv_dim
        MODULE PROCEDURE fmdiv_dzm
        MODULE PROCEDURE fmdiv_zfm
        MODULE PROCEDURE fmdiv_zim
        MODULE PROCEDURE fmdiv_zzm
        MODULE PROCEDURE fmdiv_cfm
        MODULE PROCEDURE fmdiv_cim
        MODULE PROCEDURE fmdiv_czm
        MODULE PROCEDURE fmdiv_fmi
        MODULE PROCEDURE fmdiv_fmr
        MODULE PROCEDURE fmdiv_fmd
        MODULE PROCEDURE fmdiv_fmz
        MODULE PROCEDURE fmdiv_fmc
        MODULE PROCEDURE fmdiv_fmfm
        MODULE PROCEDURE fmdiv_fmim
        MODULE PROCEDURE fmdiv_fmzm
        MODULE PROCEDURE fmdiv_imi
        MODULE PROCEDURE fmdiv_imr
        MODULE PROCEDURE fmdiv_imd
        MODULE PROCEDURE fmdiv_imz
        MODULE PROCEDURE fmdiv_imc
        MODULE PROCEDURE fmdiv_imfm
        MODULE PROCEDURE fmdiv_imim
        MODULE PROCEDURE fmdiv_imzm
        MODULE PROCEDURE fmdiv_zmi
        MODULE PROCEDURE fmdiv_zmr
        MODULE PROCEDURE fmdiv_zmd
        MODULE PROCEDURE fmdiv_zmz
        MODULE PROCEDURE fmdiv_zmc
        MODULE PROCEDURE fmdiv_zmfm
        MODULE PROCEDURE fmdiv_zmim
        MODULE PROCEDURE fmdiv_zmzm
      END INTERFACE
      INTERFACE OPERATOR (**)
        MODULE PROCEDURE fmpwr_ifm
        MODULE PROCEDURE fmpwr_iim
        MODULE PROCEDURE fmpwr_izm
        MODULE PROCEDURE fmpwr_rfm
        MODULE PROCEDURE fmpwr_rim
        MODULE PROCEDURE fmpwr_rzm
        MODULE PROCEDURE fmpwr_dfm
        MODULE PROCEDURE fmpwr_dim
        MODULE PROCEDURE fmpwr_dzm
        MODULE PROCEDURE fmpwr_zfm
        MODULE PROCEDURE fmpwr_zim
        MODULE PROCEDURE fmpwr_zzm
        MODULE PROCEDURE fmpwr_cfm
        MODULE PROCEDURE fmpwr_cim
        MODULE PROCEDURE fmpwr_czm
        MODULE PROCEDURE fmpwr_fmi
        MODULE PROCEDURE fmpwr_fmr
        MODULE PROCEDURE fmpwr_fmd
        MODULE PROCEDURE fmpwr_fmz
        MODULE PROCEDURE fmpwr_fmc
        MODULE PROCEDURE fmpwr_fmfm
        MODULE PROCEDURE fmpwr_fmim
        MODULE PROCEDURE fmpwr_fmzm
        MODULE PROCEDURE fmpwr_imi
        MODULE PROCEDURE fmpwr_imr
        MODULE PROCEDURE fmpwr_imd
        MODULE PROCEDURE fmpwr_imz
        MODULE PROCEDURE fmpwr_imc
        MODULE PROCEDURE fmpwr_imfm
        MODULE PROCEDURE fmpwr_imim
        MODULE PROCEDURE fmpwr_imzm
        MODULE PROCEDURE fmpwr_zmi
        MODULE PROCEDURE fmpwr_zmr
        MODULE PROCEDURE fmpwr_zmd
        MODULE PROCEDURE fmpwr_zmz
        MODULE PROCEDURE fmpwr_zmc
        MODULE PROCEDURE fmpwr_zmfm
        MODULE PROCEDURE fmpwr_zmim
        MODULE PROCEDURE fmpwr_zmzm
      END INTERFACE
      INTERFACE abs
        MODULE PROCEDURE fmabs_fm
        MODULE PROCEDURE fmabs_im
        MODULE PROCEDURE fmabs_zm
      END INTERFACE
      INTERFACE acos
        MODULE PROCEDURE fmacos_fm
        MODULE PROCEDURE fmacos_zm
      END INTERFACE
      INTERFACE aimag
        MODULE PROCEDURE fmaimag_zm
      END INTERFACE
      INTERFACE aint
        MODULE PROCEDURE fmaint_fm
        MODULE PROCEDURE fmaint_zm
      END INTERFACE
      INTERFACE anint
        MODULE PROCEDURE fmanint_fm
        MODULE PROCEDURE fmanint_zm
      END INTERFACE
      INTERFACE asin
        MODULE PROCEDURE fmasin_fm
        MODULE PROCEDURE fmasin_zm
      END INTERFACE
      INTERFACE atan
        MODULE PROCEDURE fmatan_fm
        MODULE PROCEDURE fmatan_zm
      END INTERFACE
      INTERFACE atan2
        MODULE PROCEDURE fmatan2_fm
      END INTERFACE
      INTERFACE btest
        MODULE PROCEDURE fmbtest_im
      END INTERFACE
      INTERFACE ceiling
        MODULE PROCEDURE fmceiling_fm
        MODULE PROCEDURE fmceiling_zm
      END INTERFACE
      INTERFACE cmplx
        MODULE PROCEDURE fmcmplx_fm
        MODULE PROCEDURE fmcmplx_im
      END INTERFACE
      INTERFACE conjg
        MODULE PROCEDURE fmconjg_zm
      END INTERFACE
      INTERFACE cos
        MODULE PROCEDURE fmcos_fm
        MODULE PROCEDURE fmcos_zm
      END INTERFACE
      INTERFACE cosh
        MODULE PROCEDURE fmcosh_fm
        MODULE PROCEDURE fmcosh_zm
      END INTERFACE
      INTERFACE dble
        MODULE PROCEDURE fmdble_fm
        MODULE PROCEDURE fmdble_im
        MODULE PROCEDURE fmdble_zm
      END INTERFACE
      INTERFACE digits
        MODULE PROCEDURE fmdigits_fm
        MODULE PROCEDURE fmdigits_im
        MODULE PROCEDURE fmdigits_zm
      END INTERFACE
      INTERFACE dim
        MODULE PROCEDURE fmdim_fm
        MODULE PROCEDURE fmdim_im
      END INTERFACE
      INTERFACE dint
        MODULE PROCEDURE fmdint_fm
        MODULE PROCEDURE fmdint_zm
      END INTERFACE
      INTERFACE dotproduct
        MODULE PROCEDURE fmdotproduct_fm
        MODULE PROCEDURE fmdotproduct_im
        MODULE PROCEDURE fmdotproduct_zm
      END INTERFACE
      INTERFACE epsilon
        MODULE PROCEDURE fmepsilon_fm
      END INTERFACE
      INTERFACE exp
        MODULE PROCEDURE fmexp_fm
        MODULE PROCEDURE fmexp_zm
      END INTERFACE
      INTERFACE exponent
        MODULE PROCEDURE fmexponent_fm
      END INTERFACE
      INTERFACE floor
        MODULE PROCEDURE fmfloor_fm
        MODULE PROCEDURE fmfloor_im
        MODULE PROCEDURE fmfloor_zm
      END INTERFACE
      INTERFACE fraction
        MODULE PROCEDURE fmfraction_fm
        MODULE PROCEDURE fmfraction_zm
      END INTERFACE
      INTERFACE huge
        MODULE PROCEDURE fmhuge_fm
        MODULE PROCEDURE fmhuge_im
        MODULE PROCEDURE fmhuge_zm
      END INTERFACE
      INTERFACE int
        MODULE PROCEDURE fmint_fm
        MODULE PROCEDURE fmint_im
        MODULE PROCEDURE fmint_zm
      END INTERFACE
      INTERFACE log
        MODULE PROCEDURE fmlog_fm
        MODULE PROCEDURE fmlog_zm
      END INTERFACE
      INTERFACE log10
        MODULE PROCEDURE fmlog10_fm
        MODULE PROCEDURE fmlog10_zm
      END INTERFACE
      INTERFACE matmul
        MODULE PROCEDURE fmmatmul_fm
        MODULE PROCEDURE fmmatmul_im
        MODULE PROCEDURE fmmatmul_zm
      END INTERFACE
      INTERFACE max
        MODULE PROCEDURE fmmax_fm
        MODULE PROCEDURE fmmax_im
      END INTERFACE
      INTERFACE maxexponent
        MODULE PROCEDURE fmmaxexponent_fm
      END INTERFACE
      INTERFACE min
        MODULE PROCEDURE fmmin_fm
        MODULE PROCEDURE fmmin_im
      END INTERFACE
      INTERFACE minexponent
        MODULE PROCEDURE fmminexponent_fm
      END INTERFACE
      INTERFACE mod
        MODULE PROCEDURE fmmod_fm
        MODULE PROCEDURE fmmod_im
      END INTERFACE
      INTERFACE modulo
        MODULE PROCEDURE fmmodulo_fm
        MODULE PROCEDURE fmmodulo_im
      END INTERFACE
      INTERFACE nearest
        MODULE PROCEDURE fmnearest_fm
      END INTERFACE
      INTERFACE nint
        MODULE PROCEDURE fmnint_fm
        MODULE PROCEDURE fmnint_im
        MODULE PROCEDURE fmnint_zm
      END INTERFACE
      INTERFACE precision
        MODULE PROCEDURE fmprecision_fm
        MODULE PROCEDURE fmprecision_zm
      END INTERFACE
      INTERFACE radix
        MODULE PROCEDURE fmradix_fm
        MODULE PROCEDURE fmradix_im
        MODULE PROCEDURE fmradix_zm
      END INTERFACE
      INTERFACE range
        MODULE PROCEDURE fmrange_fm
        MODULE PROCEDURE fmrange_im
        MODULE PROCEDURE fmrange_zm
      END INTERFACE
      INTERFACE real
        MODULE PROCEDURE fmreal_fm
        MODULE PROCEDURE fmreal_im
        MODULE PROCEDURE fmreal_zm
      END INTERFACE
      INTERFACE rrspacing
        MODULE PROCEDURE fmrrspacing_fm
      END INTERFACE
      INTERFACE scale
        MODULE PROCEDURE fmscale_fm
        MODULE PROCEDURE fmscale_zm
      END INTERFACE
      INTERFACE setexponent
        MODULE PROCEDURE fmsetexponent_fm
      END INTERFACE
      INTERFACE sign
        MODULE PROCEDURE fmsign_fm
        MODULE PROCEDURE fmsign_im
      END INTERFACE
      INTERFACE sin
        MODULE PROCEDURE fmsin_fm
        MODULE PROCEDURE fmsin_zm
      END INTERFACE
      INTERFACE sinh
        MODULE PROCEDURE fmsinh_fm
        MODULE PROCEDURE fmsinh_zm
      END INTERFACE
      INTERFACE spacing
        MODULE PROCEDURE fmspacing_fm
      END INTERFACE
      INTERFACE sqrt
        MODULE PROCEDURE fmsqrt_fm
        MODULE PROCEDURE fmsqrt_zm
      END INTERFACE
      INTERFACE tan
        MODULE PROCEDURE fmtan_fm
        MODULE PROCEDURE fmtan_zm
      END INTERFACE
      INTERFACE tanh
        MODULE PROCEDURE fmtanh_fm
        MODULE PROCEDURE fmtanh_zm
      END INTERFACE
      INTERFACE tiny
        MODULE PROCEDURE fmtiny_fm
        MODULE PROCEDURE fmtiny_im
        MODULE PROCEDURE fmtiny_zm
      END INTERFACE
      INTERFACE to_fm
        MODULE PROCEDURE fm_i
        MODULE PROCEDURE fm_r
        MODULE PROCEDURE fm_d
        MODULE PROCEDURE fm_z
        MODULE PROCEDURE fm_c
        MODULE PROCEDURE fm_fm
        MODULE PROCEDURE fm_im
        MODULE PROCEDURE fm_zm
        MODULE PROCEDURE fm_st
      END INTERFACE
      INTERFACE to_im
        MODULE PROCEDURE im_i
        MODULE PROCEDURE im_r
        MODULE PROCEDURE im_d
        MODULE PROCEDURE im_z
        MODULE PROCEDURE im_c
        MODULE PROCEDURE im_fm
        MODULE PROCEDURE im_im
        MODULE PROCEDURE im_zm
        MODULE PROCEDURE im_st
      END INTERFACE
      INTERFACE to_zm
        MODULE PROCEDURE zm_i
        MODULE PROCEDURE zm_r
        MODULE PROCEDURE zm_d
        MODULE PROCEDURE zm_z
        MODULE PROCEDURE zm_c
        MODULE PROCEDURE zm_fm
        MODULE PROCEDURE zm_im
        MODULE PROCEDURE zm_zm
        MODULE PROCEDURE zm_st
      END INTERFACE
      INTERFACE to_int
        MODULE PROCEDURE fm_2int
        MODULE PROCEDURE im_2int
        MODULE PROCEDURE zm_2int
      END INTERFACE
      INTERFACE to_sp
        MODULE PROCEDURE fm_2sp
        MODULE PROCEDURE im_2sp
        MODULE PROCEDURE zm_2sp
      END INTERFACE
      INTERFACE to_dp
        MODULE PROCEDURE fm_2dp
        MODULE PROCEDURE im_2dp
        MODULE PROCEDURE zm_2dp
      END INTERFACE
      INTERFACE to_spz
        MODULE PROCEDURE fm_2spz
        MODULE PROCEDURE im_2spz
        MODULE PROCEDURE zm_2spz
      END INTERFACE
      INTERFACE to_dpz
        MODULE PROCEDURE fm_2dpz
        MODULE PROCEDURE im_2dpz
        MODULE PROCEDURE zm_2dpz
      END INTERFACE
! ..
! .. Derived Type Declarations ..
      TYPE :: fm
        SEQUENCE
        REAL (kind(0.0D0)) :: mfm(0:lunpck)
      END TYPE fm
      TYPE :: im
        SEQUENCE
        REAL (kind(0.0D0)) :: mim(0:lunpck)
      END TYPE im
      TYPE :: zm
        SEQUENCE
        REAL (kind(0.0D0)) :: mzm(0:lunpkz)
      END TYPE zm
! ..
! .. Local Structures ..
      TYPE (fm), PRIVATE :: mtfm, mufm
      TYPE (im), PRIVATE :: mtim, muim
      TYPE (zm), PRIVATE :: mtzm, muzm
! ..
    CONTAINS


!                                                                   =

      SUBROUTINE fmeq_ifm(ival,ma)
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (INOUT) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmm2i
! ..
        CALL fmm2i(ma%mfm,ival)
      END SUBROUTINE fmeq_ifm

      SUBROUTINE fmeq_iim(ival,ma)
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (INOUT) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imm2i
! ..
        CALL imm2i(ma%mim,ival)
      END SUBROUTINE fmeq_iim

      SUBROUTINE fmeq_izm(ival,ma)
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (INOUT) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmm2i
! ..
        CALL zmm2i(ma%mzm,ival)
      END SUBROUTINE fmeq_izm

      SUBROUTINE fmeq_rfm(r,ma)
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (INOUT) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp
! ..
        CALL fmm2sp(ma%mfm,r)
      END SUBROUTINE fmeq_rfm

      SUBROUTINE fmeq_rim(r,ma)
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (INOUT) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmm2sp(mtfm%mfm,r)
      END SUBROUTINE fmeq_rim

      SUBROUTINE fmeq_rzm(r,ma)
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (INOUT) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2sp(mtfm%mfm,r)
      END SUBROUTINE fmeq_rzm

      SUBROUTINE fmeq_dfm(d,ma)
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (INOUT) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp
! ..
        CALL fmm2dp(ma%mfm,d)
      END SUBROUTINE fmeq_dfm

      SUBROUTINE fmeq_dim(d,ma)
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (INOUT) :: d
! ..
! .. External Subroutines ..
        EXTERNAL imm2dp
! ..
        CALL imm2dp(ma%mim,d)
      END SUBROUTINE fmeq_dim

      SUBROUTINE fmeq_dzm(d,ma)
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (INOUT) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,d)
      END SUBROUTINE fmeq_dzm

      SUBROUTINE fmeq_zfm(z,ma)
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (INOUT) :: z
! ..
! .. Local Scalars ..
        REAL :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp
! ..
        CALL fmm2sp(ma%mfm,r)
        z = cmplx(r,0.0)
      END SUBROUTINE fmeq_zfm

      SUBROUTINE fmeq_zim(z,ma)
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (INOUT) :: z
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL imm2dp
! ..
        CALL imm2dp(ma%mim,d)
        z = cmplx(real(d),0.0)
      END SUBROUTINE fmeq_zim

      SUBROUTINE fmeq_zzm(z,ma)
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (INOUT) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmm2z
! ..
        CALL zmm2z(ma%mzm,z)
      END SUBROUTINE fmeq_zzm

      SUBROUTINE fmeq_cfm(c,ma)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (INOUT) :: c
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp
! ..
        CALL fmm2dp(ma%mfm,d)
        c = cmplx(d,0.0D0,kind(0.0D0))
      END SUBROUTINE fmeq_cfm

      SUBROUTINE fmeq_cim(c,ma)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (INOUT) :: c
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL imm2dp
! ..
        CALL imm2dp(ma%mim,d)
        c = cmplx(d,0.0D0,kind(0.0D0))
      END SUBROUTINE fmeq_cim

      SUBROUTINE fmeq_czm(c,ma)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (INOUT) :: c
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d1, d2
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp, zmimag, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,d1)
        CALL zmimag(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,d2)
        c = cmplx(d1,d2,kind(0.0D0))
      END SUBROUTINE fmeq_czm

      SUBROUTINE fmeq_fmi(ma,ival)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,ma%mfm)
      END SUBROUTINE fmeq_fmi

      SUBROUTINE fmeq_fmr(ma,r)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,ma%mfm)
      END SUBROUTINE fmeq_fmr

      SUBROUTINE fmeq_fmd(ma,d)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,ma%mfm)
      END SUBROUTINE fmeq_fmd

      SUBROUTINE fmeq_fmz(ma,z)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        REAL :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        r = real(z)
        CALL fmsp2m(r,ma%mfm)
      END SUBROUTINE fmeq_fmz

      SUBROUTINE fmeq_fmc(ma,c)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        d = real(c,kind(0.0D0))
        CALL fmdp2m(d,ma%mfm)
      END SUBROUTINE fmeq_fmc

      SUBROUTINE fmeq_fmfm(ma,mb)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(mb%mfm,ma%mfm)
      END SUBROUTINE fmeq_fmfm

      SUBROUTINE fmeq_fmim(ma,mb)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        CALL imi2fm(mb%mim,ma%mfm)
      END SUBROUTINE fmeq_fmim

      SUBROUTINE fmeq_fmzm(ma,mb)
! .. Structure Arguments ..
        TYPE (fm), INTENT (INOUT) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(mb%mzm,ma%mfm)
      END SUBROUTINE fmeq_fmzm

      SUBROUTINE fmeq_imi(ma,ival)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,ma%mim)
      END SUBROUTINE fmeq_imi

      SUBROUTINE fmeq_imr(ma,r)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        ival = int(r)
        CALL imi2m(ival,ma%mim)
      END SUBROUTINE fmeq_imr

      SUBROUTINE fmeq_imd(ma,d)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        ival = int(d)
        CALL imi2m(ival,ma%mim)
      END SUBROUTINE fmeq_imd

      SUBROUTINE fmeq_imz(ma,z)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imfm2i
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imfm2i(mtfm%mfm,ma%mim)
      END SUBROUTINE fmeq_imz

      SUBROUTINE fmeq_imc(ma,c)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imfm2i
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imfm2i(mtfm%mfm,ma%mim)
      END SUBROUTINE fmeq_imc

      SUBROUTINE fmeq_imfm(ma,mb)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imfm2i
! ..
        CALL imfm2i(mb%mfm,ma%mim)
      END SUBROUTINE fmeq_imfm

      SUBROUTINE fmeq_imim(ma,mb)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(mb%mim,ma%mim)
      END SUBROUTINE fmeq_imim

      SUBROUTINE fmeq_imzm(ma,mb)
! .. Structure Arguments ..
        TYPE (im), INTENT (INOUT) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imfm2i, zmreal
! ..
        CALL zmreal(mb%mzm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,ma%mim)
      END SUBROUTINE fmeq_imzm

      SUBROUTINE fmeq_zmi(ma,ival)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmi2m
! ..
        CALL zmi2m(ival,ma%mzm)
      END SUBROUTINE fmeq_zmi

      SUBROUTINE fmeq_zmr(ma,r)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        COMPLEX :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmz2m
! ..
        z = cmplx(r,0.0)
        CALL zmz2m(z,ma%mzm)
      END SUBROUTINE fmeq_zmr

      SUBROUTINE fmeq_zmd(ma,d)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmdp2m(0.0D0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,ma%mzm)
      END SUBROUTINE fmeq_zmd

      SUBROUTINE fmeq_zmz(ma,z)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmz2m
! ..
        CALL zmz2m(z,ma%mzm)
      END SUBROUTINE fmeq_zmz

      SUBROUTINE fmeq_zmc(ma,c)
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx
! ..
        d = real(c,kind(0.0D0))
        CALL fmdp2m(d,mtfm%mfm)
        d = aimag(c)
        CALL fmdp2m(d,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,ma%mzm)
      END SUBROUTINE fmeq_zmc

      SUBROUTINE fmeq_zmfm(ma,mb)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,ma%mzm)
      END SUBROUTINE fmeq_zmfm

      SUBROUTINE fmeq_zmim(ma,mb)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,ma%mzm)
      END SUBROUTINE fmeq_zmim

      SUBROUTINE fmeq_zmzm(ma,mb)
! .. Structure Arguments ..
        TYPE (zm), INTENT (INOUT) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(mb%mzm,ma%mzm)
      END SUBROUTINE fmeq_zmzm

!  Reference:  The 39 Steps, John Buchan, 1915, Curtis Publishers.





!                                                                .EQ.

      FUNCTION fmleq_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmleq_ifm = fmcomp(mtfm%mfm,'EQ',ma%mfm)
      END FUNCTION fmleq_ifm

      FUNCTION fmleq_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmleq_iim = imcomp(mtim%mim,'EQ',ma%mim)
      END FUNCTION fmleq_iim

      FUNCTION fmleq_izm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmimag, zmreal
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_izm = l1 .AND. l2
      END FUNCTION fmleq_izm

      FUNCTION fmleq_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmleq_rfm = fmcomp(mtfm%mfm,'EQ',ma%mfm)
      END FUNCTION fmleq_rfm

      FUNCTION fmleq_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmleq_rim = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmleq_rim

      FUNCTION fmleq_rzm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_rzm = l1 .AND. l2
      END FUNCTION fmleq_rzm

      FUNCTION fmleq_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmleq_dfm = fmcomp(mtfm%mfm,'EQ',ma%mfm)
      END FUNCTION fmleq_dfm

      FUNCTION fmleq_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmleq_dim = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmleq_dim

      FUNCTION fmleq_dzm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmimag, zmreal
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_dzm = l1 .AND. l2
      END FUNCTION fmleq_dzm

      FUNCTION fmleq_zfm(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',ma%mfm)
        l2 = .TRUE.
        IF (aimag(z)/=0.0) l2 = .FALSE.
        fmleq_zfm = l1 .AND. l2
      END FUNCTION fmleq_zfm

      FUNCTION fmleq_zim(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        ndig = ndsave
        l2 = .TRUE.
        IF (aimag(z)/=0.0) l2 = .FALSE.
        fmleq_zim = l1 .AND. l2
      END FUNCTION fmleq_zim

      FUNCTION fmleq_zzm(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmsp2m(aimag(z),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zzm = l1 .AND. l2
      END FUNCTION fmleq_zzm

      FUNCTION fmleq_cfm(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',ma%mfm)
        l2 = .TRUE.
        IF (aimag(c)/=0.0) l2 = .FALSE.
        fmleq_cfm = l1 .AND. l2
      END FUNCTION fmleq_cfm

      FUNCTION fmleq_cim(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        ndig = ndsave
        l2 = .TRUE.
        IF (aimag(c)/=0.0) l2 = .FALSE.
        fmleq_cim = l1 .AND. l2
      END FUNCTION fmleq_cim

      FUNCTION fmleq_czm(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmleq_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmimag, zmreal
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmdp2m(aimag(c),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_czm = l1 .AND. l2
      END FUNCTION fmleq_czm

      FUNCTION fmleq_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmleq_fmi = fmcomp(ma%mfm,'EQ',mtfm%mfm)
      END FUNCTION fmleq_fmi

      FUNCTION fmleq_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmleq_fmr = fmcomp(ma%mfm,'EQ',mtfm%mfm)
      END FUNCTION fmleq_fmr

      FUNCTION fmleq_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmleq_fmd = fmcomp(ma%mfm,'EQ',mtfm%mfm)
      END FUNCTION fmleq_fmd

      FUNCTION fmleq_fmz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        l1 = fmcomp(ma%mfm,'EQ',mtfm%mfm)
        l2 = .TRUE.
        IF (aimag(z)/=0.0) l2 = .FALSE.
        fmleq_fmz = l1 .AND. l2
      END FUNCTION fmleq_fmz

      FUNCTION fmleq_fmc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        l1 = fmcomp(ma%mfm,'EQ',mtfm%mfm)
        l2 = .TRUE.
        IF (aimag(c)/=0.0) l2 = .FALSE.
        fmleq_fmc = l1 .AND. l2
      END FUNCTION fmleq_fmc

      FUNCTION fmleq_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmleq_fmfm = fmcomp(ma%mfm,'EQ',mb%mfm)
      END FUNCTION fmleq_fmfm

      FUNCTION fmleq_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm
! ..
        CALL fmint(ma%mfm,mtfm%mfm)
        IF (fmcomp(ma%mfm,'EQ',mtfm%mfm)) THEN
          CALL imi2fm(mb%mim,mtfm%mfm)
          fmleq_fmim = fmcomp(ma%mfm,'EQ',mtfm%mfm)
        ELSE
          fmleq_fmim = .FALSE.
        END IF
      END FUNCTION fmleq_fmim

      FUNCTION fmleq_fmzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(mb%mzm,mtfm%mfm)
        l1 = fmcomp(ma%mfm,'EQ',mtfm%mfm)
        l2 = .TRUE.
        IF (mb%mzm(kptimu+2)/=0) l2 = .FALSE.
        fmleq_fmzm = l1 .AND. l2
      END FUNCTION fmleq_fmzm

      FUNCTION fmleq_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmleq_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmleq_imi = imcomp(ma%mim,'EQ',mtim%mim)
      END FUNCTION fmleq_imi

      FUNCTION fmleq_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmleq_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmleq_imr = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmleq_imr

      FUNCTION fmleq_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmleq_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmleq_imd = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmleq_imd

      FUNCTION fmleq_imz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmleq_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ndig = ndsave
        l2 = .TRUE.
        IF (aimag(z)/=0.0) l2 = .FALSE.
        fmleq_imz = l1 .AND. l2
      END FUNCTION fmleq_imz

      FUNCTION fmleq_imc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmleq_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ndig = ndsave
        l2 = .TRUE.
        IF (aimag(c)/=0.0) l2 = .FALSE.
        fmleq_imc = l1 .AND. l2
      END FUNCTION fmleq_imc

      FUNCTION fmleq_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm
! ..
        CALL fmint(mb%mfm,mtfm%mfm)
        IF (fmcomp(mb%mfm,'EQ',mtfm%mfm)) THEN
          CALL imi2fm(ma%mim,mtfm%mfm)
          fmleq_imfm = fmcomp(mb%mfm,'EQ',mtfm%mfm)
        ELSE
          fmleq_imfm = .FALSE.
        END IF
      END FUNCTION fmleq_imfm

      FUNCTION fmleq_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmleq_imim = imcomp(ma%mim,'EQ',mb%mim)
      END FUNCTION fmleq_imim

      FUNCTION fmleq_imzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm, zmreal
! ..
        CALL zmreal(mb%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. mb%mzm(kptimu+2)==0) THEN
          CALL imi2fm(ma%mim,mufm%mfm)
          fmleq_imzm = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ELSE
          fmleq_imzm = .FALSE.
        END IF
      END FUNCTION fmleq_imzm

      FUNCTION fmleq_zmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmint, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. ma%mzm(kptimu+2)==0) THEN
          CALL fmi2m(ival,mufm%mfm)
          fmleq_zmi = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        ELSE
          fmleq_zmi = .FALSE.
        END IF
      END FUNCTION fmleq_zmi

      FUNCTION fmleq_zmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zmr = l1 .AND. l2
      END FUNCTION fmleq_zmr

      FUNCTION fmleq_zmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmimag, zmreal
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zmd = l1 .AND. l2
      END FUNCTION fmleq_zmd

      FUNCTION fmleq_zmz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmsp2m(aimag(z),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zmz = l1 .AND. l2
      END FUNCTION fmleq_zmz

      FUNCTION fmleq_zmc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmimag, zmreal
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL fmdp2m(aimag(c),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zmc = l1 .AND. l2
      END FUNCTION fmleq_zmc

      FUNCTION fmleq_zmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        l1 = fmcomp(mb%mfm,'EQ',mtfm%mfm)
        l2 = .TRUE.
        IF (ma%mzm(kptimu+2)/=0) l2 = .FALSE.
        fmleq_zmfm = l1 .AND. l2
      END FUNCTION fmleq_zmfm

      FUNCTION fmleq_zmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. ma%mzm(kptimu+2)==0) THEN
          CALL imi2fm(mb%mim,mufm%mfm)
          fmleq_zmim = fmcomp(mufm%mfm,'EQ',mtfm%mfm)
        ELSE
          fmleq_zmim = .FALSE.
        END IF
      END FUNCTION fmleq_zmim

      FUNCTION fmleq_zmzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmleq_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmimag, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL zmreal(mb%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        CALL zmimag(ma%mzm,mtfm%mfm)
        CALL zmimag(mb%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'EQ',mufm%mfm)
        fmleq_zmzm = l1 .AND. l2
      END FUNCTION fmleq_zmzm






!                                                                .NE.

      FUNCTION fmlne_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlne_ifm = fmcomp(mtfm%mfm,'NE',ma%mfm)
      END FUNCTION fmlne_ifm

      FUNCTION fmlne_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlne_iim = imcomp(mtim%mim,'NE',ma%mim)
      END FUNCTION fmlne_iim

      FUNCTION fmlne_izm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmimag, zmreal
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_izm = l1 .OR. l2
      END FUNCTION fmlne_izm

      FUNCTION fmlne_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlne_rfm = fmcomp(mtfm%mfm,'NE',ma%mfm)
      END FUNCTION fmlne_rfm

      FUNCTION fmlne_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlne_rim = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlne_rim

      FUNCTION fmlne_rzm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_rzm = l1 .OR. l2
      END FUNCTION fmlne_rzm

      FUNCTION fmlne_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlne_dfm = fmcomp(mtfm%mfm,'NE',ma%mfm)
      END FUNCTION fmlne_dfm

      FUNCTION fmlne_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlne_dim = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlne_dim

      FUNCTION fmlne_dzm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmimag, zmreal
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_dzm = l1 .OR. l2
      END FUNCTION fmlne_dzm

      FUNCTION fmlne_zfm(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',ma%mfm)
        l2 = .FALSE.
        IF (aimag(z)/=0.0) l2 = .TRUE.
        fmlne_zfm = l1 .OR. l2
      END FUNCTION fmlne_zfm

      FUNCTION fmlne_zim(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        ndig = ndsave
        l2 = .FALSE.
        IF (aimag(z)/=0.0) l2 = .TRUE.
        fmlne_zim = l1 .OR. l2
      END FUNCTION fmlne_zim

      FUNCTION fmlne_zzm(z,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmsp2m(aimag(z),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zzm = l1 .OR. l2
      END FUNCTION fmlne_zzm

      FUNCTION fmlne_cfm(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',ma%mfm)
        l2 = .FALSE.
        IF (aimag(c)/=0.0) l2 = .TRUE.
        fmlne_cfm = l1 .OR. l2
      END FUNCTION fmlne_cfm

      FUNCTION fmlne_cim(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        ndig = ndsave
        l2 = .FALSE.
        IF (aimag(c)/=0.0) l2 = .TRUE.
        fmlne_cim = l1 .OR. l2
      END FUNCTION fmlne_cim

      FUNCTION fmlne_czm(c,ma)
! .. Function Return Value ..
        LOGICAL :: fmlne_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmimag, zmreal
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmdp2m(aimag(c),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_czm = l1 .OR. l2
      END FUNCTION fmlne_czm

      FUNCTION fmlne_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlne_fmi = fmcomp(ma%mfm,'NE',mtfm%mfm)
      END FUNCTION fmlne_fmi

      FUNCTION fmlne_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlne_fmr = fmcomp(ma%mfm,'NE',mtfm%mfm)
      END FUNCTION fmlne_fmr

      FUNCTION fmlne_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlne_fmd = fmcomp(ma%mfm,'NE',mtfm%mfm)
      END FUNCTION fmlne_fmd

      FUNCTION fmlne_fmz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        l1 = fmcomp(ma%mfm,'NE',mtfm%mfm)
        l2 = .FALSE.
        IF (aimag(z)/=0.0) l2 = .TRUE.
        fmlne_fmz = l1 .OR. l2
      END FUNCTION fmlne_fmz

      FUNCTION fmlne_fmc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        l1 = fmcomp(ma%mfm,'NE',mtfm%mfm)
        l2 = .FALSE.
        IF (aimag(c)/=0.0) l2 = .TRUE.
        fmlne_fmc = l1 .OR. l2
      END FUNCTION fmlne_fmc

      FUNCTION fmlne_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmlne_fmfm = fmcomp(ma%mfm,'NE',mb%mfm)
      END FUNCTION fmlne_fmfm

      FUNCTION fmlne_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm
! ..
        CALL fmint(ma%mfm,mtfm%mfm)
        IF (fmcomp(ma%mfm,'EQ',mtfm%mfm)) THEN
          CALL imi2fm(mb%mim,mtfm%mfm)
          fmlne_fmim = fmcomp(ma%mfm,'NE',mtfm%mfm)
        ELSE
          fmlne_fmim = .TRUE.
        END IF
      END FUNCTION fmlne_fmim

      FUNCTION fmlne_fmzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(mb%mzm,mtfm%mfm)
        l1 = fmcomp(ma%mfm,'NE',mtfm%mfm)
        l2 = .FALSE.
        IF (mb%mzm(kptimu+2)/=0) l2 = .TRUE.
        fmlne_fmzm = l1 .OR. l2
      END FUNCTION fmlne_fmzm

      FUNCTION fmlne_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlne_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlne_imi = imcomp(ma%mim,'NE',mtim%mim)
      END FUNCTION fmlne_imi

      FUNCTION fmlne_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlne_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlne_imr = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlne_imr

      FUNCTION fmlne_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlne_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlne_imd = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlne_imd

      FUNCTION fmlne_imz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmlne_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ndig = ndsave
        l2 = .FALSE.
        IF (aimag(z)/=0.0) l2 = .TRUE.
        fmlne_imz = l1 .OR. l2
      END FUNCTION fmlne_imz

      FUNCTION fmlne_imc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmlne_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        l1 = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ndig = ndsave
        l2 = .FALSE.
        IF (aimag(c)/=0.0) l2 = .TRUE.
        fmlne_imc = l1 .OR. l2
      END FUNCTION fmlne_imc

      FUNCTION fmlne_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm
! ..
        CALL fmint(mb%mfm,mtfm%mfm)
        IF (fmcomp(mb%mfm,'EQ',mtfm%mfm)) THEN
          CALL imi2fm(ma%mim,mtfm%mfm)
          fmlne_imfm = fmcomp(mb%mfm,'NE',mtfm%mfm)
        ELSE
          fmlne_imfm = .TRUE.
        END IF
      END FUNCTION fmlne_imfm

      FUNCTION fmlne_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmlne_imim = imcomp(ma%mim,'NE',mb%mim)
      END FUNCTION fmlne_imim

      FUNCTION fmlne_imzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm, zmreal
! ..
        CALL zmreal(mb%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. mb%mzm(kptimu+2)==0) THEN
          CALL imi2fm(ma%mim,mufm%mfm)
          fmlne_imzm = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ELSE
          fmlne_imzm = .TRUE.
        END IF
      END FUNCTION fmlne_imzm

      FUNCTION fmlne_zmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmint, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. ma%mzm(kptimu+2)==0) THEN
          CALL fmi2m(ival,mufm%mfm)
          fmlne_zmi = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        ELSE
          fmlne_zmi = .TRUE.
        END IF
      END FUNCTION fmlne_zmi

      FUNCTION fmlne_zmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zmr = l1 .OR. l2
      END FUNCTION fmlne_zmr

      FUNCTION fmlne_zmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmimag, zmreal
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmi2m(0,mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zmd = l1 .OR. l2
      END FUNCTION fmlne_zmd

      FUNCTION fmlne_zmz(ma,z)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, zmimag, zmreal
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmsp2m(aimag(z),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zmz = l1 .OR. l2
      END FUNCTION fmlne_zmz

      FUNCTION fmlne_zmc(ma,c)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmimag, zmreal
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL zmreal(ma%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL fmdp2m(aimag(c),mtfm%mfm)
        CALL zmimag(ma%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zmc = l1 .OR. l2
      END FUNCTION fmlne_zmc

      FUNCTION fmlne_zmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        l1 = fmcomp(mb%mfm,'NE',mtfm%mfm)
        l2 = .FALSE.
        IF (ma%mzm(kptimu+2)/=0) l2 = .TRUE.
        fmlne_zmfm = l1 .OR. l2
      END FUNCTION fmlne_zmfm

      FUNCTION fmlne_zmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imi2fm, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mufm%mfm)
        IF (fmcomp(mufm%mfm,'EQ',mtfm%mfm) .AND. ma%mzm(kptimu+2)==0) THEN
          CALL imi2fm(mb%mim,mufm%mfm)
          fmlne_zmim = fmcomp(mufm%mfm,'NE',mtfm%mfm)
        ELSE
          fmlne_zmim = .TRUE.
        END IF
      END FUNCTION fmlne_zmim

      FUNCTION fmlne_zmzm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlne_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. Local Scalars ..
        LOGICAL :: l1, l2
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL zmimag, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL zmreal(mb%mzm,mufm%mfm)
        l1 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        CALL zmimag(ma%mzm,mtfm%mfm)
        CALL zmimag(mb%mzm,mufm%mfm)
        l2 = fmcomp(mtfm%mfm,'NE',mufm%mfm)
        fmlne_zmzm = l1 .OR. l2
      END FUNCTION fmlne_zmzm





!                                                                .GT.

      FUNCTION fmlgt_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlgt_ifm = fmcomp(mtfm%mfm,'GT',ma%mfm)
      END FUNCTION fmlgt_ifm

      FUNCTION fmlgt_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlgt_iim = imcomp(mtim%mim,'GT',ma%mim)
      END FUNCTION fmlgt_iim

      FUNCTION fmlgt_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlgt_rfm = fmcomp(mtfm%mfm,'GT',ma%mfm)
      END FUNCTION fmlgt_rfm

      FUNCTION fmlgt_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlgt_rim = fmcomp(mtfm%mfm,'GT',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_rim

      FUNCTION fmlgt_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlgt_dfm = fmcomp(mtfm%mfm,'GT',ma%mfm)
      END FUNCTION fmlgt_dfm

      FUNCTION fmlgt_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlgt_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlgt_dim = fmcomp(mtfm%mfm,'GT',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_dim

      FUNCTION fmlgt_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlgt_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlgt_fmi = fmcomp(ma%mfm,'GT',mtfm%mfm)
      END FUNCTION fmlgt_fmi

      FUNCTION fmlgt_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlgt_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlgt_fmr = fmcomp(ma%mfm,'GT',mtfm%mfm)
      END FUNCTION fmlgt_fmr

      FUNCTION fmlgt_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlgt_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlgt_fmd = fmcomp(ma%mfm,'GT',mtfm%mfm)
      END FUNCTION fmlgt_fmd

      FUNCTION fmlgt_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlgt_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmlgt_fmfm = fmcomp(ma%mfm,'GT',mb%mfm)
      END FUNCTION fmlgt_fmfm

      FUNCTION fmlgt_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlgt_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = mb%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(mb%mim,mtfm%mfm)
        fmlgt_fmim = fmcomp(ma%mfm,'GT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_fmim

      FUNCTION fmlgt_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlgt_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlgt_imi = imcomp(ma%mim,'GT',mtim%mim)
      END FUNCTION fmlgt_imi

      FUNCTION fmlgt_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlgt_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlgt_imr = fmcomp(mufm%mfm,'GT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_imr

      FUNCTION fmlgt_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlgt_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlgt_imd = fmcomp(mufm%mfm,'GT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_imd

      FUNCTION fmlgt_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlgt_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(ma%mim,mtfm%mfm)
        fmlgt_imfm = fmcomp(mtfm%mfm,'GT',mb%mfm)
        ndig = ndsave
      END FUNCTION fmlgt_imfm

      FUNCTION fmlgt_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlgt_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmlgt_imim = imcomp(ma%mim,'GT',mb%mim)
      END FUNCTION fmlgt_imim




!                                                                .GE.

      FUNCTION fmlge_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlge_ifm = fmcomp(mtfm%mfm,'GE',ma%mfm)
      END FUNCTION fmlge_ifm

      FUNCTION fmlge_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlge_iim = imcomp(mtim%mim,'GE',ma%mim)
      END FUNCTION fmlge_iim

      FUNCTION fmlge_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlge_rfm = fmcomp(mtfm%mfm,'GE',ma%mfm)
      END FUNCTION fmlge_rfm

      FUNCTION fmlge_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlge_rim = fmcomp(mtfm%mfm,'GE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlge_rim

      FUNCTION fmlge_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlge_dfm = fmcomp(mtfm%mfm,'GE',ma%mfm)
      END FUNCTION fmlge_dfm

      FUNCTION fmlge_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlge_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlge_dim = fmcomp(mtfm%mfm,'GE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlge_dim

      FUNCTION fmlge_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlge_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlge_fmi = fmcomp(ma%mfm,'GE',mtfm%mfm)
      END FUNCTION fmlge_fmi

      FUNCTION fmlge_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlge_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlge_fmr = fmcomp(ma%mfm,'GE',mtfm%mfm)
      END FUNCTION fmlge_fmr

      FUNCTION fmlge_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlge_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlge_fmd = fmcomp(ma%mfm,'GE',mtfm%mfm)
      END FUNCTION fmlge_fmd

      FUNCTION fmlge_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlge_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmlge_fmfm = fmcomp(ma%mfm,'GE',mb%mfm)
      END FUNCTION fmlge_fmfm

      FUNCTION fmlge_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlge_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = mb%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(mb%mim,mtfm%mfm)
        fmlge_fmim = fmcomp(ma%mfm,'GE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlge_fmim

      FUNCTION fmlge_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlge_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlge_imi = imcomp(ma%mim,'GE',mtim%mim)
      END FUNCTION fmlge_imi

      FUNCTION fmlge_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlge_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlge_imr = fmcomp(mufm%mfm,'GE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlge_imr

      FUNCTION fmlge_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlge_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlge_imd = fmcomp(mufm%mfm,'GE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlge_imd

      FUNCTION fmlge_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlge_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(ma%mim,mtfm%mfm)
        fmlge_imfm = fmcomp(mtfm%mfm,'GE',mb%mfm)
        ndig = ndsave
      END FUNCTION fmlge_imfm

      FUNCTION fmlge_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlge_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmlge_imim = imcomp(ma%mim,'GE',mb%mim)
      END FUNCTION fmlge_imim





!                                                                .LT.

      FUNCTION fmllt_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmllt_ifm = fmcomp(mtfm%mfm,'LT',ma%mfm)
      END FUNCTION fmllt_ifm

      FUNCTION fmllt_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmllt_iim = imcomp(mtim%mim,'LT',ma%mim)
      END FUNCTION fmllt_iim

      FUNCTION fmllt_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmllt_rfm = fmcomp(mtfm%mfm,'LT',ma%mfm)
      END FUNCTION fmllt_rfm

      FUNCTION fmllt_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmllt_rim = fmcomp(mtfm%mfm,'LT',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmllt_rim

      FUNCTION fmllt_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmllt_dfm = fmcomp(mtfm%mfm,'LT',ma%mfm)
      END FUNCTION fmllt_dfm

      FUNCTION fmllt_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmllt_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmllt_dim = fmcomp(mtfm%mfm,'LT',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmllt_dim

      FUNCTION fmllt_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmllt_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmllt_fmi = fmcomp(ma%mfm,'LT',mtfm%mfm)
      END FUNCTION fmllt_fmi

      FUNCTION fmllt_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmllt_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmllt_fmr = fmcomp(ma%mfm,'LT',mtfm%mfm)
      END FUNCTION fmllt_fmr

      FUNCTION fmllt_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmllt_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmllt_fmd = fmcomp(ma%mfm,'LT',mtfm%mfm)
      END FUNCTION fmllt_fmd

      FUNCTION fmllt_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmllt_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmllt_fmfm = fmcomp(ma%mfm,'LT',mb%mfm)
      END FUNCTION fmllt_fmfm

      FUNCTION fmllt_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmllt_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = mb%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(mb%mim,mtfm%mfm)
        fmllt_fmim = fmcomp(ma%mfm,'LT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmllt_fmim

      FUNCTION fmllt_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmllt_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmllt_imi = imcomp(ma%mim,'LT',mtim%mim)
      END FUNCTION fmllt_imi

      FUNCTION fmllt_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmllt_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmllt_imr = fmcomp(mufm%mfm,'LT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmllt_imr

      FUNCTION fmllt_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmllt_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmllt_imd = fmcomp(mufm%mfm,'LT',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmllt_imd

      FUNCTION fmllt_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmllt_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(ma%mim,mtfm%mfm)
        fmllt_imfm = fmcomp(mtfm%mfm,'LT',mb%mfm)
        ndig = ndsave
      END FUNCTION fmllt_imfm

      FUNCTION fmllt_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmllt_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmllt_imim = imcomp(ma%mim,'LT',mb%mim)
      END FUNCTION fmllt_imim





!                                                                .LE.

      FUNCTION fmlle_ifm(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlle_ifm = fmcomp(mtfm%mfm,'LE',ma%mfm)
      END FUNCTION fmlle_ifm

      FUNCTION fmlle_iim(ival,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlle_iim = imcomp(mtim%mim,'LE',ma%mim)
      END FUNCTION fmlle_iim

      FUNCTION fmlle_rfm(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlle_rfm = fmcomp(mtfm%mfm,'LE',ma%mfm)
      END FUNCTION fmlle_rfm

      FUNCTION fmlle_rim(r,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlle_rim = fmcomp(mtfm%mfm,'LE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlle_rim

      FUNCTION fmlle_dfm(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlle_dfm = fmcomp(mtfm%mfm,'LE',ma%mfm)
      END FUNCTION fmlle_dfm

      FUNCTION fmlle_dim(d,ma)
! .. Function Return Value ..
        LOGICAL :: fmlle_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlle_dim = fmcomp(mtfm%mfm,'LE',mufm%mfm)
        ndig = ndsave
      END FUNCTION fmlle_dim

      FUNCTION fmlle_fmi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlle_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        fmlle_fmi = fmcomp(ma%mfm,'LE',mtfm%mfm)
      END FUNCTION fmlle_fmi

      FUNCTION fmlle_fmr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlle_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        fmlle_fmr = fmcomp(ma%mfm,'LE',mtfm%mfm)
      END FUNCTION fmlle_fmr

      FUNCTION fmlle_fmd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlle_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        fmlle_fmd = fmcomp(ma%mfm,'LE',mtfm%mfm)
      END FUNCTION fmlle_fmd

      FUNCTION fmlle_fmfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlle_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fmlle_fmfm = fmcomp(ma%mfm,'LE',mb%mfm)
      END FUNCTION fmlle_fmfm

      FUNCTION fmlle_fmim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlle_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = mb%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(mb%mim,mtfm%mfm)
        fmlle_fmim = fmcomp(ma%mfm,'LE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlle_fmim

      FUNCTION fmlle_imi(ma,ival)
! .. Function Return Value ..
        LOGICAL :: fmlle_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        fmlle_imi = imcomp(ma%mim,'LE',mtim%mim)
      END FUNCTION fmlle_imi

      FUNCTION fmlle_imr(ma,r)
! .. Function Return Value ..
        LOGICAL :: fmlle_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlle_imr = fmcomp(mufm%mfm,'LE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlle_imr

      FUNCTION fmlle_imd(ma,d)
! .. Function Return Value ..
        LOGICAL :: fmlle_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        fmlle_imd = fmcomp(mufm%mfm,'LE',mtfm%mfm)
        ndig = ndsave
      END FUNCTION fmlle_imd

      FUNCTION fmlle_imfm(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlle_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. Local Scalars ..
        INTEGER :: ka, ndsave
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        ndsave = ndig
        ka = ma%mim(1)
        ndig = max(ka+ngrd52,ndig)
        CALL imi2fm(ma%mim,mtfm%mfm)
        fmlle_imfm = fmcomp(mtfm%mfm,'LE',mb%mfm)
        ndig = ndsave
      END FUNCTION fmlle_imfm

      FUNCTION fmlle_imim(ma,mb)
! .. Function Return Value ..
        LOGICAL :: fmlle_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        fmlle_imim = imcomp(ma%mim,'LE',mb%mim)
      END FUNCTION fmlle_imim





!                                                                   +

      FUNCTION fmadd_ifm(ival,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmadd(mtfm%mfm,ma%mfm,fmadd_ifm%mfm)
      END FUNCTION fmadd_ifm

      FUNCTION fmadd_iim(ival,ma)
! .. Function Return Value ..
        TYPE (im) :: fmadd_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imadd, imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        CALL imadd(mtim%mim,ma%mim,fmadd_iim%mim)
      END FUNCTION fmadd_iim

      FUNCTION fmadd_izm(ival,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(mtzm%mzm,ma%mzm,fmadd_izm%mzm)
      END FUNCTION fmadd_izm

      FUNCTION fmadd_rfm(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmadd(mtfm%mfm,ma%mfm,fmadd_rfm%mfm)
      END FUNCTION fmadd_rfm

      FUNCTION fmadd_rim(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmadd(mtfm%mfm,mufm%mfm,fmadd_rim%mfm)
      END FUNCTION fmadd_rim

      FUNCTION fmadd_rzm(r,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmadd, zmcmpx
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(mtzm%mzm,ma%mzm,fmadd_rzm%mzm)
      END FUNCTION fmadd_rzm

      FUNCTION fmadd_dfm(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmadd(mtfm%mfm,ma%mfm,fmadd_dfm%mfm)
      END FUNCTION fmadd_dfm

      FUNCTION fmadd_dim(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmdp2m, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmadd(mtfm%mfm,mufm%mfm,fmadd_dim%mfm)
      END FUNCTION fmadd_dim

      FUNCTION fmadd_dzm(d,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(mtzm%mzm,ma%mzm,fmadd_dzm%mzm)
      END FUNCTION fmadd_dzm

      FUNCTION fmadd_zfm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(mtzm%mzm,muzm%mzm,fmadd_zfm%mzm)
      END FUNCTION fmadd_zfm

      FUNCTION fmadd_zim(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmadd, zmcmpx, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(mtzm%mzm,muzm%mzm,fmadd_zim%mzm)
      END FUNCTION fmadd_zim

      FUNCTION fmadd_zzm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmadd, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmadd(mtzm%mzm,ma%mzm,fmadd_zzm%mzm)
      END FUNCTION fmadd_zzm

      FUNCTION fmadd_cfm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(mtzm%mzm,muzm%mzm,fmadd_cfm%mzm)
      END FUNCTION fmadd_cfm

      FUNCTION fmadd_cim(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(mtzm%mzm,muzm%mzm,fmadd_cim%mzm)
      END FUNCTION fmadd_cim

      FUNCTION fmadd_czm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(mtzm%mzm,ma%mzm,fmadd_czm%mzm)
      END FUNCTION fmadd_czm

      FUNCTION fmadd_fmi(ma,ival)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmadd(ma%mfm,mtfm%mfm,fmadd_fmi%mfm)
      END FUNCTION fmadd_fmi

      FUNCTION fmadd_fmr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmadd(ma%mfm,mtfm%mfm,fmadd_fmr%mfm)
      END FUNCTION fmadd_fmr

      FUNCTION fmadd_fmd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmadd(ma%mfm,mtfm%mfm,fmadd_fmd%mfm)
      END FUNCTION fmadd_fmd

      FUNCTION fmadd_fmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(muzm%mzm,mtzm%mzm,fmadd_fmz%mzm)
      END FUNCTION fmadd_fmz

      FUNCTION fmadd_fmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(muzm%mzm,mtzm%mzm,fmadd_fmc%mzm)
      END FUNCTION fmadd_fmc

      FUNCTION fmadd_fmfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmadd
! ..
        CALL fmadd(ma%mfm,mb%mfm,fmadd_fmfm%mfm)
      END FUNCTION fmadd_fmfm

      FUNCTION fmadd_fmim(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, imi2fm
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmadd(ma%mfm,mtfm%mfm,fmadd_fmim%mfm)
      END FUNCTION fmadd_fmim

      FUNCTION fmadd_fmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(ma%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmadd(mtzm%mzm,mb%mzm,fmadd_fmzm%mzm)
      END FUNCTION fmadd_fmzm

      FUNCTION fmadd_imi(ma,ival)
! .. Function Return Value ..
        TYPE (im) :: fmadd_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imadd, imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        CALL imadd(ma%mim,mtim%mim,fmadd_imi%mim)
      END FUNCTION fmadd_imi

      FUNCTION fmadd_imr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmadd(mufm%mfm,mtfm%mfm,fmadd_imr%mfm)
      END FUNCTION fmadd_imr

      FUNCTION fmadd_imd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmdp2m, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmadd(mufm%mfm,mtfm%mfm,fmadd_imd%mfm)
      END FUNCTION fmadd_imd

      FUNCTION fmadd_imz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmadd, zmcmpx, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(muzm%mzm,mtzm%mzm,fmadd_imz%mzm)
      END FUNCTION fmadd_imz

      FUNCTION fmadd_imc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(muzm%mzm,mtzm%mzm,fmadd_imc%mzm)
      END FUNCTION fmadd_imc

      FUNCTION fmadd_imfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmadd(mtfm%mfm,mb%mfm,fmadd_imfm%mfm)
      END FUNCTION fmadd_imfm

      FUNCTION fmadd_imim(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmadd_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imadd
! ..
        CALL imadd(ma%mim,mb%mim,fmadd_imim%mim)
      END FUNCTION fmadd_imim

      FUNCTION fmadd_imzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmadd, zmcmpx
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(muzm%mzm,mb%mzm,fmadd_imzm%mzm)
      END FUNCTION fmadd_imzm

      FUNCTION fmadd_zmi(ma,ival)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmi%mzm)
      END FUNCTION fmadd_zmi

      FUNCTION fmadd_zmr(ma,r)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmadd, zmcmpx
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmr%mzm)
      END FUNCTION fmadd_zmr

      FUNCTION fmadd_zmd(ma,d)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmd%mzm)
      END FUNCTION fmadd_zmd

      FUNCTION fmadd_zmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmadd, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmz%mzm)
      END FUNCTION fmadd_zmz

      FUNCTION fmadd_zmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmadd, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmc%mzm)
      END FUNCTION fmadd_zmc

      FUNCTION fmadd_zmfm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmadd, zmcmpx
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmadd(ma%mzm,mtzm%mzm,fmadd_zmfm%mzm)
      END FUNCTION fmadd_zmfm

      FUNCTION fmadd_zmim(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmadd, zmcmpx
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmadd(ma%mzm,muzm%mzm,fmadd_zmim%mzm)
      END FUNCTION fmadd_zmim

      FUNCTION fmadd_zmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmadd
! ..
        CALL zmadd(ma%mzm,mb%mzm,fmadd_zmzm%mzm)
      END FUNCTION fmadd_zmzm

      FUNCTION fmadd_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmadd_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,fmadd_fm%mfm)
      END FUNCTION fmadd_fm

      FUNCTION fmadd_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmadd_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,fmadd_im%mim)
      END FUNCTION fmadd_im

      FUNCTION fmadd_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmadd_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(ma%mzm,fmadd_zm%mzm)
      END FUNCTION fmadd_zm





!                                                                   -

      FUNCTION fmsub_ifm(ival,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsub
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmsub(mtfm%mfm,ma%mfm,fmsub_ifm%mfm)
      END FUNCTION fmsub_ifm

      FUNCTION fmsub_iim(ival,ma)
! .. Function Return Value ..
        TYPE (im) :: fmsub_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m, imsub
! ..
        CALL imi2m(ival,mtim%mim)
        CALL imsub(mtim%mim,ma%mim,fmsub_iim%mim)
      END FUNCTION fmsub_iim

      FUNCTION fmsub_izm(ival,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(mtzm%mzm,ma%mzm,fmsub_izm%mzm)
      END FUNCTION fmsub_izm

      FUNCTION fmsub_rfm(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, fmsub
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmsub(mtfm%mfm,ma%mfm,fmsub_rfm%mfm)
      END FUNCTION fmsub_rfm

      FUNCTION fmsub_rim(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, fmsub, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmsub(mtfm%mfm,mufm%mfm,fmsub_rim%mfm)
      END FUNCTION fmsub_rim

      FUNCTION fmsub_rzm(r,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmsub
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(mtzm%mzm,ma%mzm,fmsub_rzm%mzm)
      END FUNCTION fmsub_rzm

      FUNCTION fmsub_dfm(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmsub
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmsub(mtfm%mfm,ma%mfm,fmsub_dfm%mfm)
      END FUNCTION fmsub_dfm

      FUNCTION fmsub_dim(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmsub, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmsub(mtfm%mfm,mufm%mfm,fmsub_dim%mfm)
      END FUNCTION fmsub_dim

      FUNCTION fmsub_dzm(d,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(mtzm%mzm,ma%mzm,fmsub_dzm%mzm)
      END FUNCTION fmsub_dzm

      FUNCTION fmsub_zfm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(mtzm%mzm,muzm%mzm,fmsub_zfm%mzm)
      END FUNCTION fmsub_zfm

      FUNCTION fmsub_zim(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(mtzm%mzm,muzm%mzm,fmsub_zim%mzm)
      END FUNCTION fmsub_zim

      FUNCTION fmsub_zzm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmsub(mtzm%mzm,ma%mzm,fmsub_zzm%mzm)
      END FUNCTION fmsub_zzm

      FUNCTION fmsub_cfm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(mtzm%mzm,muzm%mzm,fmsub_cfm%mzm)
      END FUNCTION fmsub_cfm

      FUNCTION fmsub_cim(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(mtzm%mzm,muzm%mzm,fmsub_cim%mzm)
      END FUNCTION fmsub_cim

      FUNCTION fmsub_czm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(mtzm%mzm,ma%mzm,fmsub_czm%mzm)
      END FUNCTION fmsub_czm

      FUNCTION fmsub_fmi(ma,ival)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsub
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,fmsub_fmi%mfm)
      END FUNCTION fmsub_fmi

      FUNCTION fmsub_fmr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, fmsub
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,fmsub_fmr%mfm)
      END FUNCTION fmsub_fmr

      FUNCTION fmsub_fmd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmsub
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,fmsub_fmd%mfm)
      END FUNCTION fmsub_fmd

      FUNCTION fmsub_fmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(muzm%mzm,mtzm%mzm,fmsub_fmz%mzm)
      END FUNCTION fmsub_fmz

      FUNCTION fmsub_fmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(muzm%mzm,mtzm%mzm,fmsub_fmc%mzm)
      END FUNCTION fmsub_fmc

      FUNCTION fmsub_fmfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsub
! ..
        CALL fmsub(ma%mfm,mb%mfm,fmsub_fmfm%mfm)
      END FUNCTION fmsub_fmfm

      FUNCTION fmsub_fmim(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsub, imi2fm
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,fmsub_fmim%mfm)
      END FUNCTION fmsub_fmim

      FUNCTION fmsub_fmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(ma%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmsub(mtzm%mzm,mb%mzm,fmsub_fmzm%mzm)
      END FUNCTION fmsub_fmzm

      FUNCTION fmsub_imi(ma,ival)
! .. Function Return Value ..
        TYPE (im) :: fmsub_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m, imsub
! ..
        CALL imi2m(ival,mtim%mim)
        CALL imsub(ma%mim,mtim%mim,fmsub_imi%mim)
      END FUNCTION fmsub_imi

      FUNCTION fmsub_imr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, fmsub, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmsub(mufm%mfm,mtfm%mfm,fmsub_imr%mfm)
      END FUNCTION fmsub_imr

      FUNCTION fmsub_imd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmsub, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmsub(mufm%mfm,mtfm%mfm,fmsub_imd%mfm)
      END FUNCTION fmsub_imd

      FUNCTION fmsub_imz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(muzm%mzm,mtzm%mzm,fmsub_imz%mzm)
      END FUNCTION fmsub_imz

      FUNCTION fmsub_imc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(muzm%mzm,mtzm%mzm,fmsub_imc%mzm)
      END FUNCTION fmsub_imc

      FUNCTION fmsub_imfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsub, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmsub(mtfm%mfm,mb%mfm,fmsub_imfm%mfm)
      END FUNCTION fmsub_imfm

      FUNCTION fmsub_imim(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmsub_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imsub
! ..
        CALL imsub(ma%mim,mb%mim,fmsub_imim%mim)
      END FUNCTION fmsub_imim

      FUNCTION fmsub_imzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmsub
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(muzm%mzm,mb%mzm,fmsub_imzm%mzm)
      END FUNCTION fmsub_imzm

      FUNCTION fmsub_zmi(ma,ival)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmi%mzm)
      END FUNCTION fmsub_zmi

      FUNCTION fmsub_zmr(ma,r)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmsub
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmr%mzm)
      END FUNCTION fmsub_zmr

      FUNCTION fmsub_zmd(ma,d)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmd%mzm)
      END FUNCTION fmsub_zmd

      FUNCTION fmsub_zmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmsub, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmz%mzm)
      END FUNCTION fmsub_zmz

      FUNCTION fmsub_zmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmsub
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmc%mzm)
      END FUNCTION fmsub_zmc

      FUNCTION fmsub_zmfm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmsub
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmsub(ma%mzm,mtzm%mzm,fmsub_zmfm%mzm)
      END FUNCTION fmsub_zmfm

      FUNCTION fmsub_zmim(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmsub
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmsub(ma%mzm,muzm%mzm,fmsub_zmim%mzm)
      END FUNCTION fmsub_zmim

      FUNCTION fmsub_zmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmsub
! ..
        CALL zmsub(ma%mzm,mb%mzm,fmsub_zmzm%mzm)
      END FUNCTION fmsub_zmzm

      FUNCTION fmsub_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsub_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,mtfm%mfm)
        IF (mtfm%mfm(1)/=munkno) mtfm%mfm(2) = -mtfm%mfm(2)
        CALL fmeq(mtfm%mfm,fmsub_fm%mfm)
      END FUNCTION fmsub_fm

      FUNCTION fmsub_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmsub_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,mtim%mim)
        IF (mtim%mim(1)/=munkno) mtim%mim(2) = -mtim%mim(2)
        CALL imeq(mtim%mim,fmsub_im%mim)
      END FUNCTION fmsub_im

      FUNCTION fmsub_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsub_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(ma%mzm,mtzm%mzm)
        IF (mtzm%mzm(1)/=munkno) mtzm%mzm(2) = -mtzm%mzm(2)
        IF (mtzm%mzm(kptimu+1)/=munkno) THEN
          mtzm%mzm(kptimu+2) = -mtzm%mzm(kptimu+2)
        END IF
        CALL zmeq(mtzm%mzm,fmsub_zm%mzm)
      END FUNCTION fmsub_zm





!                                                                   *

      FUNCTION fmmpy_ifm(ival,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmmpyi
! ..
        CALL fmmpyi(ma%mfm,ival,fmmpy_ifm%mfm)
      END FUNCTION fmmpy_ifm

      FUNCTION fmmpy_iim(ival,ma)
! .. Function Return Value ..
        TYPE (im) :: fmmpy_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL immpyi
! ..
        CALL immpyi(ma%mim,ival,fmmpy_iim%mim)
      END FUNCTION fmmpy_iim

      FUNCTION fmmpy_izm(ival,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmmpyi
! ..
        CALL zmmpyi(ma%mzm,ival,fmmpy_izm%mzm)
      END FUNCTION fmmpy_izm

      FUNCTION fmmpy_rfm(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmmpy(mtfm%mfm,ma%mfm,fmmpy_rfm%mfm)
      END FUNCTION fmmpy_rfm

      FUNCTION fmmpy_rim(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmmpy(mtfm%mfm,mufm%mfm,fmmpy_rim%mfm)
      END FUNCTION fmmpy_rim

      FUNCTION fmmpy_rzm(r,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmmpy
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(mtzm%mzm,ma%mzm,fmmpy_rzm%mzm)
      END FUNCTION fmmpy_rzm

      FUNCTION fmmpy_dfm(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmmpy
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmmpy(mtfm%mfm,ma%mfm,fmmpy_dfm%mfm)
      END FUNCTION fmmpy_dfm

      FUNCTION fmmpy_dim(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmmpy, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmmpy(mtfm%mfm,mufm%mfm,fmmpy_dim%mfm)
      END FUNCTION fmmpy_dim

      FUNCTION fmmpy_dzm(d,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(mtzm%mzm,ma%mzm,fmmpy_dzm%mzm)
      END FUNCTION fmmpy_dzm

      FUNCTION fmmpy_zfm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(mtzm%mzm,muzm%mzm,fmmpy_zfm%mzm)
      END FUNCTION fmmpy_zfm

      FUNCTION fmmpy_zim(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(mtzm%mzm,muzm%mzm,fmmpy_zim%mzm)
      END FUNCTION fmmpy_zim

      FUNCTION fmmpy_zzm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmmpy(mtzm%mzm,ma%mzm,fmmpy_zzm%mzm)
      END FUNCTION fmmpy_zzm

      FUNCTION fmmpy_cfm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(mtzm%mzm,muzm%mzm,fmmpy_cfm%mzm)
      END FUNCTION fmmpy_cfm

      FUNCTION fmmpy_cim(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(mtzm%mzm,muzm%mzm,fmmpy_cim%mzm)
      END FUNCTION fmmpy_cim

      FUNCTION fmmpy_czm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(mtzm%mzm,ma%mzm,fmmpy_czm%mzm)
      END FUNCTION fmmpy_czm

      FUNCTION fmmpy_fmi(ma,ival)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmmpyi
! ..
        CALL fmmpyi(ma%mfm,ival,fmmpy_fmi%mfm)
      END FUNCTION fmmpy_fmi

      FUNCTION fmmpy_fmr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmmpy(ma%mfm,mtfm%mfm,fmmpy_fmr%mfm)
      END FUNCTION fmmpy_fmr

      FUNCTION fmmpy_fmd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmmpy
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmmpy(ma%mfm,mtfm%mfm,fmmpy_fmd%mfm)
      END FUNCTION fmmpy_fmd

      FUNCTION fmmpy_fmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(muzm%mzm,mtzm%mzm,fmmpy_fmz%mzm)
      END FUNCTION fmmpy_fmz

      FUNCTION fmmpy_fmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(muzm%mzm,mtzm%mzm,fmmpy_fmc%mzm)
      END FUNCTION fmmpy_fmc

      FUNCTION fmmpy_fmfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy
! ..
        CALL fmmpy(ma%mfm,mb%mfm,fmmpy_fmfm%mfm)
      END FUNCTION fmmpy_fmfm

      FUNCTION fmmpy_fmim(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, imi2fm
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmmpy(ma%mfm,mtfm%mfm,fmmpy_fmim%mfm)
      END FUNCTION fmmpy_fmim

      FUNCTION fmmpy_fmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmmpy
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(ma%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmmpy(mtzm%mzm,mb%mzm,fmmpy_fmzm%mzm)
      END FUNCTION fmmpy_fmzm

      FUNCTION fmmpy_imi(ma,ival)
! .. Function Return Value ..
        TYPE (im) :: fmmpy_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL immpyi
! ..
        CALL immpyi(ma%mim,ival,fmmpy_imi%mim)
      END FUNCTION fmmpy_imi

      FUNCTION fmmpy_imr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmmpy(mufm%mfm,mtfm%mfm,fmmpy_imr%mfm)
      END FUNCTION fmmpy_imr

      FUNCTION fmmpy_imd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmmpy, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmmpy(mufm%mfm,mtfm%mfm,fmmpy_imd%mfm)
      END FUNCTION fmmpy_imd

      FUNCTION fmmpy_imz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(muzm%mzm,mtzm%mzm,fmmpy_imz%mzm)
      END FUNCTION fmmpy_imz

      FUNCTION fmmpy_imc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(muzm%mzm,mtzm%mzm,fmmpy_imc%mzm)
      END FUNCTION fmmpy_imc

      FUNCTION fmmpy_imfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmmpy_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmmpy(mtfm%mfm,mb%mfm,fmmpy_imfm%mfm)
      END FUNCTION fmmpy_imfm

      FUNCTION fmmpy_imim(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmmpy_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL immpy
! ..
        CALL immpy(ma%mim,mb%mim,fmmpy_imim%mim)
      END FUNCTION fmmpy_imim

      FUNCTION fmmpy_imzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmmpy
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(muzm%mzm,mb%mzm,fmmpy_imzm%mzm)
      END FUNCTION fmmpy_imzm

      FUNCTION fmmpy_zmi(ma,ival)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmmpyi
! ..
        CALL zmmpyi(ma%mzm,ival,fmmpy_zmi%mzm)
      END FUNCTION fmmpy_zmi

      FUNCTION fmmpy_zmr(ma,r)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmmpy
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(ma%mzm,mtzm%mzm,fmmpy_zmr%mzm)
      END FUNCTION fmmpy_zmr

      FUNCTION fmmpy_zmd(ma,d)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(ma%mzm,mtzm%mzm,fmmpy_zmd%mzm)
      END FUNCTION fmmpy_zmd

      FUNCTION fmmpy_zmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmmpy, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmmpy(ma%mzm,mtzm%mzm,fmmpy_zmz%mzm)
      END FUNCTION fmmpy_zmz

      FUNCTION fmmpy_zmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmmpy
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmmpy(ma%mzm,mtzm%mzm,fmmpy_zmc%mzm)
      END FUNCTION fmmpy_zmc

      FUNCTION fmmpy_zmfm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmmpy
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmmpy(ma%mzm,mtzm%mzm,fmmpy_zmfm%mzm)
      END FUNCTION fmmpy_zmfm

      FUNCTION fmmpy_zmim(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmmpy
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmmpy(ma%mzm,muzm%mzm,fmmpy_zmim%mzm)
      END FUNCTION fmmpy_zmim

      FUNCTION fmmpy_zmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmmpy_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmmpy
! ..
        CALL zmmpy(ma%mzm,mb%mzm,fmmpy_zmzm%mzm)
      END FUNCTION fmmpy_zmzm





!                                                                   /

      FUNCTION fmdiv_ifm(ival,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmi2m
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmdiv(mtfm%mfm,ma%mfm,fmdiv_ifm%mfm)
      END FUNCTION fmdiv_ifm

      FUNCTION fmdiv_iim(ival,ma)
! .. Function Return Value ..
        TYPE (im) :: fmdiv_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imdiv, imi2m
! ..
        CALL imi2m(ival,mtim%mim)
        CALL imdiv(mtim%mim,ma%mim,fmdiv_iim%mim)
      END FUNCTION fmdiv_iim

      FUNCTION fmdiv_izm(ival,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmdiv
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,ma%mzm,fmdiv_izm%mzm)
      END FUNCTION fmdiv_izm

      FUNCTION fmdiv_rfm(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmdiv(mtfm%mfm,ma%mfm,fmdiv_rfm%mfm)
      END FUNCTION fmdiv_rfm

      FUNCTION fmdiv_rim(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmdiv(mtfm%mfm,mufm%mfm,fmdiv_rim%mfm)
      END FUNCTION fmdiv_rim

      FUNCTION fmdiv_rzm(r,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmdiv
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,ma%mzm,fmdiv_rzm%mzm)
      END FUNCTION fmdiv_rzm

      FUNCTION fmdiv_dfm(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmdiv(mtfm%mfm,ma%mfm,fmdiv_dfm%mfm)
      END FUNCTION fmdiv_dfm

      FUNCTION fmdiv_dim(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmdp2m, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmdiv(mtfm%mfm,mufm%mfm,fmdiv_dim%mfm)
      END FUNCTION fmdiv_dim

      FUNCTION fmdiv_dzm(d,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,ma%mzm,fmdiv_dzm%mzm)
      END FUNCTION fmdiv_dzm

      FUNCTION fmdiv_zfm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(mtzm%mzm,muzm%mzm,fmdiv_zfm%mzm)
      END FUNCTION fmdiv_zfm

      FUNCTION fmdiv_zim(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(mtzm%mzm,muzm%mzm,fmdiv_zim%mzm)
      END FUNCTION fmdiv_zim

      FUNCTION fmdiv_zzm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,ma%mzm,fmdiv_zzm%mzm)
      END FUNCTION fmdiv_zzm

      FUNCTION fmdiv_cfm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(mtzm%mzm,muzm%mzm,fmdiv_cfm%mzm)
      END FUNCTION fmdiv_cfm

      FUNCTION fmdiv_cim(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(mtzm%mzm,muzm%mzm,fmdiv_cim%mzm)
      END FUNCTION fmdiv_cim

      FUNCTION fmdiv_czm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,ma%mzm,fmdiv_czm%mzm)
      END FUNCTION fmdiv_czm

      FUNCTION fmdiv_fmi(ma,ival)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmdivi
! ..
        CALL fmdivi(ma%mfm,ival,fmdiv_fmi%mfm)
      END FUNCTION fmdiv_fmi

      FUNCTION fmdiv_fmr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmdiv(ma%mfm,mtfm%mfm,fmdiv_fmr%mfm)
      END FUNCTION fmdiv_fmr

      FUNCTION fmdiv_fmd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmdp2m
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmdiv(ma%mfm,mtfm%mfm,fmdiv_fmd%mfm)
      END FUNCTION fmdiv_fmd

      FUNCTION fmdiv_fmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(muzm%mzm,mtzm%mzm,fmdiv_fmz%mzm)
      END FUNCTION fmdiv_fmz

      FUNCTION fmdiv_fmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(muzm%mzm,mtzm%mzm,fmdiv_fmc%mzm)
      END FUNCTION fmdiv_fmc

      FUNCTION fmdiv_fmfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv
! ..
        CALL fmdiv(ma%mfm,mb%mfm,fmdiv_fmfm%mfm)
      END FUNCTION fmdiv_fmfm

      FUNCTION fmdiv_fmim(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, imi2fm
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmdiv(ma%mfm,mtfm%mfm,fmdiv_fmim%mfm)
      END FUNCTION fmdiv_fmim

      FUNCTION fmdiv_fmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmdiv
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(ma%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmdiv(mtzm%mzm,mb%mzm,fmdiv_fmzm%mzm)
      END FUNCTION fmdiv_fmzm

      FUNCTION fmdiv_imi(ma,ival)
! .. Function Return Value ..
        TYPE (im) :: fmdiv_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imdivi
! ..
        CALL imdivi(ma%mim,ival,fmdiv_imi%mim)
      END FUNCTION fmdiv_imi

      FUNCTION fmdiv_imr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmdiv(mufm%mfm,mtfm%mfm,fmdiv_imr%mfm)
      END FUNCTION fmdiv_imr

      FUNCTION fmdiv_imd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmdp2m, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmdiv(mufm%mfm,mtfm%mfm,fmdiv_imd%mfm)
      END FUNCTION fmdiv_imd

      FUNCTION fmdiv_imz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(muzm%mzm,mtzm%mzm,fmdiv_imz%mzm)
      END FUNCTION fmdiv_imz

      FUNCTION fmdiv_imc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(muzm%mzm,mtzm%mzm,fmdiv_imc%mzm)
      END FUNCTION fmdiv_imc

      FUNCTION fmdiv_imfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmdiv_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmdiv(mtfm%mfm,mb%mfm,fmdiv_imfm%mfm)
      END FUNCTION fmdiv_imfm

      FUNCTION fmdiv_imim(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmdiv_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imdiv
! ..
        CALL imdiv(ma%mim,mb%mim,fmdiv_imim%mim)
      END FUNCTION fmdiv_imim

      FUNCTION fmdiv_imzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmdiv
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(muzm%mzm,mb%mzm,fmdiv_imzm%mzm)
      END FUNCTION fmdiv_imzm

      FUNCTION fmdiv_zmi(ma,ival)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmdivi
! ..
        CALL zmdivi(ma%mzm,ival,fmdiv_zmi%mzm)
      END FUNCTION fmdiv_zmi

      FUNCTION fmdiv_zmr(ma,r)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmdiv
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(ma%mzm,mtzm%mzm,fmdiv_zmr%mzm)
      END FUNCTION fmdiv_zmr

      FUNCTION fmdiv_zmd(ma,d)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(ma%mzm,mtzm%mzm,fmdiv_zmd%mzm)
      END FUNCTION fmdiv_zmd

      FUNCTION fmdiv_zmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmdiv, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmdiv(ma%mzm,mtzm%mzm,fmdiv_zmz%mzm)
      END FUNCTION fmdiv_zmz

      FUNCTION fmdiv_zmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmdiv
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmdiv(ma%mzm,mtzm%mzm,fmdiv_zmc%mzm)
      END FUNCTION fmdiv_zmc

      FUNCTION fmdiv_zmfm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmdiv
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmdiv(ma%mzm,mtzm%mzm,fmdiv_zmfm%mzm)
      END FUNCTION fmdiv_zmfm

      FUNCTION fmdiv_zmim(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmdiv
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmdiv(ma%mzm,muzm%mzm,fmdiv_zmim%mzm)
      END FUNCTION fmdiv_zmim

      FUNCTION fmdiv_zmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdiv_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmdiv
! ..
        CALL zmdiv(ma%mzm,mb%mzm,fmdiv_zmzm%mzm)
      END FUNCTION fmdiv_zmzm





!                                                                  **

      FUNCTION fmpwr_ifm(ival,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_ifm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmpwr
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmpwr(mtfm%mfm,ma%mfm,fmpwr_ifm%mfm)
      END FUNCTION fmpwr_ifm

      FUNCTION fmpwr_iim(ival,ma)
! .. Function Return Value ..
        TYPE (im) :: fmpwr_iim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m, impwr
! ..
        CALL imi2m(ival,mtim%mim)
        CALL impwr(mtim%mim,ma%mim,fmpwr_iim%mim)
      END FUNCTION fmpwr_iim

      FUNCTION fmpwr_izm(ival,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_izm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmpwr
! ..
        CALL fmi2m(ival,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,ma%mzm,fmpwr_izm%mzm)
      END FUNCTION fmpwr_izm

      FUNCTION fmpwr_rfm(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_rfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmpwr(mtfm%mfm,ma%mfm,fmpwr_rfm%mfm)
      END FUNCTION fmpwr_rfm

      FUNCTION fmpwr_rim(r,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_rim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmpwr(mtfm%mfm,mufm%mfm,fmpwr_rim%mfm)
      END FUNCTION fmpwr_rim

      FUNCTION fmpwr_rzm(r,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_rzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmpwr
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,ma%mzm,fmpwr_rzm%mzm)
      END FUNCTION fmpwr_rzm

      FUNCTION fmpwr_dfm(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_dfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmpwr
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmpwr(mtfm%mfm,ma%mfm,fmpwr_dfm%mfm)
      END FUNCTION fmpwr_dfm

      FUNCTION fmpwr_dim(d,ma)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_dim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmpwr, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmpwr(mtfm%mfm,mufm%mfm,fmpwr_dim%mfm)
      END FUNCTION fmpwr_dim

      FUNCTION fmpwr_dzm(d,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_dzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,ma%mzm,fmpwr_dzm%mzm)
      END FUNCTION fmpwr_dzm

      FUNCTION fmpwr_zfm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(mtzm%mzm,muzm%mzm,fmpwr_zfm%mzm)
      END FUNCTION fmpwr_zfm

      FUNCTION fmpwr_zim(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(mtzm%mzm,muzm%mzm,fmpwr_zim%mzm)
      END FUNCTION fmpwr_zim

      FUNCTION fmpwr_zzm(z,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,ma%mzm,fmpwr_zzm%mzm)
      END FUNCTION fmpwr_zzm

      FUNCTION fmpwr_cfm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_cfm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(mtzm%mzm,muzm%mzm,fmpwr_cfm%mzm)
      END FUNCTION fmpwr_cfm

      FUNCTION fmpwr_cim(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_cim
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(mtzm%mzm,muzm%mzm,fmpwr_cim%mzm)
      END FUNCTION fmpwr_cim

      FUNCTION fmpwr_czm(c,ma)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_czm
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,ma%mzm,fmpwr_czm%mzm)
      END FUNCTION fmpwr_czm

      FUNCTION fmpwr_fmi(ma,ival)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_fmi
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmipwr
! ..
        CALL fmipwr(ma%mfm,ival,fmpwr_fmi%mfm)
      END FUNCTION fmpwr_fmi

      FUNCTION fmpwr_fmr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_fmr
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, fmsp2m
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmpwr(ma%mfm,mtfm%mfm,fmpwr_fmr%mfm)
      END FUNCTION fmpwr_fmr

      FUNCTION fmpwr_fmd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_fmd
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmpwr
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmpwr(ma%mfm,mtfm%mfm,fmpwr_fmd%mfm)
      END FUNCTION fmpwr_fmd

      FUNCTION fmpwr_fmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_fmz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(muzm%mzm,mtzm%mzm,fmpwr_fmz%mzm)
      END FUNCTION fmpwr_fmz

      FUNCTION fmpwr_fmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_fmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(muzm%mzm,mtzm%mzm,fmpwr_fmc%mzm)
      END FUNCTION fmpwr_fmc

      FUNCTION fmpwr_fmfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_fmfm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr
! ..
        CALL fmpwr(ma%mfm,mb%mfm,fmpwr_fmfm%mfm)
      END FUNCTION fmpwr_fmfm

      FUNCTION fmpwr_fmim(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_fmim
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, imi2fm
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmpwr(ma%mfm,mtfm%mfm,fmpwr_fmim%mfm)
      END FUNCTION fmpwr_fmim

      FUNCTION fmpwr_fmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_fmzm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmpwr
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(ma%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmpwr(mtzm%mzm,mb%mzm,fmpwr_fmzm%mzm)
      END FUNCTION fmpwr_fmzm

      FUNCTION fmpwr_imi(ma,ival)
! .. Function Return Value ..
        TYPE (im) :: fmpwr_imi
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m, impwr
! ..
        CALL imi2m(ival,mtim%mim)
        CALL impwr(ma%mim,mtim%mim,fmpwr_imi%mim)
      END FUNCTION fmpwr_imi

      FUNCTION fmpwr_imr(ma,r)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_imr
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, fmsp2m, imi2fm
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmpwr(mufm%mfm,mtfm%mfm,fmpwr_imr%mfm)
      END FUNCTION fmpwr_imr

      FUNCTION fmpwr_imd(ma,d)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_imd
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmpwr, imi2fm
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imi2fm(ma%mim,mufm%mfm)
        CALL fmpwr(mufm%mfm,mtfm%mfm,fmpwr_imd%mfm)
      END FUNCTION fmpwr_imd

      FUNCTION fmpwr_imz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_imz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(muzm%mzm,mtzm%mzm,fmpwr_imz%mzm)
      END FUNCTION fmpwr_imz

      FUNCTION fmpwr_imc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_imc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, imi2fm, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(muzm%mzm,mtzm%mzm,fmpwr_imc%mzm)
      END FUNCTION fmpwr_imc

      FUNCTION fmpwr_imfm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmpwr_imfm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmpwr(mtfm%mfm,mb%mfm,fmpwr_imfm%mfm)
      END FUNCTION fmpwr_imfm

      FUNCTION fmpwr_imim(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmpwr_imim
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL impwr
! ..
        CALL impwr(ma%mim,mb%mim,fmpwr_imim%mim)
      END FUNCTION fmpwr_imim

      FUNCTION fmpwr_imzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_imzm
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (zm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmpwr
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(muzm%mzm,mb%mzm,fmpwr_imzm%mzm)
      END FUNCTION fmpwr_imzm

      FUNCTION fmpwr_zmi(ma,ival)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmi
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmipwr
! ..
        CALL zmipwr(ma%mzm,ival,fmpwr_zmi%mzm)
      END FUNCTION fmpwr_zmi

      FUNCTION fmpwr_zmr(ma,r)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmr
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx, zmpwr
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(ma%mzm,mtzm%mzm,fmpwr_zmr%mzm)
      END FUNCTION fmpwr_zmr

      FUNCTION fmpwr_zmd(ma,d)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmd
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(ma%mzm,mtzm%mzm,fmpwr_zmd%mzm)
      END FUNCTION fmpwr_zmd

      FUNCTION fmpwr_zmz(ma,z)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmpwr, zmz2m
! ..
        CALL zmz2m(z,mtzm%mzm)
        CALL zmpwr(ma%mzm,mtzm%mzm,fmpwr_zmz%mzm)
      END FUNCTION fmpwr_zmz

      FUNCTION fmpwr_zmc(ma,c)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmc
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx, zmpwr
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,mtzm%mzm)
        CALL zmpwr(ma%mzm,mtzm%mzm,fmpwr_zmc%mzm)
      END FUNCTION fmpwr_zmc

      FUNCTION fmpwr_zmfm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmfm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (fm), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx, zmpwr
! ..
        CALL fmi2m(0,mtfm%mfm)
        CALL zmcmpx(mb%mfm,mtfm%mfm,mtzm%mzm)
        CALL zmpwr(ma%mzm,mtzm%mzm,fmpwr_zmfm%mzm)
      END FUNCTION fmpwr_zmfm

      FUNCTION fmpwr_zmim(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmim
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
        TYPE (im), INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx, zmpwr
! ..
        CALL imi2fm(mb%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,muzm%mzm)
        CALL zmpwr(ma%mzm,muzm%mzm,fmpwr_zmim%mzm)
      END FUNCTION fmpwr_zmim

      FUNCTION fmpwr_zmzm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmpwr_zmzm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmpwr
! ..
        CALL zmpwr(ma%mzm,mb%mzm,fmpwr_zmzm%mzm)
      END FUNCTION fmpwr_zmzm





!                                                                 ABS

      FUNCTION fmabs_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmabs_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmabs
! ..
        CALL fmabs(ma%mfm,fmabs_fm%mfm)
      END FUNCTION fmabs_fm

      FUNCTION fmabs_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmabs_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imabs
! ..
        CALL imabs(ma%mim,fmabs_im%mim)
      END FUNCTION fmabs_im

      FUNCTION fmabs_zm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmabs_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmabs
! ..
        CALL zmabs(ma%mzm,fmabs_zm%mfm)
      END FUNCTION fmabs_zm





!                                                                ACOS

      FUNCTION fmacos_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmacos_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmacos
! ..
        CALL fmacos(ma%mfm,fmacos_fm%mfm)
      END FUNCTION fmacos_fm

      FUNCTION fmacos_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmacos_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmacos
! ..
        CALL zmacos(ma%mzm,fmacos_zm%mzm)
      END FUNCTION fmacos_zm





!                                                               AIMAG


      FUNCTION fmaimag_zm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmaimag_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmimag
! ..
        CALL zmimag(ma%mzm,fmaimag_zm%mfm)
      END FUNCTION fmaimag_zm





!                                                                AINT

      FUNCTION fmaint_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmaint_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmint
! ..
        CALL fmint(ma%mfm,fmaint_fm%mfm)
      END FUNCTION fmaint_fm

      FUNCTION fmaint_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmaint_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmint
! ..
        CALL zmint(ma%mzm,fmaint_zm%mzm)
      END FUNCTION fmaint_zm





!                                                               ANINT

      FUNCTION fmanint_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmanint_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmnint
! ..
        CALL fmnint(ma%mfm,fmanint_fm%mfm)
      END FUNCTION fmanint_fm

      FUNCTION fmanint_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmanint_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmnint
! ..
        CALL zmnint(ma%mzm,fmanint_zm%mzm)
      END FUNCTION fmanint_zm





!                                                                ASIN

      FUNCTION fmasin_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmasin_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmasin
! ..
        CALL fmasin(ma%mfm,fmasin_fm%mfm)
      END FUNCTION fmasin_fm

      FUNCTION fmasin_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmasin_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmasin
! ..
        CALL zmasin(ma%mzm,fmasin_zm%mzm)
      END FUNCTION fmasin_zm





!                                                                ATAN

      FUNCTION fmatan_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmatan_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmatan
! ..
        CALL fmatan(ma%mfm,fmatan_fm%mfm)
      END FUNCTION fmatan_fm

      FUNCTION fmatan_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmatan_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmatan
! ..
        CALL zmatan(ma%mzm,fmatan_zm%mzm)
      END FUNCTION fmatan_zm





!                                                               ATAN2

      FUNCTION fmatan2_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmatan2_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmatn2
! ..
        CALL fmatn2(ma%mfm,mb%mfm,fmatan2_fm%mfm)
      END FUNCTION fmatan2_fm





!                                                               BTEST

      FUNCTION fmbtest_im(ma,pos)
! .. Function Return Value ..
        LOGICAL :: fmbtest_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: pos
! ..
! .. External Subroutines ..
        EXTERNAL imabs, imdiv, imi2m, immod, impwr
! ..
        CALL imi2m(2,mtim%mim)
        CALL imi2m(pos,muim%mim)
        CALL impwr(mtim%mim,muim%mim,muim%mim)
        CALL imdiv(ma%mim,muim%mim,muim%mim)
        CALL imabs(muim%mim,muim%mim)
        CALL immod(muim%mim,mtim%mim,muim%mim)
        IF (muim%mim(2)==0) THEN
          fmbtest_im = .FALSE.
        ELSE
          fmbtest_im = .TRUE.
        END IF
      END FUNCTION fmbtest_im





!                                                             CEILING

      FUNCTION fmceiling_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmceiling_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmaddi, fmeq, fmint, fmsub
! ..
        CALL fmint(ma%mfm,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,mufm%mfm)
        IF (mufm%mfm(2)==0) THEN
          CALL fmeq(ma%mfm,fmceiling_fm%mfm)
        ELSE IF (ma%mfm(2)>0) THEN
          CALL fmaddi(mtfm%mfm,1)
          CALL fmeq(mtfm%mfm,fmceiling_fm%mfm)
        ELSE
          CALL fmeq(mtfm%mfm,fmceiling_fm%mfm)
        END IF
      END FUNCTION fmceiling_fm

      FUNCTION fmceiling_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmceiling_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmcmpx
! ..
        mtzm = ceiling(real(ma))
        mufm = ceiling(aimag(ma))
        mtfm = real(mtzm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,fmceiling_zm%mzm)
      END FUNCTION fmceiling_zm





!                                                               CMPLX


      FUNCTION fmcmplx_fm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmcmplx_fm
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
        TYPE (fm), OPTIONAL, INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx
! ..
        IF (present(mb)) THEN
          CALL zmcmpx(ma%mfm,mb%mfm,fmcmplx_fm%mzm)
        ELSE
          CALL fmi2m(0,mtfm%mfm)
          CALL zmcmpx(ma%mfm,mtfm%mfm,fmcmplx_fm%mzm)
        END IF
      END FUNCTION fmcmplx_fm

      FUNCTION fmcmplx_im(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmcmplx_im
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
        TYPE (im), OPTIONAL, INTENT (IN) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx
! ..
        IF (present(mb)) THEN
          CALL imi2fm(ma%mim,mtfm%mfm)
          CALL imi2fm(mb%mim,mufm%mfm)
          CALL zmcmpx(mtfm%mfm,mufm%mfm,fmcmplx_im%mzm)
        ELSE
          CALL imi2fm(ma%mim,mtfm%mfm)
          CALL fmi2m(0,mufm%mfm)
          CALL zmcmpx(mtfm%mfm,mufm%mfm,fmcmplx_im%mzm)
        END IF
      END FUNCTION fmcmplx_im





!                                                               CONJG


      FUNCTION fmconjg_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmconjg_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmconj
! ..
        CALL zmconj(ma%mzm,fmconjg_zm%mzm)
      END FUNCTION fmconjg_zm





!                                                                 COS

      FUNCTION fmcos_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmcos_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmcos
! ..
        CALL fmcos(ma%mfm,fmcos_fm%mfm)
      END FUNCTION fmcos_fm

      FUNCTION fmcos_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmcos_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmcos
! ..
        CALL zmcos(ma%mzm,fmcos_zm%mzm)
      END FUNCTION fmcos_zm





!                                                                COSH

      FUNCTION fmcosh_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmcosh_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmcosh
! ..
        CALL fmcosh(ma%mfm,fmcosh_fm%mfm)
      END FUNCTION fmcosh_fm

      FUNCTION fmcosh_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmcosh_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmcosh
! ..
        CALL zmcosh(ma%mzm,fmcosh_zm%mzm)
      END FUNCTION fmcosh_zm





!                                                                DBLE


      FUNCTION fmdble_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdble_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,fmdble_fm%mfm)
      END FUNCTION fmdble_fm

      FUNCTION fmdble_im(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdble_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        CALL imi2fm(ma%mim,fmdble_im%mfm)
      END FUNCTION fmdble_im

      FUNCTION fmdble_zm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdble_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,fmdble_zm%mfm)
      END FUNCTION fmdble_zm





!                                                              DIGITS

      FUNCTION fmdigits_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmdigits_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmdigits_fm = ndig
      END FUNCTION fmdigits_fm

      FUNCTION fmdigits_im(ma)
! .. Function Return Value ..
        INTEGER :: fmdigits_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
        fmdigits_im = ndigmx
      END FUNCTION fmdigits_im

      FUNCTION fmdigits_zm(ma)
! .. Function Return Value ..
        INTEGER :: fmdigits_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
        fmdigits_zm = ndig
      END FUNCTION fmdigits_zm





!                                                                 DIM

      FUNCTION fmdim_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmdim_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmdim
! ..
        CALL fmdim(ma%mfm,mb%mfm,fmdim_fm%mfm)
      END FUNCTION fmdim_fm

      FUNCTION fmdim_im(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmdim_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imdim
! ..
        CALL imdim(ma%mim,mb%mim,fmdim_im%mim)
      END FUNCTION fmdim_im





!                                                                DINT

      FUNCTION fmdint_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmdint_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmint
! ..
        CALL fmint(ma%mfm,fmdint_fm%mfm)
      END FUNCTION fmdint_fm

      FUNCTION fmdint_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmdint_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmint
! ..
        CALL zmint(ma%mzm,fmdint_zm%mzm)
      END FUNCTION fmdint_zm





!                                                          DOTPRODUCT

      FUNCTION fmdotproduct_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmdotproduct_fm
! ..
! .. Intrinsic Functions ..
        INTRINSIC lbound, size
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma(:), mb(:)
! ..
! .. Local Scalars ..
        INTEGER :: j, ja, jb
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv, fmeq
! ..
        IF (size(ma)==size(mb)) THEN
          mtfm = 0
          DO j = 1, size(ma)
            ja = lbound(ma,dim=1) + j - 1
            jb = lbound(mb,dim=1) + j - 1
            mtfm = mtfm + ma(ja)*mb(jb)
          END DO
          CALL fmeq(mtfm%mfm,fmdotproduct_fm%mfm)
        ELSE
          mtfm = 1
          mufm = 0
          CALL fmdiv(mtfm%mfm,mufm%mfm,fmdotproduct_fm%mfm)
        END IF
      END FUNCTION fmdotproduct_fm

      FUNCTION fmdotproduct_im(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmdotproduct_im
! ..
! .. Intrinsic Functions ..
        INTRINSIC lbound, size
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma(:), mb(:)
! ..
! .. Local Scalars ..
        INTEGER :: j, ja, jb
! ..
! .. External Subroutines ..
        EXTERNAL imdiv, imeq
! ..
        IF (size(ma)==size(mb)) THEN
          mtim = 0
          DO j = 1, size(ma)
            ja = lbound(ma,dim=1) + j - 1
            jb = lbound(mb,dim=1) + j - 1
            mtim = mtim + ma(ja)*mb(jb)
          END DO
          CALL imeq(mtim%mim,fmdotproduct_im%mim)
        ELSE
          mtim = 1
          muim = 0
          CALL imdiv(mtim%mim,muim%mim,fmdotproduct_im%mim)
        END IF
      END FUNCTION fmdotproduct_im

      FUNCTION fmdotproduct_zm(ma,mb)
! .. Function Return Value ..
        TYPE (zm) :: fmdotproduct_zm
! ..
! .. Intrinsic Functions ..
        INTRINSIC lbound, size
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma(:), mb(:)
! ..
! .. Local Scalars ..
        INTEGER :: j, ja, jb
! ..
! .. External Subroutines ..
        EXTERNAL zmdiv, zmeq
! ..
        IF (size(ma)==size(mb)) THEN
          mtzm = 0
          DO j = 1, size(ma)
            ja = lbound(ma,dim=1) + j - 1
            jb = lbound(mb,dim=1) + j - 1
            mtzm = mtzm + ma(ja)*mb(jb)
          END DO
          CALL zmeq(mtzm%mzm,fmdotproduct_zm%mzm)
        ELSE
          mtzm = 1
          muzm = 0
          CALL zmdiv(mtzm%mzm,muzm%mzm,fmdotproduct_zm%mzm)
        END IF
      END FUNCTION fmdotproduct_zm





!                                                             EPSILON

      FUNCTION fmepsilon_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmepsilon_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmulp
! ..
        CALL fmi2m(1,mtfm%mfm)
        CALL fmulp(mtfm%mfm,fmepsilon_fm%mfm)
      END FUNCTION fmepsilon_fm





!                                                                 EXP

      FUNCTION fmexp_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmexp_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmexp
! ..
        CALL fmexp(ma%mfm,fmexp_fm%mfm)
      END FUNCTION fmexp_fm

      FUNCTION fmexp_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmexp_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmexp
! ..
        CALL zmexp(ma%mzm,fmexp_zm%mzm)
      END FUNCTION fmexp_zm





!                                                            EXPONENT

      FUNCTION fmexponent_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmexponent_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmexponent_fm = int(ma%mfm(1))
      END FUNCTION fmexponent_fm





!                                                               FLOOR

      FUNCTION fmfloor_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmfloor_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmaddi, fmeq, fmint, fmsub
! ..
        CALL fmint(ma%mfm,mtfm%mfm)
        CALL fmsub(ma%mfm,mtfm%mfm,mufm%mfm)
        IF (mufm%mfm(2)==0) THEN
          CALL fmeq(ma%mfm,fmfloor_fm%mfm)
        ELSE IF (ma%mfm(2)<0) THEN
          CALL fmaddi(mtfm%mfm,-1)
          CALL fmeq(mtfm%mfm,fmfloor_fm%mfm)
        ELSE
          CALL fmeq(mtfm%mfm,fmfloor_fm%mfm)
        END IF
      END FUNCTION fmfloor_fm

      FUNCTION fmfloor_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmfloor_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,fmfloor_im%mim)
      END FUNCTION fmfloor_im

      FUNCTION fmfloor_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmfloor_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmcmpx
! ..
        mtzm = floor(real(ma))
        mufm = floor(aimag(ma))
        mtfm = real(mtzm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,fmfloor_zm%mzm)
      END FUNCTION fmfloor_zm





!                                                            FRACTION

      FUNCTION fmfraction_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmfraction_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,mtfm%mfm)
        mtfm%mfm(1) = 0
        CALL fmeq(mtfm%mfm,fmfraction_fm%mfm)
      END FUNCTION fmfraction_fm

      FUNCTION fmfraction_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmfraction_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(ma%mzm,mtzm%mzm)
        mtzm%mzm(1) = 0
        mtzm%mzm(kptimu+1) = 0
        CALL zmeq(mtzm%mzm,fmfraction_zm%mzm)
      END FUNCTION fmfraction_zm





!                                                                HUGE

      FUNCTION fmhuge_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmhuge_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmbig
! ..
        CALL fmbig(fmhuge_fm%mfm)
      END FUNCTION fmhuge_fm

      FUNCTION fmhuge_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmhuge_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imbig
! ..
        CALL imbig(fmhuge_im%mim)
      END FUNCTION fmhuge_im

      FUNCTION fmhuge_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmhuge_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmbig, zmcmpx
! ..
        CALL fmbig(mtfm%mfm)
        CALL zmcmpx(mtfm%mfm,mtfm%mfm,fmhuge_zm%mzm)
      END FUNCTION fmhuge_zm





!                                                                 INT

      FUNCTION fmint_fm(ma)
! .. Function Return Value ..
        TYPE (im) :: fmint_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imfm2i
! ..
        CALL fmint(ma%mfm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,fmint_fm%mim)
      END FUNCTION fmint_fm

      FUNCTION fmint_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmint_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,fmint_im%mim)
      END FUNCTION fmint_im

      FUNCTION fmint_zm(ma)
! .. Function Return Value ..
        TYPE (im) :: fmint_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmint, imfm2i, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmint(mtfm%mfm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,fmint_zm%mim)
      END FUNCTION fmint_zm





!                                                                 LOG

      FUNCTION fmlog_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmlog_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmln
! ..
        CALL fmln(ma%mfm,fmlog_fm%mfm)
      END FUNCTION fmlog_fm

      FUNCTION fmlog_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmlog_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmln
! ..
        CALL zmln(ma%mzm,fmlog_zm%mzm)
      END FUNCTION fmlog_zm





!                                                               LOG10

      FUNCTION fmlog10_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmlog10_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmlg10
! ..
        CALL fmlg10(ma%mfm,fmlog10_fm%mfm)
      END FUNCTION fmlog10_fm

      FUNCTION fmlog10_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmlog10_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmlg10
! ..
        CALL zmlg10(ma%mzm,fmlog10_zm%mzm)
      END FUNCTION fmlog10_zm





!                                                              MATMUL

      FUNCTION fmmatmul_fm(ma,mb) RESULT (mc)
! .. Intrinsic Functions ..
        INTRINSIC lbound, size, ubound
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma(:,:), mb(:,:)
! ..
! .. Local Scalars ..
        INTEGER :: i, j, k
! ..
! .. Function Return Value ..
        TYPE (fm) :: mc(size(ma,dim=1),size(mb,dim=2))
! ..
        IF (size(ma,dim=2)==size(mb,dim=1)) THEN
          DO i = lbound(ma,dim=1), ubound(ma,dim=1)
            DO j = lbound(mb,dim=2), ubound(mb,dim=2)
              mtfm = 0
              DO k = lbound(ma,dim=2), ubound(ma,dim=2)
                mtfm = mtfm + ma(i,k)*mb(k,j)
              END DO
              mc(i,j) = mtfm
            END DO
          END DO
        ELSE
          mtfm = 1
          mufm = 0
          mc(1,1) = mtfm/mufm
          DO i = 1, size(ma,dim=1)
            DO j = 1, size(mb,dim=2)
              mc(i,j) = mc(1,1)
            END DO
          END DO
        END IF
      END FUNCTION fmmatmul_fm

      FUNCTION fmmatmul_im(ma,mb) RESULT (mc)
! .. Intrinsic Functions ..
        INTRINSIC lbound, size, ubound
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma(:,:), mb(:,:)
! ..
! .. Local Scalars ..
        INTEGER :: i, j, k
! ..
! .. Function Return Value ..
        TYPE (im) :: mc(size(ma,dim=1),size(mb,dim=2))
! ..
        IF (size(ma,dim=2)==size(mb,dim=1)) THEN
          DO i = lbound(ma,dim=1), ubound(ma,dim=1)
            DO j = lbound(mb,dim=2), ubound(mb,dim=2)
              mtim = 0
              DO k = lbound(ma,dim=2), ubound(ma,dim=2)
                mtim = mtim + ma(i,k)*mb(k,j)
              END DO
              mc(i,j) = mtim
            END DO
          END DO
        ELSE
          mtim = 1
          muim = 0
          mc(1,1) = mtim/muim
          DO i = 1, size(ma,dim=1)
            DO j = 1, size(mb,dim=2)
              mc(i,j) = mc(1,1)
            END DO
          END DO
        END IF
      END FUNCTION fmmatmul_im

      FUNCTION fmmatmul_zm(ma,mb) RESULT (mc)
! .. Intrinsic Functions ..
        INTRINSIC lbound, size, ubound
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma(:,:), mb(:,:)
! ..
! .. Local Scalars ..
        INTEGER :: i, j, k
! ..
! .. Function Return Value ..
        TYPE (zm) :: mc(size(ma,dim=1),size(mb,dim=2))
! ..
        IF (size(ma,dim=2)==size(mb,dim=1)) THEN
          DO i = lbound(ma,dim=1), ubound(ma,dim=1)
            DO j = lbound(mb,dim=2), ubound(mb,dim=2)
              mtzm = 0
              DO k = lbound(ma,dim=2), ubound(ma,dim=2)
                mtzm = mtzm + ma(i,k)*mb(k,j)
              END DO
              mc(i,j) = mtzm
            END DO
          END DO
        ELSE
          mtzm = 1
          muzm = 0
          mc(1,1) = mtzm/muzm
          DO i = 1, size(ma,dim=1)
            DO j = 1, size(mb,dim=2)
              mc(i,j) = mc(1,1)
            END DO
          END DO
        END IF
      END FUNCTION fmmatmul_zm





!                                                                 MAX

      FUNCTION fmmax_fm(ma,mb,mc,md,me,mf,mg,mh,mi,mj)
! .. Function Return Value ..
        TYPE (fm) :: fmmax_fm
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
        TYPE (fm), OPTIONAL, INTENT (IN) :: mc, md, me, mf, mg, mh, mi, mj
! ..
! .. External Subroutines ..
        EXTERNAL fmeq, fmmax
! ..
        CALL fmmax(ma%mfm,mb%mfm,mtfm%mfm)
        IF (present(mc)) THEN
          CALL fmmax(mtfm%mfm,mc%mfm,mtfm%mfm)
        END IF
        IF (present(md)) THEN
          CALL fmmax(mtfm%mfm,md%mfm,mtfm%mfm)
        END IF
        IF (present(me)) THEN
          CALL fmmax(mtfm%mfm,me%mfm,mtfm%mfm)
        END IF
        IF (present(mf)) THEN
          CALL fmmax(mtfm%mfm,mf%mfm,mtfm%mfm)
        END IF
        IF (present(mg)) THEN
          CALL fmmax(mtfm%mfm,mg%mfm,mtfm%mfm)
        END IF
        IF (present(mh)) THEN
          CALL fmmax(mtfm%mfm,mh%mfm,mtfm%mfm)
        END IF
        IF (present(mi)) THEN
          CALL fmmax(mtfm%mfm,mi%mfm,mtfm%mfm)
        END IF
        IF (present(mj)) THEN
          CALL fmmax(mtfm%mfm,mj%mfm,mtfm%mfm)
        END IF
        CALL fmeq(mtfm%mfm,fmmax_fm%mfm)
      END FUNCTION fmmax_fm

      FUNCTION fmmax_im(ma,mb,mc,md,me,mf,mg,mh,mi,mj)
! .. Function Return Value ..
        TYPE (im) :: fmmax_im
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
        TYPE (im), OPTIONAL, INTENT (IN) :: mc, md, me, mf, mg, mh, mi, mj
! ..
! .. External Subroutines ..
        EXTERNAL imeq, immax
! ..
        CALL immax(ma%mim,mb%mim,mtim%mim)
        IF (present(mc)) THEN
          CALL immax(mtim%mim,mc%mim,mtim%mim)
        END IF
        IF (present(md)) THEN
          CALL immax(mtim%mim,md%mim,mtim%mim)
        END IF
        IF (present(me)) THEN
          CALL immax(mtim%mim,me%mim,mtim%mim)
        END IF
        IF (present(mf)) THEN
          CALL immax(mtim%mim,mf%mim,mtim%mim)
        END IF
        IF (present(mg)) THEN
          CALL immax(mtim%mim,mg%mim,mtim%mim)
        END IF
        IF (present(mh)) THEN
          CALL immax(mtim%mim,mh%mim,mtim%mim)
        END IF
        IF (present(mi)) THEN
          CALL immax(mtim%mim,mi%mim,mtim%mim)
        END IF
        IF (present(mj)) THEN
          CALL immax(mtim%mim,mj%mim,mtim%mim)
        END IF
        CALL imeq(mtim%mim,fmmax_im%mim)
      END FUNCTION fmmax_im





!                                                         MAXEXPONENT

      FUNCTION fmmaxexponent_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmmaxexponent_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmmaxexponent_fm = int(mxexp) + 1
      END FUNCTION fmmaxexponent_fm





!                                                                 MIN

      FUNCTION fmmin_fm(ma,mb,mc,md,me,mf,mg,mh,mi,mj)
! .. Function Return Value ..
        TYPE (fm) :: fmmin_fm
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
        TYPE (fm), OPTIONAL, INTENT (IN) :: mc, md, me, mf, mg, mh, mi, mj
! ..
! .. External Subroutines ..
        EXTERNAL fmeq, fmmin
! ..
        CALL fmmin(ma%mfm,mb%mfm,mtfm%mfm)
        IF (present(mc)) THEN
          CALL fmmin(mtfm%mfm,mc%mfm,mtfm%mfm)
        END IF
        IF (present(md)) THEN
          CALL fmmin(mtfm%mfm,md%mfm,mtfm%mfm)
        END IF
        IF (present(me)) THEN
          CALL fmmin(mtfm%mfm,me%mfm,mtfm%mfm)
        END IF
        IF (present(mf)) THEN
          CALL fmmin(mtfm%mfm,mf%mfm,mtfm%mfm)
        END IF
        IF (present(mg)) THEN
          CALL fmmin(mtfm%mfm,mg%mfm,mtfm%mfm)
        END IF
        IF (present(mh)) THEN
          CALL fmmin(mtfm%mfm,mh%mfm,mtfm%mfm)
        END IF
        IF (present(mi)) THEN
          CALL fmmin(mtfm%mfm,mi%mfm,mtfm%mfm)
        END IF
        IF (present(mj)) THEN
          CALL fmmin(mtfm%mfm,mj%mfm,mtfm%mfm)
        END IF
        CALL fmeq(mtfm%mfm,fmmin_fm%mfm)
      END FUNCTION fmmin_fm

      FUNCTION fmmin_im(ma,mb,mc,md,me,mf,mg,mh,mi,mj)
! .. Function Return Value ..
        TYPE (im) :: fmmin_im
! ..
! .. Intrinsic Functions ..
        INTRINSIC present
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
        TYPE (im), OPTIONAL, INTENT (IN) :: mc, md, me, mf, mg, mh, mi, mj
! ..
! .. External Subroutines ..
        EXTERNAL imeq, immin
! ..
        CALL immin(ma%mim,mb%mim,mtim%mim)
        IF (present(mc)) THEN
          CALL immin(mtim%mim,mc%mim,mtim%mim)
        END IF
        IF (present(md)) THEN
          CALL immin(mtim%mim,md%mim,mtim%mim)
        END IF
        IF (present(me)) THEN
          CALL immin(mtim%mim,me%mim,mtim%mim)
        END IF
        IF (present(mf)) THEN
          CALL immin(mtim%mim,mf%mim,mtim%mim)
        END IF
        IF (present(mg)) THEN
          CALL immin(mtim%mim,mg%mim,mtim%mim)
        END IF
        IF (present(mh)) THEN
          CALL immin(mtim%mim,mh%mim,mtim%mim)
        END IF
        IF (present(mi)) THEN
          CALL immin(mtim%mim,mi%mim,mtim%mim)
        END IF
        IF (present(mj)) THEN
          CALL immin(mtim%mim,mj%mim,mtim%mim)
        END IF
        CALL imeq(mtim%mim,fmmin_im%mim)
      END FUNCTION fmmin_im





!                                                         MINEXPONENT

      FUNCTION fmminexponent_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmminexponent_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmminexponent_fm = -int(mxexp)
      END FUNCTION fmminexponent_fm





!                                                                 MOD

      FUNCTION fmmod_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmmod_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmmod
! ..
        CALL fmmod(ma%mfm,mb%mfm,fmmod_fm%mfm)
      END FUNCTION fmmod_fm

      FUNCTION fmmod_im(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmmod_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL immod
! ..
        CALL immod(ma%mim,mb%mim,fmmod_im%mim)
      END FUNCTION fmmod_im





!                                                              MODULO

      FUNCTION fmmodulo_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmmodulo_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmadd, fmeq, fmmod
! ..
        CALL fmmod(ma%mfm,mb%mfm,mtfm%mfm)
        IF (mtfm%mfm(2)/=0) THEN
          IF ((ma%mfm(2)>0 .AND. mb%mfm(2)<0) .OR. (ma%mfm(2)<0 .AND. mb%mfm( &
              2)>0)) THEN
            CALL fmadd(mtfm%mfm,mb%mfm,mtfm%mfm)
          END IF
        END IF
        CALL fmeq(mtfm%mfm,fmmodulo_fm%mfm)
      END FUNCTION fmmodulo_fm

      FUNCTION fmmodulo_im(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmmodulo_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imadd, imeq, immod
! ..
        CALL immod(ma%mim,mb%mim,mtim%mim)
        IF (mtim%mim(2)/=0) THEN
          IF ((ma%mim(2)>0 .AND. mb%mim(2)<0) .OR. (ma%mim(2)<0 .AND. mb%mim( &
              2)>0)) THEN
            CALL imadd(mtim%mim,mb%mim,mtim%mim)
          END IF
        END IF
        CALL imeq(mtim%mim,fmmodulo_im%mim)
      END FUNCTION fmmodulo_im





!                                                             NEAREST

      FUNCTION fmnearest_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmnearest_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
        EXTERNAL fmabs, fmadd, fmbig, fmdiv, fmi2m, fmsub, fmulp
! ..
        IF (ma%mfm(2)==0) THEN
          IF (mb%mfm(2)>=0) THEN
            CALL fmbig(mtfm%mfm)
            CALL fmi2m(1,mufm%mfm)
            CALL fmdiv(mufm%mfm,mtfm%mfm,fmnearest_fm%mfm)
          ELSE
            CALL fmbig(mtfm%mfm)
            CALL fmi2m(-1,mufm%mfm)
            CALL fmdiv(mufm%mfm,mtfm%mfm,fmnearest_fm%mfm)
          END IF
        ELSE
          IF (mb%mfm(2)>=0) THEN
            CALL fmulp(ma%mfm,mtfm%mfm)
            CALL fmabs(mtfm%mfm,mtfm%mfm)
            CALL fmadd(ma%mfm,mtfm%mfm,mufm%mfm)
            CALL fmulp(mufm%mfm,mufm%mfm)
            CALL fmabs(mufm%mfm,mufm%mfm)
            IF (fmcomp(mtfm%mfm,'LE',mufm%mfm)) THEN
              CALL fmadd(ma%mfm,mtfm%mfm,fmnearest_fm%mfm)
            ELSE
              CALL fmadd(ma%mfm,mufm%mfm,fmnearest_fm%mfm)
            END IF
          ELSE
            CALL fmulp(ma%mfm,mtfm%mfm)
            CALL fmabs(mtfm%mfm,mtfm%mfm)
            CALL fmsub(ma%mfm,mtfm%mfm,mufm%mfm)
            CALL fmulp(mufm%mfm,mufm%mfm)
            CALL fmabs(mufm%mfm,mufm%mfm)
            IF (fmcomp(mtfm%mfm,'LE',mufm%mfm)) THEN
              CALL fmsub(ma%mfm,mtfm%mfm,fmnearest_fm%mfm)
            ELSE
              CALL fmsub(ma%mfm,mufm%mfm,fmnearest_fm%mfm)
            END IF
          END IF
        END IF
      END FUNCTION fmnearest_fm





!                                                                NINT

      FUNCTION fmnint_fm(ma)
! .. Function Return Value ..
        TYPE (im) :: fmnint_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmnint, imfm2i
! ..
        CALL fmnint(ma%mfm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,fmnint_fm%mim)
      END FUNCTION fmnint_fm

      FUNCTION fmnint_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmnint_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,fmnint_im%mim)
      END FUNCTION fmnint_im

      FUNCTION fmnint_zm(ma)
! .. Function Return Value ..
        TYPE (im) :: fmnint_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmnint, imfm2i, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmnint(mtfm%mfm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,fmnint_zm%mim)
      END FUNCTION fmnint_zm





!                                                           PRECISION

      FUNCTION fmprecision_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmprecision_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmprecision_fm = int(log10(real(mbase))*(ndig-1)+1)
      END FUNCTION fmprecision_fm

      FUNCTION fmprecision_zm(ma)
! .. Function Return Value ..
        INTEGER :: fmprecision_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
        fmprecision_zm = int(log10(real(mbase))*(ndig-1)+1)
      END FUNCTION fmprecision_zm





!                                                               RADIX

      FUNCTION fmradix_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmradix_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmradix_fm = int(mbase)
      END FUNCTION fmradix_fm

      FUNCTION fmradix_im(ma)
! .. Function Return Value ..
        INTEGER :: fmradix_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
        fmradix_im = int(mbase)
      END FUNCTION fmradix_im

      FUNCTION fmradix_zm(ma)
! .. Function Return Value ..
        INTEGER :: fmradix_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
        fmradix_zm = int(mbase)
      END FUNCTION fmradix_zm





!                                                               RANGE

      FUNCTION fmrange_fm(ma)
! .. Function Return Value ..
        INTEGER :: fmrange_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
        fmrange_fm = int(mxexp*log10(real(mbase)))
      END FUNCTION fmrange_fm

      FUNCTION fmrange_im(ma)
! .. Function Return Value ..
        INTEGER :: fmrange_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
        fmrange_im = int(ndigmx*log10(real(mbase)))
      END FUNCTION fmrange_im

      FUNCTION fmrange_zm(ma)
! .. Function Return Value ..
        INTEGER :: fmrange_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
        fmrange_zm = int(mxexp*log10(real(mbase)))
      END FUNCTION fmrange_zm





!                                                                REAL


      FUNCTION fmreal_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmreal_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,fmreal_fm%mfm)
      END FUNCTION fmreal_fm

      FUNCTION fmreal_im(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmreal_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        CALL imi2fm(ma%mim,fmreal_im%mfm)
      END FUNCTION fmreal_im

      FUNCTION fmreal_zm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmreal_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,fmreal_zm%mfm)
      END FUNCTION fmreal_zm





!                                                           RRSPACING

      FUNCTION fmrrspacing_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmrrspacing_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmabs, fmeq
! ..
        CALL fmabs(ma%mfm,mtfm%mfm)
        mtfm%mfm(1) = ndig
        CALL fmeq(mtfm%mfm,fmrrspacing_fm%mfm)
      END FUNCTION fmrrspacing_fm





!                                                               SCALE


      FUNCTION fmscale_fm(ma,l)
! .. Function Return Value ..
        TYPE (fm) :: fmscale_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: l
! ..
! .. External Subroutines ..
        EXTERNAL fmeq, fmi2m, fmipwr, fmmpy
! ..
        CALL fmeq(ma%mfm,mtfm%mfm)
        IF (abs(mtfm%mfm(1)+l)<mxexp) THEN
          mtfm%mfm(1) = mtfm%mfm(1) + l
          CALL fmeq(mtfm%mfm,fmscale_fm%mfm)
        ELSE
          CALL fmi2m(int(mbase),mufm%mfm)
          CALL fmipwr(mufm%mfm,l,mufm%mfm)
          CALL fmmpy(ma%mfm,mufm%mfm,fmscale_fm%mfm)
        END IF
      END FUNCTION fmscale_fm

      FUNCTION fmscale_zm(ma,l)
! .. Function Return Value ..
        TYPE (zm) :: fmscale_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: l
! ..
! .. External Subroutines ..
        EXTERNAL zmeq, zmi2m, zmipwr, zmmpy
! ..
        CALL zmeq(ma%mzm,mtzm%mzm)
        IF (abs(mtzm%mzm(1)+l)<mxexp .AND. abs(mtzm%mzm(kptimu+ &
            1)+l)<mxexp) THEN
          mtzm%mzm(1) = mtzm%mzm(1) + l
          mtzm%mzm(kptimu+1) = mtzm%mzm(kptimu+1) + l
          CALL zmeq(mtzm%mzm,fmscale_zm%mzm)
        ELSE
          CALL zmi2m(int(mbase),muzm%mzm)
          CALL zmipwr(muzm%mzm,l,muzm%mzm)
          CALL zmmpy(ma%mzm,muzm%mzm,fmscale_zm%mzm)
        END IF
      END FUNCTION fmscale_zm





!                                                         SETEXPONENT


      FUNCTION fmsetexponent_fm(ma,l)
! .. Function Return Value ..
        TYPE (fm) :: fmsetexponent_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: l
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,mtfm%mfm)
        mtfm%mfm(1) = l
        CALL fmeq(mtfm%mfm,fmsetexponent_fm%mfm)
      END FUNCTION fmsetexponent_fm





!                                                                SIGN

      FUNCTION fmsign_fm(ma,mb)
! .. Function Return Value ..
        TYPE (fm) :: fmsign_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsign
! ..
        CALL fmsign(ma%mfm,mb%mfm,fmsign_fm%mfm)
      END FUNCTION fmsign_fm

      FUNCTION fmsign_im(ma,mb)
! .. Function Return Value ..
        TYPE (im) :: fmsign_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imsign
! ..
        CALL imsign(ma%mim,mb%mim,fmsign_im%mim)
      END FUNCTION fmsign_im





!                                                                 SIN

      FUNCTION fmsin_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsin_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmsin
! ..
        CALL fmsin(ma%mfm,fmsin_fm%mfm)
      END FUNCTION fmsin_fm

      FUNCTION fmsin_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsin_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmsin
! ..
        CALL zmsin(ma%mzm,fmsin_zm%mzm)
      END FUNCTION fmsin_zm





!                                                                SINH

      FUNCTION fmsinh_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsinh_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmsinh
! ..
        CALL fmsinh(ma%mfm,fmsinh_fm%mfm)
      END FUNCTION fmsinh_fm

      FUNCTION fmsinh_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsinh_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmsinh
! ..
        CALL zmsinh(ma%mzm,fmsinh_zm%mzm)
      END FUNCTION fmsinh_zm





!                                                             SPACING


      FUNCTION fmspacing_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmspacing_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmabs, fmulp
! ..
        CALL fmabs(ma%mfm,mtfm%mfm)
        CALL fmulp(mtfm%mfm,fmspacing_fm%mfm)
      END FUNCTION fmspacing_fm





!                                                                SQRT

      FUNCTION fmsqrt_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmsqrt_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmsqrt
! ..
        CALL fmsqrt(ma%mfm,fmsqrt_fm%mfm)
      END FUNCTION fmsqrt_fm

      FUNCTION fmsqrt_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmsqrt_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmsqrt
! ..
        CALL zmsqrt(ma%mzm,fmsqrt_zm%mzm)
      END FUNCTION fmsqrt_zm





!                                                                 TAN

      FUNCTION fmtan_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmtan_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmtan
! ..
        CALL fmtan(ma%mfm,fmtan_fm%mfm)
      END FUNCTION fmtan_fm

      FUNCTION fmtan_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmtan_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmtan
! ..
        CALL zmtan(ma%mzm,fmtan_zm%mzm)
      END FUNCTION fmtan_zm





!                                                                TANH

      FUNCTION fmtanh_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmtanh_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmtanh
! ..
        CALL fmtanh(ma%mfm,fmtanh_fm%mfm)
      END FUNCTION fmtanh_fm

      FUNCTION fmtanh_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmtanh_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmtanh
! ..
        CALL zmtanh(ma%mzm,fmtanh_zm%mzm)
      END FUNCTION fmtanh_zm





!                                                                TINY

      FUNCTION fmtiny_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fmtiny_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmbig, fmdiv, fmi2m
! ..
        CALL fmbig(mtfm%mfm)
        CALL fmi2m(1,mufm%mfm)
        CALL fmdiv(mufm%mfm,mtfm%mfm,fmtiny_fm%mfm)
      END FUNCTION fmtiny_fm

      FUNCTION fmtiny_im(ma)
! .. Function Return Value ..
        TYPE (im) :: fmtiny_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(1,fmtiny_im%mim)
      END FUNCTION fmtiny_im

      FUNCTION fmtiny_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: fmtiny_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmbig, fmdiv, fmi2m, zmcmpx
! ..
        CALL fmbig(mtfm%mfm)
        CALL fmi2m(1,mufm%mfm)
        CALL fmdiv(mufm%mfm,mtfm%mfm,mtfm%mfm)
        CALL zmcmpx(mtfm%mfm,mtfm%mfm,fmtiny_zm%mzm)
      END FUNCTION fmtiny_zm





!                                                               TO_FM

      FUNCTION fm_i(ival)
! .. Function Return Value ..
        TYPE (fm) :: fm_i
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,fm_i%mfm)
      END FUNCTION fm_i

      FUNCTION fm_r(r)
! .. Function Return Value ..
        TYPE (fm) :: fm_r
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(r,fm_r%mfm)
      END FUNCTION fm_r

      FUNCTION fm_d(d)
! .. Function Return Value ..
        TYPE (fm) :: fm_d
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(d,fm_d%mfm)
      END FUNCTION fm_d

      FUNCTION fm_z(z)
! .. Function Return Value ..
        TYPE (fm) :: fm_z
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(real(z),fm_z%mfm)
      END FUNCTION fm_z

      FUNCTION fm_c(c)
! .. Function Return Value ..
        TYPE (fm) :: fm_c
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),fm_c%mfm)
      END FUNCTION fm_c

      FUNCTION fm_fm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fm_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,fm_fm%mfm)
      END FUNCTION fm_fm

      FUNCTION fm_im(ma)
! .. Function Return Value ..
        TYPE (fm) :: fm_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        CALL imi2fm(ma%mim,fm_im%mfm)
      END FUNCTION fm_im

      FUNCTION fm_st(st)
! .. Function Return Value ..
        TYPE (fm) :: fm_st
! ..
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (IN) :: st
! ..
! .. External Subroutines ..
        EXTERNAL fmst2m
! ..
        CALL fmst2m(st,fm_st%mfm)
      END FUNCTION fm_st

      FUNCTION fm_zm(ma)
! .. Function Return Value ..
        TYPE (fm) :: fm_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,fm_zm%mfm)
      END FUNCTION fm_zm





!                                                               TO_IM

      FUNCTION im_i(ival)
! .. Function Return Value ..
        TYPE (im) :: im_i
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,im_i%mim)
      END FUNCTION im_i

      FUNCTION im_r(r)
! .. Function Return Value ..
        TYPE (im) :: im_r
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imfm2i
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,im_r%mim)
      END FUNCTION im_r

      FUNCTION im_d(d)
! .. Function Return Value ..
        TYPE (im) :: im_d
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imfm2i
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,im_d%mim)
      END FUNCTION im_d

      FUNCTION im_z(z)
! .. Function Return Value ..
        TYPE (im) :: im_z
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m, imfm2i
! ..
        CALL fmsp2m(real(z),mtfm%mfm)
        CALL imfm2i(mtfm%mfm,im_z%mim)
      END FUNCTION im_z

      FUNCTION im_c(c)
! .. Function Return Value ..
        TYPE (im) :: im_c
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, imfm2i
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL imfm2i(mtfm%mfm,im_c%mim)
      END FUNCTION im_c

      FUNCTION im_fm(ma)
! .. Function Return Value ..
        TYPE (im) :: im_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imfm2i
! ..
        CALL imfm2i(ma%mfm,im_fm%mim)
      END FUNCTION im_fm

      FUNCTION im_im(ma)
! .. Function Return Value ..
        TYPE (im) :: im_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,im_im%mim)
      END FUNCTION im_im

      FUNCTION im_st(st)
! .. Function Return Value ..
        TYPE (im) :: im_st
! ..
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (IN) :: st
! ..
! .. External Subroutines ..
        EXTERNAL imst2m
! ..
        CALL imst2m(st,im_st%mim)
      END FUNCTION im_st

      FUNCTION im_zm(ma)
! .. Function Return Value ..
        TYPE (im) :: im_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imfm2i, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL imfm2i(mtfm%mfm,im_zm%mim)
      END FUNCTION im_zm





!                                                               TO_ZM

      FUNCTION zm_i(ival)
! .. Function Return Value ..
        TYPE (zm) :: zm_i
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmi2m
! ..
        CALL zmi2m(ival,zm_i%mzm)
      END FUNCTION zm_i

      FUNCTION zm_r(r)
! .. Function Return Value ..
        TYPE (zm) :: zm_r
! ..
! .. Scalar Arguments ..
        REAL, INTENT (IN) :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, fmsp2m, zmcmpx
! ..
        CALL fmsp2m(r,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,zm_r%mzm)
      END FUNCTION zm_r

      FUNCTION zm_d(d)
! .. Function Return Value ..
        TYPE (zm) :: zm_d
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)), INTENT (IN) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, fmi2m, zmcmpx
! ..
        CALL fmdp2m(d,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,zm_d%mzm)
      END FUNCTION zm_d

      FUNCTION zm_z(z)
! .. Function Return Value ..
        TYPE (zm) :: zm_z
! ..
! .. Scalar Arguments ..
        COMPLEX, INTENT (IN) :: z
! ..
! .. External Subroutines ..
        EXTERNAL zmz2m
! ..
        CALL zmz2m(z,zm_z%mzm)
      END FUNCTION zm_z

      FUNCTION zm_c(c)
! .. Function Return Value ..
        TYPE (zm) :: zm_c
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! ..
! .. Scalar Arguments ..
        COMPLEX (kind(0.0D0)), INTENT (IN) :: c
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m, zmcmpx
! ..
        CALL fmdp2m(real(c,kind(0.0D0)),mtfm%mfm)
        CALL fmdp2m(aimag(c),mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,zm_c%mzm)
      END FUNCTION zm_c

      FUNCTION zm_fm(ma)
! .. Function Return Value ..
        TYPE (zm) :: zm_fm
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, zmcmpx
! ..
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(ma%mfm,mufm%mfm,zm_fm%mzm)
      END FUNCTION zm_fm

      FUNCTION zm_im(ma)
! .. Function Return Value ..
        TYPE (zm) :: zm_im
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m, imi2fm, zmcmpx
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmi2m(0,mufm%mfm)
        CALL zmcmpx(mtfm%mfm,mufm%mfm,zm_im%mzm)
      END FUNCTION zm_im

      FUNCTION zm_st(st)
! .. Function Return Value ..
        TYPE (zm) :: zm_st
! ..
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (IN) :: st
! ..
! .. External Subroutines ..
        EXTERNAL zmst2m
! ..
        CALL zmst2m(st,zm_st%mzm)
      END FUNCTION zm_st

      FUNCTION zm_zm(ma)
! .. Function Return Value ..
        TYPE (zm) :: zm_zm
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(ma%mzm,zm_zm%mzm)
      END FUNCTION zm_zm





!                                                              TO_INT

      FUNCTION fm_2int(ma)
! .. Function Return Value ..
        INTEGER :: fm_2int
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2i
! ..
        CALL fmm2i(ma%mfm,fm_2int)
      END FUNCTION fm_2int

      FUNCTION im_2int(ma)
! .. Function Return Value ..
        INTEGER :: im_2int
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imm2i
! ..
        CALL imm2i(ma%mim,im_2int)
      END FUNCTION im_2int

      FUNCTION zm_2int(ma)
! .. Function Return Value ..
        INTEGER :: zm_2int
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmm2i
! ..
        CALL zmm2i(ma%mzm,zm_2int)
      END FUNCTION zm_2int





!                                                               TO_SP

      FUNCTION fm_2sp(ma)
! .. Function Return Value ..
        REAL :: fm_2sp
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp
! ..
        CALL fmm2sp(ma%mfm,fm_2sp)
      END FUNCTION fm_2sp

      FUNCTION im_2sp(ma)
! .. Function Return Value ..
        REAL :: im_2sp
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmm2sp(mtfm%mfm,im_2sp)
      END FUNCTION im_2sp

      FUNCTION zm_2sp(ma)
! .. Function Return Value ..
        REAL :: zm_2sp
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2sp(mtfm%mfm,zm_2sp)
      END FUNCTION zm_2sp





!                                                               TO_DP

      FUNCTION fm_2dp(ma)
! .. Function Return Value ..
        REAL (kind(0.0D0)) :: fm_2dp
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp
! ..
        CALL fmm2dp(ma%mfm,fm_2dp)
      END FUNCTION fm_2dp

      FUNCTION im_2dp(ma)
! .. Function Return Value ..
        REAL (kind(0.0D0)) :: im_2dp
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,im_2dp)
      END FUNCTION im_2dp

      FUNCTION zm_2dp(ma)
! .. Function Return Value ..
        REAL (kind(0.0D0)) :: zm_2dp
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,zm_2dp)
      END FUNCTION zm_2dp





!                                                              TO_SPZ

      FUNCTION fm_2spz(ma)
! .. Function Return Value ..
        COMPLEX :: fm_2spz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Local Scalars ..
        REAL :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp
! ..
        CALL fmm2sp(ma%mfm,r)
        fm_2spz = cmplx(r,0.0)
      END FUNCTION fm_2spz

      FUNCTION im_2spz(ma)
! .. Function Return Value ..
        COMPLEX :: im_2spz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Local Scalars ..
        REAL :: r
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp, imi2fm
! ..
        CALL imi2fm(ma%mim,mtfm%mfm)
        CALL fmm2sp(mtfm%mfm,r)
        im_2spz = cmplx(r,0.0)
      END FUNCTION im_2spz

      FUNCTION zm_2spz(ma)
! .. Function Return Value ..
        COMPLEX :: zm_2spz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmm2z
! ..
        CALL zmm2z(ma%mzm,zm_2spz)
      END FUNCTION zm_2spz





!                                                              TO_DPZ

      FUNCTION fm_2dpz(ma)
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! .. Function Return Value ..
        COMPLEX (kind(0.0D0)) :: fm_2dpz
! ..
! .. Structure Arguments ..
        TYPE (fm), INTENT (IN) :: ma
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp
! ..
        CALL fmm2dp(ma%mfm,d)
        fm_2dpz = cmplx(d,0.0D0,kind(0.0D0))
      END FUNCTION fm_2dpz

      FUNCTION im_2dpz(ma)
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! .. Function Return Value ..
        COMPLEX (kind(0.0D0)) :: im_2dpz
! ..
! .. Structure Arguments ..
        TYPE (im), INTENT (IN) :: ma
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d
! ..
! .. External Subroutines ..
        EXTERNAL imm2dp
! ..
        CALL imm2dp(ma%mim,d)
        im_2dpz = cmplx(d,0.0D0,kind(0.0D0))
      END FUNCTION im_2dpz

      FUNCTION zm_2dpz(ma)
! ..
! .. Intrinsic Functions ..
        INTRINSIC kind
! .. Function Return Value ..
        COMPLEX (kind(0.0D0)) :: zm_2dpz
! ..
! .. Structure Arguments ..
        TYPE (zm), INTENT (IN) :: ma
! ..
! .. Local Scalars ..
        REAL (kind(0.0D0)) :: d1, d2
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp, zmimag, zmreal
! ..
        CALL zmreal(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,d1)
        CALL zmimag(ma%mzm,mtfm%mfm)
        CALL fmm2dp(mtfm%mfm,d2)
        zm_2dpz = cmplx(d1,d2,kind(0.0D0))
      END FUNCTION zm_2dpz





! Interface routines for calling with the FM, IM, and ZM derived types.

      SUBROUTINE fm_abs(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmabs
! ..
        CALL fmabs(ma%mfm,mb%mfm)
      END SUBROUTINE fm_abs

      SUBROUTINE fm_acos(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmacos
! ..
        CALL fmacos(ma%mfm,mb%mfm)
      END SUBROUTINE fm_acos

      SUBROUTINE fm_add(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmadd
! ..
        CALL fmadd(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_add

      SUBROUTINE fm_addi(ma,ival)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmaddi
! ..
        CALL fmaddi(ma%mfm,ival)
      END SUBROUTINE fm_addi

      SUBROUTINE fm_asin(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmasin
! ..
        CALL fmasin(ma%mfm,mb%mfm)
      END SUBROUTINE fm_asin

      SUBROUTINE fm_atan(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmatan
! ..
        CALL fmatan(ma%mfm,mb%mfm)
      END SUBROUTINE fm_atan

      SUBROUTINE fm_atn2(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmatn2
! ..
        CALL fmatn2(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_atn2

      SUBROUTINE fm_big(ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmbig
! ..
        CALL fmbig(ma%mfm)
      END SUBROUTINE fm_big

      SUBROUTINE fm_chsh(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmchsh
! ..
        CALL fmchsh(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_chsh

      FUNCTION fm_comp(ma,lrel,mb)
! .. Function Return Value ..
        LOGICAL :: fm_comp
! ..
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        CHARACTER (2) :: lrel
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: fmcomp
! ..
        fm_comp = fmcomp(ma%mfm,lrel,mb%mfm)
      END FUNCTION fm_comp

      SUBROUTINE fm_cos(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmcos
! ..
        CALL fmcos(ma%mfm,mb%mfm)
      END SUBROUTINE fm_cos

      SUBROUTINE fm_cosh(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmcosh
! ..
        CALL fmcosh(ma%mfm,mb%mfm)
      END SUBROUTINE fm_cosh

      SUBROUTINE fm_cssn(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmcssn
! ..
        CALL fmcssn(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_cssn

      SUBROUTINE fm_dim(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmdim
! ..
        CALL fmdim(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_dim

      SUBROUTINE fm_div(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmdiv
! ..
        CALL fmdiv(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_div

      SUBROUTINE fm_divi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmdivi
! ..
        CALL fmdivi(ma%mfm,ival,mb%mfm)
      END SUBROUTINE fm_divi

      SUBROUTINE fm_dp2m(x,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)) :: x
! ..
! .. External Subroutines ..
        EXTERNAL fmdp2m
! ..
        CALL fmdp2m(x,ma%mfm)
      END SUBROUTINE fm_dp2m

      SUBROUTINE fm_dpm(x,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)) :: x
! ..
! .. External Subroutines ..
        EXTERNAL fmdpm
! ..
        CALL fmdpm(x,ma%mfm)
      END SUBROUTINE fm_dpm

      SUBROUTINE fm_eq(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmeq
! ..
        CALL fmeq(ma%mfm,mb%mfm)
      END SUBROUTINE fm_eq

      SUBROUTINE fm_equ(ma,mb,na,nb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: na, nb
! ..
! .. External Subroutines ..
        EXTERNAL fmequ
! ..
        CALL fmequ(ma%mfm,mb%mfm,na,nb)
      END SUBROUTINE fm_equ

      SUBROUTINE fm_exp(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmexp
! ..
        CALL fmexp(ma%mfm,mb%mfm)
      END SUBROUTINE fm_exp

      SUBROUTINE fm_form(form,ma,string)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form, string
! ..
! .. External Subroutines ..
        EXTERNAL fmform
! ..
        CALL fmform(form,ma%mfm,string)
      END SUBROUTINE fm_form

      SUBROUTINE fm_fprt(form,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form
! ..
! .. External Subroutines ..
        EXTERNAL fmfprt
! ..
        CALL fmfprt(form,ma%mfm)
      END SUBROUTINE fm_fprt

      SUBROUTINE fm_i2m(ival,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmi2m
! ..
        CALL fmi2m(ival,ma%mfm)
      END SUBROUTINE fm_i2m

      SUBROUTINE fm_inp(line,ma,la,lb)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: la, lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL fminp
! ..
        CALL fminp(line,ma%mfm,la,lb)
      END SUBROUTINE fm_inp

      SUBROUTINE fm_int(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmint
! ..
        CALL fmint(ma%mfm,mb%mfm)
      END SUBROUTINE fm_int

      SUBROUTINE fm_ipwr(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmipwr
! ..
        CALL fmipwr(ma%mfm,ival,mb%mfm)
      END SUBROUTINE fm_ipwr

      SUBROUTINE fm_lg10(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmlg10
! ..
        CALL fmlg10(ma%mfm,mb%mfm)
      END SUBROUTINE fm_lg10

      SUBROUTINE fm_ln(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmln
! ..
        CALL fmln(ma%mfm,mb%mfm)
      END SUBROUTINE fm_ln

      SUBROUTINE fm_lni(ival,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmlni
! ..
        CALL fmlni(ival,ma%mfm)
      END SUBROUTINE fm_lni

      SUBROUTINE fm_m2dp(ma,x)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)) :: x
! ..
! .. External Subroutines ..
        EXTERNAL fmm2dp
! ..
        CALL fmm2dp(ma%mfm,x)
      END SUBROUTINE fm_m2dp

      SUBROUTINE fm_m2i(ma,ival)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmm2i
! ..
        CALL fmm2i(ma%mfm,ival)
      END SUBROUTINE fm_m2i

      SUBROUTINE fm_m2sp(ma,x)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        REAL :: x
! ..
! .. External Subroutines ..
        EXTERNAL fmm2sp
! ..
        CALL fmm2sp(ma%mfm,x)
      END SUBROUTINE fm_m2sp

      SUBROUTINE fm_max(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmmax
! ..
        CALL fmmax(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_max

      SUBROUTINE fm_min(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmmin
! ..
        CALL fmmin(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_min

      SUBROUTINE fm_mod(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmmod
! ..
        CALL fmmod(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_mod

      SUBROUTINE fm_mpy(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmmpy
! ..
        CALL fmmpy(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_mpy

      SUBROUTINE fm_mpyi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL fmmpyi
! ..
        CALL fmmpyi(ma%mfm,ival,mb%mfm)
      END SUBROUTINE fm_mpyi

      SUBROUTINE fm_nint(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmnint
! ..
        CALL fmnint(ma%mfm,mb%mfm)
      END SUBROUTINE fm_nint

      SUBROUTINE fm_out(ma,line,lb)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL fmout
! ..
        CALL fmout(ma%mfm,line,lb)
      END SUBROUTINE fm_out

      SUBROUTINE fm_pi(ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmpi
! ..
        CALL fmpi(ma%mfm)
      END SUBROUTINE fm_pi

      SUBROUTINE fm_prnt(ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL fmprnt
! ..
        CALL fmprnt(ma%mfm)
      END SUBROUTINE fm_prnt

      SUBROUTINE fm_pwr(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmpwr
! ..
        CALL fmpwr(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_pwr

      SUBROUTINE fm_read(kread,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kread
! ..
! .. External Subroutines ..
        EXTERNAL fmread
! ..
        CALL fmread(kread,ma%mfm)
      END SUBROUTINE fm_read

      SUBROUTINE fm_rpwr(ma,ival,jval,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival, jval
! ..
! .. External Subroutines ..
        EXTERNAL fmrpwr
! ..
        CALL fmrpwr(ma%mfm,ival,jval,mb%mfm)
      END SUBROUTINE fm_rpwr

      SUBROUTINE fm_set(nprec)
! .. Scalar Arguments ..
        INTEGER :: nprec
! ..
! .. External Subroutines ..
        EXTERNAL fmset
! ..
        CALL fmset(nprec)
      END SUBROUTINE fm_set

      SUBROUTINE fm_sign(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmsign
! ..
        CALL fmsign(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_sign

      SUBROUTINE fm_sin(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsin
! ..
        CALL fmsin(ma%mfm,mb%mfm)
      END SUBROUTINE fm_sin

      SUBROUTINE fm_sinh(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsinh
! ..
        CALL fmsinh(ma%mfm,mb%mfm)
      END SUBROUTINE fm_sinh

      SUBROUTINE fm_sp2m(x,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        REAL :: x
! ..
! .. External Subroutines ..
        EXTERNAL fmsp2m
! ..
        CALL fmsp2m(x,ma%mfm)
      END SUBROUTINE fm_sp2m

      SUBROUTINE fm_sqr(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsqr
! ..
        CALL fmsqr(ma%mfm,mb%mfm)
      END SUBROUTINE fm_sqr

      SUBROUTINE fm_sqrt(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmsqrt
! ..
        CALL fmsqrt(ma%mfm,mb%mfm)
      END SUBROUTINE fm_sqrt

      SUBROUTINE fm_st2m(string,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: string
! ..
! .. External Subroutines ..
        EXTERNAL fmst2m
! ..
        CALL fmst2m(string,ma%mfm)
      END SUBROUTINE fm_st2m

      SUBROUTINE fm_sub(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL fmsub
! ..
        CALL fmsub(ma%mfm,mb%mfm,mc%mfm)
      END SUBROUTINE fm_sub

      SUBROUTINE fm_tan(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmtan
! ..
        CALL fmtan(ma%mfm,mb%mfm)
      END SUBROUTINE fm_tan

      SUBROUTINE fm_tanh(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmtanh
! ..
        CALL fmtanh(ma%mfm,mb%mfm)
      END SUBROUTINE fm_tanh

      SUBROUTINE fm_ulp(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL fmulp
! ..
        CALL fmulp(ma%mfm,mb%mfm)
      END SUBROUTINE fm_ulp

      SUBROUTINE fm_writ(kwrite,ma)
! .. Structure Arguments ..
        TYPE (fm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kwrite
! ..
! .. External Subroutines ..
        EXTERNAL fmwrit
! ..
        CALL fmwrit(kwrite,ma%mfm)
      END SUBROUTINE fm_writ

      SUBROUTINE im_abs(ma,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imabs
! ..
        CALL imabs(ma%mim,mb%mim)
      END SUBROUTINE im_abs

      SUBROUTINE im_add(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imadd
! ..
        CALL imadd(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_add

      SUBROUTINE im_big(ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imbig
! ..
        CALL imbig(ma%mim)
      END SUBROUTINE im_big

      FUNCTION im_comp(ma,lrel,mb)
! .. Function Return Value ..
        LOGICAL :: im_comp
! ..
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. Scalar Arguments ..
        CHARACTER (2) :: lrel
! ..
! .. External Functions ..
        LOGICAL, EXTERNAL :: imcomp
! ..
        im_comp = imcomp(ma%mim,lrel,mb%mim)
      END FUNCTION im_comp

      SUBROUTINE im_dim(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imdim
! ..
        CALL imdim(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_dim

      SUBROUTINE im_div(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imdiv
! ..
        CALL imdiv(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_div

      SUBROUTINE im_divi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imdivi
! ..
        CALL imdivi(ma%mim,ival,mb%mim)
      END SUBROUTINE im_divi

      SUBROUTINE im_divr(ma,mb,mc,md)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc, md
! ..
! .. External Subroutines ..
        EXTERNAL imdivr
! ..
        CALL imdivr(ma%mim,mb%mim,mc%mim,md%mim)
      END SUBROUTINE im_divr

      SUBROUTINE im_dvir(ma,ival,mb,irem)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: irem, ival
! ..
! .. External Subroutines ..
        EXTERNAL imdvir
! ..
        CALL imdvir(ma%mim,ival,mb%mim,irem)
      END SUBROUTINE im_dvir

      SUBROUTINE im_eq(ma,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imeq
! ..
        CALL imeq(ma%mim,mb%mim)
      END SUBROUTINE im_eq

      SUBROUTINE im_fm2i(ma,mb)
! .. Structure Arguments ..
        TYPE (fm) :: ma
        TYPE (im) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imfm2i
! ..
        CALL imfm2i(ma%mfm,mb%mim)
      END SUBROUTINE im_fm2i

      SUBROUTINE im_form(form,ma,string)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form, string
! ..
! .. External Subroutines ..
        EXTERNAL imform
! ..
        CALL imform(form,ma%mim,string)
      END SUBROUTINE im_form

      SUBROUTINE im_fprt(form,ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form
! ..
! .. External Subroutines ..
        EXTERNAL imfprt
! ..
        CALL imfprt(form,ma%mim)
      END SUBROUTINE im_fprt

      SUBROUTINE im_gcd(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imgcd
! ..
        CALL imgcd(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_gcd

      SUBROUTINE im_i2fm(ma,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma
        TYPE (fm) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL imi2fm
! ..
        CALL imi2fm(ma%mim,mb%mfm)
      END SUBROUTINE im_i2fm

      SUBROUTINE im_i2m(ival,ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imi2m
! ..
        CALL imi2m(ival,ma%mim)
      END SUBROUTINE im_i2m

      SUBROUTINE im_inp(line,ma,la,lb)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: la, lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL iminp
! ..
        CALL iminp(line,ma%mim,la,lb)
      END SUBROUTINE im_inp

      SUBROUTINE im_m2dp(ma,x)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        REAL (kind(0.0D0)) :: x
! ..
! .. External Subroutines ..
        EXTERNAL imm2dp
! ..
        CALL imm2dp(ma%mim,x)
      END SUBROUTINE im_m2dp

      SUBROUTINE im_m2i(ma,ival)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL imm2i
! ..
        CALL imm2i(ma%mim,ival)
      END SUBROUTINE im_m2i

      SUBROUTINE im_max(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL immax
! ..
        CALL immax(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_max

      SUBROUTINE im_min(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL immin
! ..
        CALL immin(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_min

      SUBROUTINE im_mod(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL immod
! ..
        CALL immod(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_mod

      SUBROUTINE im_mpy(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL immpy
! ..
        CALL immpy(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_mpy

      SUBROUTINE im_mpyi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL immpyi
! ..
        CALL immpyi(ma%mim,ival,mb%mim)
      END SUBROUTINE im_mpyi

      SUBROUTINE im_mpym(ma,mb,mc,md)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc, md
! ..
! .. External Subroutines ..
        EXTERNAL immpym
! ..
        CALL immpym(ma%mim,mb%mim,mc%mim,md%mim)
      END SUBROUTINE im_mpym

      SUBROUTINE im_out(ma,line,lb)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL imout
! ..
        CALL imout(ma%mim,line,lb)
      END SUBROUTINE im_out

      SUBROUTINE im_pmod(ma,mb,mc,md)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc, md
! ..
! .. External Subroutines ..
        EXTERNAL impmod
! ..
        CALL impmod(ma%mim,mb%mim,mc%mim,md%mim)
      END SUBROUTINE im_pmod

      SUBROUTINE im_prnt(ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL imprnt
! ..
        CALL imprnt(ma%mim)
      END SUBROUTINE im_prnt

      SUBROUTINE im_pwr(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL impwr
! ..
        CALL impwr(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_pwr

      SUBROUTINE im_read(kread,ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kread
! ..
! .. External Subroutines ..
        EXTERNAL imread
! ..
        CALL imread(kread,ma%mim)
      END SUBROUTINE im_read

      SUBROUTINE im_set(nprec)
! .. Scalar Arguments ..
        INTEGER :: nprec
! ..
! .. External Subroutines ..
        EXTERNAL fmset
! ..
        CALL fmset(nprec)
      END SUBROUTINE im_set

      SUBROUTINE im_sign(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imsign
! ..
        CALL imsign(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_sign

      SUBROUTINE im_sqr(ma,mb)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL imsqr
! ..
        CALL imsqr(ma%mim,mb%mim)
      END SUBROUTINE im_sqr

      SUBROUTINE im_st2m(string,ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: string
! ..
! .. External Subroutines ..
        EXTERNAL imst2m
! ..
        CALL imst2m(string,ma%mim)
      END SUBROUTINE im_st2m

      SUBROUTINE im_sub(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (im) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL imsub
! ..
        CALL imsub(ma%mim,mb%mim,mc%mim)
      END SUBROUTINE im_sub

      SUBROUTINE im_writ(kwrite,ma)
! .. Structure Arguments ..
        TYPE (im) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kwrite
! ..
! .. External Subroutines ..
        EXTERNAL imwrit
! ..
        CALL imwrit(kwrite,ma%mim)
      END SUBROUTINE im_writ

      SUBROUTINE zm_abs(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma
        TYPE (fm) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmabs
! ..
        CALL zmabs(ma%mzm,mb%mfm)
      END SUBROUTINE zm_abs

      SUBROUTINE zm_acos(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmacos
! ..
        CALL zmacos(ma%mzm,mb%mzm)
      END SUBROUTINE zm_acos

      SUBROUTINE zm_add(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmadd
! ..
        CALL zmadd(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_add

      SUBROUTINE zm_addi(ma,ival)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmaddi
! ..
        CALL zmaddi(ma%mzm,ival)
      END SUBROUTINE zm_addi

      SUBROUTINE zm_arg(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma
        TYPE (fm) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmarg
! ..
        CALL zmarg(ma%mzm,mb%mfm)
      END SUBROUTINE zm_arg

      SUBROUTINE zm_asin(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmasin
! ..
        CALL zmasin(ma%mzm,mb%mzm)
      END SUBROUTINE zm_asin

      SUBROUTINE zm_atan(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmatan
! ..
        CALL zmatan(ma%mzm,mb%mzm)
      END SUBROUTINE zm_atan

      SUBROUTINE zm_chsh(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmchsh
! ..
        CALL zmchsh(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_chsh

      SUBROUTINE zm_cmpx(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (fm) :: ma, mb
        TYPE (zm) :: mc
! ..
! .. External Subroutines ..
        EXTERNAL zmcmpx
! ..
        CALL zmcmpx(ma%mfm,mb%mfm,mc%mzm)
      END SUBROUTINE zm_cmpx

      SUBROUTINE zm_conj(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmconj
! ..
        CALL zmconj(ma%mzm,mb%mzm)
      END SUBROUTINE zm_conj

      SUBROUTINE zm_cos(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmcos
! ..
        CALL zmcos(ma%mzm,mb%mzm)
      END SUBROUTINE zm_cos

      SUBROUTINE zm_cosh(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmcosh
! ..
        CALL zmcosh(ma%mzm,mb%mzm)
      END SUBROUTINE zm_cosh

      SUBROUTINE zm_cssn(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmcssn
! ..
        CALL zmcssn(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_cssn

      SUBROUTINE zm_div(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmdiv
! ..
        CALL zmdiv(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_div

      SUBROUTINE zm_divi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmdivi
! ..
        CALL zmdivi(ma%mzm,ival,mb%mzm)
      END SUBROUTINE zm_divi

      SUBROUTINE zm_eq(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmeq
! ..
        CALL zmeq(ma%mzm,mb%mzm)
      END SUBROUTINE zm_eq

      SUBROUTINE zm_equ(ma,mb,na,nb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: na, nb
! ..
! .. External Subroutines ..
        EXTERNAL zmequ
! ..
        CALL zmequ(ma%mzm,mb%mzm,na,nb)
      END SUBROUTINE zm_equ

      SUBROUTINE zm_exp(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmexp
! ..
        CALL zmexp(ma%mzm,mb%mzm)
      END SUBROUTINE zm_exp

      SUBROUTINE zm_form(form1,form2,ma,string)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form1, form2, string
! ..
! .. External Subroutines ..
        EXTERNAL zmform
! ..
        CALL zmform(form1,form2,ma%mzm,string)
      END SUBROUTINE zm_form

      SUBROUTINE zm_fprt(form1,form2,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: form1, form2
! ..
! .. External Subroutines ..
        EXTERNAL zmfprt
! ..
        CALL zmfprt(form1,form2,ma%mzm)
      END SUBROUTINE zm_fprt

      SUBROUTINE zm_i2m(ival,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmi2m
! ..
        CALL zmi2m(ival,ma%mzm)
      END SUBROUTINE zm_i2m

      SUBROUTINE zm_2i2m(ival1,ival2,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival1, ival2
! ..
! .. External Subroutines ..
        EXTERNAL zm2i2m
! ..
        CALL zm2i2m(ival1,ival2,ma%mzm)
      END SUBROUTINE zm_2i2m

      SUBROUTINE zm_imag(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma
        TYPE (fm) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmimag
! ..
        CALL zmimag(ma%mzm,mb%mfm)
      END SUBROUTINE zm_imag

      SUBROUTINE zm_inp(line,ma,la,lb)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: la, lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL zminp
! ..
        CALL zminp(line,ma%mzm,la,lb)
      END SUBROUTINE zm_inp

      SUBROUTINE zm_int(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmint
! ..
        CALL zmint(ma%mzm,mb%mzm)
      END SUBROUTINE zm_int

      SUBROUTINE zm_ipwr(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmipwr
! ..
        CALL zmipwr(ma%mzm,ival,mb%mzm)
      END SUBROUTINE zm_ipwr

      SUBROUTINE zm_lg10(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmlg10
! ..
        CALL zmlg10(ma%mzm,mb%mzm)
      END SUBROUTINE zm_lg10

      SUBROUTINE zm_ln(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmln
! ..
        CALL zmln(ma%mzm,mb%mzm)
      END SUBROUTINE zm_ln

      SUBROUTINE zm_m2i(ma,ival)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmm2i
! ..
        CALL zmm2i(ma%mzm,ival)
      END SUBROUTINE zm_m2i

      SUBROUTINE zm_m2z(ma,zval)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX :: zval
! ..
! .. External Subroutines ..
        EXTERNAL zmm2z
! ..
        CALL zmm2z(ma%mzm,zval)
      END SUBROUTINE zm_m2z

      SUBROUTINE zm_mpy(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmmpy
! ..
        CALL zmmpy(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_mpy

      SUBROUTINE zm_mpyi(ma,ival,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival
! ..
! .. External Subroutines ..
        EXTERNAL zmmpyi
! ..
        CALL zmmpyi(ma%mzm,ival,mb%mzm)
      END SUBROUTINE zm_mpyi

      SUBROUTINE zm_nint(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmnint
! ..
        CALL zmnint(ma%mzm,mb%mzm)
      END SUBROUTINE zm_nint

      SUBROUTINE zm_out(ma,line,lb,last1,last2)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: last1, last2, lb
! ..
! .. Array Arguments ..
        CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
        EXTERNAL zmout
! ..
        CALL zmout(ma%mzm,line,lb,last1,last2)
      END SUBROUTINE zm_out

      SUBROUTINE zm_prnt(ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. External Subroutines ..
        EXTERNAL zmprnt
! ..
        CALL zmprnt(ma%mzm)
      END SUBROUTINE zm_prnt

      SUBROUTINE zm_pwr(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmpwr
! ..
        CALL zmpwr(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_pwr

      SUBROUTINE zm_read(kread,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kread
! ..
! .. External Subroutines ..
        EXTERNAL zmread
! ..
        CALL zmread(kread,ma%mzm)
      END SUBROUTINE zm_read

      SUBROUTINE zm_real(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma
        TYPE (fm) :: mb
! ..
! .. External Subroutines ..
        EXTERNAL zmreal
! ..
        CALL zmreal(ma%mzm,mb%mfm)
      END SUBROUTINE zm_real

      SUBROUTINE zm_rpwr(ma,ival,jval,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. Scalar Arguments ..
        INTEGER :: ival, jval
! ..
! .. External Subroutines ..
        EXTERNAL zmrpwr
! ..
        CALL zmrpwr(ma%mzm,ival,jval,mb%mzm)
      END SUBROUTINE zm_rpwr

      SUBROUTINE zm_set(nprec)
! .. Scalar Arguments ..
        INTEGER :: nprec
! ..
! .. External Subroutines ..
        EXTERNAL zmset
! ..
        CALL zmset(nprec)
      END SUBROUTINE zm_set

      SUBROUTINE zm_sin(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmsin
! ..
        CALL zmsin(ma%mzm,mb%mzm)
      END SUBROUTINE zm_sin

      SUBROUTINE zm_sinh(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmsinh
! ..
        CALL zmsinh(ma%mzm,mb%mzm)
      END SUBROUTINE zm_sinh

      SUBROUTINE zm_sqr(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmsqr
! ..
        CALL zmsqr(ma%mzm,mb%mzm)
      END SUBROUTINE zm_sqr

      SUBROUTINE zm_sqrt(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmsqrt
! ..
        CALL zmsqrt(ma%mzm,mb%mzm)
      END SUBROUTINE zm_sqrt

      SUBROUTINE zm_st2m(string,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        CHARACTER (*) :: string
! ..
! .. External Subroutines ..
        EXTERNAL zmst2m
! ..
        CALL zmst2m(string,ma%mzm)
      END SUBROUTINE zm_st2m

      SUBROUTINE zm_sub(ma,mb,mc)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb, mc
! ..
! .. External Subroutines ..
        EXTERNAL zmsub
! ..
        CALL zmsub(ma%mzm,mb%mzm,mc%mzm)
      END SUBROUTINE zm_sub

      SUBROUTINE zm_tan(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmtan
! ..
        CALL zmtan(ma%mzm,mb%mzm)
      END SUBROUTINE zm_tan

      SUBROUTINE zm_tanh(ma,mb)
! .. Structure Arguments ..
        TYPE (zm) :: ma, mb
! ..
! .. External Subroutines ..
        EXTERNAL zmtanh
! ..
        CALL zmtanh(ma%mzm,mb%mzm)
      END SUBROUTINE zm_tanh

      SUBROUTINE zm_writ(kwrite,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        INTEGER :: kwrite
! ..
! .. External Subroutines ..
        EXTERNAL zmwrit
! ..
        CALL zmwrit(kwrite,ma%mzm)
      END SUBROUTINE zm_writ

      SUBROUTINE zm_z2m(zval,ma)
! .. Structure Arguments ..
        TYPE (zm) :: ma
! ..
! .. Scalar Arguments ..
        COMPLEX :: zval
! ..
! .. External Subroutines ..
        EXTERNAL zmz2m
! ..
        CALL zmz2m(zval,ma%mzm)
      END SUBROUTINE zm_z2m

    END MODULE fmzm

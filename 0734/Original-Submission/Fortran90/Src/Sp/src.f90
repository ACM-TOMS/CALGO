! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> true_false.f90
    MODULE true_false ! some convenient true false abbreviations

      IMPLICIT NONE
      PUBLIC
!-----               -------------!

      LOGICAL, PARAMETER :: true = .TRUE., false = .FALSE.

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> low.f90
    MODULE low ! precision specification for real computations

      PUBLIC :: stnd, extd

! This requests the processor to use a real implementation 'stnd'
! which provides at least 6 decimal digits of precision and an
! exponent range of at least 10 ^ +- 35.  This would be suitable for
! low accuracy computations.  It is expected that this precision will
! be available on all machines.

      INTEGER, PARAMETER :: stnd = SELECTED_REAL_KIND(6,35)
!-------------------------

! A few computations are preferably done in higher precision 'extd'. The
! numbers chosen here should be such that the underlying hardware will
! select a higher precision for kind 'extd' than for kind 'stnd', if
! this is feasible.  If a higher precision is not readily available,
! the same values may be used as are given above for 'stnd'. It is
! anticipated that on most machines this higher precision will also
! be available.

      INTEGER, PARAMETER :: extd = SELECTED_REAL_KIND(12,35)
!-------------------------

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> normal.f90
    MODULE normal ! precision specification for real computations

      PUBLIC :: stnd, extd

! This requests the processor to use a real implementation 'stnd'
! which provides at least 12 decimal digits of precision and an
! exponent range of at least 10 ^ +- 50.  It is expected that this
! precision will be available on all machines.

      INTEGER, PARAMETER :: stnd = SELECTED_REAL_KIND(12,50)
!-------------------------

! A few computations are preferably done in higher precision 'extd'. The
! numbers chosen here should be such that the underlying hardware will
! select a higher precision for kind 'extd' than for kind 'stnd', if
! this is feasible.  If a higher precision is not readily available,
! the same values may be used as are given above for 'stnd'. It is
! anticipated that on many machines this higher precision may
! not be available.

!Integer, PARAMETER :: extd = Selected_Real_Kind ( 20, 50 ) ! preferred
      INTEGER, PARAMETER :: extd = SELECTED_REAL_KIND(12,50) ! NAG f90: Sun 4
!-------------------------

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> extended.f90
    MODULE extended ! precision specification for real computations

      PUBLIC :: stnd, extd

! This requests the processor to use a real implementation 'stnd'
! which provides at least 20 decimal digits of precision and an
! exponent range of at least 10 ^ +- 80.  It is expected that this
! precision may not be available on all machines.

      INTEGER, PARAMETER :: stnd = SELECTED_REAL_KIND(20,80)
!-------------------------

! A few computations are preferably done in higher precision 'extd'. The
! numbers chosen here should be such that the underlying hardware will
! select a higher precision for kind 'extd' than for kind 'stnd', if
! this is feasible.  If a higher precision is not readily available,
! the same values may be used as are given above for 'stnd'. It is
! anticipated that on many machines, such an even higher precision may
! not be available.

      INTEGER, PARAMETER :: extd = SELECTED_REAL_KIND(30,80)
!-------------------------

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> long.f90
    MODULE long ! precision specification for real computations

      TYPE :: verylong ! not yet implemented

        PRIVATE
        INTEGER :: value
        TYPE (verylong), POINTER :: next

      END TYPE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> integers.f90
    MODULE integers ! precision specification for integer computations

      PUBLIC :: short, long

! This is provided for those machines where using short integers
! has some advantage. The range here is from -999 to +999. Note that
! 999 = 10^3 -1. No harm will be done if short integers are made
! the same as long integers.

      INTEGER, PARAMETER :: short = SELECTED_INT_KIND(3)
!-------------------------

! The range here is at least from  -9 999 999  to +9 999 999. Note
! that 10^7 - 1 = 9,999,999. This may limit the largest possible value
! of the dimension  n  of problems which can be solved. If n is to
! be larger, the 7 should be replaced with a number k so that
! n is considerably less than 10^k.

      INTEGER, PARAMETER :: long = SELECTED_INT_KIND(7)
!-------------------------

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> precision.f90
    MODULE precisions
      USE integers, ONLY : short, long
      USE normal, ONLY : stnd, extd

! This provides a convenient way of selecting the precision
! required for a computation. By simply ensuring that a leading '!'
! appears on all but exactly one of the following USE statements,
! and then recompiling all routines, the precision of an entire
! computation can be altered.

!    USE Low
!    USE Extended
!    USE Long


      PRIVATE

      PUBLIC :: stnd, extd, short, long

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> num_constants.f90
    MODULE num_constants ! Values of machine numeric characteristics

! This provides simple names for the various machine dependent
! constants with intrinsic values in Fortran 90.  All are for
! precision 'stnd'.

      USE precisions, ONLY : stnd

      IMPLICIT NONE
      PUBLIC
!-----               -------------!

      PRIVATE :: stnd, x

      REAL (stnd) :: x !dummy value
      REAL (stnd), PARAMETER :: machtol = EPSILON(x), machhuge = HUGE(x), &
        machmaxexp = MAXEXPONENT(x), machminexp = MINEXPONENT(x), &
        machdecprec = PRECISION(x), machbase = RADIX(x), &
        machdecexpr = RANGE(x), machtiny = TINY(x)

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> reals.f90
    MODULE reals ! Real constants at a given precision

! This provides names for all required real values. By only
! using real values as defined within this module, all
! problems associated with the precision of real literal
! values can be totally avoided.

! All values are at 'stnd' precision. There are values for the
! integral values to 100, the fractions in tenths from 0 to 1,
! the reciprocals of the integers to 10, some common
! fractions, and whatever miscellaneous real values are
! necessary.  A standard naming convention is used.


      USE precisions, ONLY : stnd

      IMPLICIT NONE
      PUBLIC
!-----       -------------!

      PRIVATE :: stnd

! digits
      REAL (stnd), PARAMETER :: zero = 0, one = 1, two = 2, three = 3, &
        four = 4, five = 5, six = 6, seven = 7, eight = 8, nine = 9

! tenths
      REAL (stnd), PARAMETER :: c0_1 = 0.1_stnd, c0_2 = 0.2_stnd, &
        c0_3 = 0.3_stnd, c0_4 = 0.4_stnd, c0_5 = 0.5_stnd, c0_6 = 0.6_stnd, &
        c0_7 = 0.7_stnd, c0_8 = 0.8_stnd, c0_9 = 0.9_stnd

! reciprocals
      REAL (stnd), PARAMETER :: tenth = 0.1_stnd, ninth = one/nine, &
        eighth = 0.125_stnd, seventh = one/seven, sixth = one/six, &
        fifth = 0.2_stnd, quarter = 0.25_stnd, third = one/three, &
        half = 0.5_stnd

! fractions a/b named as fa_b
      REAL (stnd), PARAMETER :: f1_29 = one/29, f2_3 = two/three, &
        f4_3 = four/three, f7_3 = seven/three

! integral values to 99
      REAL (stnd), PARAMETER :: ten = 10, c10 = 10, c40 = 40, c70 = 70, &
        c11 = 11, c41 = 41, c71 = 71, c12 = 12, c42 = 42, c72 = 72, c13 = 13, &
        c43 = 43, c73 = 73, c14 = 14, c44 = 44, c74 = 74, c15 = 15, c45 = 45, &
        c75 = 75, c16 = 16, c46 = 46, c76 = 76, c17 = 17, c47 = 47, c77 = 77, &
        c18 = 18, c48 = 48, c78 = 78, c19 = 19, c49 = 49, c79 = 79, c20 = 20, &
        c50 = 50, c80 = 80, c21 = 21, c51 = 51, c81 = 81, c22 = 22, c52 = 52, &
        c82 = 82, c23 = 23, c53 = 53, c83 = 83, c24 = 24, c54 = 54, c84 = 84, &
        c25 = 25, c55 = 55, c85 = 85, c26 = 26, c56 = 56, c86 = 86, c27 = 27, &
        c57 = 57, c87 = 87, c28 = 28, c58 = 58, c88 = 88, c29 = 29, c59 = 59, &
        c89 = 89, c30 = 30, c60 = 60, c90 = 90, c31 = 31, c61 = 61, c91 = 91, &
        c32 = 32, c62 = 62, c92 = 92, c33 = 33, c63 = 63, c93 = 93, c34 = 34, &
        c64 = 64, c94 = 94, c35 = 35, c65 = 65, c95 = 95, c36 = 36, c66 = 66, &
        c96 = 96, c37 = 37, c67 = 67, c97 = 97, c38 = 38, c68 = 68, c98 = 98, &
        c39 = 39, c69 = 69, c99 = 99

! miscellaneous integral values
      REAL (stnd), PARAMETER :: c100 = 100, c180 = 180, c200 = 200, &
        c256 = 256, c360 = 360, c400 = 400, c600 = 600, c681 = 681, &
        c991 = 991, c1162 = 1162, c2324 = 2324, c10000 = 10000, c40000 = 40000

! miscellaneous real values
! form: d.dd          named as  cd_dd
!                 nn
! form: d.dd x 10     named as  cd_ddEnn
!                -nn
! form: d.dd x 10     named as  cd_ddMnn
      REAL (stnd), PARAMETER :: c1_e6 = 1.0E6_stnd, c2_m6 = 2.0E-6_stnd, &
        c1_m13 = 1.0E-13_stnd, c1_m5 = 1.0E-5_stnd, c1_m4 = 1.0E-4_stnd, &
        c4_m2 = 4.0E-2_stnd, c1_m2 = 1.0E-2_stnd, c0_1136 = 0.1136_stnd, &
        c1_0001 = 1.0001_stnd, c1_2 = 1.2_stnd, c1_5 = 1.5_stnd, &
        c2_25 = 2.25_stnd, c2_5 = 2.5_stnd, c2_625 = 2.625_stnd, &
        c7_5 = 7.5_stnd, c10_1 = 10.1_stnd, c19_8 = 19.8_stnd, &
        c20_2 = 20.2_stnd

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> strings.f90
    MODULE strings ! Routines for simple string stuff

      USE precisions, ONLY : short

      IMPLICIT NONE
      PUBLIC
!------              -------------!

! Names for common characters
      CHARACTER (1), PARAMETER :: ampersand = '&', apostrophe = '''', &
        atsign = '@', backslash = '\\', backquote = '`', bang = '!', &
        blank = ' ', caret = '^', cbrace = '}', cbracket = ']', cparen = ')', &
        colon = ':', comma = ',', dash = '-', dollar = '$', equals = '=', &
        exclamation = '!', greaterthan = '>', hash = '#', lessthan = '<', &
        minus = '-', obrace = '{', obracket = '[', oparen = '(', &
        percent = '%', period = '.', plus = '+', quesmark = '?', quote = '"', &
        semicolon = ';', slash = '/', star = '*', tilde = '~', vertbar = '|', &
        underscore = '_'


      CHARACTER (0), PARAMETER :: null = ''

! Codes for case conversions
      INTEGER (short), PARAMETER :: toupper = 1, tolower = 2, capitalize = 3

! conversion
      INTEGER (short), PARAMETER, PRIVATE :: shift = ICHAR('a') - ICHAR('A')

! This assumes that the relation between a lower case letter and
! its corresponding upper case letter is the same for every letter.

    CONTAINS

      SUBROUTINE ascii_case_change(string,type)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        CHARACTER (*), INTENT (INOUT) :: string
        INTEGER (short), INTENT (IN) :: type

! DESCRIPTION:

!       This converts each lower case alphabetic letter in string to upper
!       case, or vice versa.  Specifically,

!          If type = ToUpper,    conversion is lower to upper
!          If type = ToLower,    conversion is upper to lower
!          If type = Capitalize, use upper for first letter; lower for rest

!       Definitions of ToUpper, ToLower and Capitalize may be obtained from
!       the host module Strings.

!       All non-alphabetic characters are left unchanged.

!   It uses the ASCII character set.

! LOCAL DECLARATIONS:

        INTEGER (short) :: i

! EXECUTION:

        DO i = 1, LEN(string)
          SELECT CASE (type)
          CASE (toupper)
            SELECT CASE (string(i:i))
            CASE ('a':'z')
              string(i:i) = ACHAR(IACHAR(string(i:i))+shift)
            END SELECT
          CASE (tolower,capitalize)
            SELECT CASE (string(i:i))
            CASE ('A':'Z')
              string(i:i) = ACHAR(IACHAR(string(i:i))-shift)
            END SELECT
          END SELECT
        END DO

        IF (type==capitalize) THEN
          SELECT CASE (string(1:1))
          CASE ('a':'z')
            string(1:1) = ACHAR(IACHAR(string(1:1))+shift)
          END SELECT
        END IF

! EXIT:
        RETURN

! FORMATS:  none.

      END SUBROUTINE


      SUBROUTINE case_change(string,type)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        CHARACTER (*), INTENT (INOUT) :: string
        INTEGER (short), INTENT (IN) :: type

! DESCRIPTION:

!       This converts each lower case alphabetic letter in string to upper
!       case, or vice versa.  Specifically,

!          If type = ToUpper,    conversion is lower to upper
!          If type = ToLower,    conversion is upper to lower
!          If type = Capitalize, use upper for first letter; lower for rest

!       Definitions of ToUpper, ToLower and Capitalize may be obtained from
!       the host module Strings.

!       All non-alphabetic characters are left unchanged.

!   It uses the underlying machine character set.

! LOCAL DECLARATIONS:

        INTEGER (short) :: i

! EXECUTION:

        DO i = 1, LEN(string)
          SELECT CASE (type)
          CASE (toupper)
            SELECT CASE (string(i:i))
            CASE ('a':'i')
              string(i:i) = CHAR(ICHAR(string(i:i))+shift)
            CASE ('j':'r')
              string(i:i) = CHAR(ICHAR(string(i:i))+shift)
            CASE ('s':'z')
              string(i:i) = CHAR(ICHAR(string(i:i))+shift)
            END SELECT
          CASE (tolower,capitalize)
            SELECT CASE (string(i:i))
            CASE ('A':'I')
              string(i:i) = CHAR(ICHAR(string(i:i))-shift)
            CASE ('J':'R')
              string(i:i) = CHAR(ICHAR(string(i:i))-shift)
            CASE ('S':'Z')
              string(i:i) = CHAR(ICHAR(string(i:i))-shift)
            END SELECT
          END SELECT
        END DO

        IF (type==capitalize) THEN
          SELECT CASE (string(1:1))
          CASE ('a':'i')
            string(1:1) = CHAR(ICHAR(string(1:1))+shift)
          CASE ('j':'r')
            string(1:1) = CHAR(ICHAR(string(1:1))+shift)
          CASE ('s':'z')
            string(1:1) = CHAR(ICHAR(string(1:1))+shift)
          END SELECT
        END IF

! EXIT:
        RETURN

! FORMATS:  none.

      END SUBROUTINE


      SUBROUTINE mid_shift(string,from,to,number)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        CHARACTER (*), INTENT (INOUT) :: string
        INTEGER (short), INTENT (IN) :: from, to, number

! DESCRIPTION:

!     This routine performs a shift of characters within string. The
!     number of characters shifted is number and they are shifted so
!     that the character in position 'from' is moved to position 'to'.
!     Characters in the to position are overwritten. Blanks replace
!     characters in the from position. Shifting may be left or right,
!     and the from and to positions may overlap.  Care is taken not to
!     alter or use any characters beyond the defined limits of the string.

! LOCAL DECLARATIONS:

        INTEGER (short) :: end, end1, n, shorten, slen
        CHARACTER (len=number) :: substring

! EXECUTION:

        slen = LEN(string)

        IF (from/=to) THEN

          end1 = from + number - 1 ! end1 is the initial position of the last
! character in the substring specified to
! be shifted.

! n is the number of characters that will be shifted and that will
! remain within the confines of string.  n may have to be reduced
! for a shift to the right to ensure that characters beyond the
! end of string are not dealt with (this is not necessary for left
! shifts since the substring moved to the left can never extend
! beyond the end of string).

          IF (end1>slen) THEN
            shorten = end1 - slen
            n = number - shorten
          ELSE
            n = number
          END IF

          IF (from<to) n = MIN(slen,to+n-1_short) - to + 1

          end = from + n - 1 ! end is the last character that will
! actually be shifted (and still remain
! within the limits of string).

! substring is a temporary fix until a compiler bug is fixed.
          substring(1:n) = string(from:end)
          string(to:MIN(slen,to+n-1_short)) = substring(1:n)

          IF (from<to) THEN ! shift to right
            string(from:MIN(to-1_short,end1)) = blank
          ELSE ! shift to left
            string(MAX(to+n,from):end) = blank
          END IF

        END IF

! EXIT:
        RETURN

! FORMATS:  none.

      END SUBROUTINE


      SUBROUTINE center(string)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        CHARACTER (*), INTENT (INOUT) :: string

! DESCRIPTION:

!     This routine shifts the nonblank characters of string so that
!     there is a balance of blanks on left and right.


! LOCAL DECLARATIONS:

        INTEGER (short) :: start, endch, clen, left

! EXECUTION:

        endch = LEN_TRIM(string) ! Find last non-blank character.

        IF (endch/=0) THEN

          start = VERIFY(string,blank) ! Find first nonblank character.

          clen = endch - start + 1 ! Compute shift and do it.
          left = 1 + (LEN(string)-clen)/2

          IF (start>left) THEN
            string(left:) = ADJUSTL(string(left:)) ! move left
          ELSE
            string(1:left+1) = ADJUSTR(string(1:left+1)) ! move right
          END IF

        END IF

! EXIT:
        RETURN

! FORMATS:  None.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> general.f90
    MODULE general ! General purpose routines
      USE strings, ONLY : blank, mid_shift
      USE precisions, ONLY : stnd, long, short
      USE reals, ONLY : one

!!!!USE   F90_UNIX    ! for use with the second version of CpuSecs below

      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: entering, leaving, indent, cpusecs, my_date_time, writemat, &
        writevec

      INTEGER (short), PARAMETER :: spaces = 2
      INTEGER (short) :: lastcount = 0

      CHARACTER (9), PARAMETER :: months(12) = (/ 'January  ', 'February ', &
        'March    ', 'April    ', 'May      ', 'June     ', 'July     ', &
        'August   ', 'September', 'October  ', 'November ', 'December '/), &
        days(0:6) = (/ 'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', &
        'Thursday ', 'Friday   ', 'Saturday '/)

      INTEGER (short), PARAMETER :: defunit = 6 ! For WriteMat and WriteVec
      INTEGER (short), PARAMETER :: defline = 80 ! Default values
      INTEGER (short), PARAMETER :: defindent = 0
      CHARACTER (1), PARAMETER :: deff = 'f'
      INTEGER (short), PARAMETER :: defw = 12
      INTEGER (short), PARAMETER :: defd = 6
      INTEGER (short), PARAMETER :: defs = 3

    CONTAINS

      FUNCTION entering(string,trace,unit,level)

! Upon entering a procedure, this function will be called.
! It will return a prefix string suitable for indenting
! output lines from the procedure. It takes the given string
! and prepends 'level' blanks, followed by a '[', and appends
! the character ']'. For example, if string were 'hi' and level
! were 7, it would return '       [hi]'.  Level is then also
! incremented by the value of Spaces.

! If trace is true, it also outputs a message that the
! routine identified by string was entered.

! ARGUMENTS:

        CHARACTER (*), INTENT (IN) :: string
        LOGICAL, INTENT (IN) :: trace
        INTEGER (short), INTENT (IN) :: unit
        INTEGER (short), INTENT (INOUT) :: level

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(string)+2+level) :: entering

! EXECUTION:

        entering(1:level) = blank
        IF (trace) WRITE (unit,'(A)') entering(1:level) // '[Into ' // &
          string // ']'
        entering(level+1:) = '[' // string // ']'

        level = level + spaces

      END FUNCTION

      SUBROUTINE leaving(string,trace,unit,level)

! This is the 'opposite' to Entering. It should be called
! just before leaving a routine. Level is reduced by
! Spaces and if trace is true, an exit message is output.

! ARGUMENTS:

        CHARACTER (*), INTENT (IN) :: string
        LOGICAL, INTENT (IN) :: trace
        INTEGER (short), INTENT (IN) :: unit
        INTEGER (short), INTENT (INOUT) :: level

! EXECUTION:

        level = level - spaces
        IF (trace) WRITE (unit,'(A)') REPEAT(blank,INT(level)) // '[Done ' // &
          string // ']'

      END SUBROUTINE

      SUBROUTINE indent(id,unit,level)

! This is also used to indent output, albeit in a manner
! different from Entering and Leaving. It simply writes
! out 'level' blanks followed by the string id in [], and
! leaves the output file marker where it is. It uses
! nonadvancing output.  If level is not present, just
! the id part is output; i.e. level is treated as zero.

! ARGUMENTS:

        CHARACTER (*), INTENT (IN) :: id
        INTEGER (short), INTENT (IN) :: unit
        INTEGER (short), INTENT (IN), OPTIONAL :: level

! LOCAL DECLARATIONS:

        INTEGER (short) :: lev

! EXECUTION:

        IF (PRESENT(level)) THEN
          lev = level
        ELSE
          lev = 0
        END IF

        WRITE (unit,'(A)',advance='no') REPEAT(blank,INT(lev)) // '[' // id // &
          ']'

      END SUBROUTINE

!!!!FUNCTION CpuSecs()

! This function obtains, from a C system routine cputime,
! the current value of the system CPU usage clock. This value (in seconds)
! is then returned as a standard precision real value.

!!!!    Real(stnd) :: CpuSecs, cputime

!!!!    CpuSecs  = cputime()

!!!!    return

!!!!end Function CpuSecs

!!!!FUNCTION CpuSecs()

! This function obtains, from the Unix module F90_UNIX,
! the current value of the system CPU usage clock. This value (in seconds)
! is then returned as a standard precision real value.
! The module F90_UNIX, provided with NAG compiler v2.0, contains the
! functions times and clock_ticks_per_second, as well as the defined type
! TMS.  The USE statement for the F90_UNIX module, located at line 6 of
! general.f90, must be uncommented when using this version of CpuSecs.

!!!!    Real(stnd) :: CpuSecs
!!!!    Type(TMS)  :: buffer

!!!!    CpuSecs  = times(buffer)
!!!!    CpuSecs  = real((buffer%utime+buffer%stime),stnd) / &
!!!!               clock_ticks_per_second()

!!!!    return

!!!!end Function CpuSecs

      FUNCTION cpusecs()

! This function obtains, from the intrinsic routine system_clock,
! the current value of the system CPU usage clock. This value
! is then converted to seconds and returned as a standard precision
! real value.

! LOCAL DECLARATIONS:

        REAL (stnd) :: cpusecs
        INTEGER :: count, count_rate, count_max
        REAL (stnd) :: secs

! EXECUTION:

        CALL SYSTEM_CLOCK(count,count_rate,count_max)

        secs = REAL(count,stnd)/REAL(count_rate,stnd)

! wraparound of clock ticker
        IF (count<lastcount) secs = secs + REAL(count_max,stnd)/REAL( &
          count_rate,stnd)

        lastcount = count

        cpusecs = secs

        RETURN

      END FUNCTION

      SUBROUTINE my_date_time(chdate)

! ARGUMENTS:

        CHARACTER (*) :: chdate

! DESCRIPTION:

!  This routine returns in chdate a 41-character date of the form given
!  in model (below). It uses the time and date as obtained from the
!  intrinsic routine Date_and_Time and converts them to the form of the
!  model given below.  The time and date are returned from date_and_time
!  as character strings, respectively, of the form:

!        time:  hhmmss.sss
!        date:  yyyymmdd

!  Note that excess blanks in the date are eliminated.
!  If chdate is more than 41 characters in length, only the
!  leftmost 41 will be altered.  If it is less than 41 in
!  length, only the leftmost characters of the date will be returned.

! LOCAL DECLARATIONS:

        INTEGER (short), PARAMETER :: pthour = 1, ptmin = 4, ptampm = 7, &
          ptmon = 24, ptday = 34, ptyear = 38, ptdayn = 13

        CHARACTER (*), PARAMETER :: model = &
          '00:00 a.m., Wednesday, September 00, 1999'

        INTEGER (short) :: kmon, to, k, modlen
        INTEGER (long) :: dayno

        CHARACTER (10) :: time
        CHARACTER (8) :: date
        CHARACTER (41) :: tdate

! EXECUTION:

        tdate = model
        modlen = LEN_TRIM(tdate)

        CALL DATE_AND_TIME(date,time)

        IF (date(7:7)=='0') date(7:7) = ' '

        tdate(ptday:ptday+1) = date(7:8)
        tdate(ptyear:ptyear+3) = date(1:4)

        READ (date(7:8),'(i2)') dayno
        READ (date(1:4),'(i4)') k

        READ (date(5:6),'(i2)') kmon
        tdate(ptmon:ptmon+8) = months(kmon)
        to = LEN_TRIM(months(kmon))

        IF (to/=9) CALL mid_shift(tdate,ptmon+9_short,ptmon+to,modlen)

        IF (kmon==1 .OR. kmon==2) THEN
          kmon = kmon + 13
          k = k - 1
        ELSE
          kmon = kmon + 1
        END IF

        dayno = dayno + INT(REAL(kmon)*30.6001)
        dayno = dayno + INT(REAL(k)*365.25)
        dayno = MOD(dayno+5,7)

        tdate(ptmin:ptmin+1) = time(3:4) ! minute
        READ (time(1:2),'(i2)') k ! hour

        IF (k>=13) THEN
          k = k - 12
          tdate(ptampm:ptampm) = 'p'
        ELSE IF (k==12) THEN
          tdate(ptampm:ptampm) = 'p'
        ELSE IF (k==0) THEN
          k = k + 12
          tdate(ptampm:ptampm) = 'a'
        ELSE
          tdate(ptampm:ptampm) = 'a'
        END IF

        WRITE (tdate(pthour:pthour+1),'(i2)') k
        tdate(ptdayn:ptdayn+8) = days(dayno)
        k = LEN_TRIM(days(dayno))

        IF (k/=9) THEN ! ==> shift over blanks.
          CALL mid_shift(tdate,ptdayn+9_short,ptdayn+k,modlen)
        END IF

        modlen = MIN(modlen,INT(LEN_TRIM(chdate),short))
        chdate(1:modlen) = tdate

        RETURN

      END SUBROUTINE

! Optional
      SUBROUTINE writemat(x,a,b,y,f,w,d,s,unit,name,indent,line)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:,:)

        REAL (stnd), INTENT (IN), OPTIONAL :: y(SIZE(x,1),SIZE(x,2)), a, b
        INTEGER (short), INTENT (IN), OPTIONAL :: unit, w, d, s, indent, line
        CHARACTER (1), INTENT (IN), OPTIONAL :: f
        CHARACTER (*), INTENT (IN), OPTIONAL :: name

! DESCRIPTION:

!   Print out a submatrix of a matrix with given format, as below.

!   Print a title for the matrix: name

!   If the value of the argument  indent  is positive, then each output
!   line is preceded by  indent  blank characters.  If  indent  is negative,
!   then only the matrix output will be indented (not the title).

!   Print the values from the matrix a*X + b*Y.
!   a and b are scalars; X and Y are matrices.

!   Print each entry in format  fw.d
!              ..f is a character 'f', 'g', 'e' or 'd'
!              ..w and d are integers
!              ..s is the number of spaces between each entry.
!              ..line is the number of characters per line.
!                  (if line=0, then line is replaced with 80)

!   All output is on the unit 'unit'.

!   Defaults are defined for all optional arguments.  See start of module.


!   LOCAL DECLARATIONS:

        REAL (stnd) :: fa, fb
        INTEGER (short) :: actunit, actw, actd, acts, actindent, actline, i, j
        CHARACTER (1) :: actf
        CHARACTER (34), PARAMETER :: defform = &
          '((??x,??(? ??.??, ?? x), ? ??.??))'
!  e.g.       ((02x,05(f 09.04, 03 x), f 09.04))
!             1234567890123456789012345678901234
        CHARACTER (34) :: form

!   EXECUTION:

!   First set appropriate default values.

        form = defform

        IF (PRESENT(unit)) THEN
          actunit = unit
        ELSE
          actunit = defunit
        END IF

        IF (PRESENT(line)) THEN
          actline = line
        ELSE
          actline = defline
        END IF

        IF (PRESENT(indent)) THEN
          actindent = indent
        ELSE
          actindent = defindent
        END IF

        IF (PRESENT(f)) THEN
          actf = f
        ELSE
          actf = deff
        END IF

        IF (PRESENT(w)) THEN
          actw = w
        ELSE
          actw = defw
        END IF

        IF (PRESENT(d)) THEN
          actd = d
        ELSE
          actd = defd
        END IF

        IF (PRESENT(s)) THEN
          acts = s
        ELSE
          acts = defs
        END IF

        IF (PRESENT(a)) THEN
          fa = a
        ELSE
          fa = one
        END IF

        IF (PRESENT(b)) THEN
          fb = b
        ELSE
          fb = one
        END IF

        IF (actindent==0) THEN
          WRITE (form(03:06),90010) blank
        ELSE
          WRITE (form(03:04),90020) ABS(actindent)
        END IF

        WRITE (form(07:08),90020) (actline-ABS(actindent)-actw)/(actw+acts)
        WRITE (form(10:10),90000) actf
        WRITE (form(12:13),90020) actw
        WRITE (form(15:16),90020) actd

        IF (acts/=0) THEN
          WRITE (form(19:20),90020) acts
        ELSE
          WRITE (form(17:22),90010) blank
        END IF

        WRITE (form(26:26),90000) actf
        WRITE (form(28:29),90020) actw
        WRITE (form(31:32),90020) actd

        DO i = 1, actindent
          WRITE (unit,90010,advance='no') blank
        END DO

        IF (PRESENT(name)) WRITE (unit,90010) name

        DO i = 1, SIZE(x,1)
          IF (PRESENT(y)) THEN
            WRITE (unit,form) (fa*x(i,j)+fb*y(i,j),j=1,SIZE(x,2))
          ELSE
            WRITE (unit,form) (fa*x(i,j),j=1,SIZE(x,2))
          END IF
        END DO

        RETURN

! FORMATS:  Also see character string form.

90000   FORMAT (A1)
90010   FORMAT (A)
90020   FORMAT (I2)

      END SUBROUTINE

! Optional
      SUBROUTINE writevec(x,a,b,y,f,w,d,s,unit,name,indent,line)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)

        REAL (stnd), INTENT (IN), OPTIONAL :: y(SIZE(x)), a, b
        INTEGER (short), INTENT (IN), OPTIONAL :: unit, w, d, s, indent, line
        CHARACTER (1), INTENT (IN), OPTIONAL :: f
        CHARACTER (*), INTENT (IN), OPTIONAL :: name

! DESCRIPTION: This functions just as WriteMat, except that it
!              prints a vector instead of a matrix.

! LOCAL DECLARATIONS:

        REAL (stnd) :: fa, fb
        INTEGER (short) :: actunit, actw, actd, acts, actindent, actline, i
        CHARACTER (1) :: actf
        CHARACTER (34), PARAMETER :: defform = &
          '((??x,??(? ??.??, ?? x), ? ??.??))'
!  e.g.       ((02x,05(f 09.04, 03 x), f 09.04))
!             1234567890123456789012345678901234
        CHARACTER (34) :: form

! EXECUTION:

!   First set appropriate default values.

        form = defform

        IF (PRESENT(unit)) THEN
          actunit = unit
        ELSE
          actunit = defunit
        END IF

        IF (PRESENT(line)) THEN
          actline = line
        ELSE
          actline = defline
        END IF

        IF (PRESENT(indent)) THEN
          actindent = indent
        ELSE
          actindent = defindent
        END IF

        IF (PRESENT(f)) THEN
          actf = f
        ELSE
          actf = deff
        END IF

        IF (PRESENT(w)) THEN
          actw = w
        ELSE
          actw = defw
        END IF

        IF (PRESENT(d)) THEN
          actd = d
        ELSE
          actd = defd
        END IF

        IF (PRESENT(s)) THEN
          acts = s
        ELSE
          acts = defs
        END IF

        IF (PRESENT(a)) THEN
          fa = a
        ELSE
          fa = one
        END IF

        IF (PRESENT(b)) THEN
          fb = b
        ELSE
          fb = one
        END IF

        IF (actindent==0) THEN
          WRITE (form(03:06),90010) blank
        ELSE
          WRITE (form(03:04),90020) ABS(actindent)
        END IF

        WRITE (form(07:08),90020) (actline-ABS(actindent)-actw)/(actw+acts)
        WRITE (form(10:10),90000) actf
        WRITE (form(12:13),90020) actw
        WRITE (form(15:16),90020) actd

        IF (acts/=0) THEN
          WRITE (form(19:20),90020) acts
        ELSE
          WRITE (form(17:22),90010) blank
        END IF

        WRITE (form(26:26),90000) actf
        WRITE (form(28:29),90020) actw
        WRITE (form(31:32),90020) actd

        DO i = 1, actindent
          WRITE (unit,90010,advance='no') blank
        END DO

        IF (PRESENT(name)) WRITE (unit,90010) name

        IF (PRESENT(y)) THEN
          WRITE (unit,form) (fa*x(i)+fb*y(i),i=1,SIZE(x))
        ELSE
          WRITE (unit,form) (fa*x(i),i=1,SIZE(x))
        END IF

        RETURN

! FORMATS:  Also see character string form.

90000   FORMAT (A1)
90010   FORMAT (A)
90020   FORMAT (I2)

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> inner_product.f90
    MODULE inner_product
      USE precisions, ONLY : stnd, long, short

! This contains routines for computing the normal inner product of
! two vectors, as given by their dot product, as well as a routine
! for doing a generalized inner product. Without modification, the
! generalized inner product is a duplicate of the normal one.

! An operator is provided for applying each inner product, namely

!       x .GIP. y       for the generalized inner product of x and y
!       x .IP.  y       for the normal inner product of x and y

! The operator .GIP. is in fact overloaded; there must be a second,
! specialized version which computes the generalized inner product
! of the vector x with the vector e[j] = (0,0,...,0,1,0,...,0), i.e.
! with the canonical jth unit vector. It is accessed as

!       x .GIP. j

! For each inner product, the norm ||.|| it induces is also available,
! as defined by

!                                  1/2
! .GNORM. x = ||x||  =  (x .GIP. x)     for the generalized norm

!                                  1/2
! .NORM.  x = ||x||  =  (x .IP.  x)     for the usual norm

! In the case of the usual norm, a special routine is used for
! computing the value with due regard to avoiding unnecessary overflow.
!----------------------------------------------------------------------

! Users who wish may replace the general inner product routine with
! one of their own. In this case, the interface should NOT be changed.
! Only the code WITHIN the functions GenInner, GenInnerJ and GenNorm
! should be modified; i.e. all changes should be restricted to these
! functions.

! It is up to the user to verify that a legal inner product is used.
! Users who do not understand this discussion may safely ignore this;
! anyone who may need this capability will undoubtedly know what this
! is all about.

! If the general inner product routine is altered, the general norm
! function must also be changed, in accordance with the definition above.
! This may be done by activating the line so marked in the routine GenNorm.

! If overflow is likely to be a problem when the generalized norm is computed
! by taking the square root of the inner product, a special general norm
! routine may be written. It is the responsibility of the user.

!----------------------------------------------------------------------


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: OPERATOR (.ip.), OPERATOR (.gip.), OPERATOR (.norm.), &
        OPERATOR (.gnorm.)

      INTERFACE OPERATOR (.gip.)
        MODULE PROCEDURE geninner, geninnerj


      END INTERFACE

      INTERFACE OPERATOR (.gnorm.)
        MODULE PROCEDURE gennorm


      END INTERFACE

      INTERFACE OPERATOR (.ip.)
        MODULE PROCEDURE inner


      END INTERFACE

      INTERFACE OPERATOR (.norm.)
        MODULE PROCEDURE norm2


      END INTERFACE

    CONTAINS

      REAL (stnd) FUNCTION geninner(x,y)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)
        REAL (stnd), INTENT (IN) :: y(SIZE(x))

! DESCRIPTION:

!   This is made available for users wishing to redefine the
!   underlying metric of the space involved by redefining the
!   inner product used for the updating of the quasi-Newton
!   matrices.

! EXECUTION:

! Replace code below this line if required---------------

        geninner = DOT_PRODUCT(x,y)

! Replace code above this line if required---------------

      END FUNCTION


      REAL (stnd) FUNCTION geninnerj(x,j)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)
        INTEGER (long), INTENT (IN) :: j

! DESCRIPTION:

!   This is made available for users wishing to redefine the
!   underlying metric of the space involved by redefining the
!   inner product used for the updating of the quasi-Newton
!   matrices.  Specifically, this routine computes the general
!   inner product of the vector x and the jth unit coordinate
!   vector  [0,...,0,1,0,...,0], where the 1 is in the jth
!   position.  In general, this may *not* just equal  x[j].

! EXECUTION:

! Replace code below this line if required---------------

        geninnerj = x(j)

! Replace code above this line if required---------------

      END FUNCTION


      REAL (stnd) FUNCTION norm2(v)
        USE reals, ONLY : zero, one
        USE num_constants, ONLY : machtiny, machhuge, machtol


        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: v(:)

! DESCRIPTION:

!     This computes the 2-norm (i.e. the Euclidean norm) of the vector  v
!     of length n, with due regard to avoiding overflow and underflow.

!     The routine is based on snrm2 from the blas (in linpack), but this
!     version is written in Fortran 90. It is machine independent.

!     The machine constants MachTiny (the smallest magnitude), MachHuge(the
!     largest magnitude), and MachTol (epsilon) are used to calculate the
!     constants cutlo and cuthi.  Three different cases must be considered
!     when calculating the norm:

!        (1)  All components of v are below cutlo.

!               To avoid underflow, each component is divided by sqrt(min)/n
!               and then the regular Euclidean norm of this modified vector
!               is calculated.  This result is then multiplied by
!               sqrt(min)/n  in order to get the correct value for the norm.

!        (2)  One or more components are greater than cuthi.

!                 To avoid overflow, the same method as in case (1) is used
!                 with a scaling factor of   sqrt(max)*n .

!        (3)  All components are less than cuthi, with at least one
!             component greater than cutlo.

!                 The regular formula for the Euclidean norm is used.

! PARAMETERS:

        INTEGER (short), PARAMETER :: null = 0, small = 1, normal = 2, &
          large = 3

! LOCAL DECLARATIONS:

        INTEGER (long) :: i, n
        INTEGER (short) :: range

        REAL (stnd) :: cutlo, cuthi, summ, xmax

! EXECUTION:

        n = SIZE(v)

        IF (n<=0) THEN
          norm2 = zero
          RETURN
        END IF

        cutlo = SQRT(machtiny/machtol)
        cuthi = SQRT(machhuge)/n

        summ = zero
        range = null

! Evaluate the norm by accumulating a scaled sum of squares and
! adjusting the scaling as numbers of increasingly large magnitude
! are found.

        DO i = 1, n

          SELECT CASE (range)

          CASE (normal)
            IF (ABS(v(i))<cuthi) THEN
              summ = summ + v(i)**2
            ELSE
              range = large
              xmax = ABS(v(i))
              summ = one + (summ/v(i))/v(i)
            END IF

          CASE (small)
            IF (ABS(v(i))<=cutlo) THEN
              IF (ABS(v(i))<=xmax) THEN
                summ = summ + (v(i)/xmax)**2
              ELSE
                summ = one + (xmax/v(i))**2
                xmax = ABS(v(i))
              END IF
            ELSE IF (ABS(v(i))>=cuthi) THEN
              range = large
              xmax = ABS(v(i))
              summ = one + (summ/v(i))/v(i)
            ELSE
              range = normal
              summ = (summ*xmax)*xmax + v(i)**2
            END IF

          CASE (large)
            IF (ABS(v(i))<=xmax) THEN
              summ = summ + (v(i)/xmax)**2
            ELSE
              summ = one + summ*(xmax/v(i))**2
              xmax = ABS(v(i))
            END IF

          CASE (null)
            IF (ABS(v(i))==zero) THEN
!                                 just fall through...
            ELSE IF (ABS(v(i))<=cutlo) THEN
              range = small
              xmax = ABS(v(i))
              summ = one
            ELSE IF (ABS(v(i))>=cuthi) THEN
              range = large
              xmax = ABS(v(i))
              summ = one
            ELSE
              range = normal
              summ = v(i)**2
            END IF

          END SELECT

        END DO

        SELECT CASE (range)
        CASE (normal,null)
          norm2 = SQRT(summ)
        CASE DEFAULT
          norm2 = xmax*SQRT(summ)
        END SELECT

! EXIT:
        RETURN

! FORMATS:  none are defined.

      END FUNCTION

      REAL (stnd) FUNCTION gennorm(x)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)

! DESCRIPTION:

!   This is made available for users wishing to redefine the
!   underlying metric of the space involved by redefining the
!   inner product used for the updating of the quasi-Newton
!   matrices.  It is used in conjunction with GenInner and indeed
!   this routine computes the norm of x associated with the
!   redefined inner product.

! EXECUTION:

! GenNorm = sqrt ( x .GIP. x )  ! Activate this line if desired.

! The above line *defines* the generalized norm; the user may wish
! however to evaluate it in some other fashion, possibly to prevent
! unwanted overflows, as in the routine Norm2.

        gennorm = norm2(x) ! Remove this line if inner product changed.

      END FUNCTION


      REAL (stnd) FUNCTION inner(x,y)

        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)
        REAL (stnd), INTENT (IN) :: y(SIZE(x))

! DESCRIPTION:

!   This is provided so that an operator may be defined for
!   computation of the standard inner product.

! EXECUTION:

        inner = DOT_PRODUCT(x,y)

      END FUNCTION


    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> supp_codes.f90
    MODULE supp_codes
      USE precisions, ONLY : short

! See Module Min_Codes for an explanation of how these codes work.


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: str_evaluation, str_derivative, str_returnstate

!=================
! For Evaluate_f |
!=================

!Str_Evaluation
! Function evaluation
! just evaluate f at x
! evaluate both f and g at x
! just evaluate g at x
      INTEGER (short), PARAMETER, PUBLIC :: justf = 0, both = 1, justg = 2, &
        noop = 3 ! Evaluate_f just calls the
!   user routine and then returns.


!Str_Derivative
! Derivative calculation modes
! Compute from analytic formulae
! Compute both and compare
! Like CompareTest, just iteration 1
      INTEGER (short), PARAMETER, PUBLIC :: analytic = 1, comparetest = 2, &
        firsttest = 3, differences = 4 ! Compute using finite differences

!Str_ReturnState
! Return from Evaluation
! f and/or g succesfully computed.
! User requested abort.
! Function count limit exceeded.
! f could not be computed.
! g could not be computed.
      INTEGER (short), PARAMETER, PUBLIC :: ok = 0, abort = -1, limit = -2, &
        nof = -3, nog = -4, noforg = -5 ! Neither f nor g could be computed.

!================
! For Test_Done |
!================

! type of norm
! use absolute sum norm
! use sqrt(sum of squares) norm
! use max absolute term norm
      INTEGER (short), PARAMETER, PUBLIC :: l1 = 1, l2 = 2, linf = 3, g2 = 4 ! use sqrt(x .GIP. x) norm


    CONTAINS


      FUNCTION str_evaluation(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (5), PARAMETER :: str(justf:noop) = (/ 'justF', 'both ', &
          'justG', 'noOp '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_evaluation

! DESCRIPTION: This function returns the specified value for the evaluation
!              to be done.

! EXECUTION:
        str_evaluation = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_derivative(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (11), PARAMETER :: str(analytic:differences) = (/ &
          'Analytic   ', 'CompareTest', 'FirstTest  ', 'Differences'/)

        CHARACTER (LEN_TRIM(str(i))) :: str_derivative

! DESCRIPTION: This function returns the specified value for the derivative
!              method.

! EXECUTION:
        str_derivative = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_returnstate(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (6), PARAMETER :: str(noforg:ok) = (/ 'NoForG', 'NoG   ', &
          'NoF   ', 'Limit ', 'Abort ', 'OK    '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_returnstate

! DESCRIPTION: This function returns the specified value for the return
!              condition from the function evaluation.

! EXECUTION:
        str_returnstate = str(i)
        RETURN

      END FUNCTION


    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> supp_states.f90
    MODULE supp_states

      USE precisions, ONLY : stnd, long, short

      IMPLICIT NONE
      PUBLIC
!------              -------------!

!=================
! For Evaluate_f |
!=================

      TYPE :: evalstate
        INTEGER (short) :: scalef, derivatives, expense, evallimit, &
          evaltraceunit
        LOGICAL :: tracef, traceg, tracedervtest
      END TYPE

      TYPE :: evalcounts
        INTEGER (long) :: fevals, gevals
        REAL (stnd) :: time
      END TYPE

      TYPE :: evalerrors
        REAL (stnd) :: worst, average, accumsum
        INTEGER (long) :: gradcnt, index, accumnumber
      END TYPE

!====================
! For Print_Iterate |
!====================

      TYPE :: printstate
        INTEGER (long) :: frequency
        INTEGER (short) :: printunit
        LOGICAL :: printx, printgrad
      END TYPE

      TYPE :: lastprint
        INTEGER (short) :: unit
        INTEGER (long) :: iter
        LOGICAL :: prx
        LOGICAL :: prg
      END TYPE

      TYPE :: printvalues
        REAL (stnd) :: time
        INTEGER (long) :: nextpoint
        TYPE (lastprint) :: last
      END TYPE

!================
! For Test_Done |
!================

      TYPE :: termstate
        INTEGER (short) :: thenorm
        LOGICAL :: traceterm
        INTEGER (short) :: termtraceunit
        LOGICAL :: usegrad, usestep, useshanno, usefunc
        REAL (stnd) :: fatx0, normgatx0
      END TYPE

      TYPE :: termvalues
        REAL (stnd) :: normgsq, normxsq, diffsq
      END TYPE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> supp_defs.f90
    MODULE supp_defs ! Default values for support routines
      USE reals, ONLY : zero, one
      USE true_false, ONLY : false, true
      USE supp_codes, ONLY : analytic, l2
      USE supp_states, ONLY : evalstate, lastprint, printvalues, evalerrors, &
        evalcounts, termstate, printstate
      USE precisions, ONLY : short, stnd, long


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: defaultevalstate, evalstate, defaultprintstate, printstate, &
        defaulttermstate, termstate, initevalcts, evalcounts, initevalers, &
        evalerrors, initprvals, printvalues

!=================
! For Evaluate_f |
!=================

! Default values
      INTEGER (short), PARAMETER, PRIVATE :: scalef = 0, &
        derivatives = analytic, expense = 1, evaltraceunit = 6

      INTEGER (long), PARAMETER, PRIVATE :: evallimit = 300

      LOGICAL, PARAMETER, PRIVATE :: tracef = false, traceg = false, &
        tracedervtest = false

      TYPE (evalstate) :: defaultevalstate = evalstate(scalef,derivatives, &
        expense,evallimit,evaltraceunit,tracef,traceg,tracedervtest)


! Initial values
      INTEGER (long), PARAMETER, PRIVATE :: fevals = 0, gevals = 0

      REAL (stnd), PARAMETER, PRIVATE :: etime = zero

      TYPE (evalcounts) :: initevalcts = evalcounts(fevals,gevals,etime)


! Initial values for derivative est.
      REAL (stnd), PARAMETER, PRIVATE :: worst = zero, average = zero, &
        accumsum = zero

      INTEGER (long), PARAMETER, PRIVATE :: iteration = 0, index = 0, &
        accumnumber = 0

      TYPE (evalerrors) :: initevalers = evalerrors(worst,average,accumsum, &
        iteration,index,accumnumber)

!====================
! For Print_Iterate |
!====================

! Default Values
      INTEGER (long), PARAMETER, PRIVATE :: frequency = 10000

! Default Values
      INTEGER (short), PARAMETER, PRIVATE :: printunit = 6

      LOGICAL, PARAMETER, PRIVATE :: printx = false, printgrad = false

      TYPE (printstate) :: defaultprintstate = printstate(frequency,printunit, &
        printx,printgrad)


! Initial Values
      REAL (stnd), PARAMETER, PRIVATE :: time = zero

      INTEGER (long), PARAMETER, PRIVATE :: nextpoint = 0

      TYPE (printvalues) :: initprvals = printvalues(time,nextpoint, &
        lastprint(printunit,-1,false,false))

!================
! For Test_Done |
!================

! Default Values
      INTEGER (short), PARAMETER, PRIVATE :: thenorm = l2, termtraceunit = 6

      LOGICAL, PARAMETER, PRIVATE :: traceterm = false, usegrad = true, &
        usestep = false, useshanno = true, usefunc = false

      REAL (stnd), PARAMETER, PRIVATE :: fatx0 = one, normgatx0 = one

      TYPE (termstate) :: defaulttermstate = termstate(thenorm,traceterm, &
        termtraceunit,usegrad,usestep,useshanno,usefunc,fatx0,normgatx0)

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> fscale.f90
    MODULE fscale ! 'USE'ed in modules Support and EvaluateF
      USE reals, ONLY : one, three, two
      USE supp_codes, ONLY : justf, justg, both
      USE true_false, ONLY : true
      USE precisions, ONLY : stnd, short


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: f_scale

    CONTAINS

      SUBROUTINE f_scale(ft,scale,fscaled,gscale,which)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: ft
        INTEGER (short), INTENT (IN) :: scale
        REAL (stnd), INTENT (OUT) :: fscaled, gscale

        INTEGER (short), INTENT (IN), OPTIONAL :: which

! DESCRIPTION:

!   This subroutine applies one of several scalings (linear or nonlinear)
!   to a function, thereby modifying the function and/or gradient value.

!   On Entry:

!           Ft    - the present function value

!           Scale - the code for the type of scale desired, where the scaling
!                   function is one of the following:

!                     1:   ff(z) = 1 + z
!                     2:   ff(z) = z*z
!                     3:   ff(z) = -1 / (1 + z*z)
!                     4:   ff(z) =  sqrt(1 + z*z)
!                     5:   ff(z) = z*z*z

!           which - See Evaluate_F. If it is omitted, the value both is assumed.

!   On Exit:
!           FScaled - the scaled function value.
!           GScale  - gradient scaling factor.

! LOCAL DECLARATIONS:

        LOGICAL :: dof, dog

! EXECUTION:

        IF ( .NOT. PRESENT(which)) THEN
          dof = true
          dog = true
        ELSE
          dof = which == justf .OR. which == both
          dog = which == justg .OR. which == both
        END IF

        SELECT CASE (INT(scale))

        CASE (1) ! ff(z) = 1 + f(z)

          IF (dof) fscaled = ft + one
          IF (dog) gscale = one

        CASE (2) ! ff(z) = z*z

          IF (dof) fscaled = ft*ft
          IF (dog) gscale = two*ft

        CASE (3) ! ff(z) = -1/(1+z**2)

          IF (dof) fscaled = -one/(one+ft**2)
          IF (dog) gscale = two*ft*(one/(one+ft**2))**2

        CASE (4) ! ff(z) = sqrt(1+z**2)

          IF (dof) fscaled = SQRT(one+ft**2)
          IF (dog) gscale = ft/SQRT(one+ft**2)

        CASE (5) ! ff(z) = z*z*z

          IF (dof) fscaled = ft*ft*ft
          IF (dog) gscale = three*ft*ft

        CASE DEFAULT ! ff(z) = z

          IF (dof) fscaled = ft
          IF (dog) gscale = one

        END SELECT

        RETURN

! FORMATS:  none.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> evaluatef.f90
    MODULE evaluatef ! 'USE'ed in module Support
      USE general, ONLY : indent, writevec, cpusecs
      USE fscale, ONLY : f_scale
      USE reals, ONLY : zero, one, c100
      USE num_constants, ONLY : machtol
      USE supp_codes, ONLY : justf, noforg, nog, nof, abort, limit, ok, &
        analytic, comparetest, firsttest, both, justg
      USE supp_defs, ONLY : initevalcts, initevalers
      USE true_false, ONLY : false
      USE precisions, ONLY : stnd, long, short
      USE supp_states, ONLY : evalstate, evalerrors, evalcounts



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: evaluate_f, evalstate, evalcounts, evalerrors

    CONTAINS

! optional
      RECURSIVE SUBROUTINE evaluate_f(fname,x,fx,g,job,ev,cts,first,test, &
          level)


        INTERFACE
          SUBROUTINE fname(x,fx,g,job)
            USE precisions, ONLY : stnd, short
            REAL (stnd), INTENT (IN) :: x(:)
            REAL (stnd), INTENT (OUT) :: fx, g(SIZE(x))
            INTEGER (short), INTENT (INOUT) :: job
          END SUBROUTINE
        END INTERFACE

! ARGUMENTS:

        REAL (stnd), INTENT (INOUT) :: x(:), fx
        REAL (stnd), INTENT (OUT) :: g(SIZE(x))

        INTEGER (short), INTENT (INOUT) :: job

        TYPE (evalstate), INTENT (IN) :: ev
        TYPE (evalcounts), INTENT (INOUT) :: cts

        LOGICAL, INTENT (IN), OPTIONAL :: first
        TYPE (evalerrors), INTENT (INOUT), OPTIONAL :: test
        INTEGER (short), INTENT (IN), OPTIONAL :: level

! DESCRIPTION:

!   This subroutine evaluates a test function at the given point x.  It returns
!   the value of the function and/or the value of the gradient at x.  It allows
!   the application of a nonlinear scaling to the function if desired (see
!   ScaleF below). It also allows the use of finite differences (see Derivatives
!   below). It can also act as a noop, i.e. as a do nothing routine (see
!   job below).

!  On Entry:
!  --------

!   First  This must be present and true each time a new problem is begun.
!          It causes initialization of counts described below to be done.

!   FName  is the name of the function to evaluate.  There must be a subroutine
!          provided of the form

!                    Subroutine FName(x,fx,g,job)

!          (x, fx, g and job mean the same as in this subroutine.)
!          The precise interface definition for FName appears above.

!   x      contains the value of the n-coordinates x(1),...,x(n) at which to
!          evaluate the function.

!   job =  justf   only evaluate the function.
!       =  both    evaluate both.
!       =  justg   only evaluate the gradient.
!       =  noop    actually, if job has any value other than one of the first
!                  three, then just call FName with this same code for job;
!                  i.e. Evaluate_f should do nothing. This is intended for
!                  the possible convenience of the writer of FName.

!   Ev     This record contains a number of control fields. These are explained
!          below.

!   level  If present, indent trace output by this number of spaces.

!   Cts    See below.

!   Test   See Finite Difference Error Estimates below.

!  On Exit:
!  --------

!   fx     contains the function value (with the scaling applied if required).

!   g      contains the gradient value (with the scaling applied if required).

!   job = OK      The request made on the call was completed satisfactorily.
!                 fx and/or g are available as requested.
!         Abort   The minimization routine which called Evaluate_f is hereby
!                 requested to exit immediately to the routine which called
!                 it.  This can be used by the routine FName to trigger
!                 premature termination due to circumstances of which the
!                 minimization routine may not be aware.
!         Limit   Terminate the minimization; the preset limit on the number
!                 of function evaluations allowed has been exceeded.
!         NoF     The function value could not be determined.
!         NoG     The gradient value could not be determined.
!         NoForG  Neither fx nor g could be evaluated.

!   Cts     See below.

!  Control values in Ev:   Ev contains the following fields (amongst others).
!  -------------------

!   ScaleF      controls the nonlinear scaling of FName.

!          = 0     no effect.

!          = k > 0 This routine computes and returns  ff( FName(x) ), where ff
!                  is the  k-th  of the nonlinear functions of one variable
!                  defined in the routine F_Scale.

!                Note that for certain scalings, the function value is needed to
!                scale the gradient vector. So, in this case,  if you call
!                Evaluate_f just for a gradient value, it is necessary to supply
!                the function value in fx as well, in order to do the scaling.

!   Derivatives This specifies the method by which derivatives are to be
!               computed, when requested. The choice is between:

!      Analytic Use analytic formulae which must be coded and available in the
!               user routine FName.

!      Differences Use finite difference approximations. In this case, the user
!               routine FName may ignore calls with job /= justf, and need only
!               be able to compute function values.  Further comments appear in
!               the discussion of finite difference computations (below).

!      CompareTest  In this case both analytic and finite differences are
!               computed.  They are then compared and a record is kept to see to
!               what extent they disagree. A record of the level of agreement is
!               available through the optional argument Test.  A more complete
!               description is also given below.

!      FirstTest   This case is precisely the same as for CompareTest, with
!               the sole exception that the testing only takes place on the
!               first call to Evaluate_f.

!   EvalLimit   The maximum value allowed for the count Cts%FEvals.

!         = 0      On entry specifies no maximum, i.e. any number of function
!                  evaluations may be done. This is inadvisable.

!         = k > 0  Specifies the maximum number of times that FName may be
!                called.  If the function evaluation count in FEvals is greater
!                than or equal to  EvalLimit on entry to Evaluate_f, then the
!                function is not evaluated and the return code job is set as
!                above. Note that the count in FEvals does NOT include function
!                evaluations used for computing finite difference gradients.

!   Expense     This may be used to artificially increase the expense of a
!               function evaluation. If Expense is >1, then the function is
!               simply evaluated Expense times. Of course, Expense-1 of these
!               are totally redundant. This can be used for testing purposes.

!   It is possible to trace parts of the execution of Evaluate_f.

!   TraceF        if true, the function value is to be printed
!   TraceG        if true, the gradient value is to be printed
!   TraceDervTest if true, the function values are to be printed during testing
!                          mode

!   EvalTraceUnit  the unit number for output of computed values

!   Note that an error message is printed when the maximum number of function
!   evaluations is exceeded, provided either TraceF or TraceG is true.

!  Counts available on exit.
!  ------------------------

!   The following counts are available in the record Cts on exit.  As noted
!   above, Evaluate_f must be called once with First = true to ensure that these
!   are correctly initialized.

!   FEvals  the number of calls to evaluate the function, i.e. calls with job
!           = justf or both. Note that calls to evaluate the function for the
!           purpose of doing finite differences are not included in this count.

!   GEvals  the number of calls to evaluate the gradient, i.e. calls with
!           job = justg or both.

!   time    the amount of cpu time spent in  Evaluate_f.  The time used in
!           the final scaling is included in the timing.  Timing commences on
!           entry to Evaluate_f, and ends just before return from Evaluate_f.

!  Finite Difference Error Estimates
!  ---------------------------------

!   If the derivative mode is CompareTest, then the accuracy of finite
!   difference/analytic derivative calculations is monitored. On return, if
!   Test is present, one may find in the fields of Test:

!      Worst      The worst error which occurred.
!      Average    Estimate of the average number of figures of agreement.
!      GradCnt    Count of gradient evaluations where the worst error occurred.
!      Index      Component of the gradient where the worst error occurred.

!   To be specific, when in test mode, each component of the analytic derivative
!   is computed, and these are returned in g as the gradient.  As well, for each
!   component, a finite difference approximation is computed (as described
!   below) and the relative difference between that and the analytic component
!   is determined.  This quantity is monitored, and the largest such value is
!   recorded. In addition, Index records in which component of that gradient the
!   error occurred, and GradCnt tells which gradient evaluation was in
!   progress when the error occurred; i.e. GradCnt just records the current
!   value of Cts%GEvals.  Note that Index and GradCnt only refer to the point
!   at which the largest error occurred.

!   If the function and gradient evaluations are correct, one would normally
!   expect the relative error to be of the order of 10**-(t/2), where  t  is the
!   number of figures of relative accuracy of the machine in use.  However, as
!   the minimum is approached and the gradient components generally become very
!   small, this relative accuracy may be much worse than expected. Therefore we
!   also maintain an estimate of the average agreement. Here, for each component
!   of each gradient computation, we compute the base 10 log of the relative
!   accuracy; this is roughly the number of significant figures of agreement
!   between the two values. This quantity is monitored and  Average is returned
!   as the average value of the number of significant figures of agreement.

!   When function and gradient computations are correct, Worst will generally be
!   at least as small as  10**(-t/2), although it can be more like 10**(-t/4).
!   Gross blunders will usually give the worst error a value very near to 1, but
!   not always.  If all is well, Average will usually be about t/2; blunders
!   will often result in Average being near 0 or 1.

!   Finite difference computations

!   For first derivatives, simple forward differences are used.
!   To estimate the i-th component of the gradient of f, we compute

!                 ( f(x + h*e[i]) - f(x) ) / h,

!   where h = eps * abs(x[i]).  When x[i] = 0, we just choose h = eps.
!   Here eps is sqrt(MachTol), where MachTol is the machine accuracy. This is
!   used when Derivatives = Differences, CompareTest or FirstTest.

!   When Derivatives = CompareTest, more information is required; thus we also
!   compute f(x + sqrt(h)*e[i]).  This means that when in test mode, twice as
!   many function evaluations are needed.  This is required to eliminate
!   scaling effects in the estimate of figures of agreement.

! SUBROUTINES:
!               FName    ...the user routine.
!               F_Scale  ...performs scaling.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'Eval'

! LOCAL DECLARATIONS:

        LOGICAL :: firsttime, docompare, dof, dog

        INTEGER (short) :: rememberbad, dervs, job_return, calls
        INTEGER (long) :: count, kk, zdcnt, zindex, zgcnt

        REAL (stnd) :: ft, fv, tt, scale, rh, zerr, zserr, fval, fval2, h, &
          terr

! EXECUTION:

        IF (ev%tracef .OR. ev%traceg) THEN
          CALL indent(id,ev%evaltraceunit,level)
          WRITE (ev%evaltraceunit,90000) job, justf, justg, both
        END IF

        IF (job/=justf .AND. job/=justg .AND. job/=both) THEN ! noop call.
          CALL fname(x,fx,g,job)
          GO TO 10
        END IF

        dervs = ev%derivatives

        IF (PRESENT(first)) THEN
          firsttime = first
        ELSE
          firsttime = false
        END IF

        IF (firsttime) THEN
          cts = initevalcts
          IF (PRESENT(test)) test = initevalers
          IF (ev%derivatives==firsttest) dervs = comparetest
        ELSE
          IF (ev%derivatives==firsttest) dervs = analytic

          IF (ev%evallimit>0 .AND. cts%fevals>=ev%evallimit) THEN
            GO TO 20
          END IF
        END IF

        dof = job == justf .OR. job == both
        dog = job == justg .OR. job == both

        cts%time = cts%time - cpusecs()

        IF (dof) cts%fevals = cts%fevals + 1
        IF (dog) cts%gevals = cts%gevals + 1

!-----First compute required function and/or gradient values.

        docompare = PRESENT(test) .AND. dervs == comparetest

        IF (docompare) THEN
          zerr = test%worst
          zserr = test%accumsum
          zdcnt = test%accumnumber
          zindex = test%index
          zgcnt = test%gradcnt
        END IF
        rememberbad = ok

EXPENSE: DO kk = 1, ev%expense ! Repeat if required to simulate expensive call.

          IF (dervs==analytic .OR. .NOT. dog) THEN
            calls = 0
          ELSE
            calls = SIZE(x) ! extra calls to FName needed.
          END IF

          job_return = job ! First compute  f(x)  and/or  g(x).
          CALL fname(x,fval,g,job_return)

          SELECT CASE (job_return)
          CASE (limit,abort,nof,nog,noforg)
            rememberbad = job_return
            CYCLE
          END SELECT

          IF (job==justg) THEN
            IF (scale/=0) ft = fx
          ELSE
            ft = fval
          END IF

NCALLS:   DO count = 1, calls !Do extra calls, if required: function values only.

            tt = x(count) ! save x

            h = SQRT(machtol)
            IF (tt/=zero) THEN ! finite difference step
              h = h*ABS(tt)
            END IF

            x(count) = tt + h ! Compute  f( x + h * e[count] )

            job_return = justf
            CALL fname(x,fval,g,job_return)

            x(count) = tt ! restore x

            SELECT CASE (job_return)
            CASE (limit,abort,nof,nog,noforg)
              rememberbad = job_return
              CYCLE EXPENSE
            END SELECT

DOCOMPAR:   IF (docompare) THEN

              IF (ev%tracedervtest) THEN
                CALL indent(id,ev%evaltraceunit,level)
                WRITE (ev%evaltraceunit,90030) g(count), count, (fval-ft)/h
              END IF

! Estimate error, and leave computed analytic gradients in g.
! Use f at x + a * e[count], for a = h and sqrt(h).

              rh = SQRT(h)
              x(count) = tt + rh

              job_return = justf
              CALL fname(x,fval2,g,job_return)
              x(count) = tt

              SELECT CASE (job_return)
              CASE (limit,abort,nof,nog,noforg)
                rememberbad = job_return
                CYCLE EXPENSE
              END SELECT

              IF (ABS(fval2-ft)>c100*machtol*ABS(ft)) THEN
                terr = (fval-ft-h*g(count))/(fval2-ft-rh*g(count))
                IF (tt>one) terr = terr/tt

                terr = MAX(MIN(one,ABS(terr)),machtol) ! in [MachTol,1]

                zserr = zserr - LOG10(terr) ! Number of figures agreement.
                zdcnt = zdcnt + 1

                IF (ev%tracedervtest) THEN
                  CALL indent(id,ev%evaltraceunit,level)
                  WRITE (ev%evaltraceunit,90020) terr, -LOG10(terr)
                END IF
                IF (terr>ABS(zerr)) THEN
                  zindex = count
                  zgcnt = cts%gevals
                  zerr = SIGN(terr,zerr)
                END IF
              ELSE
                zerr = -ABS(zerr) ! Flag case with excessive cancellation.
                IF (ev%tracedervtest) THEN
                  CALL indent(id,ev%evaltraceunit,level)
                  WRITE (ev%evaltraceunit,90010)
                END IF
              END IF
            ELSE ! DOCOMPAR  Estimate gradients with forward finite difference
              g(count) = (fval-ft)/h
            END IF DOCOMPAR

          END DO NCALLS

          IF (ev%scalef/=0) THEN ! Do scaling: define Fv and scale.
            CALL f_scale(ft,ev%scalef,fv,scale,job)
          ELSE
            fv = ft
            scale = one
          END IF

          IF (dof) THEN ! Revise function and gradient as necessary.
            fx = fv
          END IF

          IF (dog .AND. scale/=one) THEN
            g = scale*g
          END IF

        END DO EXPENSE

        IF (docompare) THEN
          test%worst = zerr
          test%index = zindex
          test%gradcnt = zgcnt
          test%accumnumber = zdcnt
          test%accumsum = zserr
          test%average = zserr/zdcnt
        END IF

        job = rememberbad

! EXIT:
        cts%time = cts%time + cpusecs()

        IF (ev%tracef .AND. job==ok) THEN
          CALL indent(id,ev%evaltraceunit,level)
          WRITE (ev%evaltraceunit,90040) fx
        END IF

        IF (ev%traceg .AND. job==ok) THEN
          CALL writevec(g,f='g',w=15_short,d=8_short,unit=ev%evaltraceunit, &
            name='gradient',indent=level)
        END IF

10      RETURN !  NoOp

20      IF (ev%tracef .OR. ev%traceg) THEN ! max. number of func. evals exceeded.
          CALL indent(id,ev%evaltraceunit,level)
          WRITE (ev%evaltraceunit,90050)
        END IF
        job = limit
        RETURN

! FORMATS:

90000   FORMAT (' job (f,g,fg)=',4I3)
90010   FORMAT (' Excessive error in gradient estimation.')
90020   FORMAT (' Error estimate in gradient estimation: ', &
          G15.7/' Estimated figures of agreement:        ',G9.2)
90030   FORMAT (' Analytic gradient   ',G22.15,' (component ',I3, &
          ')'/' Estimated derivative',G22.15)
90040   FORMAT (' Function = ',G26.16)
90050   FORMAT (/ &
          ' The number of function evaluations allowed has been exceeded.')

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> print.f90
    MODULE print ! 'USE'ed in module Support
      USE general, ONLY : cpusecs, writevec, indent
      USE true_false, ONLY : false
      USE precisions, ONLY : long, short, stnd
      USE supp_states, ONLY : evalcounts, printvalues, printstate


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: print_iterate, evalcounts, printvalues, printstate

    CONTAINS

! OPTIONAL
      SUBROUTINE print_iterate(n,iteration,fx,nrmg,counts,c,v,x,g,force,first, &
          copy,copyunit,level)


! ARGUMENTS:

        INTEGER (long), INTENT (IN) :: n, iteration

        REAL (stnd), INTENT (IN) :: fx, nrmg

        TYPE (evalcounts), INTENT (IN) :: counts
        TYPE (printstate), INTENT (IN) :: c
        TYPE (printvalues), INTENT (INOUT) :: v


        REAL (stnd), INTENT (IN), OPTIONAL :: x(n), g(n)
        LOGICAL, INTENT (IN), OPTIONAL :: force, first

        LOGICAL, INTENT (IN), OPTIONAL :: copy
        INTEGER (short), INTENT (IN), OPTIONAL :: copyunit
        INTEGER (short), INTENT (IN), OPTIONAL :: level

! DESCRIPTION:

!   This routine prints the value Fx of some function, and the norm nrmg of
!   its gradient, at some point x, along with (optionally) the point x, and
!   (optionally) the value of the gradient g at the point x.

!   Arguments.

!   n           the dimension of the problem
!   iteration   current iteration number
!   Fx          the function value at x
!   nrmg        the norm of the gradient at x
!   Counts      function and gradient evaluation counts (from Evaluate_f)
!   C           the control record
!                   PrintUnit   unit for output
!                   frequency   interval between output
!                   PrintX      true to include point x in output
!                   PrintGrad   true to include gradient g in output
!   V           current data; needed for next call
!                   time        cumulative time herein
!                   NextPoint   next point to print at
!                   last..iter  iteration no of last point printed
!                         unit  unit used for last print
!                         prx   value of printx from last print
!                         prg   value of printgrad from last print

!   Optional Arguments (see below).

!   x        the current point
!   g        the gradient value at x
!   force    force output
!   first    this is first call
!   copy     duplicate output wanted
!   copyunit duplicate unit
!   level    if present, indent trace output by this number of spaces


!  1. Control is under C%frequency.

!     frequency  <= 0   there is no output.

!     frequency  >  0   print every frequency-th iteration:

!               the iteration number in iteration
!               the function value in Fx
!               the no. of function evaluations in Counts@FEvals
!               the norm of the gradient in nrmg
!               the no. of gradient evaluations in Counts@GEvals
!               the function/gradient evaluation time
!               the point x (if C%PrintX  is true and x is present)
!               the gradient g (if C%PrintGrad is true and g is present)

!  2. V%NextPoint records the number of the next iteration at which to print.
!     If, on entry, iteration = NextPoint, generate the output.
!     If on entry the value of iteration is already greater than that of
!     nextpoint, nextpoint is repeatedly incremented by frequency until
!     that is no longer true. It is initialized to 0 on a call with first
!     present and true.

!  3. Setting force to true will cause the current point to be printed, even if
!     iteration does not match the value of NextPoint.  This is useful for
!     forcing printing of the final point reached.

!  4. The routine is careful not to repeat a printing request.  If the output
!     unit or the status of C%PrintX  or C%PrintGrad or the iteration count is
!     different, then the printing is done; otherwise it is considered a
!     repeat of a previous request and it is ignored.  Information about the
!     previous iterate is stored in V%last and should not be altered by
!     the calling routine.

!  5. V%time is used for accumulating the time spent in the print routine.
!     It is started from zero on a call with first present and true,
!     and each call increments time by the time spent in the routine.

!  6. If copy if present and true and copyunit is present,
!     all output is duplicated on the designated unit.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'Prnt'

! LOCAL DECLARATIONS:

        LOGICAL :: isfirst, forceit, docopy

        REAL (stnd) :: secs

! EXECUTION:

        secs = cpusecs()

        IF (PRESENT(first)) THEN
          isfirst = first
        ELSE
          isfirst = false
        END IF

        IF (PRESENT(force)) THEN
          forceit = force
        ELSE
          forceit = false
        END IF

        IF (PRESENT(copy)) THEN
          docopy = copy .AND. PRESENT(copyunit)
        ELSE
          docopy = false
        END IF

        IF (isfirst) THEN
          v%nextpoint = 0
          v%time = -secs
        ELSE
          v%time = v%time - secs
        END IF

! Need to print?
PRINT:  IF (isfirst .OR. (iteration/=v%last%iter) .OR. (c%printunit/=v%last% &
            unit) .OR. (c%printx .NEQV. v%last%prx) .OR. &
            (c%printgrad .NEQV. v%last%prg)) THEN

          IF (c%frequency>0) THEN
            DO WHILE (v%nextpoint<iteration)
              v%nextpoint = v%nextpoint + c%frequency
            END DO
          END IF

          IF ((c%printunit/=0) .AND. (c%frequency>0) .AND. (forceit .OR. ( &
              iteration==v%nextpoint))) THEN

!        Save information defining this print request.
            v%last%iter = iteration
            v%last%unit = c%printunit
            v%last%prx = c%printx
            v%last%prg = c%printgrad

! Print iteration number, function value, norm of g, and
!       number of function/gradient evaluations.

            CALL indent(id,c%printunit,level)
            WRITE (c%printunit,90000) iteration, fx, counts%fevals, nrmg, &
              counts%gevals, counts%time
            IF (docopy) THEN
              CALL indent(id,copyunit,level)
              WRITE (copyunit,90000) iteration, fx, counts%fevals, nrmg, &
                counts%gevals, counts%time
            END IF

! if appropriate, also print x and g.

            IF (PRESENT(x) .AND. c%printx) THEN
              CALL writevec(x,f='g',unit=c%printunit,name='['//id//'] x->', &
                indent=level)
              IF (docopy) THEN
                CALL writevec(x,f='g',unit=copyunit,name='['//id//'] x->', &
                  indent=level)
              END IF
            END IF

            IF (PRESENT(g) .AND. c%printgrad) THEN
              CALL writevec(g,f='g',unit=c%printunit,name='['//id//'] g->', &
                indent=level)
              IF (docopy) THEN
                CALL writevec(g,f='g',unit=copyunit,name='['//id//'] g->', &
                  indent=level)
              END IF
            END IF
          END IF

! Update counter.
          IF (iteration==v%nextpoint .AND. c%frequency>0) &
            v%nextpoint = v%nextpoint + c%frequency

        END IF PRINT

! EXIT:

        v%time = v%time + cpusecs()

! FORMATS:

!999 format(' ',' ...pt ',i3,'; f=',g23.16,'(#',i3,') ||g||=',        &
!                           e7.2, '(#',i3,'); '/,17x,f8.3,' secs'    )
90000   FORMAT (' ','pt ',I3,'; f=',G23.16,'(#',I3,') ||g||=',E7.2,'(#',I3, &
          '); ',F8.3,' s')

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> testdone.f90
    MODULE testdone ! 'USE'ed in module Support
      USE reals, ONLY : zero, one
      USE general, ONLY : indent
      USE true_false, ONLY : true, false
      USE supp_codes, ONLY : l1, g2, linf, l2
      USE precisions, ONLY : stnd, short
      USE supp_states, ONLY : termstate, termvalues
      USE inner_product, ONLY : OPERATOR (.ip.), OPERATOR (.gip.)


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: test_done, termstate, termvalues

    CONTAINS

      LOGICAL FUNCTION test_done(eps,c,fx,g,xi,xim1,v,level)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: eps
        TYPE (termstate), INTENT (IN) :: c
        REAL (stnd), INTENT (IN), OPTIONAL :: fx, g(:), xi(:), xim1(:)
        TYPE (termvalues), INTENT (OUT), OPTIONAL, TARGET :: v
        INTEGER (short), INTENT (IN), OPTIONAL :: level

! DESCRIPTION:

!   This routine is used to test whether or not to terminate a minimization
!   routine.  It provides a means of using uniform criteria for different
!   routines.  A choice of criteria is provided, according to values which are
!   passed in the the state structure C.

!   Note that in one case, one test cannot be applied; i.e. on the first point,
!   since a pair of successive points is required, one can not use the test on
!   the stepsize.  In this case, the argument xim1 must NOT be passed to the
!   routine, and the function value is always returned as false if one applies
!   the test which looks at more than one point.

!   Note that Test_Done will always return a true value if the norm of the
!   gradient computes to 0, since in this case, no further steps can be taken
!   anyway.

!   Here ||u|| denotes the appropriate norm of the vector u.

!   n     is the length of the vectors.

!   Some of the actions of the routine are determined by values set in the
!   state structure C. Appropriate integer codes are defined in Supp_Codes.
!   The fields of C have the following effects:

!     TheNorm  = L1    use the   L1      (abssum) norm of vectors.
!              = L2    use the   L2   (Euclidean) norm of vectors.
!              = Linf  use the maximum (infinity) norm of vectors.
!              = G2    use the generalized        norm of vectors.

!     TraceTerm  If the trace argument is set to true, then the result of
!            each test will be printed on the trace unit.  Note that when
!            several tests are being applied, the trace will show each of the
!            possible test values, whereas normally an exit occurs as soon as
!            any test fails.

!   Tests selected: there are currently four. They are applied independently,
!   and all tests chosen must be passed. The arguments required, depending
!   on the tests in use, must be present.

!   The tests may be scaled relative to certain values, normally the
!   value of the function and gradient at the initial point. See below.

!      UseGrad     Test if the appropriate norm of g is <= eps.  This first
!                  type of test is most commonly used to see if the gradient is
!                  sufficiently small.  Thus the test applied is

!                           ||g|| <= eps * NormGatX0

!      UseStep     Test if the appropriate norm of the difference between xi and
!                  xim1 is <= eps.  The test is absolute if the norm of xi is
!                  less than 1, and relative otherwise.  This type of test is
!                  normally used to test the distance between successive points.
!                  Thus the test is

!                       ||xi-xim1|| <= eps * max(1,||xi||)

!      UseShanno   Use a test appearing in Shanno's CONMIN using x and g.
!                  Termination is indicated when

!                              ||g||
!                           ------------  <=   eps * NormGatX0
!                           max(1,||xi||)

!      UseFunc     Terminate if the function value is sufficiently small.  This
!                  test would normally only be used in a relative manner.
!                  Thus the test is

!                         |Fx| <= eps * FatX0

!   It is often desired to make termination tests relative to the function
!   and/or gradient values at the initial point.  In the tests above, the values
!   FatX0  and  NormGatX0 are used; these may be thought of as the function
!   value at X0, along with the norm of the gradient at that point.  The values
!   for FatX0 and NormGatX0 must be set by the user before calling Test_Done.
!   If relative tests are not desired, these values should be set to 1. These
!   values are fields in C.

!   Other points to note are:

!   Some tests are actually done by comparing the squares of the norms against
!   eps**2.  Thus it is possible that this version of this routine might
!   generate an unwanted overflow or underflow.

!   Note:  Neither  g  nor xi nor xim1 is altered by this routine.  Only those
!          vectors used in the test are actually referenced. For example, if the
!          only test to be applied is UseGrad, then neither xi nor xim1 is
!          referenced, and they can be omitted from the call.  Vectors which
!          should be present but are not are treated as zero.

!   On Return:  If the desired tests are *all* passed, then Test_Done is set
!               to true; otherwise it is set to false.

!               In addition, the vector norms computed during application of the
!               tests are available.  Of course, only those which were actually
!               computed in applying the desired tests will be defined. These
!               are returned in the fields of V, if V is present.

!               We have specifically:

!        gsq    norm squared of  g,       gsq    = ||g||**2
!        xsq    norm squared of  xi,      xsq    = ||xi||**2
!        diffsq norm squared of  xi-xim1, diffsq = ||xi-xim1||**2

!   level    if present, indent trace output by this number of spaces

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'Term'

! LOCAL DECLARATIONS:

        LOGICAL :: pass, passthis

        REAL (stnd), POINTER :: gsq, xsq, diffsq
        REAL (stnd), TARGET :: locgsq, locxsq, locdiffsq

! EXECUTION:

        IF (PRESENT(v)) THEN ! return values if V present in call
          gsq => v%normgsq
          xsq => v%normxsq
          diffsq => v%diffsq
        ELSE
          gsq => locgsq
          xsq => locxsq
          diffsq => locdiffsq
        END IF

        IF (c%usegrad .OR. c%useshanno) THEN ! (norm of g) ^2

          IF (PRESENT(g)) THEN
            SELECT CASE (c%thenorm)
            CASE (l1)
              gsq = SUM(ABS(g))**2
            CASE (l2)
              gsq = g .ip. g
            CASE (linf)
              gsq = MAXVAL(ABS(g))**2
            CASE (g2)
              gsq = g .gip. g
            END SELECT
          ELSE
            gsq = zero
          END IF

          IF (gsq==zero) THEN
            test_done = true
            GO TO 10
          END IF
        END IF

        IF (c%usestep .OR. c%useshanno) THEN ! (norm of xi) ^2

          IF (PRESENT(xi)) THEN
            SELECT CASE (c%thenorm)
            CASE (l1)
              xsq = SUM(ABS(xi))**2
            CASE (l2)
              xsq = xi .ip. xi
            CASE (linf)
              xsq = MAXVAL(ABS(xi))**2
            CASE (g2)
              xsq = xi .gip. xi
            END SELECT
          ELSE
            xsq = zero
          END IF
        END IF

        IF (c%usestep) THEN ! ( norm of xi-xim1 ) ^2

          IF (PRESENT(xim1) .AND. PRESENT(xi)) THEN
            SELECT CASE (c%thenorm)
            CASE (l1)
              diffsq = SUM(ABS(xi-xim1))**2
            CASE (l2)
              diffsq = (xi-xim1) .ip. (xi-xim1)
            CASE (linf)
              diffsq = MAXVAL(ABS(xi-xim1))**2
            CASE (g2)
              diffsq = (xi-xim1) .gip. (xi-xim1)
            END SELECT
          ELSE
            test_done = false
            GO TO 10
          END IF
        END IF

        IF (c%traceterm) THEN
          CALL indent(id,c%termtraceunit,level)
          WRITE (c%termtraceunit,90010) eps
        END IF

        pass = true
        IF (c%usegrad) THEN
          passthis = gsq <= (eps*c%normgatx0)**2
          pass = pass .AND. passthis
          IF (c%traceterm) THEN
            CALL indent(id,c%termtraceunit,level)
            WRITE (c%termtraceunit,90000) passthis, '(grad) gsq = ', gsq, &
              ' NormGatX0= ', c%normgatx0
          END IF
        END IF

        IF (c%usestep .AND. (pass .OR. c%traceterm)) THEN
          passthis = diffsq <= eps**2*MAX(one,xsq)
          pass = pass .AND. passthis
          IF (c%traceterm) THEN
            CALL indent(id,c%termtraceunit,level)
            WRITE (c%termtraceunit,90000) passthis, '(step) diffsq= ', diffsq, &
              '  xsq= ', xsq
          END IF
        END IF

        IF (c%useshanno .AND. (pass .OR. c%traceterm)) THEN
          passthis = gsq <= (eps*c%normgatx0)**2*MAX(one,xsq)
          pass = pass .AND. passthis
          IF (c%traceterm) THEN
            CALL indent(id,c%termtraceunit,level)
            WRITE (c%termtraceunit,90000) passthis, '(shxg) g = ', gsq, &
              '  xsq= ', xsq, '  NormGatX0= ', c%normgatx0
          END IF
        END IF

        IF (c%usefunc .AND. (pass .OR. c%traceterm)) THEN
          passthis = ABS(fx) <= eps*ABS(c%fatx0)
          pass = pass .AND. passthis
          IF (c%traceterm) THEN
            CALL indent(id,c%termtraceunit,level)
            WRITE (c%termtraceunit,90000) passthis, '(func) fx = ', fx, &
              '  fx0= ', c%fatx0
          END IF
        END IF

        test_done = pass

10      RETURN

! FORMATS:

90000   FORMAT (' Pass= ',L1,3(A,G8.3))
90010   FORMAT (' eps = ',G8.3)

      END FUNCTION

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> support.f90
    MODULE support ! Support routines for minimization -'USE'ed in module MinimizeF
      USE fscale, ONLY : f_scale
      USE testdone, ONLY : test_done
      USE evaluatef, ONLY : evaluate_f
      USE print, ONLY : print_iterate


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: print_iterate, evaluate_f, test_done, f_scale

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> min_codes.f90
    MODULE min_codes ! Simulation of Pascal enumerated types
      USE precisions, ONLY : short

!---
! This provides a facility for selecting multiple cases, when there
! are more than 2, i.e. when logical variables are not suitable.

! The codes may be accessed with a statement

!   USE Min_Codes

! A specific code, e.g. Resume, may be accessed as in

!   USE Min_Codes, only : Resume

! For each group of codes, which correspond to the cases possible
! for some situation , there is a subroutine Str_*(i) which will return
! a string corresponding to the current value of i. For example,
! Str_Entrystate(Resume), which is in fact Str_Entrystate(7), will
! return the string 'Resume'.

! The subroutine Control_Codes will print the entry status and the values
! of the minimization control codes that have been set.

! The subroutine Error_Codes will print the completion status if an error
! is detected.  If bad arguments have been provided, the routine will also
! specify the arguments that are in error.

!------------------------------------------------------------------------


      IMPLICIT NONE
      PRIVATE
!-----               -------------!

      PUBLIC :: str_entrystate, str_exitstate, str_errorcodes, str_method, &
        str_dointerp, str_startalpha, str_startstep, str_updateform, &
        str_scalegamma, str_htest, str_seth0

!------------------------------------------------------------------------
!Str_EntryState

! Entry Status of minimization.
! Usual case, no RC, no F or G
! No RC, but F, G at x0 given
! Start from scratch with RC
! Reentry for RC: continue
! Reentry for RC: F not available
! Reentry for RC: G not available
! Reentry for RC: F & G not available
      INTEGER (short), PARAMETER, PUBLIC :: normal = 0, normalwithfg = 1, &
        revcommstart = 2, revcommrestart = 3, revcommnof = 4, revcommnog = 5, &
        revcommnoforg = 6, resume = 7 ! Resume a checkpointed run

!------------------------------------------------------------------------
!Str_ExitState

! Completion Status of minimization
! Normal Completion
! RC: compute F
! RC: compute G
! RC: compute F and G
! Quit: storage request error
! Quit: initial point is min
! Quit: f undefd at initial point
! Quit: bad argument; see Error
! Quit: line search failed
! Quit: non descent direction
! Quit: too many f evaluations
! Quit: m=0 invalid with ProductForm
      INTEGER (short), PARAMETER, PUBLIC :: done = 0, revcommf = 1, &
        revcommg = 2, revcommfandg = 3, storageerror = 4, initialmin = 5, &
        initialundefd = 6, badinput = 7, lsfail = 8, nodescent = 9, &
        excessfevals = 10, invalidm = 11, abortmin = 12 ! Abort requested by Evaluate_f

!------------------------------------------------------------------------
!Str_ErrorCodes

! Error codes regarding input
! n must be 2 or more
! accuracy must be positive
! Unacceptable state
! memory must be nonnegative
! limit must be nonnegative
! print frequency must be nonnegative
! chkpt frequency must be nonnegative
! error opening/reading checkpnt file
! unit must be non-negative integer
! no such method defined
! productform not used with SD or CG
! must be non-negative integer
! no such derivative mode
! test memory must be non-negative
! scale must be positive integer
! expense must be positive integer
! unit must be non-negative integer
! unit must be non-negative integer
! no such kind of norm
! unit must be non-negative integer
! bad interpolation strategy
! bad method for first alpha
! bad gamma strategy
! bad restart test type
! bad update choice
! bad step length for CG alpha
! rho strictly positive
! beta between zero and one, inclusive
      INTEGER (short), PARAMETER, PUBLIC :: badn = 0, badacc = 1, &
        badstate = 2, badmemory = 3, badevallim = 4, badfreq = 5, &
        badcheckpoint = 6, badcheckopen = 7, badcheckunit = 8, badmethod = 9, &
        badcombination = 10, badterms = 11, badderiv = 12, badsysmemory = 13, &
        badscalef = 14, badexpense = 15, badevalunit = 16, badprintunit = 17, &
        badnorm = 18, badtermunit = 19, badinterp = 20, badalpha = 21, &
        badgamma = 22, badhtest = 23, badupdate = 24, badstep = 25, &
        badrho = 26, badbeta = 27, badseth0 = 28 ! no such kind of initial H0 setting

!------------------------------------------------------------------------
!Str_Method

! Choice of method
! Use steepest descent steps
! Use H0 and conjugate gradient
! Use Shanno's 2-step CONMIN
! Use given number of terms
! Variable number of terms
! Alg. according to storage available
! Use dynamic strategy based on time
      INTEGER (short), PARAMETER, PUBLIC :: sd = 0, cg = 1, conmin = 2, &
        fixterms = 3, variable = 4, available = 5, dynamic = 6, qn = 7 ! Use QN update matrix

!------------------------------------------------------------------------
!Str_DoInterp

! Interpolation control
! on every step, force quad. interpolation
! just on d[m+1] and later, force quad. int.
! just on d[m+2] and later, force quad. int.
! just on d[m+3] and later, force quad. int.
      INTEGER (short), PARAMETER, PUBLIC :: every = 0, on1 = 1, on2 = 2, &
        on3 = 3, never = 4 ! never force quadratic interpolation

!------------------------------------------------------------------------
!Str_StartAlpha

! Initial alpha strategy
! never use alpha=1 strategy
! start with alpha=1 BEFORE pt x[m+1] found
! start with alpha=1 BEFORE pt x[m+2] found
! start with alpha=1 BEFORE pt x[m+3] found
      INTEGER (short), PARAMETER, PUBLIC :: neverqn = 0, before1qn = 1, &
        before2qn = 2, before3qn = 3, alwaysqn = 4 ! always start with alpha=1

!------------------------------------------------------------------------
!Str_StartStep

! When initial alpha size /= 1
! Fletcher's choice for CG step
      INTEGER (short), PARAMETER, PUBLIC :: fletcher = 0, powell = 1 ! Powell's choice for CG step

!------------------------------------------------------------------------
!Str_UpdateForm

! Form of maintaining H
! Buckley LeNir method
! Nocedal's approach
      INTEGER (short), PARAMETER, PUBLIC :: sumform = 0, productform = 1, &
        factoredform = 2 ! Powell's factored ZZ^T

!------------------------------------------------------------------------
!Str_ScaleGamma

! Use gamma for step scaling
! don't use this strategy
! use it on first QN step only
      INTEGER (short), PARAMETER, PUBLIC :: nogammascale = 0, gammafirst = 1, &
        gammaall = 2 ! use it on all QN steps

!------------------------------------------------------------------------
!Str_HTest

! Restart test
! No restarts implemented
! Use the original Powell test
      INTEGER (short), PARAMETER, PUBLIC :: norestart = 0, usei = 1, useh = 2 ! Use the test with metric H

!------------------------------------------------------------------------
!Str_SetH0

! Form of setting H
! Identity matrix
! Diagonal matrix
! First computed matrix    See H0v_Multiply
! Second computed matrix   See H0v_Multiply
      INTEGER (short), PARAMETER, PUBLIC :: ident = 0, diagonal = 1, &
        computed1 = 2, computed2 = 3, computed3 = 4 ! Third computed matrix    See H0v_Multiply

!------------------------------------------------------------------------

! Form of allocating space
! Allocate a term of H
! Allocate Hi
! Allocate H0
! Allocate Temp
      INTEGER (short), PARAMETER, PUBLIC :: allocnext = 1, allochi = 2, &
        alloch0 = 3, alloctemp = 4, allocu = 5 ! Allocate u

!------------------------------------------------------------------------


    CONTAINS


      FUNCTION str_entrystate(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (14), PARAMETER :: str(normal:resume) = (/ 'Normal        ', &
          'NormalWithFG  ', 'RevCommStart  ', 'RevCommRestart', &
          'RevCommNoF    ', 'RevCommNoG    ', 'RevCommNoForG ', &
          'Resume        '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_entrystate

! DESCRIPTION: This function returns the specified value of the entry state.

! EXECUTION:
        str_entrystate = str(i)
        RETURN

      END FUNCTION

      FUNCTION str_exitstate(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (13), PARAMETER :: str(done:abortmin) = (/ 'Done         ', &
          'RevCommF     ', 'RevCommG     ', 'RevCommFandG ', 'StorageError ', &
          'InitialMin   ', 'Initialundefd', 'BadInput     ', 'LSFail       ', &
          'NoDescent    ', 'ExcessFEvals ', 'InvalidM     ', 'AbortMin     '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_exitstate

! DESCRIPTION: This function returns the specified value of the exit state.

! EXECUTION:
        str_exitstate = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_errorcodes(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (14), PARAMETER :: str(badn:badseth0) = (/ 'BadN          ', &
          'BadAcc        ', 'BadState      ', 'BadMemory     ', &
          'BadEvalLim    ', 'BadFreq       ', 'BadCheckPoint ', &
          'BadCheckOpen  ', 'BadCheckUnit  ', 'BadMethod     ', &
          'BadCombination', 'BadTerms      ', 'BadDeriv      ', &
          'BadSysMemory  ', 'BadScaleF     ', 'BadExpense    ', &
          'BadEvalUnit   ', 'BadPrintUnit  ', 'BadNorm       ', &
          'BadTermUnit   ', 'BadInterp     ', 'BadAlpha      ', &
          'BadGamma      ', 'BadHTest      ', 'BadUpdate     ', &
          'BadStep       ', 'BadRho        ', 'BadBeta       ', &
          'BadSetH0      '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_errorcodes

! DESCRIPTION: This function returns the specified value of the error code.

! EXECUTION:
        str_errorcodes = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_method(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (9), PARAMETER :: str(sd:qn) = (/ 'SD       ', 'CG       ', &
          'ConMin   ', 'FixTerms ', 'Variable ', 'Available', 'Dynamic  ', &
          'QN       '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_method

! DESCRIPTION: This function returns the specified value of Us%Method.

! EXECUTION:
        str_method = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_dointerp(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (5), PARAMETER :: str(every:never) = (/ 'Every', 'On1  ', &
          'On2  ', 'On3  ', 'Never'/)

        CHARACTER (LEN_TRIM(str(i))) :: str_dointerp

! DESCRIPTION: This function returns the specified value of Us%DoInterpolation.

! EXECUTION:
        str_dointerp = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_startalpha(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (9), PARAMETER :: str(neverqn:alwaysqn) = (/ 'NeverQN  ', &
          'Before1QN', 'Before2QN', 'Before3QN', 'AlwaysQN '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_startalpha

! DESCRIPTION: This function returns the specified value of Us%StartAlpha.

! EXECUTION:
        str_startalpha = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_startstep(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (8), PARAMETER :: str(fletcher:powell) = (/ 'Fletcher', &
          'Powell  '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_startstep

! DESCRIPTION: This function returns the specified value of Us%StartStep.

! EXECUTION:
        str_startstep = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_updateform(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (12), PARAMETER :: str(sumform:factoredform) = (/ &
          'SumForm     ', 'ProductForm ', 'FactoredForm'/)

        CHARACTER (LEN_TRIM(str(i))) :: str_updateform

! DESCRIPTION: This function returns the specified value of Us%UpdateForm.

! EXECUTION:
        str_updateform = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_scalegamma(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (12), PARAMETER :: str(nogammascale:gammaall) = (/ &
          'NoGammaScale', 'GammaFirst  ', 'GammaAll    '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_scalegamma

! DESCRIPTION: This function returns the specified value of Us%ScaleGamma.

! EXECUTION:
        str_scalegamma = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_htest(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (9), PARAMETER :: str(norestart:useh) = (/ 'NoRestart', &
          'UseI     ', 'UseH     '/)

        CHARACTER (LEN_TRIM(str(i))) :: str_htest

! DESCRIPTION: This function returns the specified value of Us%HTest.

! EXECUTION:
        str_htest = str(i)
        RETURN

      END FUNCTION


      FUNCTION str_seth0(i)

! ARGUMENTS:
        INTEGER (short) :: i

! LOCAL DECLARATIONS:

        CHARACTER (9), PARAMETER :: str(ident:computed3) = (/ 'Ident    ', &
          'Diagonal ', 'Computed1', 'Computed2', 'Computed3'/)

        CHARACTER (LEN_TRIM(str(i))) :: str_seth0

! DESCRIPTION: This function returns the specified value of Us%SetH0.

! EXECUTION:
        str_seth0 = str(i)
        RETURN

      END FUNCTION


    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> min_states.f90
    MODULE min_states ! the type definitions for control information.
      USE precisions, ONLY : long, short, stnd


      IMPLICIT NONE
      PRIVATE
!-----               -------------!

      PUBLIC :: hterm, hmatrix, userdefined, sharedvars, tracelist, &
        searchvals, minlocal, numerrs

      INTEGER (long), PARAMETER :: numerrs = 10

      TYPE :: hterm
        REAL (stnd), POINTER :: data(:)
        TYPE (hterm), POINTER :: next, prev
      END TYPE

      TYPE :: hmatrix ! structure of the QN matrix

        TYPE (hterm), POINTER :: firstterm, lastterm, nextterm
! used in product and updateh
        REAL (stnd), POINTER :: h0(:), hi(:), temp(:), u(:) ! used in factored module
        INTEGER (long) :: nterms, mterms

! Note that the matrix H may be:
!  o  a symmetric n x n matrix, stored as a vector of its upper half
!     in row order, and pointed to by Hi;
!  o  a linked list of terms which define H as an update of H0.  Then
!     H0 is stored as a vector which defines a diagonal matrix; it is
!     given as a pointer.  FirstTerm points at the first term of the
!     list, once it has been allocated.  NextTerm points at the space
!     allocated for the next term in the list until it is linked into
!     the list.  LastTerm always points to the last term in the
!     sequence.  Each term of the list contains a one-dimensional array
!     pointed at by Data.

      END TYPE

      TYPE :: userdefined

        INTEGER (short) :: method, dointerpolation, startalpha, scalegamma, &
          htest, updateform, startstep, seth0

        INTEGER (long) :: checkpoint, systemmemory

        REAL (stnd) :: decreaseinf, rho, beta, parh0

        LOGICAL :: scalecolumns, relativetof0, relativetog0, &
          quadinterpolation, ignoreinterval, countfromrestart
      END TYPE

      TYPE :: sharedvars
        INTEGER (short) :: methodinuse, error(0:numerrs)
        INTEGER (long) :: memory, memoryused, restartct, restartpt, forcect, &
          ncalls, ct, it, ixstart, ixnu, ixeta, ixs, ixu, ixy, ixc, ixend, &
          ixleng
        LOGICAL :: qnpart, dorestart, doforcerestart, steepdstep, limmemory, &
          doneinterpolation, lsfinished, validf
        REAL (stnd) :: accuracy, eta, iterationtime
      END TYPE

      TYPE :: tracelist ! for trace control. Any change: see also internal
! subroutine in Minimize_f.

! unit
        INTEGER (short) :: u, level
        LOGICAL :: input, flow, steptypes, lsalpha, lsreal, lsflow, update, &
          values, vectors, xandd, cubic
      END TYPE

      TYPE :: searchvals ! for Line_Search

        LOGICAL :: goodpt, alphalbdd, alphaubdd

        REAL (stnd) :: currmg, lwbd, flbd, dtglbd, upbd
      END TYPE

      TYPE :: minlocal

        INTEGER (long) :: checknext, memoryneeded, n

        REAL (stnd) :: eps, nrmg, nrmd, alpha, fp, dgp, ap, fmin, flast, &
          glast, dglast, dg0, nrmx

        LOGICAL :: firsteval, cold, checkflag, unknownm, all_allocated, &
          noprint
      END TYPE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> min_defaults.f90
    MODULE min_defaults
      USE min_states, ONLY : userdefined, tracelist
      USE true_false, ONLY : false, true
      USE reals, ONLY : one, zero, fifth
      USE precisions, ONLY : short, stnd, long
      USE min_codes, ONLY : dynamic, ident, powell, sumform, usei, gammafirst, &
        alwaysqn, never


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: defaultuserdefined, defaulttraces, checkfile, checkfilelen, &
        checkunit, maxupdates, userdefined, tracelist

      CHARACTER (*), PARAMETER :: checkfile = 'ChkPt'

      INTEGER (short), PARAMETER :: checkfilelen = 100, checkunit = 9

      INTEGER (long), PARAMETER :: maxupdates = 7

      INTEGER (short), PARAMETER :: method = dynamic, dointerpolation = never, &
        startalpha = alwaysqn, scalegamma = gammafirst, htest = usei, &
        updateform = sumform, startstep = powell, seth0 = ident

      INTEGER (long), PARAMETER :: checkpoint = 0, systemmemory = 0

      REAL (stnd), PARAMETER :: decreaseinf = -one, rho = fifth, beta = one, &
        parh0 = zero

      LOGICAL, PARAMETER :: scalecolumns = false, relativetof0 = false, &
        relativetog0 = false, quadinterpolation = true, &
        ignoreinterval = false, countfromrestart = false

      TYPE (userdefined) :: defaultuserdefined = userdefined(method, &
        dointerpolation,startalpha,scalegamma,htest,updateform,startstep, &
        seth0,checkpoint,systemmemory,decreaseinf,rho,beta,parh0,scalecolumns, &
        relativetof0,relativetog0,quadinterpolation,ignoreinterval, &
        countfromrestart)

      INTEGER (short), PARAMETER :: trunit = 6, trlevel = 0

      LOGICAL, PARAMETER :: trinput = false, trflow = false, &
        trsteptypes = false, trlsalpha = false, trlsreal = false, &
        trlsflow = false, trupdate = false, trvalues = false, &
        trvectors = false, trxandd = false, trcubic = false


! for trace control
      TYPE (tracelist) :: defaulttraces = tracelist(trunit,trlevel,trinput, &
        trflow,trsteptypes,trlsalpha,trlsreal,trlsflow,trupdate,trvalues, &
        trvectors,trxandd,trcubic)

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> systemstate.f90
    MODULE systemstate ! 'USE'ed in modules Minimize, MinimizeF, CheckPoint,
      USE supp_states, ONLY : evalstate, stnd, termvalues, termstate, &
        printvalues, printstate, evalerrors, evalcounts
      USE min_states, ONLY : userdefined, hmatrix, tracelist, minlocal, &
        searchvals, sharedvars
!                    and Restart.


      IMPLICIT NONE
      PRIVATE
!-----               -------------!

      PUBLIC :: minimizestate, userdefined, sharedvars, searchvals, minlocal, &
        evalstate, evalcounts, evalerrors, printstate, printvalues, termstate, &
        termvalues, tracelist, hmatrix


      TYPE :: minimizestate

        TYPE (userdefined) :: user
        TYPE (sharedvars) :: shared
        TYPE (searchvals) :: lssave
        TYPE (minlocal) :: local

        TYPE (evalstate) :: eval
        TYPE (evalcounts) :: evalcts
        TYPE (evalerrors) :: evalers
        TYPE (printstate) :: print
        TYPE (printvalues) :: prvals
        TYPE (termstate) :: term
        TYPE (termvalues) :: termvals

        TYPE (tracelist) :: trace

        REAL (stnd), POINTER :: d(:), xx(:), gg(:)
        TYPE (hmatrix) :: h

      END TYPE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> checkpoint.f90
    MODULE checkpoint ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving
      USE true_false, ONLY : true, false
      USE min_defaults, ONLY : checkfilelen, checkfile, checkunit
      USE min_states, ONLY : hterm
      USE systemstate, ONLY : minimizestate
      USE precisions, ONLY : stnd, long, short


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: check_point, minimizestate

    CONTAINS

! Optional arguments
      SUBROUTINE check_point(x,fx,g,state,c,iteration,file,outunit)

! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:), fx, g(SIZE(x))
        INTEGER (short), INTENT (IN) :: state
        INTEGER (long), INTENT (IN) :: iteration
        TYPE (minimizestate), INTENT (INOUT) :: c

        CHARACTER (*), INTENT (IN), OPTIONAL :: file
        INTEGER (short), INTENT (IN), OPTIONAL :: outunit

! DESCRIPTION:

!     Given the arguments and state vector for routine Minimize_f
!     (except FName), this routine writes these values to a file to guard
!     against losing information in the event of a system crash.  Note that
!     part of the output consists of the data in a linked list.

!     It accepts as optional arguments a file name and a unit number.

!     In order to avoid the possibility of a crash during checkpointing, it
!     uses two output files, and writes to them alternately. The first is
!     indicated by 'a' being appended to the file name, which is done when
!     C%Local%CheckFlag is true, and the second by an appended 'b', when
!     C%Local%CheckFlag is false.  Note that the routine automatically appends
!     the letters 'a' or 'b', and so the user should not include them in the
!     file name if it is specified.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'ChkPt'

! LOCAL DECLARATIONS:

        INTEGER (short) :: unit
        INTEGER (long) :: lent

        CHARACTER (len=LEN(id)+2+c%trace%level) :: pre

        TYPE (hterm), POINTER :: term

        CHARACTER ((checkfilelen)) :: outfile

        CHARACTER (1) :: suffix

! EXECUTION:

        pre = entering(id,c%trace%flow,c%trace%u,c%trace%level)

        IF (iteration==c%local%checknext) THEN ! Checkpoint this iteration.

! Open appropriate file for output.

          IF (PRESENT(outunit)) THEN
            unit = outunit
          ELSE
            unit = checkunit
          END IF

          IF (c%local%checkflag) THEN
            suffix = 'a'
          ELSE
            suffix = 'b'
          END IF

          IF (PRESENT(file)) THEN
            lent = LEN_TRIM(file)
            outfile = file(1:lent) // suffix ! User-supplied file name
          ELSE
            lent = LEN_TRIM(checkfile)
            outfile = checkfile(1:lent) // suffix ! Default file name
          END IF

! access appropriate file
          OPEN (unit,file=outfile,form='unformatted',status='unknown', &
            position='rewind')

! Write values.

          WRITE (unit) x, fx
          WRITE (unit) g, state

          WRITE (unit) c%trace
          WRITE (unit) c%user
          WRITE (unit) c%shared
          WRITE (unit) c%lssave
          WRITE (unit) c%local
          WRITE (unit) c%eval
          WRITE (unit) c%evalcts
          WRITE (unit) c%evalers
          WRITE (unit) c%print
          WRITE (unit) c%prvals
          WRITE (unit) c%term
          WRITE (unit) c%termvals
          WRITE (unit) c%d
          WRITE (unit) c%xx
          WRITE (unit) c%gg

          WRITE (unit) c%h%nterms
          WRITE (unit) c%h%mterms

          IF (ASSOCIATED(c%h%h0)) THEN
            WRITE (unit) true
            WRITE (unit) c%h%h0
          ELSE
            WRITE (unit) false
          END IF

          IF (ASSOCIATED(c%h%hi)) THEN
            WRITE (unit) true
            WRITE (unit) c%h%hi
          ELSE
            WRITE (unit) false
          END IF

          IF (ASSOCIATED(c%h%temp)) THEN
            WRITE (unit) true
            WRITE (unit) c%h%temp
          ELSE
            WRITE (unit) false
          END IF

          IF (ASSOCIATED(c%h%u)) THEN
            WRITE (unit) true
            WRITE (unit) c%h%u
          ELSE
            WRITE (unit) false
          END IF

          IF (ASSOCIATED(c%h%firstterm)) THEN

            term => c%h%firstterm
            WRITE (unit) term%data
            term => term%next

! linear list
            DO WHILE (ASSOCIATED(term) .AND. .NOT. ASSOCIATED(term,c%h% &
                firstterm)) ! circular list
              WRITE (unit) term%data
              term => term%next
            END DO

          END IF

          CLOSE (unit)

! Prepare variables for the next checkpoint.

          c%local%checkflag = .NOT. c%local%checkflag

          c%local%checknext = c%local%checknext + c%user%checkpoint

        END IF

        CALL leaving(id,c%trace%flow,c%trace%u,c%trace%level)

! FORMATS:  none.

! EXIT

        RETURN

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> myallocate.f90
    MODULE myallocate ! 'USE'ed in modules MinimizeF, Initialize
      USE min_codes, ONLY : allocnext, allocu, alloctemp, alloch0, allochi
      USE min_states, ONLY : sharedvars, userdefined, hmatrix, tracelist
      USE precisions, ONLY : long, short
      USE general, ONLY : entering, leaving


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: my_allocate, my_deallocate, sharedvars, userdefined, &
                hmatrix, tracelist

    CONTAINS

      SUBROUTINE my_allocate(size,type,enough,us,sh,tr,h)

!   ARGUMENTS:

        INTEGER (long), INTENT (IN) :: size
        INTEGER (short), INTENT (IN) :: type
        INTEGER (short), INTENT (INOUT) :: enough
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr
        TYPE (hmatrix), INTENT (INOUT) :: h

!   DESCRIPTION:

!       This routine explicitly allocates required space for the
!       matrix H.  The specific allocation depends on the value of 'type'.
!       The amount to be allocated is specified by either the argument 'size'
!       or by the shared variables ixstart and ixend.

!       The variable SystemMemory provides a way for the user to set a limit
!       on the amount of memory to be used.  If SystemMemory is 0 (the default),
!       then the allocation will automatically be attempted and will succeed or
!       fail depending on how much memory the system has available.
!       If SystemMemory is positive, allocation will only be done if the
!       total memory used (including that needed for the current allocation)
!       is within the value specified by SystemMemory.  This is true even if
!       the available memory exceeds the value of SystemMemory.
!       The variable MemoryUsed keeps track of the amount of memory allocated.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'MyAlloc'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

!   EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (us%systemmemory/=0) sh%memoryused = sh%memoryused + size

        IF (us%systemmemory==0 .OR. sh%memoryused<=us%systemmemory) THEN
          SELECT CASE (type)
          CASE (allocnext) ! Allocate update term
            ALLOCATE (h%nextterm)
            ALLOCATE (h%nextterm%data(sh%ixstart:sh%ixend),STAT=enough)
          CASE (allochi) ! Allocate QN matrix
            ALLOCATE (h%hi(1:size),STAT=enough)
          CASE (alloch0) ! Allocate diagonal scaling matrix
            ALLOCATE (h%h0(1:size),STAT=enough)
          CASE (alloctemp) ! Allocate temporary vector
            ALLOCATE (h%temp(1:size),STAT=enough)
          CASE (allocu) ! Allocate vector u for factored updates
            ALLOCATE (h%u(1:size),STAT=enough)
          END SELECT
        ELSE
          sh%memoryused = sh%memoryused - size
          enough = 1
        END IF

! EXIT:
        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

      END SUBROUTINE

      SUBROUTINE my_deallocate(sh,tr,h)

!   ARGUMENTS:

        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr
        TYPE (hmatrix), INTENT (INOUT) :: h

!   DESCRIPTION:

!       This routine explicitly deallocates space for the matrix H.


! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'MyDealloc'

! LOCAL DECLARATIONS:

        INTEGER (long) :: i
        CHARACTER (len=LEN(id)+2+tr%level) :: pre

!   EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        sh%memoryused = 0

        IF (sh%limmemory) THEN
          DO i = 1, h%nterms
            IF (ASSOCIATED(h%firstterm)) THEN
              h%nextterm => h%firstterm%next
              DEALLOCATE (h%firstterm%data)
              DEALLOCATE (h%firstterm)
              h%firstterm => h%nextterm
            END IF
          END DO

          IF (ASSOCIATED(h%h0)) DEALLOCATE (h%h0)

        ELSE ! QN case
          IF (ASSOCIATED(h%hi)) DEALLOCATE (h%hi)
        END IF

        IF (ASSOCIATED(h%temp)) DEALLOCATE (h%temp)
        IF (ASSOCIATED(h%u)) DEALLOCATE (h%u)

! EXIT:
        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> control.f90
    MODULE control !  'USE'ed in module MinimizeF
      USE min_codes, ONLY : str_entrystate, str_seth0, str_htest, &
        str_scalegamma, str_updateform, str_startstep, str_startalpha, &
        str_dointerp, str_method
      USE min_states, ONLY : userdefined, tracelist
      USE precisions, ONLY : short


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: control_codes, userdefined, tracelist

    CONTAINS

      SUBROUTINE control_codes(entrystate,us,tr)

! ARGUMENTS:

        INTEGER (short), INTENT (IN) :: entrystate
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (tracelist), INTENT (IN) :: tr

! DESCRIPTION:

!     This routine will write out the meanings of the variable EntryState
!     and the integer control codes of Us.

! EXECUTION:

        WRITE (tr%u,90000) ' Status on entry is ' // &
          str_entrystate(entrystate)

        WRITE (tr%u,90000) ' Integer control values:'

        WRITE (tr%u,90000) ' Method          = ' // str_method(us%method)
        WRITE (tr%u,90000) ' DoInterpolation = ' // str_dointerp(us% &
          dointerpolation)
        WRITE (tr%u,90000) ' StartAlpha      = ' // str_startalpha(us% &
          startalpha)
        WRITE (tr%u,90000) ' StartStep       = ' // str_startstep(us%startstep &
          )
        WRITE (tr%u,90000) ' UpdateForm      = ' // str_updateform(us% &
          updateform)
        WRITE (tr%u,90000) ' ScaleGamma      = ' // str_scalegamma(us% &
          scalegamma)
        WRITE (tr%u,90000) ' HTest           = ' // str_htest(us%htest)
        WRITE (tr%u,90000) ' SetH0           = ' // str_seth0(us%seth0)

! FORMATS:

90000   FORMAT (A)

! EXIT:
        RETURN

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> cubic.f90
    MODULE cubic ! 'USE'ed in module LineSearch
      USE general, ONLY : entering, leaving
      USE true_false, ONLY : false, true
      USE reals, ONLY : one, two, zero, three
      USE num_constants, ONLY : machhuge
      USE precisions, ONLY : stnd, extd
      USE min_states, ONLY : tracelist


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: cubic_interpolate, tracelist

    CONTAINS

      SUBROUTINE cubic_interpolate(a,fa,fpa,b,fb,fpb,left,right,x,exactint,tr, &
          ignoreinterval)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: a, fa, fpa, b, fb, fpb, left, right
        REAL (stnd), INTENT (OUT) :: x
        LOGICAL, INTENT (OUT) :: exactint

        TYPE (tracelist), INTENT (INOUT) :: tr

        LOGICAL, INTENT (IN), OPTIONAL :: ignoreinterval

! DESCRIPTION:

!        Given the points a and b, along with the function values
!     fa and fb and slopes fpa and fpb at each point, this routine
!     finds the point x  at which the cubic fitted to the data
!     has its minimum.
!        The values left and right define an interval.  If there is no minimum
!     or if it lies outside the interval, x is returned as one of the end
!     points, as appropriate.  ExactInt is returned as true if the value x
!     returned is equal to that obtained from the cubic interpolation.
!        If IgnoreInterval is present and true, the interval [left,right] is
!     ignored; it can be thought of as being from 0 to + infinity, except
!     that if the interpolated value is negative, the left boundary of the
!     interval is used, regardless of the value of IgnoreInterval.  In the
!     one case where ignoring the interval is appropriate, i.e. for a
!     quadratic function, this should never happen.
!        The interpolation is computed following details given by Lemarechal.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'Cubic'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        LOGICAL :: extreme, order, bbigger, pbigger, ignore

        REAL (stnd) :: p, disc, aleft, aright
        REAL (stnd) :: sgn, apr, bpr, num, xc, biggest

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (PRESENT(ignoreinterval)) THEN
          ignore = ignoreinterval
        ELSE
          ignore = false
        END IF

        aleft = MIN(left,right)
        aright = MAX(left,right)

        IF (tr%cubic) THEN
          WRITE (tr%u,90010) pre // ' a= ', a, '  fa= ', fa, '  fpa= ', fpa
          WRITE (tr%u,90010) pre // ' b= ', b, '  fb= ', fb, '  fpb= ', fpb
          WRITE (tr%u,90020) pre // ' interval [', aleft, ',', aright, ']'
        END IF

        extreme = false
        order = left <= right .EQV. a <= b
        sgn = SIGN(one,b-a)

        IF (tr%cubic) WRITE (tr%u,90030) pre // ' order=', order, '  sgn=', &
          sgn

EQUAL:  IF (a==b) THEN
          IF (tr%cubic) WRITE (tr%u,90000) pre // ' points equal.'
          x = a
          exactint = false
        ELSE

          p = REAL(fpa,extd) + REAL(fpb,extd) - REAL(three,extd)*REAL(fb-fa, &
            extd)/REAL(b-a,extd)

          IF (SIGN(one,fpb)/=SIGN(one,fpa)) THEN
            IF (p/=zero) THEN
              disc = REAL(one,extd) - (REAL(fpa,extd)/p)*(REAL(fpb/p,extd))
              disc = ABS(p)*SQRT(disc)
            ELSE
              disc = SQRT(ABS(REAL(fpa,extd)))*SQRT(ABS(REAL(fpb,extd)))
            END IF
          ELSE
            IF (tr%cubic) WRITE (tr%u,90000) pre // ' sign(fpa)=sign(fpb).'
            biggest = MAX(ABS(fpa),ABS(fpb),ABS(p))
            bbigger = biggest == ABS(fpb)
            pbigger = biggest == ABS(p)
            IF (tr%cubic) WRITE (tr%u,90020) pre // ' p= ', p, '  biggest= ', &
              biggest
            IF (tr%cubic) WRITE (tr%u,90040) pre // ' MachRtHuge= ', &
              SQRT(machhuge)
            IF (biggest<=SQRT(machhuge)) THEN
              disc = REAL(p**2,extd) - REAL(fpa,extd)*REAL(fpb,extd)
              IF (tr%cubic) WRITE (tr%u,90020) pre // ' p = ', p, '  disc=', &
                disc
            ELSE IF (pbigger) THEN
              disc = REAL(p,extd) - (REAL(fpb,extd)/REAL(p,extd))*REAL(fpa, &
                extd)
            ELSE IF (bbigger) THEN
              disc = (REAL(p,extd)/REAL(fpb,extd))*REAL(p,extd) - &
                REAL(fpa,extd)
            ELSE
              disc = (REAL(p,extd)/REAL(fpa,extd))*REAL(p,extd) - &
                REAL(fpb,extd)
            END IF
            IF (tr%cubic) WRITE (tr%u,90040) pre // ' disc=', disc
            IF (disc>=zero) THEN
              IF (biggest<=SQRT(machhuge)) THEN
                disc = SQRT(disc)
              ELSE
                disc = SQRT(disc)*SQRT(biggest)
              END IF
              IF (tr%cubic) WRITE (tr%u,90040) pre // ' disc=', disc
            ELSE
              exactint = false
              IF (fpa<zero) THEN
                x = aright
              ELSE
                x = aleft
              END IF
              IF (tr%cubic) WRITE (tr%u,90000) pre // ' no minimum!'
              GO TO 10
            END IF

          END IF

          disc = sgn*disc
          IF (tr%cubic) WRITE (tr%u,90040) pre // ' disc=', disc

          apr = REAL(fpa,extd) + REAL(fpb,extd) + REAL(two*p,extd)
          bpr = REAL(fpa,extd) + REAL(p,extd)
          IF (tr%cubic) WRITE (tr%u,90020) pre // ' apr= ', apr, '  bpr= ', &
            bpr

          IF (sgn*bpr<zero) THEN
            IF (tr%cubic) WRITE (tr%u,90000) pre // ' using regular form.'
            x = a + fpa*(b-a)/(bpr-disc)
            IF (tr%cubic) WRITE (tr%u,90040) pre // ' predict x=', x
          ELSE
            num = disc + bpr
            IF (tr%cubic) WRITE (tr%u,90000) pre // ' using alternate form.'
            IF (tr%cubic) WRITE (tr%u,90040) pre // ' num=', num
            IF (ABS((a-b)*num)>=(aright-aleft)*ABS(apr)) THEN
              x = aright
              extreme = true
              IF (tr%cubic) WRITE (tr%u,90040) pre // ' cut off to x=', x
            ELSE
              x = a + num*(b-a)/apr
              IF (tr%cubic) WRITE (tr%u,90040) pre // ' predict x=', x
            END IF
          END IF

          xc = x
          IF (x<zero) THEN
            x = aleft ! use left bound if interpolated value negative
          ELSE IF ( .NOT. ignore) THEN
            x = MAX(x,aleft)
            x = MIN(x,aright)
          END IF
          exactint = .NOT. extreme .AND. xc == x

          IF (tr%cubic) WRITE (tr%u,90050) pre // ' x= ', x, '  xc= ', xc, &
            '  ExactInt= ', exactint, '  extreme= ', extreme
        END IF EQUAL

10      CONTINUE
        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:

90000   FORMAT (A)
90010   FORMAT (3(A,G15.7))
90020   FORMAT (2(A,G15.7),A)
90030   FORMAT (A,L2,A,G15.7)
90040   FORMAT (A,G15.7)
90050   FORMAT (2(A,G15.7),2(A,L1))

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> errorcodes.f90
    MODULE errorcodes !  'USE'ed in module MinimizeF
      USE min_codes, ONLY : storageerror, badseth0, badbeta, badrho, badstep, &
        badupdate, badhtest, badgamma, badalpha, badinterp, badtermunit, &
        badnorm, badprintunit, badevalunit, badexpense, badscalef, &
        badsysmemory, badderiv, badterms, badcombination, badmethod, &
        badcheckunit, badcheckopen, badcheckpoint, badfreq, badmemory, &
        badstate, badacc, badn, abortmin, invalidm, excessfevals, nodescent, &
        lsfail, badinput, initialundefd, initialmin
      USE precisions, ONLY : short, long
      USE min_states, ONLY : sharedvars, tracelist


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: error_codes, sharedvars, tracelist

    CONTAINS

      SUBROUTINE error_codes(finalstate,sh,tr)

! ARGUMENTS:

        INTEGER (short), INTENT (IN) :: finalstate
        TYPE (sharedvars), INTENT (IN) :: sh
        TYPE (tracelist), INTENT (IN) :: tr

! DESCRIPTION:

!   This routine will write out the meaning of the code FinalState.

! LOCAL DECLARATIONS:

        INTEGER (long) :: i

! EXECUTION:

        SELECT CASE (finalstate) ! write out the specific error condition
        CASE (storageerror)
          WRITE (tr%u,90000) 'Insufficient storage for attempted allocation'
        CASE (initialmin)
          WRITE (tr%u,90000) 'Initial point is a critical point'
        CASE (initialundefd)
          WRITE (tr%u,90000) 'Function, gradient, or both are undefined at ' &
            // 'the initial point'
        CASE (badinput)
          WRITE (tr%u,90000) 'Unacceptable values supplied for the following ' &
            // 'arguments to Minimize_f:'
        CASE (lsfail)
          WRITE (tr%u,90000) 'Failure in line search: abnormally small ' // &
            'steepest descent step from a cold start'
        CASE (nodescent)
          WRITE (tr%u,90000) 'Current direction is a non-descent direction'
        CASE (excessfevals)
          WRITE (tr%u,90000) 'The preset limit on the allowed number of ' // &
            'function evaluations has been exceeded'
        CASE (invalidm)
          WRITE (tr%u,90010) 'No update terms have been stored: this is', &
            'invalid for the ProductForm (Nocedal) method'
        CASE (abortmin)
          WRITE (tr%u,90010) 'Evaluate_f has detected an abort request', &
            'from the function evaluation routine'
        END SELECT

        IF (finalstate==badinput) THEN

          DO i = 1, sh%error(0)
            SELECT CASE (sh%error(i))
            CASE (badn)
              WRITE (tr%u,90000) '  N: must be >= 2'
            CASE (badacc)
              WRITE (tr%u,90000) '  Accuracy: must be > 0'
            CASE (badstate)
              WRITE (tr%u,90000) '  State: value undefined'
            CASE (badmemory)
              WRITE (tr%u,90000) '  Memory: must be >= 0'
            CASE (badfreq)
              WRITE (tr%u,90000) '  Frequency: must be >= 0'
            CASE (badcheckpoint)
              WRITE (tr%u,90000) '  CheckPoint: must be >= 0'
            CASE (badcheckopen)
              WRITE (tr%u,90000) '  CheckFile: cannot be opened or read'
            CASE (badcheckunit)
              WRITE (tr%u,90000) '  CheckUnit: must be non-negative integer'
            CASE (badmethod)
              WRITE (tr%u,90000) '  Method: value undefined'
            CASE (badcombination)
              WRITE (tr%u,90000) '  Cannot use ProductForm with SD or CG'
            CASE (badterms)
              WRITE (tr%u,90000) '  Terms: must be non-negative integer'
            CASE (badderiv)
              WRITE (tr%u,90000) '  Derivatives: value undefined'
            CASE (badsysmemory)
              WRITE (tr%u,90000) '  SystemMemory: must be >= 0'
            CASE (badscalef)
              WRITE (tr%u,90000) '  ScaleF: must be positive integer'
            CASE (badexpense)
              WRITE (tr%u,90000) '  Expense: must be positive integer'
            CASE (badevalunit)
              WRITE (tr%u,90000) '  EvalTraceUnit: must be integer >= 0'
            CASE (badprintunit)
              WRITE (tr%u,90000) '  PrintUnit: must be integer >= 0'
            CASE (badnorm)
              WRITE (tr%u,90000) '  TheNorm: value undefined'
            CASE (badtermunit)
              WRITE (tr%u,90000) '  TermTraceUnit: must be integer >= 0'
            CASE (badinterp)
              WRITE (tr%u,90000) '  DoInterpolation: value undefined'
            CASE (badalpha)
              WRITE (tr%u,90000) '  StartAlpha: value undefined'
            CASE (badgamma)
              WRITE (tr%u,90000) '  ScaleGamma: value undefined'
            CASE (badhtest)
              WRITE (tr%u,90000) '  HTest: value undefined'
            CASE (badupdate)
              WRITE (tr%u,90000) '  UpdateForm: value undefined'
            CASE (badstep)
              WRITE (tr%u,90000) '  StartStep: value undefined'
            CASE (badrho)
              WRITE (tr%u,90000) '  rho: must be > 0'
            CASE (badbeta)
              WRITE (tr%u,90000) '  beta: must be between 0 and 1 inclusive'
            CASE (badseth0)
              WRITE (tr%u,90000) '  SetH0: value undefined'
            END SELECT
          END DO
        END IF

! FORMATS:

90000   FORMAT (A)
90010   FORMAT (A/A)

! EXIT:
        RETURN

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> h0v.f90
    MODULE h0v ! 'USE'ed in modules MinimizeF, Sum, Product, Factored
      USE general, ONLY : entering, leaving
      USE min_codes, ONLY : ident, computed3, computed2, computed1, diagonal
      USE reals, ONLY : zero, two, one
      USE precisions, ONLY : stnd, long
      USE min_states, ONLY : hmatrix, tracelist, userdefined


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: h0v_multiply, hmatrix, tracelist, userdefined

    CONTAINS

      SUBROUTINE h0v_multiply(h,v,h0v,us,tr,u,ka)


! ARGUMENTS:

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h

        REAL (stnd), INTENT (IN) :: v(:)
        REAL (stnd), INTENT (OUT) :: h0v(SIZE(v))

        TYPE (userdefined), INTENT (IN) :: us

        TYPE (tracelist), INTENT (INOUT) :: tr

        REAL (stnd), INTENT (INOUT), OPTIONAL :: u(SIZE(v))

        INTEGER (long), INTENT (IN), OPTIONAL :: ka

! DESCRIPTION:

!     The basic purpose of this routine is to compute the value of the product
!     of the matrix  H0  and the vector  v.  Note that H0 may be defined
!     as Z0 * Z0^T, where Z0 is given.

!     On entry the following values are required:

!       H%H0    Vector of length n that stores a matrix, where n = size(v).
!               H0 is only required on entry if SetH0 = Diagonal (see below).
!       v       Vector of length  n  that is multiplied by  H0.
!       u       Vector of length  n.  See ka below.
!       k=ka    If k=0,  compute H0 * v  (the default if ka not present).
!               If k=-1, compute H0 * v as for k=0.  Also return u = Z0^T * v.
!               If k=-2, compute Z0*Z0^T and return it in the array H%Hi.
!                  Note that H0v is not altered in this case. Also, Hi should
!                  be returned in symmetric, upper half form, stored by row.
!                  Hi must be of an appropriate size and must have been
!                  allocated.
!               If k=1,...,n,  compute Z0[k] - (v^T * Z0[k])*u, where Z0[k] is
!                  the kth column of Z0, and return it in H0v. U must be passed
!                  in.
!               It is required that -2 <= k <= n; otherwise the routine will
!               exit with no calculations done.

!      Codes (found in Us):

!       SetH0   Flag to indicate the form that H0 takes:
!               Ident:     H0 identity ( Z0 = I ).
!               Diagonal:  H0 diagonal ( Z0 = SQRT(H0) ) (for upwards
!                             compatability).
!               Computed1: H0 computed in this routine according to the formula:
!                             H0 = Z0 * Z0^T where Z0(I,J) = (I-J)**2
!               Computed2: H0 computed in this routine according to the formula:
!                             H0 = Z0 * Z0^T where Z0(I,J) = 1 for I >= J
!                                                          = 2 for I <  J
!               Computed3: H0 computed in this routine according to the formula:
!                             H0 = Z0 * Z0^T where Z0(I,J) = 1/(I+J+ParH0)

!       ParH0   Used in the specification of H0 if SetH0 = Computed3.

!     On exit from the routine the following values are returned:
!       H0v    Specified product, for k >= -1.
!       u      Also contains Z0^T*v if k = -1.
!       Hi     Contains Z0*Z0^T if k = -2.


! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'H0v'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: i, j, jj, m, n, k

        REAL (stnd) :: zv, t1, t2

        REAL (stnd), POINTER :: h0(:), hi(:)

! EXECUTION

        pre = entering(id,tr%flow,tr%u,tr%level)

        n = SIZE(v)

        h0 => h%h0
        hi => h%hi

        IF (PRESENT(ka)) THEN
          k = ka
        ELSE
          k = 0
        END IF

        IF (k>n .OR. k<-2) RETURN ! Must have k between -2 and n.

        IF (us%seth0==ident .OR. us%seth0==diagonal .AND. .NOT. ASSOCIATED(h0) &
            ) THEN

! H0 is the identity matrix.  Z0 = I.

          SELECT CASE (k)
          CASE (0,-1) !    Compute H0 * v.

            h0v = v
            IF (k==-1) u = v

          CASE (-2) !    Compute Hi and store upper half of it by row.

            hi = zero
            j = 1
            DO i = 1, n
              hi(j) = one
              j = j + n - i + 1
            END DO

          CASE (1:) !    Compute Z0[k] - (v^T*Z0[k])u.

            h0v = -v(k)*u
            h0v(k) = h0v(k) + one
          END SELECT

        ELSE IF (us%seth0==diagonal) THEN

! H0 is a diagonal matrix stored in vector H%H0.  Z0 = SQRT( H0 ).

          SELECT CASE (k)
          CASE (0,-1) !   Compute H0 * v.
            h0v = h0*v

            IF (k==-1) u = SQRT(h0)*v

          CASE (-2) !   Compute H0 and store upper half of it by row.
            hi = zero
            j = 1
            DO i = 1, n
              hi(j) = h0(i)
              j = j + n - i + 1
            END DO

          CASE (1:) !   Compute Z0[k] - (v^T*Z0[k])u.
            h0v = (-v(k)*SQRT(h0(k)))*u
            h0v(k) = h0v(k) + SQRT(h0(k))
          END SELECT

        ELSE IF (us%seth0==computed1) THEN

! Z0(i,j) = (i-j)**2.  Note that Z0 = Z0^T.

          SELECT CASE (k)
          CASE (0,-1) !    Compute H0 * v.

            h0v = zero
            DO j = 1, n ! zv = Z0(I,J)^T * v(I)
              zv = zero
              DO i = 1, n
                zv = zv + (i-j)**2*v(i)
              END DO

              IF (k==-1) u(j) = zv

! Z0(I,J) * zv
              DO i = 1, n
                h0v(i) = h0v(i) + (i-j)**2*zv
              END DO
            END DO

          CASE (-2) !    Compute H0 and store upper half of it by row.

            m = 0
            hi = zero
            DO i = 1, n
              DO j = i, n
                m = m + 1
                DO jj = 1, n
                  hi(m) = hi(m) + (i-jj)**2*(jj-j)**2
                END DO
              END DO
            END DO

          CASE (1:) !    Compute Z0[k] - (v^T*Z0[k])u.

            zv = zero
            DO i = 1, n
              zv = zv + (i-k)**2*v(i) ! v^T*Z0[k]
            END DO

            DO i = 1, n
              h0v(i) = (i-k)**2 - zv*u(i)
            END DO

          END SELECT

        ELSE IF (us%seth0==computed2) THEN

! Z0(i,j) is 1 for i >= j and 2 for i < j.

          SELECT CASE (k)
          CASE (0,-1) !       Compute H0 * v.

            h0v = zero
            DO j = 1, n
              zv = zero
              DO i = 1, j - 1
                zv = zv + two*v(i)
              END DO
              DO i = j, n
                zv = zv + v(i)
              END DO

              IF (k==-1) u(j) = zv

              DO i = 1, j - 1
                h0v(i) = h0v(i) + two*zv
              END DO
              DO i = j, n
                h0v(i) = h0v(i) + zv
              END DO
            END DO

          CASE (-2) !    Compute H0 and store upper half of it by row.

            m = 0
            DO i = 1, n
              DO j = i, n
                m = m + 1
                hi(m) = zero
                DO jj = 1, n
                  IF (i>=jj) THEN
                    t1 = one
                  ELSE
                    t1 = two
                  END IF

                  IF (jj<=j) THEN
                    t2 = one
                  ELSE
                    t2 = two
                  END IF

                  hi(m) = hi(m) + t1*t2
                END DO
              END DO
            END DO

          CASE (1:) !     Compute Z0[k] - (v^T*Z0[k])u.

            zv = zero
            DO i = 1, k - 1
              zv = zv + two*v(i)
            END DO
            DO i = k, n
              zv = zv + v(i)
            END DO

            DO i = 1, k - 1
              h0v(i) = two - zv*u(i)
            END DO
            DO i = k, n
              h0v(i) = one - zv*u(i)
            END DO

          END SELECT

        ELSE IF (us%seth0==computed3) THEN

! Z0(i,j) = 1/(i+j+ParH0).  Note that Z0 = Z0^T.

          SELECT CASE (k)
          CASE (0,-1) !    Compute H0 * v.

            h0v = zero
            DO j = 1, n
              zv = zero
              DO i = 1, n
                zv = zv + v(i)/(i+j+us%parh0)
              END DO

              IF (k==-1) u(j) = zv

              DO i = 1, n
                h0v(i) = h0v(i) + zv/(i+j+us%parh0)
              END DO
            END DO

          CASE (-2) !    Compute H0 and store upper half of it by row.

            m = 0
            DO i = 1, n
              DO j = i, n
                m = m + 1
                hi(m) = zero
                DO jj = 1, n
                  hi(m) = hi(m) + one/((i+jj+us%parh0)*(jj+j+us%parh0))
                END DO
              END DO
            END DO

          CASE (1:) !    Compute Z0[k] - (v^T*Z0[k])u.

            zv = zero
            DO i = 1, n
              zv = zv + v(i)/(i+k+us%parh0)
            END DO

            DO i = 1, n
              h0v(i) = one/(i+k+us%parh0) - zv*u(i)
            END DO

          END SELECT

        END IF

! EXIT:

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:  None.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> factored.f90
    MODULE factored ! 'USE'ed in module MinimizeF
      USE general, ONLY : writevec, leaving, entering
      USE reals, ONLY : zero
      USE h0v, ONLY : hmatrix, h0v_multiply
      USE min_codes, ONLY : diagonal
      USE min_states, ONLY : userdefined, hterm, tracelist, sharedvars
      USE precisions, ONLY : stnd, long
      USE inner_product, ONLY : OPERATOR (.gip.)



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: factored_multiply, hmatrix, userdefined, hterm, &
                tracelist, sharedvars

    CONTAINS

! ndelta, scalmu, deltai, gammai, etai,&
      SUBROUTINE factored_multiply(h,g,hg,u,us,sh,tr)

! ARGUMENTS:

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h
        REAL (stnd), INTENT (IN) :: g(:)
        REAL (stnd), INTENT (INOUT) :: u(SIZE(g))
!   Real (stnd),          intent(INOUT)         :: g(:), u(size(g))!, scalmu
!   Real (stnd),          intent(IN)            :: deltai(size(g)), &
!                                                  gammai(size(g)), &
!                                                  ndelta,          &
!                                                  etai
        REAL (stnd), INTENT (OUT) :: hg(SIZE(g))
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr

! DESCRIPTION:

!     Given the quasi-Newton update matrix  H (in ZZ' form) and
!     given the vector g, this routine computes

!               hg = H * g  .

!     It also returns the intermediate value u = Z' * g.

!     It assumes here that Z represents the matrix Z[i] at the point x[i].

!     Each iteration defines a "term" of the update, where each
!     term requires 3n or 4n entries of H, namely for

!                   u[i], sh[i] and yh[i], and, optionally, scale[i].

!     Here    s = x[i] - x[i-1] = alpha * d[i+1]
!             y = g[i] - g[i-1]
!             r = sqrt (s' * y)

!       and this defines the values

!             sh = s / r
!             yh = y / r

!       The vector u defines the orthogonal rotations.

!       If  ScaleColumns  is true, then the scaling vector scale[i] is also
!       stored for each i.  Let Z[i]*k denote the kth column of Z[i].  Then,
!       by definition, scale[i]_1 contains the current value of
!                                 sigma = ||Z[i]_1||,
!                      scale[i]_k, k > 1, contains the value by which Z[i]*k
!                                 should be multiplied.
!       Note that scale[i]_1 = ||Z[i]*1|| = sigma  and
!                 scale[i]_k * Z[i]*k >= sigma.
!       Columns Z[i]*k of length exceeding sigma are left unchanged.

!       On entry, scale[j], j < i, must be defined, along with scale[i]_1.
!       The values for scale[i]_k, k > 1, are determined and stored when the
!       columns of Z[i] are found.

!       If i > m, the following are also required:
!          deltai - step[i]
!          gammai - dgrad[i]
!          ndelt  - norm of delt at step i
!          etai   - deltai' * gammai

!       mu is used for recording the current scaling factor ndelt.

!       Note that at most m values are stored for step and dgrad,
!       but that there is always one more value of s stored, so the
!       maximum may be m+1. Also note, and this is related, that
!       deltai, gammai and etai are only used when i > m. Otherwise,
!       they are ignored.

!       Calculation of Hg also requires use of several temporary
!       areas, namely
!               v       an n by m+1 matrix
!               w       an n   vector
!               z       an n   vector
!               phi     an m+1 vector

!      NDELTA - norm of delta-hat at level i.

!     The tracing vector Tr is accessed.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'FacMult'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: k, j, lowrow, hirow, n, i
        INTEGER (long), POINTER :: m

! AUTOMATIC
        REAL (stnd) :: w(1:SIZE(g)), z(1:SIZE(g)), phi(0:h%mterms), &
          v(1:SIZE(g),0:h%mterms)
        REAL (stnd), POINTER :: uj(:), shj(:), yhj(:), sh1mk(:)

        TYPE (hterm), POINTER :: term, startterm

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        m => h%mterms
        i = sh%ct
        n = SIZE(g)

        IF (tr%vectors) THEN
          CALL writevec(g,f='g',unit=tr%u,name=pre//' g=',indent=-tr%level)
          CALL writevec(u,f='g',unit=tr%u,name=pre//' u=',indent=-tr%level)
        END IF
        IF (tr%values) WRITE (tr%u,90000) pre // ' i = ', i, ' m =', m

        IF (us%scalecolumns .AND. us%seth0==diagonal) THEN
          IF (i==1) THEN
! appscl = ndelta
          ELSE
! appscl = min( scalmu, ndelta )
          END IF
        ELSE
!appscl = ONE
        END IF

! IF ( i <= m ) scalmu = appscl

        IF (i==0) THEN
          CALL h0v_multiply(h,g,hg,us,tr,u=u,ka=-1)
        ELSE

          hg = zero

          startterm => h%firstterm ! to initialize sh1mk to something

PICKDIAGONAL: DO k = n + 1, 2 - i, -1
! k >= 2: k-1 is starting column along bottom row of table
! k <= 1: 2-k is starting row along left hand side of table

            lowrow = MAX(0,2-k)
            hirow = MIN(n-k,i-1) ! old: , m)

            term => startterm
            CALL getnextterm(0)
            IF (k<=1) THEN
              sh1mk => shj
              CALL getnextterm(1)
              startterm => term
            END IF

EACHROW:    DO j = lowrow, hirow ! Do rows up one diagonal...

              IF (tr%flow) WRITE (tr%u,*) pre // 'k, j = ', k, j

! j     is current row in table, so that
! j+k-1 is current column and
! j+k   is next column to right
! we are computing z for row j+1, column j+k

              IF (j==0) THEN ! starting from row 0.

                CALL h0v_multiply(h,yhj,w,us,tr,u=shj,ka=k-1)

              ELSE IF (j+k==2) THEN ! in column 1 since j+k-1=1

                IF (tr%vectors) CALL writevec(sh1mk,f='g',unit=tr%u, &
                  name=pre//' j+k=2: sh1mk = ',indent=-tr%level)
                IF (tr%vectors) CALL writevec(shj,f='g',unit=tr%u, &
                  name=pre//' shj=',indent=-tr%level)
                w = sh1mk - (yhj.gip.sh1mk)*shj

              ELSE !   if ( j <= m-1 ) then

                w = z - (yhj.gip.z)*shj

!else    ! must be past row m.

!tmp = gammai .GIP. z
!w   = z - tmp/etai*deltai

              END IF

              IF (tr%vectors) CALL writevec(w,f='g',unit=tr%u,name=pre//' w=', &
                indent=-tr%level)
              IF (tr%flow) WRITE (tr%u,*) pre // 'Computing z for col j+k=', &
                j + k, ' row j+1=', j + 1

              z = (-w+uj(j+k-1)/phi(j)*v(:,j))*SQRT(phi(j)/(phi(j)+uj(j+ &
                k-1)**2))
              CALL writevec(z,f='g',unit=tr%u,name=pre//'z=',indent=-tr%level)

!IF ( Us%ScaleColumns ) z = z * scalej(j+k-1)

              IF (j+k>2) THEN ! not column number 1
                IF (tr%flow) WRITE (tr%u,*) pre // 'Updating v, phi', j, &
                  ' with u at col j+k-1=', j + k - 1
                v(:,j) = v(:,j) + uj(j+k-1)*w
                phi(j) = phi(j) + uj(j+k-1)**2
                IF (tr%vectors) CALL writevec(v(:,j),f='g',unit=tr%u, &
                  name=pre//' v=',indent=-tr%level)
                IF (tr%values) WRITE (tr%u,*) pre // ' phi[', j, ']=', phi(j)
              END IF

              IF (j<hirow) THEN
                CALL getnextterm(1)
              END IF

            END DO EACHROW

! Have z for top row i, col k+i-1,
!     or for rhs col n, row n-k+1

            j = hirow + 1
!            IF ( Us%ScaleColumns .and.  j <= m ) cscale(k+j-1,j) = appscl

ONRHS:      IF (k>n-i+1) THEN ! j = n-k+1
              IF (tr%flow) WRITE (tr%u,*) pre // 'Initializing v, phi ', j

              IF (k==n+1) THEN
                CALL h0v_multiply(h,yhj,v(:,j),us,tr,u=shj,ka=n)
                v(:,j) = uj(n)*v(:,j)
              ELSE
                v(:,j) = uj(n)*(z-(yhj.gip.z)*shj)
              END IF
              IF (tr%vectors) CALL writevec(v(:,j),f='g',unit=tr%u, &
                name=pre//' v=',indent=-tr%level)
              phi(j) = uj(n)**2
              IF (tr%vectors) WRITE (tr%u,*) pre // 'phi[', j, ']=', phi(j)
            ELSE ONRHS ! Along top, j = i
              IF (tr%flow) WRITE (tr%u,*) pre // 'Top Row'

              IF (us%scalecolumns) THEN
!scalej(k+j-1) = max( ONE, scalej(1)/.NORM. z )
!z = scalej(k+j-1) * z
              END IF

              IF (k+j==2) THEN
                u(1) = sh1mk .gip. g
              ELSE
                u(k+j-1) = z .gip. g
              END IF
              hg = hg + u(k+j-1)*z
              IF (tr%values) WRITE (tr%u,*) pre // 'u[', k + j - 1, ']=', &
                u(k+j-1)
              IF (tr%vectors) CALL writevec(hg,f='g',unit=tr%u, &
                name=pre//' hg=',indent=-tr%level)

!               if(tr) then
!                  if ( trv ) then
!                      call zzwmat(0,' z',one,z,zero,z,n,1,n,1,1,
!    -                    'f',9,4,1,80,tru)
!                      if ( k-1+j .le. i .and. k-1+j .ge. 2) then
!                          if (tr%vectors) write(tru,*) pre//'  checking...j,k,k-1+j,'
!    -                                   ,'1-k=' ,j,k,k-1+j,1-k
!                          tmp = z(1)/shj(1,1-k)
!                          call zzwmat(0,'z-delt1-k',one,z,-tmp,
!    -                      shj(1,1-k), n,1,n,1,1,'f',9,4,1,80,tru)
!                      end if
!                  end if
!                  TMP = INNER(N, yhj(1,J), Z, NONORM,IW,RW,DW)
!                  if (tr%vectors) write(tru,*) pre//' zT*yhj[',j,']=',TMP
!               end if

            END IF ONRHS
!if (tr%vectors) write(Tr%u,*) pre//'hg = ', hg
!if (tr%vectors) write(Tr%u,*) pre//'phi(j)   = ', phi(j)

          END DO PICKDIAGONAL

          IF (tr%flow) WRITE (tr%u,*) pre // 'Top Row: computing u[1]', &
            ' with s hat ', i - 1
!       IF ( i <= m ) then
!           sh1mk => shj   ! help - test
!    if (tr%vectors) call WriteVec( sh1mk, &
!       f='g', unit=Tr%u, name=pre//' j+k=2: sh1mk = ',indent=-Tr%level)
!    u(1) = sh1mk .GIP. g
!           if (tr%vectors) write(Tr%u,*) pre//'u(1)  = ', u(1)    ! help
!           hg   = hg  +  u(1) * sh1mk
!       else
!           u(1) = ( deltai .GIP. v ) / sqrt(etai)
!            hg   = hg  +  u(1)/sqrt(etai) * deltai
!       end if

        END IF

        IF (tr%vectors) CALL writevec(u,f='g',unit=tr%u,name=pre//' u=', &
          indent=-tr%level)
        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:

90000   FORMAT (2(A,I5))

      CONTAINS

        SUBROUTINE getnextterm(nterms)

          INTEGER, INTENT (IN) :: nterms

          INTEGER :: i

          DO i = 1, nterms
            IF (ASSOCIATED(term%next)) term => term%next
          END DO

          uj => term%data(sh%ixu+1:sh%ixu+n)
          shj => term%data(sh%ixs+1:sh%ixs+n)
          yhj => term%data(sh%ixy+1:sh%ixy+n)
!IF ( Us%ScaleColumns ) scalej  => Term%Data(Sh%ixc+1 : Sh%ixc+n )

          IF (tr%vectors .AND. tr%values) THEN
            CALL writevec(uj,f='g',unit=tr%u,name=pre//' uj->', &
              indent=-tr%level)
            CALL writevec(shj,f='g',unit=tr%u,name=pre//' shj->', &
              indent=-tr%level)
            CALL writevec(yhj,f='g',unit=tr%u,name=pre//' yhj->', &
              indent=-tr%level)
!call WriteVec(scalej, f='g', unit=Tr%u, name=pre//' scalej->', &
!indent=-Tr%level)
          END IF

        END SUBROUTINE

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> initialize.f90
    MODULE initialize ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving
      USE reals, ONLY : zero, one
      USE min_codes, ONLY : diagonal, factoredform, storageerror, alloch0
      USE myallocate, ONLY : my_allocate
      USE precisions, ONLY : stnd, long, short
      USE min_states, ONLY : hmatrix, tracelist, userdefined, sharedvars



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: initialize_h, hmatrix, tracelist, userdefined, sharedvars

    CONTAINS

      SUBROUTINE initialize_h(x,g,h,finalstate,sh,us,tr)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: x(:)
        REAL (stnd), INTENT (IN) :: g(SIZE(x))

        INTEGER (short), INTENT (INOUT) :: finalstate

        TYPE (hmatrix), INTENT (INOUT) :: h
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (tracelist), INTENT (INOUT) :: tr

! DESCRIPTION:

!     The purpose of this routine is to define the initial
!     matrix H0 = H[0].  This is stored in n locations of one
!     of the fields of H and defines a diagonal or identity matrix.

!     To be more specific, a diagonal matrix  H is defined with
!     elements  H(1,1), H(2,2), ... , H(n,n), but for storage
!     convenience, H is actually defined as a vector of n elements
!     and these n values are stored in H%H0(1),...,H%H0(n).

!     Note that, if Us%SetH0 is  Ident  on entry, then
!     H is by default the identity. No storage is allocated for H%H0.

!     On entry, the current point  x  and the gradient  g  at  x
!     must be defined.

!     Both  x  and  g  are used to compute the diagonal scaling
!     entries of  H.  The scaling used is quite primitive and not
!     particularly to be recommended.  The main point is that the
!     facility is available, and anyone so desiring can easily implement
!     their own scaling.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'InitH'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre
        INTEGER (long) :: k, j, n
        INTEGER (short) :: enough


! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (sh%limmemory .AND. us%seth0==diagonal) THEN

          CALL my_allocate(SIZE(x),alloch0,enough,us,sh,tr,h)
! Define diagonal
! scaling matrix
          IF (enough/=0_short) THEN
            finalstate = storageerror
            GO TO 10
          END IF

          WHERE (g/=zero)
            h%h0 = ABS(x/g)
          ELSEWHERE
            h%h0 = ABS(x)
          END WHERE

          IF (us%updateform==factoredform) h%h0 = SQRT(h%h0)

        ELSE IF ( .NOT. sh%limmemory) THEN
          n = SIZE(x)
          j = 1

          h%hi = zero

          DO k = 1, n
            h%hi(j) = one ! Note: only half of H.
            j = j + n - k + 1
          END DO

        END IF

10      CALL leaving(id,tr%flow,tr%u,tr%level)

! FORMATS:  none.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> linesearch.f90
    MODULE linesearch ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving
      USE cubic, ONLY : tracelist, cubic_interpolate
      USE reals, ONLY : c0_9, one, zero, c0_3, three, c1_m2, ten, c1_m4, tenth
      USE min_codes, ONLY : never, on3, on2, on1
      USE true_false, ONLY : false, true
      USE precisions, ONLY : stnd, long
      USE min_states, ONLY : searchvals, hmatrix, sharedvars, userdefined



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: line_search, searchvals, hmatrix, sharedvars, &
                userdefined, tracelist

    CONTAINS

      SUBROUTINE line_search(alpha,f,dg,f0,dg0,ap,fp,dgp,width,ls,us,sh,h,tr)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: f, dg, f0, dg0
        REAL (stnd), INTENT (OUT) :: dgp, width
        REAL (stnd), INTENT (INOUT) :: alpha, fp, ap

        TYPE (searchvals), INTENT (INOUT), TARGET :: ls
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (hmatrix), INTENT (IN), TARGET :: h
        TYPE (tracelist), INTENT (INOUT) :: tr

! DESCRIPTION:

!     This routine performs one internal iteration of the line search.

!     First, note that the execution of this routine is very much influenced
!     by a number of variables which appear in the calling routine Minimize_f.
!     These variables are obtained from the record Sh of shared variables,
!     and from the record Us of user-defined values.

!     Assume that the current search is along a direction  d  from a starting
!     point x-beg, and that the current point along that line is  x.  Assume
!     that the previous point considered along this line was  x-prev; thus,
!     on the first call for a line search along a given direction  d  from a
!     point x-beg, x-prev is just x-beg.  Then, on entry to line_search:

!        alpha  is the step length to x (so x  is x-beg + alpha*d).
!        f      is the function value at x.
!        dg     is the inner product of d and the gradient at  x.

!        f0     is the function value at  x-beg.
!        dg0    is the inner product of d and the gradient at x-beg.

!        ap     is alpha at x-prev.
!        fp     is the function value at the previous point x-prev.
!        dgp    is the inner product of d and the gradient at x-prev.

!       The following are obtained from the record Sh:

!        Validf is true if f and dg are defined at alpha.

!        ncalls is a count of how many times the function has been evaluated
!               along this direction  d, including the evaluation at  x, but
!               not including the evaluation at  x-beg.

!        DoneInterpolation is initially false, but it is set to true when a
!               point is computed via interpolation and accepted as the next
!               trial point. This is used to prevent termination without having
!               done an interpolation.

!        ct     is the iteration number of the current direction  d
!               and of the point to be reached, namely  x.

!---on exit from Line_Search:

!    LSFinished will be returned as true if the value alpha input to Line_Search
!               defines a point at which the line search can be terminated.
!               Otherwise it should be returned as false and a new trial value
!               for alpha determined.

!        width  is the width of the interval bounding an acceptable value of
!               alpha.  If no upper bound is known, width is the distance
!               between the current alpha and the lower bound.

!        alpha  if LSFinished is false, this contains the next value of alpha to
!               be considered. In this case, the values for  ap, dgp and fp
!               will have been updated.

!     ap, dgp, fp  if LSFinished is false, and a new value is defined in alpha,
!              then the "previous" point becomes the point just calculated, so
!              fp, dgp and ap should be redefined as the values  f, dg and alpha
!              input to this routine.

!              Note that these values are *not* updated if there were no valid
!              function or gradient values at the previous point.

!       The following are obtained from the record Us:

!        DoInterpolation, QuadInterpolation, IgnoreInterval

!       See file min.doc for an explanation of these values.


! SUBROUTINES:   Cubic_Interpolation

! PARAMETERS:

        REAL (stnd), PARAMETER :: nearly1 = c0_9, bitsmall = tenth, &
          small = c1_m4, extrapolate = ten, initmargin = c1_m2, &
          xpndmargin = three, maxmargin = c0_3

        CHARACTER (*), PARAMETER :: id = 'LS'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        LOGICAL :: accept, forceinterp, nointerpneeded, interpdone, intpt, &
          test1, first, interpt, test2

        REAL (stnd) :: tp0, at, left, right, slice

        INTEGER (long), POINTER :: m
        REAL (stnd), POINTER :: lwbd, flbd, dtglbd, upbd, currmargin
        LOGICAL, POINTER :: goodpt, alphalbdd, alphaubdd


! EXECUTION:

        m => h%mterms
        currmargin => ls%currmg
        goodpt => ls%goodpt
        lwbd => ls%lwbd
        flbd => ls%flbd
        dtglbd => ls%dtglbd
        alphalbdd => ls%alphalbdd
        upbd => ls%upbd
        alphaubdd => ls%alphaubdd

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (tr%lsflow) THEN
          WRITE (tr%u,90020) pre // &
            ' m = 0    LSFinished DoneInterpolation ct ncalls', m == 0, &
            sh%lsfinished, sh%doneinterpolation, sh%ct, sh%ncalls
          WRITE (tr%u,90030) pre // ' values:' // &
            ' m  DoInterpolation LimMemory  IgnoreInterval QuadInterpolation', &
            h%mterms, us%dointerpolation, sh%limmemory, us%ignoreinterval, &
            us%quadinterpolation
          WRITE (tr%u,90090) sh%validf, 0.0, f0, dg0, ap, fp, dgp, alpha, f, &
            dg
        END IF

        first = sh%ncalls == 1

        IF (first) THEN

          alphalbdd = false
          lwbd = zero
          flbd = f0
          dtglbd = dg0

          alphaubdd = false
          upbd = zero

          goodpt = true
          currmargin = initmargin
        END IF

! Test whether the steplength criteria have been met.

        tp0 = f0 + small*alpha*dg0
        test1 = f < tp0
        test2 = dg >= nearly1*dg0

        IF (tr%lsreal) WRITE (tr%u,90010) pre // ' tp0 = ', tp0

        IF (sh%validf) THEN
          accept = test1 .AND. test2
        ELSE
          accept = false
        END IF

ACCEPTABLE: IF (accept) THEN

          IF (tr%lsflow) WRITE (tr%u,90000) pre // ' acceptable.'

! The basic acceptance test has been passed.  We must test whether the
! point may be immediately accepted, or if it is necessary to force
! another step because a required interpolation step has not yet
! been done.

! See if quadratic interpolation to be forced.

          nointerpneeded = us%dointerpolation == never .OR. &
            (us%dointerpolation==on1 .AND. sh%ct<m+1) .OR. &
            (us%dointerpolation==on2 .AND. sh%ct<m+2) .OR. &
            (us%dointerpolation==on3 .AND. sh%ct<m+3)

          forceinterp = .NOT. nointerpneeded .AND. us%dointerpolation /= never

! See if line search is done. First test if an interpolation has been
! done. Use the appropriate meaning of an "interpolation", i.e.
! according to Us%QuadInterpolation, either actually check for a formal
! interpolation, or else just do as Shanno and make sure at least
! 2 points have been considered.

          interpdone = (us%quadinterpolation .AND. sh%doneinterpolation) .OR. &
            ( .NOT. us%quadinterpolation .AND. .NOT. first)

          sh%lsfinished = interpdone .OR. .NOT. forceinterp .OR. dg == zero

          IF ( .NOT. sh%lsfinished) THEN
            IF (dg>zero) THEN
              upbd = alpha
              alphaubdd = true
            ELSE
              lwbd = alpha
              alphalbdd = true
              flbd = f
              dtglbd = dg
            END IF
          END IF

        ELSE ACCEPTABLE
          IF (tr%lsflow .AND. sh%validf) THEN
            WRITE (tr%u,90040) pre // ' not accepted; f<f0-abit = ', test1, &
              ' slope test = ', test2, ' AlphaUBdd = ', alphaubdd
            WRITE (tr%u,90100) pre // ' req''d reduction, f0-f, slope' // &
              ' limit = ', f0 - tp0, f0 - f, nearly1*dg0
          END IF

          sh%lsfinished = false

          IF ( .NOT. sh%validf) THEN
            upbd = alpha
            alphaubdd = false
          ELSE IF (f>=tp0) THEN
            upbd = alpha
            alphaubdd = true
          ELSE
            lwbd = alpha
            flbd = f
            dtglbd = dg
            alphalbdd = true
          END IF

        END IF ACCEPTABLE

        IF (tr%lsflow) WRITE (tr%u,90050) pre // ' done? ' // &
          'accept LSFinished ForceInterp InterpDone NoInterpNeeded', accept, &
          sh%lsfinished, forceinterp, interpdone, nointerpneeded

DONESEARCH: IF ( .NOT. sh%lsfinished) THEN

! Line search not done. A new point must be tried. Use cubic
! interpolation to find the trial point  at.

          IF (tr%lsreal) WRITE (tr%u,90060) pre // ' LwBd: ', lwbd, &
            ' AlphaLBdd: ', alphalbdd, ' UpBd: ', upbd, ' AlphaUBdd: ', &
            alphaubdd
          IF (upbd/=zero) THEN

            IF ( .NOT. alphaubdd .OR. .NOT. goodpt) THEN
              at = lwbd + bitsmall*(upbd-lwbd)
              IF (tr%lsreal) WRITE (tr%u,90010) pre // &
                ' taking midinterval' // ' alpha->', at
              interpt = false
            ELSE
              interpt = true
              IF (ap>upbd .AND. alphalbdd) THEN
                ap = lwbd
                fp = flbd
                dgp = dtglbd
              END IF
            END IF

          ELSE

            interpt = false
            left = alpha*(one+initmargin)
            right = extrapolate*alpha

            CALL cubic_interpolate(alpha,f,dg,ap,fp,dgp,left,right,at,intpt, &
              tr,us%ignoreinterval)
            sh%doneinterpolation = intpt

            IF (tr%lsreal) WRITE (tr%u,90070) pre // ' extrapolated in [', &
              left, ',', right, ']', '          to get alpha->', at, &
              ' with exact interpolate-> ', intpt
          END IF

          IF (interpt) THEN

            IF (goodpt) THEN

              slice = currmargin*(upbd-lwbd)
              left = lwbd + slice
              right = upbd - slice

              CALL cubic_interpolate(alpha,f,dg,ap,fp,dgp,left,right,at,intpt, &
                tr,us%ignoreinterval)
              sh%doneinterpolation = intpt

              IF (tr%lsreal) WRITE (tr%u,90070) pre // ' interpolating in [', &
                left, ',', right, ']', '          to get alpha->', at, &
                ' with exact interpolate-> ', intpt

              IF (intpt) THEN
                currmargin = initmargin
              ELSE
                currmargin = MIN(maxmargin,currmargin*xpndmargin)
              END IF

            ELSE
              at = lwbd + bitsmall*(upbd-lwbd)
              IF (tr%lsreal) WRITE (tr%u,90010) pre // &
                ' taking midinterval' // ' alpha->', alpha
            END IF

          END IF

          IF (sh%validf) THEN
            ap = alpha
            fp = f
            dgp = dg

            alpha = at
            goodpt = sh%validf
          ELSE
            alpha = at
            goodpt = false
          END IF

          IF (upbd/=0) THEN
            width = upbd - lwbd
          ELSE
            width = alpha - lwbd
          END IF

          IF (tr%lsreal) WRITE (tr%u,90010) pre // ' exit with alpha->', alpha
          IF (tr%lsflow) WRITE (tr%u,90080) pre // ' exit with GoodPt: ', &
            goodpt, '  DoneInterpolation: ', sh%doneinterpolation
          IF (tr%lsreal) WRITE (tr%u,90010) pre // ' exit with width->', width
        END IF DONESEARCH

! EXIT:

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:

90000   FORMAT (A)
90010   FORMAT (A,G15.7)
90020   FORMAT (A/3L11,10X,2I4)
90030   FORMAT (A/5X,2I11,3L13)
90040   FORMAT (A,L1,A,L1,A,L1)
90050   FORMAT (A/7X,5L10)
90060   FORMAT (A,G15.7,A,L1,A,G15.7,A,L1)
90070   FORMAT (2(A,G15.7),A/A,G15.7,A,L1)
90080   FORMAT (A,L1,A,L1)

90090   FORMAT ('     (valid data = ',L1,')      alpha       ', &
          '        f          dir''l derivative'/'      first   point ', &
          3G19.11/'      last    point ',3G19.11/'      current point ', &
          3G19.11)
90100   FORMAT (A,3G11.3)

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> product.f90
    MODULE product ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving, writevec
      USE min_codes, ONLY : gammafirst
      USE h0v, ONLY : hmatrix, h0v_multiply
      USE min_states, ONLY : userdefined, hterm, tracelist, sharedvars
      USE precisions, ONLY : stnd, long
      USE inner_product, ONLY : OPERATOR (.gip.)



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: product_multiply, hmatrix, userdefined, hterm, &
                tracelist, sharedvars

    CONTAINS

      SUBROUTINE product_multiply(h,v,hv,us,sh,tr)


! ARGUMENTS:

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h
        REAL (stnd), INTENT (IN) :: v(:)
        REAL (stnd), INTENT (OUT) :: hv(SIZE(v))
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr

! AUTOMATIC ARRAYS:

        REAL (stnd) :: aro(h%nterms) ! temporary vector

! DESCRIPTION:

!     Given the quasi-Newton update matrix  H (in product form) and
!     given the vector v, this routine computes

!               hv = H * v  .

!     If no terms have been allocated for H, then H is just H%H0.
!     If this is also unallocated, H is just taken as the identity I.

!     Each update "term" of H requires 2n+1 memory locations. The order is

!                   eta[i],  s[i] and  y[i].

!     Here    n    = the dimension of the problem
!             s    = x[i] - x[i-1] = alpha * d
!             y    = g[i] - g[i-1]
!             eta  = s' * y

!     These terms are stored in circular, i.e. wraparound, fashion
!     in a linked list.

!     The following variables are defined in the control record Us for
!     the minimization process, and are described in detail elsewhere:

!       beta          It is the parameter defining the Broyden family of
!                     updates. Beta must be 1.0 for product updates.
!       GammaScaling  The value GammaAll is not allowed with product form
!                     updates, although the value GammaFirst is permitted.

!     The tracing vector Tr is accessed as well.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'PrMult'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: k, n

        REAL (stnd) :: sv

        REAL (stnd), POINTER :: s(:), y(:), hi(:)
        REAL (stnd), POINTER :: eta

        TYPE (hterm), POINTER :: term

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

ISEMPTY: IF (ASSOCIATED(h%firstterm)) THEN

          term => h%firstterm%prev ! First half of Nocedal's recursion.
          n = SIZE(v)

          CALL assignhi

          sv = s .gip. v
          k = h%nterms
          aro(k) = sv/eta
          hv = v - aro(k)*y

REVERSE:  DO WHILE ( .NOT. ASSOCIATED(term,h%firstterm)) !circular list

            term => term%prev
            CALL assignhi ! remaining iterations of the first half.

            sv = hv .gip. s
            k = k - 1
            aro(k) = sv/eta
            hv = hv - aro(k)*y

          END DO REVERSE

! Set hv = H0 * hv. H0 is the initial positive definite matrix.
! Array Temp is used here to store the product of H*hv, which is
! then assigned to hv.

          CALL h0v_multiply(h,hv,h%temp,us,tr)
          hv = h%temp

          IF (us%scalegamma==gammafirst) THEN
            hv = (eta/(y.gip.y))*hv
          END IF

FORWARD:  DO
            CALL assignhi
            hv = hv + (aro(k)-(y.gip.hv)/eta)*s
            k = k + 1
            term => term%next ! terms of the second half of the product.
            IF (ASSOCIATED(term,h%firstterm)) EXIT
          END DO FORWARD

        ELSE

          CALL h0v_multiply(h,v,hv,us,tr)

        END IF ISEMPTY

        IF (tr%vectors .AND. tr%values) CALL writevec(hv,f='g',unit=tr%u, &
          name=pre//' H*v->',indent=-tr%level)

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:  none.

      CONTAINS

        SUBROUTINE assignhi ! 'break up'  term into its parts

          hi => term%data
          eta => hi(sh%ixeta)
          s => hi(sh%ixs+1:sh%ixs+n)
          y => hi(sh%ixy+1:sh%ixy+n)

        END SUBROUTINE

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> qnewton.f90
    MODULE qnewton ! 'USE'ed in modules MinimizeF, Dynamic, Update
      USE general, ONLY : entering, leaving, writevec
      USE precisions, ONLY : stnd, long
      USE min_states, ONLY : hmatrix, tracelist
      USE inner_product, ONLY : OPERATOR (.ip.)


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: qnewton_multiply, hmatrix, tracelist

    CONTAINS

      SUBROUTINE qnewton_multiply(h,v,hv,tr)


! ARGUMENTS:

        TYPE (hmatrix), INTENT (IN) :: h
        REAL (stnd), INTENT (IN) :: v(:)
        REAL (stnd), INTENT (OUT) :: hv(SIZE(v))
        TYPE (tracelist), INTENT (INOUT) :: tr

! AUTOMATIC ARRAYS:

        INTEGER (long) :: p(SIZE(v))

! DESCRIPTION:

!     Given the quasi-Newton update matrix  H (stored in row order
!     as the upper half of a symmetric matrix) and given the vector v,
!     this routine computes

!               hv = H * v  .

!     If no storage has been allocated for H%Hi, then H is just H%H0,
!     and, if this is also unallocated, H is just taken as the identity I.

!     The tracing vector Tr is also accessed.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'QNMult'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: i, j, k, n

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (ASSOCIATED(h%hi)) THEN
          n = SIZE(v)

          j = 1
          DO k = 1, n
            p(1:k-1) = p(1:k-1) + 1
            p(k:n) = (/ (i,i=j,j+n-k) /)
            hv(k) = h%hi(p) .ip. v
            j = j + n - k + 1
          END DO
        ELSE IF (ASSOCIATED(h%h0)) THEN
          hv = h%h0*v
        ELSE
          hv = v
        END IF

        IF (tr%vectors .AND. tr%values) CALL writevec(hv,f='g',unit=tr%u, &
          name=pre//' H*v->',indent=-tr%level)

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:  none.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> dynamic.f90
    MODULE dynamic_m ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving
      USE qnewton, ONLY : hmatrix, qnewton_multiply
      USE min_codes, ONLY : sumform, productform
      USE reals, ONLY : zero, one
      USE precisions, ONLY : stnd, long
      USE min_states, ONLY : userdefined, tracelist, sharedvars
      USE inner_product, ONLY : OPERATOR (.gip.)



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: dynamic_update_h, hmatrix, userdefined, tracelist, &
                sharedvars

    CONTAINS

      SUBROUTINE dynamic_update_h(s,y,u,h,us,sh,tr)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: s(:), y(SIZE(s))

        REAL (stnd), INTENT (INOUT) :: u(SIZE(s))

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h

        TYPE (userdefined), INTENT (INOUT) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr

! DESCRIPTION:

!     The basic purpose of this routine is to compute the value of the
!     update matrix  H  at the latest point when switching to the QN method.

! SUBROUTINES:  QNewton_Multiply
!               .GIP.

        CHARACTER (*), PARAMETER :: id = 'Dynamic_Update'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: k, j, n
        REAL (stnd) :: nu, eta
        REAL (stnd), POINTER :: term(:), s1(:), y1(:)

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        n = SIZE(s)

!   QN   case: symmetric half update stored.

        h%hi = zero ! Initialize H.
        j = 1
        DO k = 1, n
          h%hi(j) = one ! Note: only half of H.
          j = j + n - k + 1
        END DO

! *** First iteration ***

        term => h%firstterm%data

        SELECT CASE (us%updateform)
        CASE (sumform)
          y1 => term(sh%ixu+1:sh%ixu+n)
          s1 => term(sh%ixs+1:sh%ixs+n)
          nu = term(sh%ixnu)
        CASE (productform)
          s1 => term(sh%ixs+1:sh%ixs+n)
          y1 => term(sh%ixy+1:sh%ixy+n)
          nu = y1 .gip. y1
        END SELECT

! Scale initial quasi-Newton matrix.
! Store the initial Hessian, which is  H = (s'y/y'y)*I =
! (eta/nu)*I.  Then we need to recalculate the initial nu =
! y'*H*y = (eta/nu)*(nu above) = eta, and to find u = H*y = nu*Iy = nu *yy.

        eta = term(sh%ixeta)
        h%hi = (eta/nu)*h%hi

! First calculate u = H0*y and nu = y'*H0*y.  Remember that only
! the symmetric upper half of H is stored (in row order).

        CALL qnewton_multiply(h,y1,u,tr)
        nu = y1 .gip. u ! Calculate  nu = y'*H*y

! Now calculate the updated approximate Hessian H^.
! Use the BFGS update: nu, eta and H*y are known.

        y1 = -u + (one+nu/eta)*s1 ! Overwrite u or y.

        j = 1
        DO k = 1, n
          h%hi(j:j+n-k) = h%hi(j:j+n-k) + ((s1.gip.k)/eta)*y1(k:n) - &
            ((u.gip.k)/eta)*s1(k:n)
          j = j + n - k + 1
        END DO

! *** Second iteration ***

        term => h%firstterm%next%data

        s1 => term(sh%ixs+1:sh%ixs+n)

        CALL qnewton_multiply(h,y,u,tr) ! Calculate u = H1*y
        nu = y .gip. u ! Calculate nu = y'*H1*y

! Now calculate the updated approximate Hessian H^ using the BFGS update.

        y1 = -u + (one+nu/sh%eta)*s1

        j = 1
        DO k = 1, n
          h%hi(j:j+n-k) = h%hi(j:j+n-k) + ((s1.gip.k)/sh%eta)*y1(k:n) - &
            ((u.gip.k)/sh%eta)*s1(k:n)
          j = j + n - k + 1
        END DO

! EXIT:

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:  None.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> restart.f90
    MODULE restart ! 'USE'ed in module MinimizeF
      USE min_codes, ONLY : sumform, badcheckopen, productform, factoredform
      USE general, ONLY : entering, leaving
      USE true_false, ONLY : false, true
      USE min_defaults, ONLY : checkunit
      USE min_states, ONLY : hterm
      USE precisions, ONLY : stnd, long, short
      USE systemstate, ONLY : minimizestate


      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: restart_run, minimizestate

    CONTAINS

      SUBROUTINE restart_run(x,fx,g,state,c,file,finalstate,inunit) ! Optional argument


! ARGUMENTS:

        REAL (stnd), INTENT (OUT) :: x(:), fx, g(SIZE(x))
        INTEGER (short), INTENT (OUT) :: state
        TYPE (minimizestate), INTENT (OUT) :: c
        CHARACTER (*), INTENT (IN) :: file
        INTEGER (short), INTENT (INOUT) :: finalstate

        INTEGER (short), INTENT (IN), OPTIONAL :: inunit

! DESCRIPTION:

!     This routine reads checkpointed data from the specified checkpoint
!     file and prepares to continue the computation from that point.

!     The arguments and the state vector for routine Minimize_f (except FName)
!     are restored by reading checkpointed data from
!     the named File.  File *must* be provided; the user must specify
!     exactly what file is to be used for resuming the computation.  Note that
!     two files were used for checkpointing: one with an 'a' at the end of
!     its name and one with a 'b' at the end.  The argument File must include
!     the 'a' or the 'b'.  If the specified File cannot be successfully opened,
!     the routine will attempt to open the other checkpoint file.
!     If neither checkpoint file can be opened, or if any error occurs
!     in reading the data from an open file, an error code
!     is set in FinalState and the main routine will immediately
!     terminate without any attempt to continue the computation.

!     A unit number for reading the file may be given in the optional
!     argument InUnit.
!     It is assumed that the memory is available as per the initial run.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'Restart'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2) :: pre

        INTEGER (long) :: i, lent, end
        INTEGER (short) :: unit
        LOGICAL :: used_a, isthere

        TYPE (hterm), POINTER :: term, term0

! EXECUTION:

        IF (PRESENT(inunit)) THEN
          unit = inunit
        ELSE
          unit = checkunit
        END IF

! Open specified file.

        lent = LEN_TRIM(file)
        used_a = file(lent:lent) == 'a'

        OPEN (unit,file=file(1:lent),form='unformatted',status='old',err=10, &
          position='rewind')

        GO TO 20

! If the specified checkpoint file cannot be opened, try opening the
! alternate file.

10      IF (used_a) THEN
          OPEN (unit,file=file(1:lent-1)//'b',form='unformatted',status='old', &
            err=30,position='rewind')
          used_a = false
        ELSE
          OPEN (unit,file=file(1:lent-1)//'a',form='unformatted',status='old', &
            err=30,position='rewind')
          used_a = true
        END IF

! Read saved values.

20      READ (unit,err=30) x, fx
        READ (unit,err=30) g, state
        READ (unit,err=30) c%trace

        c%trace%level = 2 ! Indent one level
        pre = entering(id,c%trace%flow,c%trace%u,c%trace%level)

        READ (unit,err=30) c%user
        READ (unit,err=30) c%shared
        READ (unit,err=30) c%lssave
        READ (unit,err=30) c%local
        READ (unit,err=30) c%eval
        READ (unit,err=30) c%evalcts
        READ (unit,err=30) c%evalers
        READ (unit,err=30) c%print
        READ (unit,err=30) c%prvals
        READ (unit,err=30) c%term
        READ (unit,err=30) c%termvals

        ALLOCATE (c%d(1:SIZE(x)))
        ALLOCATE (c%xx(1:SIZE(x)))
        ALLOCATE (c%gg(1:SIZE(x)))
        READ (unit,err=30) c%d
        READ (unit,err=30) c%xx
        READ (unit,err=30) c%gg

        READ (unit,err=30) c%h%nterms
        READ (unit,err=30) c%h%mterms

        READ (unit,err=30) isthere
        IF (isthere) THEN
          ALLOCATE (c%h%h0(1:c%local%n))
          READ (unit,err=30) c%h%h0
        ELSE
          NULLIFY (c%h%h0)
        END IF

        READ (unit,err=30) isthere
        IF (isthere) THEN
          ALLOCATE (c%h%hi(1:(c%local%n*(c%local%n+1))/2))
          READ (unit,err=30) c%h%hi
        ELSE
          NULLIFY (c%h%hi)
        END IF

        READ (unit,err=30) isthere
        IF (isthere) THEN
          ALLOCATE (c%h%temp(1:c%local%n))
          READ (unit,err=30) c%h%temp
        ELSE
          NULLIFY (c%h%temp)
        END IF

        READ (unit,err=30) isthere
        IF (isthere) THEN
          ALLOCATE (c%h%u(1:c%local%n))
          READ (unit,err=30) c%h%u
        ELSE
          NULLIFY (c%h%u)
        END IF

        IF (c%h%nterms/=0) THEN
          ALLOCATE (c%h%firstterm)
          term => c%h%firstterm

          IF (c%local%all_allocated) THEN
            end = c%h%mterms
          ELSE
            end = c%h%nterms
          END IF

          DO i = 1, end

            ALLOCATE (term%data(c%shared%ixstart:c%shared%ixend))
            READ (unit,err=30) term%data

            SELECT CASE (c%user%updateform)
            CASE (sumform,factoredform) ! linear list - always
              NULLIFY (term%prev) ! nullify prev pointer.
            CASE (productform) ! circular list - set
              IF (i>=2) term%prev => term0 ! prev pointer.
            END SELECT

            term0 => term

            IF (i==c%h%nterms) c%h%lastterm => term ! Only used in SumForm
! and FactoredForm.

            IF (i/=end) THEN
              ALLOCATE (term%next)
              term => term%next
            END IF

          END DO

          SELECT CASE (c%user%updateform)
          CASE (sumform,factoredform) ! linear list
            NULLIFY (term%next) ! set next pointer for last term

          CASE (productform) ! circular list
            c%h%firstterm%prev => term ! set prev pointer for first term
            term%next => c%h%firstterm ! set next pointer for last term
! set pointer for next term
            IF (c%h%nterms==c%h%mterms) c%h%nextterm => c%h%firstterm%prev
          END SELECT

        ELSE
          NULLIFY (c%h%firstterm)
        END IF

        CLOSE (unit)

        c%local%checkflag = .NOT. used_a

        CALL leaving(id,c%trace%flow,c%trace%u,c%trace%level)
        c%trace%level = 0 ! This is correct value for return to Minimize_f

! EXIT
        RETURN

30      finalstate = badcheckopen
        RETURN

! FORMATS:  none.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> sum.f90
    MODULE sum ! 'USE'ed in modules MinimizeF and Update
      USE general, ONLY : entering, leaving, writevec
      USE true_false, ONLY : true, false
      USE reals, ONLY : one, two
      USE min_codes, ONLY : gammaall, gammafirst
      USE h0v, ONLY : hmatrix, h0v_multiply
      USE min_states, ONLY : userdefined, hterm, tracelist, sharedvars
      USE precisions, ONLY : stnd, long
      USE inner_product, ONLY : OPERATOR (.gip.)



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: sum_multiply, hmatrix, userdefined, hterm, &
                tracelist, sharedvars

    CONTAINS

      SUBROUTINE sum_multiply(h,v,hv,us,sh,tr)


! ARGUMENTS:

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h
        REAL (stnd), INTENT (IN) :: v(:)
        REAL (stnd), INTENT (OUT) :: hv(SIZE(v))
        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (IN) :: sh
        TYPE (tracelist), INTENT (INOUT) :: tr


!   DESCRIPTION:

!     Given the quasi-Newton update matrix  H (in sum form) and
!     given the vector v, this routine computes

!               hv = H * v  .

!     Each set of entries is called a "term" of the update,
!     and each update "term" of H requires 2n+2 memory locations

!                   nu[i], eta[i],  u[i] and  s[i].

!     Here     n    = the dimension of the problem
!              s    = x[i] - x[i-1] = alpha * d
!              y    = g[i] - g[i-1]
!              u    = H * y
!              nu   = y' * H * y
!              eta  = s' * y

!     If no storage has been allocated for the terms of H, then H is just H%H0.
!     If this is also unallocated, then H is just taken as the identity I.

!     The updates are stored in a linear linked list with at most m entries.

!     The following variables are defined in the control record Us for
!     the minimization process, and are described in detail elsewhere:

!       beta          It is the parameter defining the Broyden family of
!                     updates. The form used at each point is
!                            H^ = H(DFP) + beta * nu * w'w
!                     so that beta = 1 gives the BFGS update.
!       GammaScaling  controls the Shanno, Oren, Spedicato scaling of each term.

!     The tracing vector Tr is accessed as well.

! PARAMETERS:

        CHARACTER (*), PARAMETER :: id = 'SUMult'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long) :: n, i
        LOGICAL :: first

        REAL (stnd) :: sv, uv, beta, gamma, sigma, mu

        REAL (stnd), POINTER :: s(:), u(:), hi(:)
        REAL (stnd), POINTER :: eta, nu

        TYPE (hterm), POINTER :: term

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        CALL h0v_multiply(h,v,hv,us,tr) ! Set hv = H0 * v.

        beta = us%beta

!   Compute the terms of the product.

        n = SIZE(v)
        term => h%firstterm
        first = true

        DO i = 1, h%nterms ! linear linked list

          hi => term%data ! Compute next term and add into hv.

          nu => hi(sh%ixnu)
          eta => hi(sh%ixeta)
          u => hi(sh%ixu+1:sh%ixu+n)
          s => hi(sh%ixs+1:sh%ixs+n)

          IF (tr%values) WRITE (tr%u,90000) pre // ' nu= ', nu, '  eta=', eta

          uv = u .gip. v
          sv = s .gip. v

          IF (tr%values) WRITE (tr%u,90000) pre // ' uv= ', uv, '   sv= ', sv

! If gamma scaling is required, set gamma = eta/nu, and use the
! modified update formula which can be derived from Shanno's
! work. This only applies to the BFGS update.

          IF ((beta==one) .AND. ((us%scalegamma==gammaall) .OR. (us%scalegamma &
              ==gammafirst .AND. first))) THEN

            gamma = eta/nu
            IF (tr%values) WRITE (tr%u,90010) pre // ' gamma= ', gamma

            hv = gamma*hv
            mu = -sv/nu
            sigma = (two*sv/eta) - (uv/nu)

          ELSE IF (beta==one) THEN
            mu = -sv/eta
            sigma = -(one+nu/eta)*mu - uv/eta
          ELSE
            mu = ((beta-one)*uv/nu) - (beta*sv/eta)
            sigma = sv*(eta+beta*nu)/(eta*eta) - (beta*uv/eta)
          END IF

          IF (tr%values) WRITE (tr%u,90000) pre // ' mu= ', mu, '  sigma= ', &
            sigma

          hv = hv + mu*u + sigma*s

          IF (tr%vectors .AND. tr%values) CALL writevec(hv,f='g',unit=tr%u, &
            name=pre//' h*v->',indent=-tr%level)

          term => term%next
          first = false
        END DO

        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:

90000   FORMAT (2(A,G15.7))
90010   FORMAT (A,G15.7)

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> updateh.f90
    MODULE update ! 'USE'ed in module MinimizeF
      USE general, ONLY : entering, leaving
      USE reals, ONLY : one
      USE qnewton, ONLY : qnewton_multiply
      USE min_codes, ONLY : sumform, productform, factoredform
      USE sum, ONLY : hmatrix, sum_multiply
      USE true_false, ONLY : false, true
      USE precisions, ONLY : stnd, long
      USE min_states, ONLY : userdefined, tracelist, minlocal, sharedvars
      USE inner_product, ONLY : OPERATOR (.gip.)



      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: update_h, hmatrix, userdefined, tracelist, &
                minlocal, sharedvars

    CONTAINS

      SUBROUTINE update_h(s,y,u,h,us,sh,lc,tr)


! ARGUMENTS:

        REAL (stnd), INTENT (IN) :: s(:)

        REAL (stnd), INTENT (INOUT) :: u(SIZE(s))

        REAL (stnd), INTENT (IN) :: y(SIZE(s))

        TYPE (hmatrix), INTENT (INOUT), TARGET :: h

        TYPE (userdefined), INTENT (IN) :: us
        TYPE (sharedvars), INTENT (INOUT) :: sh
        TYPE (minlocal), INTENT (INOUT) :: lc
        TYPE (tracelist), INTENT (INOUT) :: tr

! DESCRIPTION:

!     The basic purpose of this routine is to compute the value of
!     the update matrix  H  at the new point.

!     In this description, "H" will denote the update matrix defined
!     when the current point is reached; "H^" will denote the update
!     matrix to be computed and used in forming the next search
!     direction.  On exit, the matrix H must have been updated to H^.

!     This routine assumes that memory for the update has been allocated
!     as required as it updates the required pointers (i.e. no allocation
!     is done in this routine).

!     On entry the following values are required.

!       s      The step taken on the iteration just completed.
!       y      The change in gradient from the previous point.
!              This may also be used as a scratch vector, in the QN case only.
!       H      The current matrix H
!       u      = Z^T * g, for the factored updates only.

!     On exit, in the sum update case, the value u is also returned, where

!       u = H*y

!     Depending on the current strategy, some of the following values must
!     be defined in Sh or Us on entry. These are explained more fully elsewhere.

!     Sh: RestartCt  The number of restarts done.
!         alpha      Step length; needed for factored updates.
!         SteepDStep True if the last search direction was steepest descent.
!         DoRestart  A flag which is true when this is a restart point.
!         QNpart     True when we are in the quasi-Newton part of the code
!         eta        = s'* y.
!     Us: CountFromRestart How many points until the next restart is forced.

!     In the quasi-Newton cases, the update will have been done in
!     place, i.e. the new matrix H^ will just have overwritten the old.

!     In the limited memory cases, another term will have been stored,
!     unless the space for updates has been used. In any case, in the
!     event of a restart, H will have redefined by a single update term.
!     There are no restarts in the product form case, and when the memory
!     limit is reached, earlier update terms are simply overwritten in a
!     circular fashion.

! SUBROUTINES:  Sum_Multiply
!               QNewton_Multiply
!               .GIP.

        CHARACTER (*), PARAMETER :: id = 'Update'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2+tr%level) :: pre

        INTEGER (long), POINTER :: m
        REAL (stnd), POINTER :: term(:)

        INTEGER (long) :: k, j, n
        REAL (stnd) :: nu

! EXECUTION:

        pre = entering(id,tr%flow,tr%u,tr%level)

        m => h%mterms
        n = SIZE(s)

STRATEGY: IF (sh%limmemory) THEN

ISRESTART: IF (sh%dorestart) THEN

            IF (tr%update) WRITE (tr%u,90000) pre // ' restart! H%nTerms->1'
            sh%restartct = sh%restartct + 1 ! Count number of restarts

! Set point to force the next restart.
            IF (us%countfromrestart) THEN
              sh%restartpt = 1 + n
            ELSE
              sh%restartpt = m + 1 + n
            END IF

            NULLIFY (h%h0) ! After restart, diagonal scaling always I.

            IF (m==0) THEN ! can't save updates; next step must be
              sh%ct = 0 ! steepest descent. Either SD or straight CG.
              h%nterms = 0
              sh%qnpart = false
            ELSE ! construct update from last two points
              sh%ct = 1
              h%nterms = 1
              h%lastterm => h%firstterm
              term => h%firstterm%data
              sh%qnpart = true

! Save the current s and u = H*y = I*y = y (Beale restart)
! vectors and save  nu = y'*H*y = y'y and eta = s'y in the
! first term of H.

              IF (us%updateform==sumform) THEN
                term(sh%ixnu) = y .gip. y
                term(sh%ixeta) = sh%eta
                term(sh%ixu+1:sh%ixu+n) = y
                term(sh%ixs+1:sh%ixs+n) = s

                IF (tr%update) WRITE (tr%u,90010) pre // ' saved: nu= ', &
                  term(sh%ixnu), '  eta= ', term(sh%ixeta)
              ELSE IF (us%updateform==factoredform) THEN
                term(sh%ixu+1:sh%ixu+n) = -lc%alpha*u
                term(sh%ixs+1:sh%ixs+n) = s/SQRT(sh%eta)
                term(sh%ixy+1:sh%ixy+n) = y/SQRT(sh%eta)
!IF ( Us%ScaleColumns) &
!term(Sh%ixc+1:Sh%ixc+1) = Lc%alpha*Lc%nrmd/sqrt(Sh%eta)
              END IF

            END IF

          ELSE ISRESTART ! not a restart
            SELECT CASE (us%updateform)

            CASE (sumform,factoredform)

! Compute u = H*y, Note that
! the computation is the same for the CG or QN parts.

              IF (us%updateform==sumform) THEN
                IF (tr%update) WRITE (tr%u,90000) pre // ' SUM form update'

                CALL sum_multiply(h,y,u,us,sh,tr)

                nu = y .gip. u
              ELSE IF (us%updateform==factoredform) THEN
                IF (tr%update) WRITE (tr%u,90000) pre // &
                  ' FACTORED form update'
              END IF

SAVE:         IF (sh%qnpart) THEN ! Save nu,eta,u and s in array H.

                h%nterms = h%nterms + 1

                IF (lc%all_allocated) THEN
                  h%lastterm => h%lastterm%next
                ELSE
                  IF ( .NOT. ASSOCIATED(h%firstterm)) THEN
                    h%firstterm => h%nextterm
                  ELSE
                    h%lastterm%next => h%nextterm
                  END IF
                  h%lastterm => h%nextterm
                  NULLIFY (h%nextterm%next)
                  NULLIFY (h%nextterm%prev)
                END IF
                term => h%lastterm%data

                IF (us%updateform==sumform) THEN
                  term(sh%ixnu) = nu
                  term(sh%ixeta) = sh%eta
                  term(sh%ixu+1:sh%ixu+n) = u
                  term(sh%ixs+1:sh%ixs+n) = s
                ELSE IF (us%updateform==factoredform) THEN
                  term(sh%ixu+1:sh%ixu+n) = -lc%alpha*u
                  term(sh%ixs+1:sh%ixs+n) = s/SQRT(sh%eta)
                  term(sh%ixy+1:sh%ixy+n) = y/SQRT(sh%eta)
!IF ( Us%ScaleColumns) &
!term(Sh%ixc+1:Sh%ixc+1) = Lc%alpha*Lc%nrmd/sqrt(Sh%eta)
                END IF

                lc%all_allocated = .NOT. lc%unknownm .AND. h%nterms >= m

                IF (us%updateform==sumform .AND. tr%update) THEN
                  WRITE (tr%u,90020) pre // ' no restart; H%nTerms =', &
                    h%nterms
                  WRITE (tr%u,90010) pre // ' saved: nu= ', term(sh%ixnu), &
                    '  eta= ', term(sh%ixeta)
                ELSE IF (us%updateform==factoredform .AND. tr%update) THEN
                  WRITE (tr%u,90000) pre // ' saved Powell update term'
                END IF
              ELSE SAVE
                IF (tr%update) WRITE (tr%u,90020) pre // &
                  ' no restart; H%nTerms stays at', h%nterms
              END IF SAVE

            CASE (productform) ! circular linked list; no restarts

              IF (tr%update) WRITE (tr%u,90000) pre // ' PRODUCT form update'

              IF (h%nterms==m) THEN
                h%firstterm => h%firstterm%next
                h%nextterm => h%nextterm%next
              ELSE IF ( .NOT. ASSOCIATED(h%firstterm)) THEN
                h%firstterm => h%nextterm
                h%firstterm%prev => h%firstterm
                h%firstterm%next => h%firstterm
                h%nterms = 1
              ELSE
                h%firstterm%prev%next => h%nextterm
                h%nextterm%next => h%firstterm
                h%nextterm%prev => h%firstterm%prev
                h%firstterm%prev => h%nextterm
                h%nterms = h%nterms + 1
              END IF
              term => h%nextterm%data
              term(sh%ixeta) = sh%eta
              term(sh%ixs+1:sh%ixs+n) = s
              term(sh%ixy+1:sh%ixy+n) = y

              lc%all_allocated = .NOT. lc%unknownm .AND. h%nterms >= m

              IF (tr%update) WRITE (tr%u,90000) pre // &
                ' saved Nocedal update term'

            END SELECT
          END IF ISRESTART

        ELSE STRATEGY !   QN   case: symmetric half update stored.

          IF (tr%flow) WRITE (tr%u,90030) pre // ' SteepDStep= ', &
            sh%steepdstep

          IF (sh%steepdstep) THEN !  scale initial quasi-Newton matrix

! Store the initial Hessian, which is  H = (s'y/y'y)*I =
! (eta/nu)*I.  Then we need to recalculate the initial nu :=
! y'*H*y = (eta/nu)*(nu) = eta, and to find u = H*y =
! nu*Iy = nu *y.

            nu = y .gip. y ! Calculate  nu = y'*H*y = y'*I*y = y'*y
            h%hi = (sh%eta/nu)*h%hi !help - remove to compare QN to Factored
          END IF

! First calculate u = Hi*y and nu = y'*Hi*y.  Remember that only
! the symmetric upper half of H is stored (in row order).

          CALL qnewton_multiply(h,y,u,tr)
          nu = y .gip. u ! Calculate  nu = y'*H*y = y'*u

! Now calculate the updated approximate Hessian H^.
! Use the BFGS update: nu, eta and H*y are known.

          h%temp = -u + (one+nu/sh%eta)*s

          j = 1
          DO k = 1, n
            h%hi(j:j+n-k) = h%hi(j:j+n-k) + ((s.gip.k)/sh%eta)*h%temp(k:n) - &
              ((u.gip.k)/sh%eta)*s(k:n)
            j = j + n - k + 1
          END DO

        END IF STRATEGY
        CALL leaving(id,tr%flow,tr%u,tr%level)
        RETURN

! FORMATS:

90000   FORMAT (A)
90010   FORMAT (2(A,G15.7))
90020   FORMAT (A,I8,A)
90030   FORMAT (A,G15.7)

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> minimizef.f90
    MODULE minimizef ! 'USE'ed in module Minimize
      USE general, ONLY : entering, leaving, writevec, cpusecs
      USE myallocate, ONLY : userdefined, my_deallocate, my_allocate, &
        sharedvars
      USE min_codes, ONLY : excessfevals, revcommg, revcommf, nodescent, &
        invalidm, allocnext, revcommfandg, lsfail, initialmin, initialundefd, &
        allocu, alloctemp, storageerror, allochi, badinput, revcommnoforg, &
        revcommnog, revcommnof, revcommrestart, resume, badseth0, computed3, &
        computed2, computed1, diagonal, ident, badbeta, badrho, badstep, &
        powell, fletcher, badupdate, factoredform, sumform, badhtest, useh, &
        usei, norestart, badgamma, gammaall, gammafirst, nogammascale, &
        badalpha, alwaysqn, before3qn, before2qn, before1qn, neverqn, &
        badinterp, never, on3, on2, on1, every, badtermunit, badnorm, &
        badprintunit, badevalunit, badexpense, badscalef, badsysmemory, &
        badderiv, badterms, badmethod, badcombination, productform, qn, &
        dynamic, available, variable, fixterms, conmin, cg, sd, badcheckunit, &
        badcheckpoint, badfreq, badevallim, badmemory, badstate, revcommstart, &
        normalwithfg, normal, badacc, badn, abortmin
      USE errorcodes, ONLY : error_codes
      USE checkpoint, ONLY : check_point
      USE dynamic_m, ONLY : dynamic_update_h
      USE qnewton, ONLY : qnewton_multiply
      USE factored, ONLY : factored_multiply
      USE product, ONLY : product_multiply
      USE update, ONLY : hmatrix, update_h
      USE sum, ONLY : sum_multiply
      USE min_defaults, ONLY : defaultuserdefined, maxupdates, defaulttraces
      USE linesearch, ONLY : searchvals, line_search
      USE reals, ONLY : zero, two, five, c0_9, c180, one
      USE h0v, ONLY : h0v_multiply
      USE initialize, ONLY : initialize_h
      USE support, ONLY : evaluate_f, print_iterate, test_done
      USE control, ONLY : control_codes
      USE supp_codes, ONLY : both, noforg, nog, nof, g2, linf, l2, l1, &
        differences, firsttest, comparetest, analytic, abort, limit, ok
      USE num_constants, ONLY : machtol
      USE restart, ONLY : minimizestate, restart_run
      USE systemstate, ONLY : tracelist, termvalues, termstate, printstate, &
        evalstate, minlocal
      USE precisions, ONLY : stnd, long, short
      USE min_states, ONLY : numerrs
      USE supp_defs, ONLY : defaultevalstate, defaulttermstate, initprvals, &
        defaultprintstate
      USE true_false, ONLY : true, false
      USE inner_product, ONLY : OPERATOR (.gnorm.), OPERATOR (.ip.), &
        OPERATOR (.gip.)

      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: minimize_f, minimizestate, tracelist

    CONTAINS

! Optional Arguments
! Commonly changed
! For Evaluate_f
! For Print_Iterate
! For Test_Done
! For Minimization
      RECURSIVE SUBROUTINE minimize_f(f,x,fx,g,accuracy,state,memory,c, &
          evallimit,frequency,checkpoint,checkfile,checkunit,method,terms, &
          derivatives,decreaseinf,systemmemory,exactls,scalef,expense,tracef, &
          traceg,tracedervtest,evaltraceunit,printunit,printgrad,printx, &
          usegrad,usestep,useshanno,usefunc,thenorm,traceterm,termtraceunit, &
          relativetof0,relativetog0,dointerpolation,quadinterpolation, &
          startalpha,scalegamma,htest,updateform,startstep,ignoreinterval, &
          countfromrestart,rho,beta,seth0,parh0,scalecolumns,traces)


! ARGUMENTS:

        INTERFACE
          SUBROUTINE f(x,fx,g,job)
            USE precisions, ONLY : stnd, short
            REAL (stnd), INTENT (IN) :: x(:)
            REAL (stnd), INTENT (OUT) :: fx, g(SIZE(x))
            INTEGER (short), INTENT (INOUT) :: job
          END SUBROUTINE
        END INTERFACE

        REAL (stnd), INTENT (INOUT) :: fx, x(:), g(SIZE(x))

        REAL (stnd), INTENT (IN) :: accuracy

        INTEGER (short), INTENT (INOUT) :: state
        INTEGER (long), INTENT (IN) :: memory
        TYPE (minimizestate), INTENT (INOUT), TARGET :: c

        CHARACTER (*), INTENT (IN), OPTIONAL :: checkfile

        INTEGER (long), INTENT (IN), OPTIONAL :: evallimit, frequency, &
          checkpoint, terms, systemmemory

        INTEGER (short), INTENT (IN), OPTIONAL :: checkunit, method, &
          derivatives, expense, thenorm, scalef, evaltraceunit, printunit, &
          termtraceunit, dointerpolation, startalpha, scalegamma, htest, &
          updateform, startstep, seth0

        TYPE (tracelist), INTENT (IN), OPTIONAL :: traces

        LOGICAL, INTENT (IN), OPTIONAL :: exactls, tracef, traceg, &
          tracedervtest, printgrad, printx, usegrad, usestep, useshanno, &
          usefunc, traceterm, quadinterpolation, ignoreinterval, &
          countfromrestart, relativetof0, relativetog0, scalecolumns

        REAL (stnd), INTENT (IN), OPTIONAL :: decreaseinf, rho, beta, parh0


! DESCRIPTION:
!                 See      "min.doc"

! SUBROUTINES:

!     F                    external procedure passed as argument.

!     Initialize_H         initial diagonal matrix
!     Cubic_Interpolation  cubic interpolation
!     Line_Search          line search loop
!     Sum_Multiply         matrix vector multiplication with sums
!     Product_Multiply     matrix vector multiplication with products
!     Factored_Multiply    matrix vector multiplication with factors
!     Update_H             update H

!     Evaluate_f           See Support
!     Print_Iterate        See Support
!     Test_Done            See Support

! LOCAL DECLARATIONS:

!-----General Declarations.

        INTEGER (short) :: entrystate, finalstate, enough, which

        INTEGER (long) :: i, qnstorage

        LOGICAL :: less, toosmall, qnstep, force1, forcerestart, &
          dorestarttest, lots

        REAL (stnd) :: dgal, tp0, tp1, tp2, width, stg, utg, nu, sigma, gamma, &
          mu

        TYPE (userdefined), POINTER :: us
        TYPE (sharedvars), POINTER :: sh
        TYPE (searchvals), POINTER :: ls
        TYPE (minlocal), POINTER :: lc
        TYPE (evalstate), POINTER :: ev
        TYPE (printstate), POINTER :: pr
        TYPE (termstate), POINTER :: tm
        TYPE (termvalues), POINTER :: tv
        TYPE (tracelist), POINTER :: tr
        TYPE (hmatrix), POINTER :: h
        INTEGER (long), POINTER :: m

! PARAMETERS:

        REAL (stnd), PARAMETER :: nearly_1 = c0_9

        CHARACTER (*), PARAMETER :: id = 'Min'

! LOCAL DECLARATIONS:

        CHARACTER (len=LEN(id)+2) :: pre

! EXECUTION:

!>>>>>>>>>> PHASE  0:  Describe Phases.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!  The code is in "phases".  The flow is forward to the end in each
!  phase. All phases are exited only at the end of the phase and
!  flow proceeds to the start of another phase, or it exits  the
!  algorithm to statement  900.

!  An exception is a jump to 920 and a return in phase VII
!  if reverse communication is being used, along with a reentry
!  from the top of phase I back to continue from the point of exit
!  at 215.

!  A second exception is in the resumption of a run from a checkpointed
!  file.  A jump to the start of phase III occurs in this case.

!>>>>>>>>>> PHASE  I:  Initial Set Up.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        us => c%user
        sh => c%shared
        ls => c%lssave
        lc => c%local
        ev => c%eval
        pr => c%print
        tm => c%term
        tv => c%termvals
        tr => c%trace
        h => c%h
        m => c%h%mterms

        IF (PRESENT(traces)) THEN
          c%trace = traces
        ELSE
          c%trace = defaulttraces
        END IF

        entrystate = state
        finalstate = ok

        IF (state==resume) THEN ! Resume run from a checkpointed file
          CALL restart_run(x,fx,g,state,c,checkfile,finalstate, &
            inunit=checkunit)

          IF (PRESENT(traces)) c%trace = traces ! allows tracing to be done

          entrystate = state
          IF (finalstate==ok) THEN
            pre = entering(id,tr%flow,tr%u,tr%level)
            GO TO 20 ! go to point just after routine Check_Point called
          ELSE
            GO TO 70 ! exit with error
          END IF
        END IF

        pre = entering(id,tr%flow,tr%u,tr%level)

        IF (anytr()) WRITE (tr%u,90000) tr

        lc%noprint = false

        IF (entrystate==revcommrestart .OR. entrystate==revcommnof .OR. &
            entrystate==revcommnog .OR. entrystate==revcommnoforg) THEN

! Supplementary call with reverse communication.

          sh%validf = entrystate == revcommrestart
          GO TO 40
        END IF

        lc%eps = five*machtol

        lc%checkflag = true

        lc%n = SIZE(x)
        CALL checkdata ! Check all parameter settings valid.

        IF (sh%error(0)/=0) THEN
          finalstate = badinput
          lc%noprint = true
        END IF

        IF (finalstate==ok) THEN
          CALL initialize_min

          lc%firsteval = true

          sh%restartct = 0 ! Initialize counts
          sh%forcect = 0

          lc%all_allocated = false
          lc%unknownm = false

          sh%memoryused = 0
          h%nterms = 0

          enough = 0

          SELECT CASE (us%updateform)
          CASE (sumform)
            sh%ixleng = 2*lc%n + 2
            sh%ixstart = -1
            sh%ixnu = -1
            sh%ixeta = 0
            sh%ixu = sh%ixeta
            sh%ixs = sh%ixu + lc%n
          CASE (productform)
            sh%ixleng = 2*lc%n + 1
            sh%ixstart = 0
            sh%ixeta = 0
            sh%ixs = sh%ixeta
            sh%ixy = sh%ixs + lc%n
          CASE (factoredform)
            sh%ixleng = 3*lc%n
            sh%ixstart = 1
            sh%ixs = sh%ixstart - 1
            sh%ixu = sh%ixs + lc%n
            sh%ixy = sh%ixu + lc%n
!IF ( Us%ScaleColumns ) then
!    Sh%ixleng = Sh%ixleng + Lc%n
!    Sh%ixc    = Sh%ixy    + Lc%n
!end if
          END SELECT

          sh%ixend = sh%ixstart + sh%ixleng - 1

          sh%methodinuse = us%method
          sh%iterationtime = zero
          sh%restartpt = zero

          SELECT CASE (sh%methodinuse)
!      ----     Allocate Memory

          CASE (sd,cg)
            sh%limmemory = true
            m = 0
            lc%all_allocated = true
          CASE (conmin)
            sh%limmemory = true
            c%h%mterms = 1
            us%dointerpolation = every
            us%startalpha = before2qn
            us%quadinterpolation = false
            us%ignoreinterval = false
          CASE (fixterms) ! m is set in routine Initialize_Min.
            sh%limmemory = true
          CASE (variable)
            sh%limmemory = true
            lc%unknownm = true
          CASE (available)
            qnstorage = (lc%n*(lc%n+1))/2 + lc%n
            IF (sh%memory>0) THEN
              lots = sh%memory >= qnstorage
            ELSE
              lots = true
            END IF
            IF (lots) THEN
              CALL my_allocate(qnstorage-lc%n,allochi,enough,us,sh,tr,h)
            END IF
            IF (lots .AND. enough==0) THEN
              sh%limmemory = false
              us%dointerpolation = never
              us%startalpha = alwaysqn
              us%quadinterpolation = false
              us%ignoreinterval = false
            ELSE
              sh%limmemory = true
              lc%unknownm = true
              us%dointerpolation = on1
              us%startalpha = before2qn
              us%quadinterpolation = true
              us%ignoreinterval = false
            END IF
          CASE (dynamic)
            sh%limmemory = true
            m = 5
            us%dointerpolation = on1
            us%startalpha = before2qn
            us%quadinterpolation = true
            us%ignoreinterval = false
          CASE (qn)
            sh%limmemory = false
            qnstorage = (lc%n*(lc%n+1))/2 + lc%n
            CALL my_allocate(qnstorage-lc%n,allochi,enough,us,sh,tr,h)
            IF (enough/=0 .OR. sh%memory>0 .AND. sh%memory<qnstorage) THEN
              lc%noprint = true
              finalstate = storageerror
            END IF

          END SELECT

          IF (PRESENT(exactls)) THEN
            IF (exactls) THEN ! Set arguments to do an exact line search
              us%dointerpolation = every
              us%startalpha = alwaysqn
              us%quadinterpolation = true
              us%ignoreinterval = true
              tm%usegrad = true
              tm%usestep = false
              tm%useshanno = false
              tm%usefunc = false
            END IF
          END IF

! Check that there is sufficient memory for limited memory cases.
! The value of MemoryNeeded only includes memory that is explicitly
! allocated.  m is 0 unless it is set otherwise.

          IF (sh%limmemory .AND. sh%memory>0) THEN
            SELECT CASE (us%updateform)
            CASE (sumform)
              lc%memoryneeded = m*sh%ixleng + 3*lc%n
            CASE (productform)
              lc%memoryneeded = m*sh%ixleng + 4*lc%n
            CASE (factoredform)
              lc%memoryneeded = m*sh%ixleng + 4*lc%n
!IF ( Us%ScaleColumns ) &
!    Lc%MemoryNeeded = Lc%MemoryNeeded + Lc%n
            END SELECT
            IF (us%seth0==diagonal) lc%memoryneeded = lc%memoryneeded + lc%n

            IF (sh%memory<lc%memoryneeded) THEN
              lc%noprint = true
              finalstate = storageerror
            END IF
          END IF

!  Allocate Temp in the QN or Limited Memory/ProductForm cases.
!  Allocate u in the Limited Memory/FactoredForm case.
!  Temp and u are never both allocated.

          IF ( .NOT. sh%limmemory .OR. us%updateform==productform) THEN
            CALL my_allocate(lc%n,alloctemp,enough,us,sh,tr,h)
          ELSE IF (us%updateform==factoredform) THEN
            CALL my_allocate(lc%n,allocu,enough,us,sh,tr,h)
          END IF
          IF (enough/=0) THEN
            lc%noprint = true
            finalstate = storageerror
          END IF

        END IF

        IF (finalstate==ok .AND. entrystate==normal) THEN

          CALL doevalf ! Get initial function value
          SELECT CASE (which)
          CASE (nof,nog,noforg)
            finalstate = initialundefd
          CASE DEFAULT
          END SELECT
          lc%firsteval = false
        END IF

        IF (finalstate==ok) THEN

          lc%nrmg = .gnorm. g

          IF (tr%values .AND. .NOT. tr%input) THEN
            WRITE (tr%u,90010) pre // ' fx        = ', fx
            WRITE (tr%u,90010) pre // ' norm of g = ', lc%nrmg
          END IF

          IF (us%relativetof0) THEN ! Initialize the termination tests.
            tm%fatx0 = fx
          ELSE
            tm%fatx0 = one
          END IF

          IF (us%relativetog0) THEN
            tm%normgatx0 = lc%nrmg
          ELSE
            tm%normgatx0 = one
          END IF

          IF (tr%input) THEN
            WRITE (tr%u,90110) lc%n, sh%memory, sh%accuracy
            IF (us%decreaseinf==zero) THEN
              WRITE (tr%u,90120) fx
            ELSE IF (us%decreaseinf<zero) THEN
              WRITE (tr%u,90130) us%decreaseinf
            ELSE
              WRITE (tr%u,90140) us%decreaseinf
            END IF

            CALL control_codes(entrystate,us,tr)
            WRITE (tr%u,90150) us%rho, us%beta, us%quadinterpolation, &
              us%countfromrestart, us%ignoreinterval, us%relativetof0, &
              us%relativetog0, sh%limmemory, lc%eps, tm%fatx0, tm%normgatx0
            IF (sh%limmemory) THEN
              WRITE (tr%u,90170) m
            ELSE
              WRITE (tr%u,90160) qnstorage
            END IF
          END IF

! Test if the initial point is the minimizer.

          IF (test_done(sh%accuracy,tm,fx,g,x,v=tv,level=tr%level)) &
            finalstate = initialmin
        END IF

        IF (finalstate/=ok) THEN
          GO TO 70
        END IF

!>>>>>>>>>>PHASE  II:  "Cold Start" With Steepest Descent.<<<<<<<<<

!     Calculate the initial search direction: dg0 is the current
!     directional derivative of f along d, while nrmg is the norm of g.

!     Initialize ct, which is used to determine whether a Beale
!     restart should be done. i.e. a restart must be forced after
!     n steps without one. Initialize SteepDStep, which indicates
!     that the current search direction is a negative gradient direction.
!     The current point is x[0].

        sh%it = -1 ! Initialize iteration counter
10      sh%steepdstep = true
        lc%cold = true
        sh%ct = 0

        sh%doforcerestart = true
        IF (sh%limmemory) THEN
          sh%restartpt = lc%n
          h%nterms = 0
          sh%qnpart = sh%methodinuse /= sd .AND. sh%methodinuse /= cg
          sh%doforcerestart = .NOT. sh%qnpart
        END IF

        sh%it = sh%it + 1
        CALL print_iterate(lc%n,sh%it,fx,lc%nrmg,c%evalcts,c%print,c%prvals,x, &
          g,first=true,level=tr%level)

        CALL initialize_h(x,g,h,finalstate,sh,us,tr)
        IF (finalstate/=ok) THEN
          lc%noprint = true
          GO TO 70
        END IF

        IF (sh%limmemory .AND. us%updateform==factoredform) THEN
          CALL h0v_multiply(h,g,c%d,us,tr,u=h%u,ka=-1)
        ELSE
          CALL h0v_multiply(h,g,c%d,us,tr)
        END IF

        c%d = -c%d

        lc%nrmd = .gnorm. c%d
        lc%dg0 = c%d .gip. g

!>>>>>>>>>> PHASE  III: Start Iteration along d[ct].<<<<<<<<<<<<<<<<

!     Begin the major iteration loop. Ncalls is used to guarantee that
!     at least two points have been tried when Method=LimMemory (see
!     DoInterpolation).  Fmin is the current function value. Force a
!     restart after n steps. Output (if desired) at end of each iteration.

20      lc%fmin = fx
        sh%ncalls = 0
        lc%nrmx = MAX(one,.gnorm.x)
        sh%doneinterpolation = false

        sh%iterationtime = sh%iterationtime - cpusecs()

        IF (tr%values) WRITE (tr%u,90010) pre // ' norm of x = ', lc%nrmx
        IF (tr%xandd) CALL writevec(x,f='g',unit=tr%u,name=pre//' x->', &
          indent=-tr%level)
        IF (tr%values) WRITE (tr%u,90010) pre // ' norm of d = ', lc%nrmd
        IF (tr%xandd) CALL writevec(c%d,f='g',unit=tr%u,name=pre//' d->', &
          indent=-tr%level)

        sh%ct = sh%ct + 1 ! So ct = index of point to which the search will lead.
!       = the index of the current search direction.

!>>>>>>>>>> PHASE  IV: Initialize Alpha for Line Search.<<<<<<<<<<<<

        sh%lsfinished = false

        IF (tr%lsflow) WRITE (tr%u,90020) pre // ' start ls :>'
        IF (tr%lsflow) WRITE (tr%u,90030) pre // &
          ' ct    QNpart  SteepDStep  cold  IgnoreInterval', sh%ct, sh%qnpart, &
          sh%steepdstep, lc%cold, us%ignoreinterval

        IF (tr%lsreal) THEN
          IF (ABS(lc%dg0)<=lc%nrmg*lc%nrmd) THEN
            WRITE (tr%u,90040) pre // ' angle(d,-g)=', &
              angle(-lc%dg0/(lc%nrmg*lc%nrmd)), ' deg.'
          ELSE
            WRITE (tr%u,90050) pre // &
              ' warning: NO ANGLE(d,-g): dg0>nrmg*nrmd', lc%dg0, &
              lc%nrmg*lc%nrmd
          END IF
        END IF

        IF (lc%cold) THEN
          IF (tr%lsflow) WRITE (tr%u,90020) pre // ' first case alpha.'

! First iteration. Scale step to one. Use estimate DecreaseInF.

          IF (us%decreaseinf==zero) THEN
            tp1 = two*ABS(fx)/lc%nrmg
          ELSE IF (us%decreaseinf>zero) THEN
            tp1 = two*us%decreaseinf/lc%nrmg
          ELSE
            tp1 = one
          END IF

          IF (sh%limmemory .AND. us%seth0==diagonal) THEN
            lc%alpha = tp1
          ELSE
            lc%alpha = tp1/lc%nrmg
          END IF

        ELSE

          qnstep = m /= 0 .AND. (us%startalpha/=neverqn) .AND. &
            (us%startalpha==before1qn .AND. sh%ct<m+1 .OR. &
            us%startalpha==before2qn .AND. sh%ct<m+2 .OR. &
            us%startalpha==before3qn .AND. sh%ct<m+3)

          force1 = qnstep .OR. (us%startalpha==alwaysqn)

          IF (force1) THEN

            lc%alpha = one
            IF (tr%lsflow) WRITE (tr%u,90020) pre // ' force alpha to 1.'

          ELSE
            IF (us%startstep==fletcher) THEN
              IF (tr%lsflow) WRITE (tr%u,90020) pre // &
                ' Fletcher scale alpha.'
              lc%alpha = lc%alpha*two*(fx-lc%flast)/(lc%dg0)
            ELSE IF (us%startstep==powell) THEN
              IF (tr%lsflow) WRITE (tr%u,90020) pre // &
                ' Sh/Powell scale alpha'
              lc%alpha = lc%alpha*(lc%dglast/lc%dg0)
            END IF
          END IF

        END IF
        IF (tr%lsflow) WRITE (tr%u,90010) pre // ' end of phase IV, alpha = ', &
          lc%alpha

!>>>>>>>>>> PHASE  V: Initialize Line Search.<<<<<<<<<<<<<<<<<<<<<<<

!   The line search fits a cubic to fx and dgal, the function and its
!   derivative at alpha, and to fp and dgp, the function and its deri-
!   vative at the previous trial point ap, where the derivatives are
!   along d.  Initialize ap, fp and dgp.

        lc%ap = zero
        lc%fp = lc%fmin
        lc%dgp = lc%dg0

!   Save the current derivative along d and the function value to
!   scale the initial step along the next search vector.

        lc%dglast = lc%dg0
        lc%flast = lc%fmin

        c%xx = x !   Store the current x and g.
        c%gg = g

!   This next little loop avoids the possibility of a
!   ridiculously small value for alpha.

        DO WHILE (fx+lc%alpha*lc%dg0>=fx+nearly_1*lc%alpha*lc%dg0)
          lc%alpha = two*lc%alpha
        END DO
        width = lc%alpha

!>>>>>>>>>> PHASE  VI: Test for Line Search Failure.<<<<<<<<<<<<<<<<

30      CONTINUE

        IF (tr%lsalpha) WRITE (tr%u,90010) pre // ' ls alpha->', lc%alpha

        IF (tr%lsreal) THEN
          WRITE (tr%u,90070) pre // ' values:   ap= ', lc%ap, '  fp    = ', &
            lc%fp
          WRITE (tr%u,90070) pre // '          dgp= ', lc%dgp, '  dglast= ', &
            lc%dglast
          WRITE (tr%u,90070) pre // '          dg0= ', lc%dg0, '  flast = ', &
            lc%flast
          WRITE (tr%u,90070) pre // '         fmin= ', lc%fmin, '  nrmd  = ', &
            lc%nrmd
        END IF

        toosmall = width*lc%nrmd <= lc%eps*lc%nrmx

        IF (toosmall) THEN

!       This is an abnormally small step. Test if the direction
!       is a gradient direction. If not, try one before aborting
!       the run; i.e. do a total restart from scratch unless this
!       step is already a steepest descent step from a cold start.

          IF (tr%lsflow) WRITE (tr%u,90020) pre // ' alpha too small.'
          IF (tr%values) WRITE (tr%u,90070) pre // ' eps= ', lc%eps, &
            ' width= ', width

          IF (lc%cold) THEN
            finalstate = lsfail
            GO TO 70
          ELSE
            GO TO 10
          END IF

        END IF

!>>>>>>>>>> PHASE  VII: Line Search Loop.<<<<<<<<<<<<<<<<<<<<<<<<<<<

!   LSFinished is set to true when the line search is deemed complete.
!   Each loop determines a new value for alpha and returns to 2000
!   unless the search has been deemed complete.

        x = c%xx + lc%alpha*c%d !   Compute the new trial point.

!   Evaluate the function at the trial point.

        IF (entrystate==revcommstart .OR. entrystate==revcommrestart .OR. &
            entrystate==revcommnoforg) THEN

          lc%noprint = true ! Exit for reverse communication.
          finalstate = revcommfandg !  (re-entry will be to 215)
          GO TO 70
        ELSE
          sh%iterationtime = sh%iterationtime + cpusecs()
          CALL doevalf
          sh%iterationtime = sh%iterationtime - cpusecs()

! If NormalWithFG, increase counts by 1 since F and G have been
! evaluated once prior to entry to Minimize_f.

          IF (entrystate==normalwithfg .AND. lc%firsteval) THEN
            c%evalcts%fevals = c%evalcts%fevals + 1
            c%evalcts%gevals = c%evalcts%gevals + 1
            lc%firsteval = false
          END IF

        END IF

40      IF (finalstate==ok) THEN

          sh%ncalls = sh%ncalls + 1

          dgal = c%d .gip. g !Compute dir'l derivative of f along d at alpha.
          lc%nrmg = .gnorm. g

          IF (tr%lsreal) THEN
            WRITE (tr%u,90010) pre // ' norm of g->', lc%nrmg

            IF (ABS(dgal)<=lc%nrmg*lc%nrmd .AND. lc%nrmg/=zero) THEN
              WRITE (tr%u,90040) pre // ' angle of d to -g->', &
                angle(-dgal/(lc%nrmg*lc%nrmd)), ' degrees'
            ELSE
              WRITE (tr%u,90070) pre // ' warning on angle of d to -g' // &
                ' dgal= ', dgal, ' > nrmg*nrmd= ', lc%nrmg*lc%nrmd
            END IF

            WRITE (tr%u,90060) pre // ' search: alpha       nrmd         ' // &
              'eps          fx', lc%alpha, lc%nrmd, lc%eps, fx
            IF (tr%vectors) CALL writevec(x,f='g',unit=tr%u,name=pre//' x-> ', &
              indent=-tr%level)
          END IF

          CALL line_search(lc%alpha,fx,dgal,lc%fmin,lc%dglast,lc%ap,lc%fp, &
            lc%dgp,width,ls,us,sh,h,tr)

          IF ( .NOT. sh%lsfinished) THEN
            DO i = 1, lc%n ! Check points not actually identical from roundoff.
              tp0 = c%xx(i) + lc%alpha*c%d(i)
              IF (tp0/=c%xx(i) .AND. tp0/=x(i)) GO TO 50
            END DO
            width = zero ! If identical, force termination with error.
50          GO TO 30
          END IF
        ELSE
          GO TO 70
        END IF

!   Flow continues to PHASE VIII if the line search is done
!   or returns to 2000 if not.

!>>>>>>>>>> PHASE  VIII: Termination Test.<<<<<<<<<<<<<<<<<<<<<<<<<<

        less = test_done(sh%accuracy,tm,fx,g,x,c%xx,tv,tr%level)

        IF (tr%update) WRITE (tr%u,90080) pre // ' term? less->', less

        IF ( .NOT. less) THEN
          sh%it = sh%it + 1
          sh%iterationtime = sh%iterationtime + cpusecs()
          CALL print_iterate(lc%n,sh%it,fx,lc%nrmg,c%evalcts,c%print,c%prvals, &
            x,g,first=false,level=tr%level)
          sh%iterationtime = sh%iterationtime - cpusecs()
        ELSE
          GO TO 70
        END IF

!>>>>>>>>>> PHASE  IX: Test if Restart Needed.<<<<<<<<<<<<<<<<<<<<<<

! Search continues.

        c%d = lc%alpha*c%d ! d[ct]=alpha*d[ct], so full step vector s is in  d.

! Allocate space, if necessary and if possible.  Lc%All_Allocated is
! also reset in routine UpDate_H.

        IF (sh%limmemory .AND. .NOT. lc%all_allocated) THEN

          IF (sh%methodinuse==variable .OR. sh%methodinuse==available) THEN

! Cease allocation if there is insufficient memory for another
! term or if the maximum number of updates is reached.

            lc%memoryneeded = lc%memoryneeded + sh%ixleng

            IF ((sh%memory>0 .AND. lc%memoryneeded>sh%memory) .OR. &
                h%nterms==maxupdates) THEN
              enough =  1_short
              GO TO 60
            END IF
          END IF

          CALL my_allocate(sh%ixleng,allocnext,enough,us,sh,tr,h)
          IF (enough==0 .AND. lc%unknownm) m = m + 1 ! Increase m

60        IF (enough/=0) THEN
            lc%all_allocated = true
            m = h%nterms
            lc%unknownm = false
          END IF
        END IF

        IF (sh%limmemory .AND. us%updateform==productform .AND. m==0) THEN
          finalstate = invalidm ! m = 0 invalid using ProductForm
          GO TO 70
        END IF

!Check if a restart is to be forced
        forcerestart = sh%limmemory .AND. sh%doforcerestart .AND. &
          (us%updateform==sumform .OR. us%updateform==factoredform) .AND. &
          (sh%ct>sh%restartpt .OR. sh%methodinuse==sd)

        IF (tr%update) WRITE (tr%u,90090) pre // &
          ' ForceRestart    ct    RestartPt    SteepDStep', forcerestart, &
          sh%ct, sh%restartpt, sh%steepdstep

        IF (sh%limmemory) THEN
          SELECT CASE (us%updateform)

          CASE (sumform,factoredform) ! Determine part of algorithm we are in
! for next step.

            sh%qnpart = (forcerestart .AND. m/=0) .OR. &
              (sh%qnpart .AND. (lc%unknownm .OR. sh%ct<=m))

            dorestarttest = .NOT. sh%qnpart .AND. sh%ct > m + 1

          CASE (productform)
            sh%qnpart = true

          END SELECT
        END IF

        IF (forcerestart) THEN

          sh%dorestart = true

        ELSE IF (sh%limmemory .AND. (us%updateform==sumform .OR. us%updateform &
            ==factoredform) .AND. dorestarttest .AND. us%htest/=norestart) &
            THEN

          IF (tr%update) WRITE (tr%u,90020) pre // ' LimMemory part: restart?'

! Must be in   LimMemory sequence, so must check if
! restart is needed according to Powell criterion.  Can apply
! in metric defined by  H  or by  I; i.e.  using  g'*H*g, or
! g'*g.  Compute values for restart test.

          IF (us%htest==useh) THEN

! Powell's test with H as currently defined.
! Use xx as temporary storage for H*g.

            CALL sum_multiply(h,g,c%xx,us,sh,tr)

            tp1 = c%xx .gip. c%gg
            tp2 = c%xx .gip. g
          ELSE
            tp1 = g .gip. c%gg ! Ordinary test: Powell's test with H = I .
            tp2 = lc%nrmg**2
          END IF

          IF (tr%update) WRITE (tr%u,90100) pre // ' restart if tp1(', tp1, &
            ') > rho*tp2 (', us%rho*tp2, ')'

! Set restart flag if tau[ct] > rho; note that tau = tp1/tp2
! but the test is done without the divide.

          sh%dorestart = ABS(tp1) > ABS(us%rho*tp2)

          IF (sh%dorestart) sh%forcect = sh%forcect + 1

        ELSE

          IF (tr%update) WRITE (tr%u,90020) pre // ' no restart test.'
          sh%dorestart = false

        END IF

!>>>>>>>>>> PHASE  X: Update for Next Step.<<<<<<<<<<<<<<<<<<<<<<<<<

!   We now call a routine to update H from its value at the last point to
!   its value at the point which we have just reached at the end of this
!   line search. The details of the updating are in Update_H.

        c%gg = g - c%gg
        sh%eta = c%d .gip. c%gg

        IF (sh%limmemory .AND. us%updateform==factoredform) c%xx = h%u
        CALL update_h(c%d,c%gg,c%xx,h,us,sh,lc,tr)

!>>>>>>>>>> PHASE   XI:  Compute New Direction.<<<<<<<<<<<<<<<<<<<<<

COMPUTE_D: IF (sh%limmemory) THEN

! First get the negative of the new search direction into d.

          SELECT CASE (us%updateform)

          CASE (sumform)

            IF (sh%qnpart .OR. sh%dorestart) THEN
              CALL sum_multiply(h,g,c%d,us,sh,tr)
            ELSE
! Note that u = H*y is in xx from update.
! Compute H^*g by putting  H*g into gg and then working in
! the new update term. This must be done separately since
! the new update may not have been saved. Note that this is
! not an initial step, so we only do the gamma scaling if
! ScaleGamma=GammaAll.

              stg = c%d .ip. g
              utg = c%xx .ip. g
              nu = c%gg .ip. c%xx

              CALL sum_multiply(h,g,c%gg,us,sh,tr)

              IF (us%scalegamma==gammaall) THEN
                sigma = (two*stg/sh%eta) - (utg/nu)
                mu = -stg/nu
                gamma = sh%eta/nu
              ELSE
                sigma = ((one+nu/sh%eta)*stg-utg)/sh%eta
                mu = -stg/sh%eta
                gamma = one
              END IF
! Compute  H^*g into d.
              c%d = mu*c%xx + sigma*c%d + gamma*c%gg
            END IF

          CASE (productform)
            CALL product_multiply(h,g,c%d,us,sh,tr)
          CASE (factoredform)
            CALL factored_multiply(h,g,c%d,h%u,us,sh,tr)
          END SELECT

          IF (tr%update) WRITE (tr%u,90020) pre // ' new d using LimMemory.'

        ELSE COMPUTE_D ! QN case:  Calculate search direction d(ct+1) = -H^*g
!           H^ is in H.

          CALL qnewton_multiply(h,g,c%d,tr)

        END IF COMPUTE_D

        c%d = -c%d
        lc%nrmd = .gnorm. c%d
        lc%dg0 = c%d .gip. g

!   Test for a downhill direction.

        IF (lc%dg0>=zero) THEN
          IF (anytr()) THEN
            WRITE (tr%u,90020) pre // &
              ' ***failing*** on nondownhill direction'
            WRITE (tr%u,90040) pre // ' ***dg0->', lc%dg0, '***'
          END IF
          finalstate = nodescent
        ELSE
          sh%steepdstep = m == 0 .AND. sh%dorestart
        END IF

        sh%iterationtime = sh%iterationtime + cpusecs()

! If Dynamic method is used, check if switch to QN method should be done.
! C%EvalCts%time is the function evaluation time.  Sh%IterationTime is the
! time for the iteration (it does not include the print or function times).

        IF (sh%methodinuse==dynamic .AND. m>=2 .AND. sh%it==2 .AND. &
            c%evalcts%time*two/(lc%n-1)>sh%iterationtime) THEN

          IF (tr%flow) WRITE (tr%u,90020) pre // &
            ' Attempting to change method to QN'

          qnstorage = (lc%n*(lc%n+1))/2 + lc%n
          IF (sh%memory>0) THEN
            lots = sh%memory >= qnstorage
          ELSE
            lots = true
          END IF
          IF (lots) THEN
            CALL my_allocate(qnstorage-lc%n,allochi,enough,us,sh,tr,h)

! Temp has already been allocated if doing ProductForm updates.
            IF (us%updateform/=productform) CALL my_allocate(lc%n,alloctemp, &
              enough,us,sh,tr,h)
          END IF
          IF (lots .AND. enough==0) THEN
            sh%methodinuse = qn
            sh%limmemory = false
            us%dointerpolation = never
            us%startalpha = alwaysqn
            us%quadinterpolation = false
            us%ignoreinterval = false
            m = 0

            CALL dynamic_update_h(c%d,c%gg,c%xx,h,us,sh,tr)

            DEALLOCATE (h%firstterm%next%data) ! Deallocate the 2 Limited-
            DEALLOCATE (h%firstterm%next) ! Memory terms.
            DEALLOCATE (h%firstterm%data)
            DEALLOCATE (h%firstterm)
            IF (us%updateform==factoredform) DEALLOCATE (h%u)
            IF (us%seth0==diagonal) DEALLOCATE (h%h0)

            IF (tr%flow) WRITE (tr%u,90020) pre // ' Method changed to QN'
          ELSE
            IF (tr%flow) WRITE (tr%u,90020) pre // ' Method not changed to QN'
          END IF

        END IF

        IF (finalstate/=ok) THEN
          GO TO 70
        ELSE
          lc%cold = false
          IF (us%checkpoint/=0) CALL check_point(x,fx,g,state,c,sh%it, &
            file=checkfile,outunit=checkunit)
          GO TO 20
        END IF

! EXIT

70      IF ( .NOT. lc%noprint) THEN
          sh%it = sh%it + 1
          CALL print_iterate(lc%n,sh%it,fx,lc%nrmg,c%evalcts,c%print,c%prvals, &
            x,g,force=true,level=tr%level)
        END IF

        IF (anytr()) CALL error_codes(finalstate,sh,tr)

        state = finalstate
! Deallocate memory if
! not doing reverse
! communication.
! No memory allocated if input bad
        IF (state/=revcommf .AND. state/=revcommg .AND. &
          state/=revcommfandg .AND. state/=badinput) CALL my_deallocate(sh,tr, &
          h)

        CALL leaving(id,tr%flow,tr%u,tr%level)

        RETURN

! FORMATS:

90000   FORMAT ( &
          ' Trace flags: unit  level  input   flow   step  lsalpha lsreal'/ &
          13X,I3,4X,I3,5L7/ &
          '             lsflow update values vectors XandD cubic'/9X,6L7)
90010   FORMAT (A,G15.7)
90020   FORMAT (A)
90030   FORMAT (A/I8,4L9)
90040   FORMAT (A,G15.7,A)
90050   FORMAT (A,2G15.7,A)
90060   FORMAT (A/11X,4(G12.6,1X))
90070   FORMAT (A,G15.7,A,G15.7)
90080   FORMAT (A,L1)
90090   FORMAT (A/4X,L9,2I11,2L12)
90100   FORMAT (A,G15.7,A,G15.7,A)

90110   FORMAT (' **** Minimize_f entered and initialization complete ****'// &
          ' dimension  = ',I5,T40,' memory available is ', &
          I7/' accuracy requested = ',G15.7)
90120   FORMAT (' Expected reduction in f equals initial function', &
          ' value of ',G15.7)
90130   FORMAT (' Expected reduction in f is unknown(DecreaseInF=',G8.1,')')
90140   FORMAT (' Expected reduction in f is ',G15.7)
90150   FORMAT (' Real control values  rho = ',G15.7,' beta = ', &
          G15.7/ &
          ' Logical control values   QuadInterpolation  CountFromRestart ', &
          ' IgnoreInterval'/17X,3L18/ &
          '                             RelativeToF0      RelativeToG0   '/ &
          17X,2L18//' The following has been set during initialization: ', &
          ' LimMemory (',L1,')'/' machine relative accuracy eps = ', &
          E8.2/' termination relative to ',G14.7,'(f); ',G14.7,'(g)')
90160   FORMAT (' storage of ',I6,' sufficient; using qn algorithm.'/)
90170   FORMAT (' storage limited; using ',I3,' updates.'/)

!----------!

      CONTAINS
!----------!

        REAL (stnd) FUNCTION angle(x)

! Return angle, in degrees, whose cos is x.

          IMPLICIT NONE
!-------------!

! ARGUMENTS:

          REAL (stnd), INTENT (IN) :: x

! EXECUTION:

          angle = c180/ACOS(-one)*ACOS(x)

        END FUNCTION


        LOGICAL FUNCTION anytr()

! returns true if any trace flag is set

          anytr = tr%input .OR. tr%flow .OR. tr%steptypes .OR. tr%lsalpha .OR. &
            tr%lsreal .OR. tr%lsflow .OR. tr%update .OR. tr%values .OR. &
            tr%vectors .OR. tr%xandd .OR. tr%cubic

        END FUNCTION


        SUBROUTINE checkdata

! Check that all input values are correct. Set a sequence of
! up to NumErrs error flags, one for each error found.  These
! are set in the array Sh%errors, with Sh%errors(0) containing
! the number of errors found.

          i = 0

          IF (lc%n<=1) THEN
            i = MIN(i+1,numerrs)
            sh%error(i) = badn
          END IF

          IF (accuracy<=zero) THEN
            i = MIN(i+1,numerrs)
            sh%error(i) = badacc
          END IF

          IF (entrystate/=normal .AND. entrystate/=normalwithfg .AND. &
              entrystate/=revcommstart) THEN
            i = MIN(i+1,numerrs)
            sh%error(i) = badstate
          END IF

          IF (memory<0) THEN
            i = MIN(i+1,numerrs)
            sh%error(i) = badmemory
          END IF

          IF (PRESENT(evallimit)) THEN
            IF (evallimit<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badevallim
            END IF
          END IF

          IF (PRESENT(frequency)) THEN
            IF (frequency<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badfreq
            END IF
          END IF

          IF (PRESENT(checkpoint)) THEN
            IF (checkpoint<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badcheckpoint
            END IF
          END IF

          IF (PRESENT(checkunit)) THEN
            IF (checkunit<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badcheckunit
            END IF
          END IF

          IF (PRESENT(method)) THEN
            SELECT CASE (method)
            CASE (sd,cg,conmin,fixterms,variable,available,dynamic,qn)
              IF ((method==sd .OR. method==cg) .AND. PRESENT(updateform)) THEN
                IF (updateform==productform) THEN
                  i = MIN(i+1,numerrs)
                  sh%error(i) = badcombination
                END IF
              END IF
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badmethod
            END SELECT
          END IF

          IF (PRESENT(terms)) THEN
            IF (terms<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badterms
            END IF
          END IF

          IF (PRESENT(derivatives)) THEN
            SELECT CASE (derivatives)
            CASE (analytic,comparetest,firsttest,differences)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badderiv
            END SELECT
          END IF

          IF (PRESENT(systemmemory)) THEN
            IF (systemmemory<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badsysmemory
            END IF
          END IF

          IF (PRESENT(scalef)) THEN
            IF (scalef<=0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badscalef
            END IF
          END IF

          IF (PRESENT(expense)) THEN
            IF (expense<1) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badexpense
            END IF
          END IF

          IF (PRESENT(evaltraceunit)) THEN
            IF (evaltraceunit<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badevalunit
            END IF
          END IF

          IF (PRESENT(printunit)) THEN
            IF (printunit<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badprintunit
            END IF
          END IF

          IF (PRESENT(thenorm)) THEN
            SELECT CASE (thenorm)
            CASE (l1,l2,linf,g2)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badnorm
            END SELECT
          END IF

          IF (PRESENT(termtraceunit)) THEN
            IF (termtraceunit<0) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badtermunit
            END IF
          END IF

          IF (PRESENT(dointerpolation)) THEN
            SELECT CASE (dointerpolation)
            CASE (every,on1,on2,on3,never)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badinterp
            END SELECT
          END IF

          IF (PRESENT(startalpha)) THEN
            SELECT CASE (startalpha)
            CASE (neverqn,before1qn,before2qn,before3qn,alwaysqn)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badalpha
            END SELECT
          END IF

          IF (PRESENT(scalegamma)) THEN
            SELECT CASE (scalegamma)
            CASE (nogammascale,gammafirst,gammaall)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badgamma
            END SELECT
          END IF

          IF (PRESENT(htest)) THEN
            SELECT CASE (htest)
            CASE (norestart,usei,useh)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badhtest
            END SELECT
          END IF

          IF (PRESENT(updateform)) THEN
            SELECT CASE (updateform)
!Case(SumForm, ProductForm, FactoredForm)
            CASE (sumform,productform)
! OK
            CASE (factoredform) ! FactoredForm is not yet implemented
              i = MIN(i+1,numerrs)
              sh%error(i) = badupdate
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badupdate
            END SELECT
          END IF

          IF (PRESENT(startstep)) THEN
            SELECT CASE (startstep)
            CASE (fletcher,powell)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badstep
            END SELECT
          END IF

          IF (PRESENT(rho)) THEN
            IF (rho<=zero) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badrho
            END IF
          END IF

          IF (PRESENT(beta)) THEN
            IF (beta<zero .OR. beta>one) THEN
              i = MIN(i+1,numerrs)
              sh%error(i) = badbeta
            END IF
          END IF

          IF (PRESENT(seth0)) THEN
            SELECT CASE (seth0)
            CASE (ident,diagonal,computed1,computed2,computed3)
! OK
            CASE DEFAULT
              i = MIN(i+1,numerrs)
              sh%error(i) = badseth0
            END SELECT
          END IF

          sh%error(0) = i

          RETURN

        END SUBROUTINE


        SUBROUTINE initialize_min

!  This routine is separated out from the main body of the code just
!  as a convenience. Its purpose is to initialize many of the quantities
!  used by the code.  In particular, it sets initial values first to
!  default values, and then resets any which for which optional
!  arguments are present.

          IMPLICIT NONE
!-------------!

          us = defaultuserdefined

          ev = defaultevalstate

          pr = defaultprintstate
          c%prvals = initprvals

          tm = defaulttermstate

          ALLOCATE (c%d(1:lc%n),c%xx(1:lc%n),c%gg(1:lc%n))

          NULLIFY (c%h%firstterm)
          NULLIFY (c%h%lastterm)
          NULLIFY (c%h%nextterm)
          NULLIFY (c%h%hi)
          NULLIFY (c%h%h0)
          NULLIFY (c%h%temp) ! used in routines Update_H and Product_Multiply
          NULLIFY (c%h%u) ! used in routines Minimize_f and Factored_Multiply

          sh%accuracy = accuracy
          sh%memory = memory

          IF (PRESENT(evallimit)) ev%evallimit = evallimit
          IF (PRESENT(frequency)) pr%frequency = frequency
          IF (PRESENT(checkpoint)) us%checkpoint = checkpoint
          IF (PRESENT(method)) us%method = method
          IF (PRESENT(terms)) THEN
            h%mterms = terms
          ELSE
            h%mterms = 0
          END IF
          IF (PRESENT(derivatives)) ev%derivatives = derivatives
          IF (PRESENT(decreaseinf)) us%decreaseinf = decreaseinf
          IF (PRESENT(systemmemory)) us%systemmemory = systemmemory
          IF (PRESENT(scalef)) ev%scalef = scalef
          IF (PRESENT(expense)) ev%expense = expense
          IF (PRESENT(tracef)) ev%tracef = tracef
          IF (PRESENT(traceg)) ev%traceg = traceg
          IF (PRESENT(tracedervtest)) ev%tracedervtest = tracedervtest
          IF (PRESENT(evaltraceunit)) ev%evaltraceunit = evaltraceunit
          IF (PRESENT(printunit)) pr%printunit = printunit
          IF (PRESENT(printgrad)) pr%printgrad = printgrad
          IF (PRESENT(printx)) pr%printx = printx
          IF (PRESENT(usegrad)) tm%usegrad = usegrad
          IF (PRESENT(usestep)) tm%usestep = usestep
          IF (PRESENT(useshanno)) tm%useshanno = useshanno
          IF (PRESENT(usefunc)) tm%usefunc = usefunc
          IF (PRESENT(thenorm)) tm%thenorm = thenorm
          IF (PRESENT(traceterm)) tm%traceterm = traceterm
          IF (PRESENT(termtraceunit)) tm%termtraceunit = termtraceunit
          IF (PRESENT(relativetof0)) us%relativetof0 = relativetof0
          IF (PRESENT(relativetog0)) us%relativetog0 = relativetog0
          IF (PRESENT(dointerpolation)) us%dointerpolation = dointerpolation
          IF (PRESENT(quadinterpolation)) us%quadinterpolation = &
            quadinterpolation
          IF (PRESENT(startalpha)) us%startalpha = startalpha
          IF (PRESENT(scalegamma)) us%scalegamma = scalegamma
          IF (PRESENT(htest)) us%htest = htest
          IF (PRESENT(updateform)) us%updateform = updateform
          IF (PRESENT(startstep)) us%startstep = startstep
          IF (PRESENT(ignoreinterval)) us%ignoreinterval = ignoreinterval
          IF (PRESENT(countfromrestart)) us%countfromrestart = &
            countfromrestart
          IF (PRESENT(rho)) us%rho = rho
          IF (PRESENT(beta)) us%beta = beta
          IF (PRESENT(seth0)) us%seth0 = seth0
          IF (PRESENT(parh0)) us%parh0 = parh0
          IF (PRESENT(scalecolumns)) us%scalecolumns = scalecolumns

          lc%checknext = us%checkpoint

        END SUBROUTINE


        SUBROUTINE doevalf

          which = both
          CALL evaluate_f(f,x,fx,g,which,ev,c%evalcts,first=lc%firsteval, &
            test=c%evalers,level=tr%level)
          lc%noprint = true
          finalstate = ok
          sh%validf = false
          SELECT CASE (which)
          CASE (limit)
            finalstate = excessfevals
            lc%noprint = false
          CASE (abort)
            finalstate = abortmin
          CASE DEFAULT
            lc%noprint = false
            sh%validf = true
          END SELECT
        END SUBROUTINE

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> minimize.f90
    MODULE minimize ! The Variable Storage Minimization Package
      USE precisions, ONLY : stnd, extd, short, long
      USE min_codes, ONLY : normal, done
      USE minimizef, ONLY : minimize_f
      USE systemstate, ONLY : minimizestate

! For detailed documentation, see the file min.doc and
! the references contained therein.





      IMPLICIT NONE
      PRIVATE
!------              -------------!

      PUBLIC :: minimizestate, minimize_f, normal, done, stnd, long, short, &
        extd

    END

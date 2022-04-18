


    FUNCTION d1mach(i)
!***PURPOSE  RETURNS DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS
!***DESCRIPTION

!     D1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBPROGRAM WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE

!          D = D1MACH(I)

!     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF D ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.

!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  D1MACH

!***FIRST EXECUTABLE STATEMENT  D1MACH
! .. Function Return Value ..
      REAL (kind(0.0D0)) :: d1mach
! ..
! .. Parameters ..
      REAL (kind(0.0D0)), PARAMETER :: base = 2.0D0, zero = 0.0D0
! ..
! .. Scalar Arguments ..
      INTEGER :: i
! ..
! .. Intrinsic Functions ..
      INTRINSIC epsilon, huge, kind, log10, tiny
! ..
      SELECT CASE (i)
      CASE (1)
        d1mach = tiny(zero)
      CASE (2)
        d1mach = huge(zero)
      CASE (3)
        d1mach = epsilon(zero)
      CASE (4)
        d1mach = base*epsilon(zero)
      CASE (5)
        d1mach = log10(base)
      END SELECT

    END FUNCTION d1mach


    FUNCTION i1mach(i)
!***PURPOSE  RETURN INTEGER MACHINE DEPENDENT CONSTANTS.
!***DESCRIPTION

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   THESE MACHINE CONSTANT ROUTINES MUST BE ACTIVATED FOR
!   A PARTICULAR ENVIRONMENT.
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     I1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE

!          K = I1MACH(I)

!     WHERE I=1,...,16.  THE (OUTPUT) VALUE OF K ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.

!  I/O UNIT NUMBERS.
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.

!  WORDS.
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.

!  INTEGERS.
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM

!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!    I1MACH( 7) = A, THE BASE.
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.

!  FLOATING-POINT NUMBERS.
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )

!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, THE BASE.

!  SINGLE-PRECISION
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
! .. Function Return Value ..
      INTEGER :: i1mach
! ..
! .. Scalar Arguments ..
      INTEGER :: i
! ..
! .. Local Scalars ..
      REAL (KIND(0.0d0)) :: d
      REAL :: r
! ..
! .. Intrinsic Functions ..
      INTRINSIC digits, huge, maxexponent, minexponent, radix
! ..
      SELECT CASE (i)
      CASE (1)
! Standard input                                   
        i1mach = 5
      CASE (2)
! Standard output                                  
        i1mach = 6
      CASE (3)
! Standard punch :-)                               
        i1mach = 6
      CASE (4)
! Standard error                                   
        i1mach = 6
      CASE (5)
!  Number of bits /integer (+1 for the s
        i1mach = digits(i) + 1
      CASE (6)
! Number of characters / integer :-)               
        i1mach = 4
      CASE (7)
        i1mach = radix(i) ! base of integers                           
      CASE (8)
        i1mach = digits(i) ! number of base radix digits in integer    
      CASE (9)
        i1mach = huge(i) ! Maximum integer                             
      CASE (10)
        i1mach = radix(r) ! base of floating point                     
      CASE (11)
        i1mach = digits(r) ! number of base radix digits in sp         
      CASE (12)
        i1mach = minexponent(r) ! minimun sp exponent                  
      CASE (13)
        i1mach = maxexponent(r) ! maximum sp exponent                  
      CASE (14)
        i1mach = digits(d) ! number of base radix digits in dp         
      CASE (15)
        i1mach = minexponent(d) ! minimun dp exponent                  
      CASE (16)
        i1mach = maxexponent(d) ! maximum dp exponent                  
      END SELECT
    END FUNCTION i1mach

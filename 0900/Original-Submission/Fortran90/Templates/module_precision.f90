!*********************************************************************************************************************************
!*
!* MODULE: PRECISION
!*
!* PURPOSE: sets precision in all variables
!*
!* DEPENDENCIES:
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables
!*          - HIGH (integer): number of bytes for double precision variables
!*          - PLO (integer): precision for logical variables
!*          - PCH (integer): precision for character variables
!*          - PIN (integer): precision for integer variables
!*          - PRE (integer): precision for real variables
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module precision
  !*
  implicit none
  !*
  logical , private :: logicalvar !* a logical variable
  character , private :: charactervar !* a character variable
  integer , private :: integervar !* an integer variable
  !*
  integer , parameter , public :: low = kind( 0.0e0 ) !* a real single precision variable
  integer , parameter , public :: high = kind( 0.0d0 ) !* a real double precision variable
  !*
  integer , parameter , public :: plo = kind( logicalvar ) !* kind for logical variables
  integer , parameter , public :: pch = kind( charactervar ) !* kind for character variables
  integer , parameter , public :: pin = kind( integervar ) !* kind for integer variables
  integer , parameter , public :: pre = high !* kind for real variables
  !*
end module precision

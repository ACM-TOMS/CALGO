!*********************************************************************************************************************************
!*
!* MODULE: PARAMETERS
!*
!* PURPOSE: set parameters regarding output files
!*
!* DEPENDENCIES:
!*               - PRECISION
!*
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*
!*          - UNITTRUTH (integer): file unit for truth
!*          - UNITMODEL (integer): file unit for model
!*          - UNITKF (integer): file unit for ekf
!*          - UNITRRSQRTKF (integer): file unit for rrsqrtkf
!*          - UNITENKF (integer): file unit for enkf
!*          - UNITRRSQRTENKF (integer): file unit for rrsqrtenkf
!*          - UNITOBSER (integer): file unit for observations
!*          - FILETRUTH (character*100): file name for truth
!*          - FILEMODEL (character*100): file name for model
!*          - FILEEKF (character*100): file name for ekf
!*          - FILERRSQRTKF (character*100): file name for rrsqrtkf
!*          - FILEENKF (character*100): file name for enkf
!*          - FILERRSQRTENKF (character*100): file name for rrsqrtenkf
!*          - FILEOBSER (character*100): file name for observations
!*
!* SCOPE: this module belongs to the example
!*
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
!*
module parameters
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) :: unittruth
  integer( kind = pin ) :: unitmodel
  integer( kind = pin ) :: unitekf
  integer( kind = pin ) :: unitrrsqrtkf
  integer( kind = pin ) :: unitenkf
  integer( kind = pin ) :: unitrrsqrtenkf
  integer( kind = pin ) :: unitobser
  character( kind = pch , len = 100_pin ) :: filetruth
  character( kind = pch , len = 100_pin ) :: filemodel
  character( kind = pch , len = 100_pin ) :: fileekf
  character( kind = pch , len = 100_pin ) :: filerrsqrtkf
  character( kind = pch , len = 100_pin ) :: fileenkf
  character( kind = pch , len = 100_pin ) :: filerrsqrtenkf
  character( kind = pch , len = 100_pin ) :: fileobser
  !*
contains
  !*
  !* Set parameters 
  subroutine parameters_parametersup()
    implicit none
    unittruth = 10_pin
    unitmodel = 20_pin
    unitekf = 30_pin
    unitrrsqrtkf = 40_pin
    unitenkf = 50_pin
    unitrrsqrtenkf = 60_pin
    unitobser = 70_pin
    filetruth = 'truth.txt'
    filemodel = 'model.txt'
    fileekf = 'ekf.txt'
    filerrsqrtkf = 'rrsqrtkf.txt'
    fileenkf = 'enkf.txt'
    filerrsqrtenkf = 'rrsqrtenkf.txt'
    fileobser = 'obser.txt'
  end subroutine parameters_parametersup
  !*
end module parameters

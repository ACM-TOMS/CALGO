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
!*          - UNITASSIMEKF (integer): file unit for ekf  
!*          - UNITASSIMRRSQRTKF (integer): file unit for rrsqrtkf  
!*          - UNITASSIMENKF (integer): file unit for enkf  
!*          - UNITASSIMRRSQRTENKF (integer): file unit for rrsqrtenkf  
!*          - UNITOBSER (integer): file unit for observations 
!*          - FILETRUTH (character*100): file name for truth       
!*          - FILEMODEL (character*100): file name for model  
!*          - FILEASSIMEKF (character*100): file name for ekf
!*          - FILEASSIMRRSQRTKF (character*100): file name for rrsqrtkf
!*          - FILEASSIMENKF (character*100): file name for enkf  
!*          - FILEASSIMRRSQRTENKF (character*100): file name for rrsqrtenkf  
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
module parameters
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: unittruth !* file unit for truth    
  integer( kind = pin ) , public :: unitmodel !* file unit for model  
  integer( kind = pin ) , public :: unitassimekf !* file unit for ekf  
  integer( kind = pin ) , public :: unitassimrrsqrtkf !* file unit for rrsqrtkf  
  integer( kind = pin ) , public :: unitassimenkf !* file unit for enkf  
  integer( kind = pin ) , public :: unitassimrrsqrtenkf !* file unit for rrsqrtenkf  
  integer( kind = pin ) , public :: unitobser !* file unit for observations 
  character( kind = pch , len = 100_pin ) , public :: filetruth !* file name for truth       
  character( kind = pch , len = 100_pin ) , public :: filemodel !* file name for model  
  character( kind = pch , len = 100_pin ) , public :: fileassimekf !* file name for ekf
  character( kind = pch , len = 100_pin ) , public :: fileassimrrsqrtkf !* file name for rrsqrtkf
  character( kind = pch , len = 100_pin ) , public :: fileassimenkf !* file name for enkf  
  character( kind = pch , len = 100_pin ) , public :: fileassimrrsqrtenkf !* file name for rrsqrtenkf  
  character( kind = pch , len = 100_pin ) , public :: fileobser !* file name for observations  
  !*
contains
  !*
  !* setting parameters
  subroutine parameters_parametersup()
    implicit none
    unittruth = 10_pin  
    unitmodel = 20_pin
    unitassimekf = 30_pin
    unitassimrrsqrtkf = 40_pin
    unitassimenkf = 50_pin
    unitassimrrsqrtenkf = 60_pin
    unitobser = 70_pin
    filetruth = 'truth.txt'
    filemodel = 'model.txt'
    fileassimekf = 'ekf.txt'
    fileassimrrsqrtkf = 'rrsqrtkf.txt'
    fileassimenkf = 'enkf.txt'
    fileassimrrsqrtenkf = 'rrsqrtenkf.txt'
    fileobser = 'obser.txt'
  end subroutine parameters_parametersup
  !*
end module parameters

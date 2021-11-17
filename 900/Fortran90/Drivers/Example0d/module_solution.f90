!*********************************************************************************************************************************
!*
!* MODULE: SOLUTION
!*
!* PURPOSE: sets the exact solution (useful to compute errors)
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
!* SCOPE: This module belongs to the example
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
module solution
  !*
  use precision
  !*
  implicit none
  !*
contains
  !*
  !* exact solution - useful to compare results
  function solution_solution( t )
    implicit none
    real( kind = pre ) , intent( in ) :: t
    real( kind = pre ) :: solution_solution
    solution_solution = 2.0_pre + sin( t )
  end function solution_solution
  !*
end module solution

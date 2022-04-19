!*********************************************************************************************************************************
!*
!* MODULE: TOOLS
!*
!* PURPOSE: provides several tools
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
!* SCOPE: this module belongs to the library
!* 
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module tools
  !*
  use precision
  !*
  implicit none
  !*
contains
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: TOOLS_IJ2S
  !*
  !* PURPOSE: maps indices of a matrix of size NX x NY to indices of a vector of size NX * NY using fortran style
  !*
  !* INPUTS:
  !*         - NX (integer): first dimension of a matrix
  !*         - I (integer): first coordinate of a matrix element ( 1 <= I <= NX )
  !*         - J (integer): second coordinate of a matrix element ( 1 <= J <= NY )
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*         - S (integer): index of the vector element ( 1 <= S <= NX * NY ) corresponding to the ( I , J ) matrix element accor&
  !*                        &ding to the fortran style
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine tools_ij2s( nx , i , j , s )
    implicit none
    integer( kind = pin ) , intent( in ) :: nx
    integer( kind = pin ) , intent( in ) :: i
    integer( kind = pin ) , intent( in ) :: j
    integer( kind = pin ) , intent( out ) :: s
    s = ( j - 1_pin ) * nx + i
  end subroutine tools_ij2s
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: TOOLS_S2IJ
  !*
  !* PURPOSE: maps indices of a vector of size NX * NY to indices of a matrix of size NX x NY using fortran style
  !*
  !* INPUTS:
  !*         - NX (integer): first dimension of a matrix
  !*         - S (integer): index of the vector element ( 1 <= S <= NX * NY )
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - I (integer): first coordinate of the matrix element corresponding to the S vector element according to the fortra&
  !*                         &n style 
  !*          - J (integer): second coordinate of the matrix element corresponding to the S vector element according to the fortr&
  !*                         &an style 
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine tools_s2ij( nx , s , i , j )
    implicit none
    integer( kind = pin ) , intent( in ) :: nx
    integer( kind = pin ) , intent( in ) :: s
    integer( kind = pin ) , intent( out ) :: i
    integer( kind = pin ) , intent( out ) :: j
    i = mod( s , nx )
    if ( i == 0_pin ) i = nx
    j = ( s - i ) / nx + 1_pin
  end subroutine tools_s2ij
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: TOOLS_IJK2S
  !*
  !* PURPOSE: maps indices of a matrix of size NX x NY x NZ to indices of a vector of size NX * NY * NZ using fortran style
  !*
  !* INPUTS:
  !*         - NX (integer): first dimension of a matrix
  !*         - NY (integer): second dimension of a matrix
  !*         - I (integer): first coordinate of a matrix element ( 1 <= I <= NX )
  !*         - J (integer): second coordinate of a matrix element ( 1 <= J <= NY )
  !*         - K (integer): third coordinate of a matrix element ( 1 <= K <= NZ )
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS: 
  !*          - S (integer): index of the vector element ( 1 <= S <= NX * NY * NZ ) corresponding to the ( I , J , K ) matrix ele&
  !*                         &ment according to the fortran style
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine tools_ijk2s( nx , ny , i , j , k , s )
    implicit none
    integer( kind = pin ) , intent( in ) :: nx
    integer( kind = pin ) , intent( in ) :: ny
    integer( kind = pin ) , intent( in ) :: i
    integer( kind = pin ) , intent( in ) :: j
    integer( kind = pin ) , intent( in ) :: k
    integer( kind = pin ) , intent( out ) :: s
    s = ( k - 1 ) * ny * nx + ( j - 1_pin ) * nx + i
  end subroutine tools_ijk2s
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: TOOLS_S2IJK
  !*
  !* PURPOSE: maps indices of a vector of size NX * NY * NZ to indices of a matrix of size NX x NY x NZ using fortran style
  !*
  !* INPUTS:
  !*         - NX (integer): first dimension of a matrix
  !*         - NY (integer): second dimension of a matrix
  !*         - S (integer): index of a vector component ( 1 <= S <= NX * NY * NZ )
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - I (integer): first coordinate of the matrix element corresponding to the S vector element according to the fortra&
  !*                         &n style 
  !*          - J (integer): second coordinate of the matrix element corresponding to the S vector element according to the fortr&
  !*                         &an style 
  !*          - K (integer): third coordinate of the matrix element corresponding to the S vector element according to the fortra&
  !*                         &n style 
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine tools_s2ijk( nx , ny , s , i , j , k )
    implicit none
    integer( kind = pin ) , intent( in ) :: nx
    integer( kind = pin ) , intent( in ) :: ny
    integer( kind = pin ) , intent( in ) :: s
    integer( kind = pin ) , intent( out ) :: i
    integer( kind = pin ) , intent( out ) :: j
    integer( kind = pin ) , intent( out ) :: k
    integer( kind = pin ) :: nxny
    integer( kind = pin ) :: ij
    nxny = nx * ny
    ij = mod( s , nxny )
    if ( ij == 0_pin ) ij = nxny
    k = ( s - ij ) / nxny + 1_pin
    i = mod( ij , nx )
    if ( i == 0_pin ) i = nx
    j = ( ij - i ) / nx + 1_pin
  end subroutine tools_s2ijk
  !*
end module tools

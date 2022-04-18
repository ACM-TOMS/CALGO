module mod_INSERTING
! **********************************************************************
!     Author : M.YOUAN 
!
!     Date of last modification : 24.4.02
!      
!     Description :this module is based one two  chained list ( one for
!     collection of matrix and a another for elements of each matrix) .
!     Subroutines are used to create,accede to,delete components of these 
!     lists  
! **********************************************************************

use blas_sparse_namedconstants
use properties
implicit none
interface access_element
  module procedure iaccess_element
  module procedure saccess_element
  module procedure daccess_element
  module procedure caccess_element
  module procedure zaccess_element
end interface

interface access_matrix
  module procedure iaccess_matrix
  module procedure saccess_matrix
  module procedure daccess_matrix
  module procedure caccess_matrix
  module procedure zaccess_matrix
end interface


!****************************************
type i_inpnt1
integer::row_ind,col_ind
integer::value
end type i_inpnt1

type i_inblock
integer ::row_block_ind,col_block_ind
integer,dimension(:,:),pointer::value
end type i_inblock

type i_invblock
   integer ::row_vblock_ind,col_vblock_ind
   integer,dimension(:,:),pointer::value
end type i_invblock

type i_inelement
type(i_inblock)::blin
type(i_inpnt1)::pntin
type(i_invblock)::vblin
end type i_inelement

type i_element
integer::number
type(i_inelement)::contents
type(i_element),pointer::pntr
end type i_element

type i_matrix
integer,dimension(6)::DIM
integer::property,number,new
character*11::format
integer,dimension(:),pointer::sub_rows,sub_cols,trb,tre
type(i_element),pointer::i_element_start
type(i_matrix),pointer::pntr
end type  i_matrix

!****************************************

type d_inpnt1
integer::row_ind,col_ind
real(kind=dp)::value
end type d_inpnt1

type d_inblock
integer ::row_block_ind,col_block_ind
real(kind=dp),dimension(:,:),pointer::value
end type d_inblock

type d_invblock
integer ::row_vblock_ind,col_vblock_ind
real(kind=dp),dimension(:,:),pointer::value
end type d_invblock

type d_inelement
type(d_inblock)::blin
type(d_inpnt1)::pntin
type(d_invblock)::vblin
end type d_inelement

type d_element
integer::number
type(d_inelement)::contents
type(d_element),pointer::pntr
end type d_element

type d_matrix
integer,dimension(6)::DIM
integer::property,number,new
character*11::format
integer,dimension(:),pointer::sub_rows,sub_cols,trb,tre
type(d_element),pointer::d_element_start
type(d_matrix),pointer::pntr
end type  d_matrix

!*****************************************

type s_inpnt1
integer::row_ind,col_ind
real(kind=sp)::value
end type s_inpnt1

type s_inblock
integer ::row_block_ind,col_block_ind
real(kind=sp),dimension(:,:),pointer::value
end type s_inblock

type s_invblock
integer ::row_vblock_ind,col_vblock_ind
real(kind=sp),dimension(:,:),pointer::value
end type s_invblock

type s_inelement
type(s_inblock)::blin
type(s_inpnt1)::pntin
type(s_invblock)::vblin
end type s_inelement

type s_element
integer::number
type(s_inelement)::contents
type(s_element),pointer::pntr
end type s_element

type s_matrix
integer,dimension(6)::DIM
integer::property,number,new
character*11::format
integer,dimension(:),pointer::sub_rows,sub_cols,trb,tre
type(s_element),pointer::s_element_start
type(s_matrix),pointer::pntr
end type  s_matrix

!****************************************

type c_inpnt1
integer::row_ind,col_ind
complex(kind=sp)::value
end type c_inpnt1

type c_inblock
integer ::row_block_ind,col_block_ind
complex(kind=sp),dimension(:,:),pointer::value
end type c_inblock

type c_invblock
integer ::row_vblock_ind,col_vblock_ind
complex(kind=sp),dimension(:,:),pointer::value
end type c_invblock

type c_inelement
type(c_inblock)::blin
type(c_inpnt1)::pntin
type(c_invblock)::vblin
end type c_inelement

type c_element
integer::number
type(c_inelement)::contents
type(c_element),pointer::pntr
end type c_element

type c_matrix
integer,dimension(6)::DIM
integer::property,number,new
character*11::format
integer,dimension(:),pointer::sub_rows,sub_cols,trb,tre
type(c_element),pointer::c_element_start
type(c_matrix),pointer::pntr
end type  c_matrix

!****************************************

type z_inpnt1
integer::row_ind,col_ind
complex(kind=dp)::value
end type z_inpnt1

type z_inblock
integer ::row_block_ind,col_block_ind
complex(kind=dp),dimension(:,:),pointer::value
end type z_inblock

type z_invblock
integer ::row_vblock_ind,col_vblock_ind
complex(kind=dp),dimension(:,:),pointer::value
end type z_invblock

type z_inelement
type(z_inblock)::blin
type(z_inpnt1)::pntin
type(z_invblock)::vblin
end type z_inelement

type z_element
integer::number
type(z_inelement)::contents
type(z_element),pointer::pntr
end type z_element

type z_matrix
integer,dimension(6)::DIM
integer::property,number,new
character*11::format
integer,dimension(:),pointer::sub_rows,sub_cols,trb,tre
type(z_element),pointer::z_element_start
type(z_matrix),pointer::pntr
end type  z_matrix

!*****************************************

type(i_matrix), pointer,SAVE :: i_matrix_start
type(d_matrix), pointer,SAVE :: d_matrix_start
type(s_matrix), pointer,SAVE :: s_matrix_start
type(c_matrix), pointer,SAVE :: c_matrix_start
type(z_matrix), pointer,SAVE :: z_matrix_start

logical, SAVE, PRIVATE :: iins_init = .FALSE.
logical, SAVE, PRIVATE :: dins_init = .FALSE.
logical, SAVE, PRIVATE :: sins_init = .FALSE.
logical, SAVE, PRIVATE :: cins_init = .FALSE.
logical, SAVE, PRIVATE :: zins_init = .FALSE.

contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "integerINSERTING_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "singlePrecisionINSERTING_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "doublePrecisionINSERTING_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "singleComplexINSERTING_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "doubleComplexINSERTING_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_INSERTING

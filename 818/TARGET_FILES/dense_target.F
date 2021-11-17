      module mod_dense_mat_algos
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : DENSE MATRIX ALGORITHMS FOR BLOCK SPARSE MATRICES
! **********************************************************************

      use properties

      implicit none

      interface block_mult_vec
        module procedure iblock_mult_vec
        module procedure sblock_mult_vec
        module procedure dblock_mult_vec
        module procedure cblock_mult_vec
        module procedure zblock_mult_vec
      end interface

      interface block_Z_mult_vec
        module procedure iblock_Z_mult_vec
        module procedure sblock_Z_mult_vec
        module procedure dblock_Z_mult_vec
        module procedure cblock_Z_mult_vec
        module procedure zblock_Z_mult_vec
      end interface

      interface block_T_mult_vec
        module procedure iblock_T_mult_vec
        module procedure sblock_T_mult_vec
        module procedure dblock_T_mult_vec
        module procedure cblock_T_mult_vec
        module procedure zblock_T_mult_vec
      end interface

      interface block_H_mult_vec
        module procedure iblock_H_mult_vec
        module procedure sblock_H_mult_vec
        module procedure dblock_H_mult_vec
        module procedure cblock_H_mult_vec
        module procedure zblock_H_mult_vec
      end interface

      interface invert_left_lower
        module procedure iinvert_left_lower
        module procedure sinvert_left_lower
        module procedure dinvert_left_lower
        module procedure cinvert_left_lower
        module procedure zinvert_left_lower
      end interface

      interface invert_T_left_lower
        module procedure iinvert_T_left_lower
        module procedure sinvert_T_left_lower
        module procedure dinvert_T_left_lower
        module procedure cinvert_T_left_lower
        module procedure zinvert_T_left_lower
      end interface

      interface invert_right_upper
        module procedure iinvert_right_upper
        module procedure sinvert_right_upper
        module procedure dinvert_right_upper
        module procedure cinvert_right_upper
        module procedure zinvert_right_upper
      end interface
 
      interface invert_T_right_upper
        module procedure iinvert_T_right_upper
        module procedure sinvert_T_right_upper
        module procedure dinvert_T_right_upper
        module procedure cinvert_T_right_upper
        module procedure zinvert_T_right_upper
      end interface


      contains
! **********************************************************************
!
! **********************************************************************
#include "integerdense_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singlePrecisiondense_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doublePrecisiondense_source.F"
! **********************************************************************
!
! **********************************************************************
#include "singleComplexdense_source.F"
! **********************************************************************
!
! **********************************************************************
#include "doubleComplexdense_source.F"
! **********************************************************************
!
! **********************************************************************
      end module mod_dense_mat_algos

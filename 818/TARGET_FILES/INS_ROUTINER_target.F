      module mod_INS_ROUTINER
      use mod_INSERTING
      use SparseBLAS1
      use properties

      interface INS_entry
      module procedure iINS_entry
      module procedure sINS_entry
      module procedure dINS_entry
      module procedure cINS_entry
      module procedure zINS_entry
      end interface
      
      interface INS_block
      module procedure iINS_block
      module procedure sINS_block
      module procedure dINS_block
      module procedure cINS_block
      module procedure zINS_block
      end interface
      
      
      interface INS_bl_entr
      module procedure iINS_bl_entr
      module procedure sINS_bl_entr
      module procedure dINS_bl_entr
      module procedure cINS_bl_entr
      module procedure zINS_bl_entr
      end interface

      interface INS_varblock
      module procedure iINS_varblock
      module procedure sINS_varblock
      module procedure dINS_varblock
      module procedure cINS_varblock
      module procedure zINS_varblock
      end interface

      interface INS_varbl_entr
      module procedure iINS_varbl_entr
      module procedure sINS_varbl_entr
      module procedure dINS_varbl_entr
      module procedure cINS_varbl_entr
      module procedure zINS_varbl_entr
      end interface


      contains

! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "integerINS_ROUTINER_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "singlePrecisionINS_ROUTINER_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "doublePrecisionINS_ROUTINER_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "singleComplexINS_ROUTINER_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
#include "doubleComplexINS_ROUTINER_source.F"
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      end module mod_INS_ROUTINER

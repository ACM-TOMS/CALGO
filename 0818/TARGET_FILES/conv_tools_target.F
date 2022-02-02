      module mod_conv_tools
      use blas_sparse_namedconstants

      interface b_up_order
      module procedure ib_up_order
      module procedure sb_up_order
      module procedure db_up_order
      module procedure cb_up_order
      module procedure zb_up_order
      end interface

      interface A_row_col
      module procedure iA_row_col
      module procedure sA_row_col
      module procedure dA_row_col
      module procedure cA_row_col
      module procedure zA_row_col
      end interface

      interface detect_diag
      module procedure idetect_diag
      module procedure sdetect_diag
      module procedure ddetect_diag
      module procedure cdetect_diag
      module procedure zdetect_diag
      end interface

      interface Ab_row_col
      module procedure iAb_row_col
      module procedure sAb_row_col
      module procedure dAb_row_col
      module procedure cAb_row_col
      module procedure zAb_row_col
      end interface

      interface detect_bdiag
      module procedure idetect_bdiag
      module procedure sdetect_bdiag
      module procedure ddetect_bdiag
      module procedure cdetect_bdiag
      module procedure zdetect_bdiag
      end interface

      interface pre_usconv_coo2csr
      module procedure ipre_usconv_coo2csr
      module procedure spre_usconv_coo2csr
      module procedure dpre_usconv_coo2csr
      module procedure cpre_usconv_coo2csr
      module procedure zpre_usconv_coo2csr
      end interface

      interface pre_usconv_coo2csc
      module procedure ipre_usconv_coo2csc
      module procedure spre_usconv_coo2csc
      module procedure dpre_usconv_coo2csc
      module procedure cpre_usconv_coo2csc
      module procedure zpre_usconv_coo2csc
      end interface

      interface pre_usconv_bco2bsc
      module procedure ipre_usconv_bco2bsc
      module procedure spre_usconv_bco2bsc
      module procedure dpre_usconv_bco2bsc
      module procedure cpre_usconv_bco2bsc
      module procedure zpre_usconv_bco2bsc
      end interface

      interface pre_usconv_bco2bsr
      module procedure ipre_usconv_bco2bsr
      module procedure spre_usconv_bco2bsr
      module procedure dpre_usconv_bco2bsr
      module procedure cpre_usconv_bco2bsr
      module procedure zpre_usconv_bco2bsr
      end interface

      interface pre_usconv_coo2dia
      module procedure ipre_usconv_coo2dia
      module procedure spre_usconv_coo2dia
      module procedure dpre_usconv_coo2dia
      module procedure cpre_usconv_coo2dia
      module procedure zpre_usconv_coo2dia
      end interface

      interface pre_usconv_dia2coo
      module procedure ipre_usconv_dia2coo
      module procedure spre_usconv_dia2coo
      module procedure dpre_usconv_dia2coo
      module procedure cpre_usconv_dia2coo
      module procedure zpre_usconv_dia2coo
      end interface

      interface pre_usconv_bco2bdi
      module procedure ipre_usconv_bco2bdi
      module procedure spre_usconv_bco2bdi
      module procedure dpre_usconv_bco2bdi
      module procedure cpre_usconv_bco2bdi
      module procedure zpre_usconv_bco2bdi
      end interface

      interface pre_usconv_bdi2bco
      module procedure ipre_usconv_bdi2bco
      module procedure spre_usconv_bdi2bco
      module procedure dpre_usconv_bdi2bco
      module procedure cpre_usconv_bdi2bco
      module procedure zpre_usconv_bdi2bco
      end interface


      contains

      subroutine  up_order(INDX,RES_INDX)
      implicit none
      integer,pointer,dimension(:) ::INDX
      integer,dimension(:),allocatable ::tes
      integer,pointer,dimension(:) ::RES_INDX 
      integer,dimension(1)::c
      integer ::i,s
      integer :: dummy
      intrinsic maxval
      intrinsic minloc
      
      s=size(INDX)
      allocate(tes(s))
      tes=INDX
      dummy = maxval(tes)+1
      do i=1,s
         c=minloc(tes)
         RES_INDX(i)=c(1)
         tes(c(1))=dummy
      end do
      deallocate(tes)
      end subroutine up_order

      function counter(INDX,value)
      implicit none
      integer ,pointer,dimension(:)::INDX
      integer ,intent(in)::value
      integer ::counter,s,j,k
      s=size(INDX)
      k=0
      do j=1,s
         if(INDX(j)==value) then
            k=k+1
         end if
      end do
      counter=k
      end function  counter

      subroutine PNTR(PNTRB,PNTRE,M_K,INDX)
      implicit none
      integer ,pointer,dimension(:)::PNTRB,PNTRE
      integer ,pointer,dimension(:) :: INDX
      integer ,intent(in) :: M_K
      integer ::j,s 
      
      s=size(INDX)
      PNTRB(1)=1
      PNTRE(M_K)=s+1
      do j=2,M_K
         PNTRB(j)=PNTRB(j-1)+counter(INDX,j-1)
         PNTRE(j-1)=PNTRB(j)
      end do
      end subroutine PNTR

      subroutine final_order(JNDX,final_indx,row_subdv)
      implicit none
      integer,pointer,dimension(:)::JNDX,row_subdv
      integer,pointer,dimension(:)::final_indx
      integer,pointer,dimension(:) :: test_int,test_ind
      integer ::d,k,s,i
      d=1
      s=size(row_subdv)
      do i=1,s
         if(row_subdv(i)>0) then
            allocate(test_int(row_subdv(i)))
            allocate(test_ind(row_subdv(i)))
            test_int=JNDX((/(i,i=d,d+row_subdv(i)-1,1)/))
            call up_order(test_int,test_ind)
            do k=1,row_subdv(i)
               final_indx(d+k-1)=test_ind(k)+d-1
            end do
            deallocate(test_int)
            deallocate(test_ind)
         end if 
         d=d+row_subdv(i) 
      end do
      end subroutine final_order

      subroutine PNTR_INV(PNTRE,INDX)
      implicit none
      integer,pointer ,dimension(:)::PNTRE
      integer,pointer ,dimension(:)::INDX
      integer :: i,j,s
   
      s=size(PNTRE)
      do j=1,PNTRE(1)-1
         INDX(j)=1
      end do
      do i=1,s-1
         if(PNTRE(i).ne.PNTRE(i+1)) then
            do j=PNTRE(i),PNTRE(i+1)-1
               INDX(j)=i+1
            end do
         end if
      end do  
      end subroutine PNTR_INV

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "integerconv_tools_source.F"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#include "singlePrecisionconv_tools_source.F"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#include "doublePrecisionconv_tools_source.F"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#include "singleComplexconv_tools_source.F"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
#include "doubleComplexconv_tools_source.F"
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      end module mod_conv_tools
